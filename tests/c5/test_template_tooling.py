#!/usr/bin/env python3
"""Fail-closed tests for the C5 target-independent template tooling.

CPU only. Run with:  CUDA_VISIBLE_DEVICES="" /usr/bin/python3 tests/c5/test_template_tooling.py
"""
import json
import pathlib
import sys
import tempfile

ROOT = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "c5_prep"))
import template_lib as T  # noqa: E402

FIXTURE = ROOT / "data/6x18.cif"
RESULTS = []


def check(name, fn):
    try:
        fn()
        RESULTS.append((name, True, ""))
    except AssertionError as exc:
        RESULTS.append((name, False, f"assertion: {exc}"))
    except Exception as exc:  # noqa: BLE001
        RESULTS.append((name, False, f"{type(exc).__name__}: {exc}"))


def expect_fail(fn, label):
    try:
        fn()
    except T.TemplateError:
        return
    raise AssertionError(f"{label} did not fail closed")


def t_fixture_parses():
    atoms, meta = T.parse_structure(FIXTURE)
    assert len(atoms) > 1000, len(atoms)
    assert meta["_exptl.method"].upper() == "ELECTRON MICROSCOPY", meta["_exptl.method"]
    ch = T.chains(atoms)
    # 6X18 has six polymer chains (A-F); G-K are waters/ions.
    assert "F" in ch and len(ch) == 6, sorted(ch)


def t_manifest_has_required_fields():
    m = T.build_manifest([FIXTURE], snapshot={"source": "local fixture 6x18.cif"})
    assert not m["rejected"], m["rejected"]
    assert len(m["templates"]) == 6, len(m["templates"])
    t = next(x for x in m["templates"] if x["chain"] == "F")
    for k in ("pdb_id", "chain", "method", "resolution_a", "release_date",
              "source_sha256", "coverage", "sequence_sha256"):
        assert k in t, k
    assert t["pdb_id"] == "6X18", t["pdb_id"]
    assert abs(t["resolution_a"] - 2.1) < 1e-6, t["resolution_a"]
    assert t["method"] == "ELECTRON MICROSCOPY"
    assert len(t["sequence"]) == t["length"]
    assert len(t["source_sha256"]) == 64


def t_predicted_filename_rejected():
    with tempfile.TemporaryDirectory() as td:
        p = pathlib.Path(td) / "AF-P02786-F1-model_v4.cif"
        p.write_text(FIXTURE.read_text())
        expect_fail(lambda: T.build_manifest([p]), "predicted filename")


def t_missing_method_rejected():
    with tempfile.TemporaryDirectory() as td:
        p = pathlib.Path(td) / "model_1.pdb"
        p.write_text("ATOM      1  CA  ALA A   1      1.000   2.000   3.000  1.00  0.00           C\n")
        expect_fail(lambda: T.build_manifest([p]), "prediction markers")


def t_plddt_atoms_rejected():
    with tempfile.TemporaryDirectory() as td:
        p = pathlib.Path(td) / "exp.pdb"
        p.write_text(
            "ATOM      1  N   ALA A   1      1.000   2.000   3.000  1.00 90.00           N\n"
            "ATOM      2  CA  ALA A   1      2.000   2.000   3.000  1.00 90.00           C\n"
            "ATOM      3  C   ALA A   1      3.000   2.000   3.000  1.00 90.00           C\n"
            "ATOM      4  O   ALA A   1      4.000   2.000   3.000  1.00 90.00           O\n"
            "ATOM      5 pLDDT ALA A   1      2.000   2.000   3.000  1.00 90.00           X\n")
        expect_fail(lambda: T.build_manifest([p]), "pLDDT atoms")


def t_unsupported_extension_rejected():
    with tempfile.TemporaryDirectory() as td:
        p = pathlib.Path(td) / "thing.xyz"
        p.write_text("x")
        expect_fail(lambda: T.build_manifest([p]), "unsupported extension")


def t_missing_file_rejected():
    expect_fail(lambda: T.build_manifest([ROOT / "data/does-not-exist.cif"]), "missing file")


def t_target_validator_accepts_valid():
    rec = T.validate_target_input({"name": "demo", "sequence": "MKTAYIAKQRQI"})
    assert rec["length"] == 12 and rec["numbering_end"] == 12
    rec2 = T.validate_target_input({"name": "n", "sequence": "mktayiakqrqi",
                                    "numbering_start": 100})
    assert rec2["numbering_start"] == 100 and rec2["numbering_end"] == 111
    assert rec2["sequence"] == "MKTAYIAKQRQI"


def t_target_validator_rejects():
    for bad, label in (
        ({}, "missing name"),
        ({"name": "x"}, "missing sequence"),
        ({"name": "x", "sequence": ""}, "empty sequence"),
        ({"name": "x", "sequence": "MKTA*YIAK"}, "non-AA"),
        ({"name": "x", "sequence": "MKT"}, "too short"),
        ({"name": "x", "sequence": "M" * 300}, "too long"),
        ({"name": "x", "sequence": "MKTAYIAKQRQI", "numbering_start": -5}, "bad numbering"),
    ):
        expect_fail(lambda b=bad: T.validate_target_input(b), label)


def t_alignment_identical():
    seq = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQ"
    m = T.align_sequences(seq, seq)
    assert m["identity"] == 1.0, m
    assert m["aligned_pairs"] == len(seq)
    assert m["query_coverage"] == 1.0 and m["template_coverage"] == 1.0


def t_alignment_mutation_detected():
    a = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQ"
    b = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVA"
    m = T.align_sequences(a, b)
    assert m["identity"] < 1.0, m
    assert m["aligned_pairs"] == len(a), m


def t_alignment_handles_indel():
    a = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQ"
    b = "MKTAYIAKQ" + "GG" + "RQISFVKSHFSRQLEERLGLIEVQ"
    m = T.align_sequences(a, b)
    assert m["aligned_pairs"] >= len(a) - 1, m
    assert m["identity"] > 0.9, m


def t_superposition_perfect_and_noisy():
    q = [(0, 0, 0), (3.8, 0, 0), (3.8, 3.8, 0), (0, 3.8, 0), (0, 0, 3.8)]
    assert T.superimpose_ca(q, q)["ca_rmsd"] == 0.0
    shifted = [(x + 10, y - 4, z + 2) for x, y, z in q]
    assert T.superimpose_ca(q, shifted)["ca_rmsd"] == 0.0
    noisy = [(x + (0.4 if i == 0 else 0), y, z) for i, (x, y, z) in enumerate(q)]
    assert T.superimpose_ca(q, noisy)["ca_rmsd"] > 0.1
    expect_fail(lambda: T.superimpose_ca(q, q[:4]), "unequal CA lists")


def t_ranking_is_deterministic_and_sensitive():
    entries = [
        {"pdb_id": "AAAA", "chain": "A", "identity": 0.9, "similarity": 0.9,
         "query_coverage": 0.9, "template_coverage": 0.9, "resolution_a": 1.5},
        {"pdb_id": "BBBB", "chain": "B", "identity": 0.5, "similarity": 0.6,
         "query_coverage": 0.9, "template_coverage": 0.9, "resolution_a": 1.0},
        {"pdb_id": "CCCC", "chain": "C", "identity": 0.9, "similarity": 0.9,
         "query_coverage": 0.9, "template_coverage": 0.9, "resolution_a": None},
    ]
    r1 = T.rank_templates(entries)
    r2 = T.rank_templates(list(reversed(entries)))
    assert [x["pdb_id"] for x in r1] == [x["pdb_id"] for x in r2], "order unstable"
    assert r1[0]["pdb_id"] == "AAAA", [x["pdb_id"] for x in r1]
    assert [x["rank"] for x in r1] == [1, 2, 3]
    # Sequence agreement dominates; missing resolution is a penalty, not a veto.
    assert r1[-1]["pdb_id"] == "BBBB", [x["pdb_id"] for x in r1]


def t_ranking_missing_resolution_last_among_ties():
    base = {"similarity": 0.9, "query_coverage": 0.9, "template_coverage": 0.9,
            "identity": 0.8}
    entries = [dict(base, pdb_id="ZZZZ", chain="Z", resolution_a=None),
               dict(base, pdb_id="YYYY", chain="Y", resolution_a=2.9)]
    order = [x["pdb_id"] for x in T.rank_templates(entries)]
    assert order == ["YYYY", "ZZZZ"], order


def t_ranking_rejects_unknown_weight():
    expect_fail(lambda: T.rank_templates([], weights={"magic": 1.0}), "unknown weight")


def t_manifest_writer_validates():
    m = T.build_manifest([FIXTURE])
    with tempfile.TemporaryDirectory() as td:
        p = T.write_manifest(m, pathlib.Path(td) / "m.json")
        assert json.loads(p.read_text())["templates"]
        broken = json.loads(json.dumps(m))
        broken["templates"][0]["sequence_sha256"] = "0" * 64
        expect_fail(lambda: T.write_manifest(broken, pathlib.Path(td) / "x.json"),
                    "hash mismatch")
        for key in ("generated_at", "tool", "templates", "rejected", "snapshot"):
            b = json.loads(json.dumps(m))
            b.pop(key)
            expect_fail(lambda bb=b: T.write_manifest(bb, pathlib.Path(td) / "y.json"),
                        f"missing {key}")


def t_extract_chain_cli():
    import subprocess
    with tempfile.TemporaryDirectory() as td:
        out = pathlib.Path(td) / "F.pdb"
        r = subprocess.run([sys.executable, str(ROOT / "c5_prep/extract_chain.py"),
                            str(FIXTURE), "--chain", "F", "--out", str(out)],
                           capture_output=True, text=True)
        assert r.returncode == 0, r.stderr
        atoms, _ = T.parse_structure(out)
        ch = T.chains(atoms)
        assert list(ch) == ["F"], list(ch)
        assert len(ch["F"]) == 384, len(ch["F"])
        # predicted input must be refused without the escape hatch
        bad = pathlib.Path(td) / "AF-demo.cif"
        bad.write_text(FIXTURE.read_text())
        r2 = subprocess.run([sys.executable, str(ROOT / "c5_prep/extract_chain.py"),
                             str(bad), "--chain", "F", "--out", str(pathlib.Path(td) / "b.pdb")],
                            capture_output=True, text=True)
        assert r2.returncode != 0, "predicted input was accepted"


TESTS = [(k[2:], v) for k, v in sorted(globals().items()) if k.startswith("t_") and callable(v)]
for name, fn in TESTS:
    check(name, fn)

ok = sum(1 for _, p, _ in RESULTS if p)
for name, passed, err in RESULTS:
    print(f"{'PASS' if passed else 'FAIL'} {name}{'' if passed else '  <- ' + err}")
print(f"\n{ok}/{len(RESULTS)} passed")
sys.exit(0 if ok == len(RESULTS) else 1)
