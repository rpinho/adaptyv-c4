#!/usr/bin/env python3
"""Tests for the apical L6 recomputation: chain parsing, fail-closed behavior, sign convention."""
import importlib.util
import pathlib
import sys
import tempfile

ROOT = pathlib.Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "apical", ROOT / "scripts/c4_linkage_apical_recompute_20261025.py")
M = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(M)


def write_pdb(path, chains):
    """chains = {'A': n_ca, 'B': n_ca}; writes ATOM CA lines only."""
    with open(path, "w") as fh:
        for chain, n in chains.items():
            for i in range(1, n + 1):
                fh.write(f"ATOM  {i:5d}  CA  ALA {chain}{i:4d}    "
                         f"{float(i):8.3f}{0.0:8.3f}{0.0:8.3f}  1.00  0.00           C\n")
        fh.write("END\n")


def test_chain_counts_basic():
    with tempfile.TemporaryDirectory() as td:
        p = pathlib.Path(td) / "x.pdb"
        write_pdb(p, {"A": 195, "B": 70})
        assert M.chain_ca_counts(p) == {"A": 195, "B": 70}


def test_chain_counts_ignores_duplicate_atoms():
    with tempfile.TemporaryDirectory() as td:
        p = pathlib.Path(td) / "x.pdb"
        write_pdb(p, {"A": 3, "B": 2})
        text = p.read_text()
        p.write_text(text.replace("END\n", "") + text.splitlines()[0] + "\nEND\n")
        assert M.chain_ca_counts(p) == {"A": 3, "B": 2}


def test_chain_counts_empty():
    with tempfile.TemporaryDirectory() as td:
        p = pathlib.Path(td) / "x.pdb"
        p.write_text("END\n")
        assert M.chain_ca_counts(p) == {}


def test_sha256_stable():
    with tempfile.TemporaryDirectory() as td:
        p = pathlib.Path(td) / "x.pdb"
        p.write_text("hello\n")
        assert M.sha256(p) == M.sha256(p)
        assert len(M.sha256(p)) == 64


def test_expected_constants():
    assert M.EXPECTED_RECEPTOR_CA == 195
    assert M.RECEPTOR_CHAIN == "B"
    assert M.BINDER_CHAIN == "A"


def test_sign_controls_pass():
    ctl = M.A.controls()
    assert ctl["all_ok"], ctl
    assert ctl["acid_on_L"] > 0
    assert ctl["acid_off_L"] < 0


def test_real_apical_pdb_is_195():
    p = M.APICAL_PDBS / "tfr1_c4_acid_site2_1_complex.pdb"
    assert p.exists(), p
    counts = M.chain_ca_counts(p)
    assert counts == {"A": 70, "B": 195}, counts


def test_real_surrogate_pdb_is_154():
    p = ROOT / "results/c4_rebuild_20261024/remote_pdbs/tfr1_c4_acid_site2_1_complex.pdb"
    assert p.exists(), p
    counts = M.chain_ca_counts(p)
    assert counts == {"A": 70, "B": 154}, counts


def test_artifact_is_342_and_has_no_site1_41():
    j = __import__("json").loads(M.ARTIFACT.read_text())
    assert j["n"] == 342, j["n"]
    names = {r["name"] for r in j["records"]}
    assert "tfr1_c4_acid_site1_41" not in names


def test_verify_mode_reproduces_sample():
    """Verify a 3-row sample against the saved artifact."""
    import subprocess
    r = subprocess.run([sys.executable, str(ROOT / "scripts/c4_linkage_apical_recompute_20261025.py"),
                        "--limit", "3"], capture_output=True, text=True, timeout=600)
    assert r.returncode == 0, r.stderr[-2000:]
    assert "VERIFY OK" in r.stdout, r.stdout[-2000:]


if __name__ == "__main__":
    failed = 0
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and callable(fn):
            try:
                fn()
                print(f"PASS {name}")
            except Exception as e:
                failed += 1
                print(f"FAIL {name}: {e}")
    print(f"{'FAILED' if failed else 'OK'} ({failed} failures)")
    sys.exit(1 if failed else 0)
