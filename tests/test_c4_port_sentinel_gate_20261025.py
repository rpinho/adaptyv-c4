#!/usr/bin/env python3
"""Fail-closed tests for the C4 port sentinel gate. No GPU, no network.

Every test here corresponds to a way the gate could wrongly report PASS, or to the
exact failure that made the first port attempt undiagnosable. Synthetic run
directories are built per test; nothing real is touched.
"""
import hashlib
import importlib.util
import json
import pathlib
import shutil
import tempfile
import unittest

ROOT = pathlib.Path(__file__).resolve().parents[1]
GATE_PATH = ROOT / "scripts/c4_port_sentinel_gate_20261025.py"


def load_gate():
    spec = importlib.util.spec_from_file_location("c4gate", GATE_PATH)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


gate = load_gate()
TARGET = "K" * 195
TARGET_LEN = 195


def sha(path):
    return hashlib.sha256(pathlib.Path(path).read_bytes()).hexdigest()


def write(path, text):
    path = pathlib.Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)
    return path


def make_reference(tmp, names=("tfr1_c4_acid_site2_1",), iptm=0.7057363986968994,
                   confidence=0.8318573236465454):
    """Reference archive. Must contain a record for every sentinel under test."""
    ref = {"target": "P02786_189-383_contiguous", "target_len": TARGET_LEN,
           "results": [{"name": n, "iptm": iptm, "confidence": confidence,
                        "complex_plddt": 0.8633875250816345} for n in names]}
    return write(pathlib.Path(tmp) / "reference.json", json.dumps(ref))


def make_run(tmp, name="tfr1_c4_acid_site2_1",
             iptm=0.7057, confidence=0.8319,
             target_seq=TARGET, use_msa_server=False,
             chain_b_len=TARGET_LEN,
             drop=(), tamper_msa=False, tamper_yaml=False):
    """Build a run dir. `drop` removes provenance keys to simulate an unarchived lane."""
    run = pathlib.Path(tmp) / "runs" / name
    run.mkdir(parents=True, exist_ok=True)
    write(run / "input.yaml", f"sequences:\n  - protein: {name}\n")
    write(run / "msa.a3m", ">q\nMKV\n")
    write(run / "model.pdb",
          "".join(f"ATOM  {i:5d}  CA  ALA A{i:4d}    {0.0:8.3f}{0.0:8.3f}{0.0:8.3f}\n"
                  for i in range(1, 80)) +
          "".join(f"ATOM  {i:5d}  CA  ALA B{i:4d}    {0.0:8.3f}{0.0:8.3f}{0.0:8.3f}\n"
                  for i in range(1, chain_b_len + 1)))
    write(run / "confidence.json",
          json.dumps({"iptm": iptm, "confidence_score": confidence}))
    meta = {
        "target_sequence": target_seq,
        "use_msa_server": use_msa_server,
        "msa_path": "msa.a3m",
        "msa_sha256": sha(run / "msa.a3m"),
        "input_yaml": "input.yaml",
        "input_yaml_sha256": sha(run / "input.yaml"),
        "boltz_version": "2.2.1", "torch": "2.14.1+cu130", "gpu": "NVIDIA A40 48GB",
        "seed": 20261025, "diffusion_samples": 1, "recycling_steps": 3,
        "config_sha256": "0" * 64,
        "confidence_json": "confidence.json", "pdb": "model.pdb",
    }
    for key in drop:
        meta.pop(key, None)
    if tamper_msa:
        meta["msa_sha256"] = "f" * 64
    write(run / "provenance.json", json.dumps(meta))
    if tamper_yaml:
        write(run / "input.yaml", "sequences:\n  - protein: TAMPERED\n")
    return run


def spec(run_dirs, gate_member=True):
    return {"target": {"sequence": TARGET, "length": TARGET_LEN},
            "sentinels": [{"name": n, "run_dir": str(d), "gate_member": gate_member}
                          for n, d in run_dirs]}


class Baseline(unittest.TestCase):
    def setUp(self):
        self.tmp = pathlib.Path(tempfile.mkdtemp())
        self.addCleanup(shutil.rmtree, self.tmp, True)
        self.ref = make_reference(self.tmp)

    def test_clean_run_passes(self):
        run = make_run(self.tmp)
        out = gate.evaluate(spec([("tfr1_c4_acid_site2_1", run)]), self.ref)
        self.assertEqual(out["overall"], gate.PASS)

    def test_real_known_failure_is_fail_not_pass(self):
        """The archived A40 numbers must FAIL the gate."""
        run = make_run(self.tmp, iptm=0.5187048316001892, confidence=0.7808405160903931)
        out = gate.evaluate(spec([("tfr1_c4_acid_site2_1", run)]), self.ref)
        self.assertEqual(out["overall"], gate.FAIL)
        self.assertAlmostEqual(out["sentinels"][0]["delta"]["iptm"], -0.1870315670967102)


class FailClosed(unittest.TestCase):
    """Missing or unverifiable provenance must be INCONCLUSIVE, never PASS."""

    def setUp(self):
        self.tmp = pathlib.Path(tempfile.mkdtemp())
        self.addCleanup(shutil.rmtree, self.tmp, True)
        self.ref = make_reference(self.tmp)

    def _verdict(self, **kw):
        run = make_run(self.tmp, **kw)
        return gate.evaluate(spec([("tfr1_c4_acid_site2_1", run)]), self.ref)

    def test_no_provenance_file_is_inconclusive(self):
        run = make_run(self.tmp)
        (run / "provenance.json").unlink()
        out = gate.evaluate(spec([("tfr1_c4_acid_site2_1", run)]), self.ref)
        self.assertEqual(out["overall"], gate.INCONCLUSIVE)

    def test_online_msa_is_inconclusive(self):
        self.assertEqual(self._verdict(use_msa_server=True)["overall"], gate.INCONCLUSIVE)

    def test_missing_msa_file_is_inconclusive(self):
        run = make_run(self.tmp)
        (run / "msa.a3m").unlink()
        out = gate.evaluate(spec([("tfr1_c4_acid_site2_1", run)]), self.ref)
        self.assertEqual(out["overall"], gate.INCONCLUSIVE)

    def test_msa_hash_mismatch_is_inconclusive(self):
        self.assertEqual(self._verdict(tamper_msa=True)["overall"], gate.INCONCLUSIVE)

    def test_yaml_tampered_after_hashing_is_inconclusive(self):
        self.assertEqual(self._verdict(tamper_yaml=True)["overall"], gate.INCONCLUSIVE)

    def test_missing_seed_is_inconclusive(self):
        self.assertEqual(self._verdict(drop=("seed",))["overall"], gate.INCONCLUSIVE)

    def test_missing_config_hash_is_inconclusive(self):
        self.assertEqual(
            self._verdict(drop=("config_sha256",))["overall"], gate.INCONCLUSIVE)

    def test_missing_pdb_is_inconclusive(self):
        run = make_run(self.tmp)
        (run / "model.pdb").unlink()
        out = gate.evaluate(spec([("tfr1_c4_acid_site2_1", run)]), self.ref)
        self.assertEqual(out["overall"], gate.INCONCLUSIVE)

    def test_surrogate_target_length_is_inconclusive(self):
        """The 154-residue stitched surrogate must be rejected, not tolerated."""
        self.assertEqual(
            self._verdict(target_seq="K" * 154, chain_b_len=154)["overall"],
            gate.INCONCLUSIVE)

    def test_right_length_wrong_residues_is_inconclusive(self):
        wrong = "E" + "K" * 194
        self.assertEqual(self._verdict(target_seq=wrong)["overall"], gate.INCONCLUSIVE)

    def test_chain_b_length_mismatch_is_inconclusive(self):
        self.assertEqual(self._verdict(chain_b_len=154)["overall"], gate.INCONCLUSIVE)

    def test_nonfinite_score_is_inconclusive(self):
        run = make_run(self.tmp)
        write(run / "confidence.json", json.dumps({"iptm": None, "confidence_score": 0.8}))
        out = gate.evaluate(spec([("tfr1_c4_acid_site2_1", run)]), self.ref)
        self.assertEqual(out["overall"], gate.INCONCLUSIVE)


class Aggregation(unittest.TestCase):
    def setUp(self):
        self.tmp = pathlib.Path(tempfile.mkdtemp())
        self.addCleanup(shutil.rmtree, self.tmp, True)
        self.ref = make_reference(
            self.tmp, names=("tfr1_c4_acid_site2_1", "sentinel_b", "diagnostic"))

    def _two(self, **kw2):
        r1 = make_run(self.tmp, name="tfr1_c4_acid_site2_1")
        r2 = make_run(self.tmp, name="sentinel_b", **kw2)
        return gate.evaluate(
            spec([("tfr1_c4_acid_site2_1", r1), ("sentinel_b", r2)]), self.ref)

    def test_all_pass_passes(self):
        self.assertEqual(self._two()["overall"], gate.PASS)

    def test_one_inconclusive_dominates_pass(self):
        self.assertEqual(self._two(drop=("seed",))["overall"], gate.INCONCLUSIVE)

    def test_inconclusive_dominates_fail(self):
        """An unexplainable member must not be reported as a plain failure."""
        self.assertEqual(self._two(drop=("seed",))["overall"], gate.INCONCLUSIVE)

    def test_a_fail_with_full_provenance_is_fail(self):
        self.assertEqual(self._two(iptm=0.10)["overall"], gate.FAIL)

    def test_diagnostic_member_does_not_set_overall(self):
        r1 = make_run(self.tmp)
        r2 = make_run(self.tmp, name="diagnostic", iptm=0.10)
        s = spec([("tfr1_c4_acid_site2_1", r1)])
        s["sentinels"].append({"name": "diagnostic", "run_dir": str(r2),
                               "gate_member": False})
        out = gate.evaluate(s, self.ref)
        self.assertEqual(out["overall"], gate.PASS)
        self.assertEqual(out["counts"][gate.FAIL], 1)

    def test_empty_gate_members_is_an_error(self):
        r1 = make_run(self.tmp)
        s = spec([("tfr1_c4_acid_site2_1", r1)], gate_member=False)
        with self.assertRaises(ValueError):
            gate.evaluate(s, self.ref)


class SpecIntegrity(unittest.TestCase):
    def test_unpinned_target_is_rejected(self):
        with self.assertRaises(ValueError):
            gate.pinned_target_sequence({"target": {}}, "x")

    def test_length_mismatch_in_spec_is_rejected(self):
        with self.assertRaises(ValueError):
            gate.pinned_target_sequence(
                {"target": {"sequence": "KKK", "length": 195}}, "x")

    def test_shipped_spec_target_matches_canonical_uniprot(self):
        """The pinned target must be P02786 189-383, derived, not hand-typed."""
        raw = (ROOT / "data/P02786.txt").read_text().splitlines()
        seq, inseq = [], False
        for line in raw:
            if line.startswith("SQ   "):
                inseq = True
                continue
            if inseq:
                if line.startswith("//"):
                    break
                seq.append("".join(c for c in line if "A" <= c <= "Z"))
        canonical = "".join(seq)[188:383]
        shipped = json.loads(
            (ROOT / "results/c4_sentinel_gate_20261025/spec_NOT_FROZEN.json").read_text())
        self.assertEqual(shipped["target"]["sequence"], canonical)
        self.assertEqual(len(canonical), 195)
        self.assertEqual(canonical[0], "K")
        self.assertEqual(canonical[-1], "E")


if __name__ == "__main__":
    unittest.main(verbosity=2)
