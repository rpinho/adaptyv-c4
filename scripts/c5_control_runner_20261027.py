#!/usr/bin/env python3
"""
C5 POSITIVE-CONTROL RUNNER — validates the scoring pipeline discriminates G12D from WT.

Runs ON THE POD. Three predictions:
  1. TCR + G12D pMHC  (should bind — KD 63 µM experimentally)
  2. TCR + WT pMHC    (should NOT bind — NB experimentally)
  3. TCR + peptide-free HLA+b2m (should NOT bind)

Pass criterion: ipTM(G12D) > ipTM(WT) AND ipTM(G12D) > ipTM(pepfree).
If this fails, the pipeline cannot discriminate and must not be used for design ranking.

Based on the C4 rescore runner but simplified for 3 fixed inputs.
"""

import argparse
import glob
import hashlib
import json
import os
import pathlib
import re
import subprocess
import sys
import time


SETTINGS = {
    "predictor": "colabfold_alphafold2_multimer_v3",
    "model_type": "alphafold2_multimer_v3",
    "num_recycle": 3,
    "num_models": 5,
    "rank": "multimer",
    "msa_mode": "mmseqs2_uniref_env",
}

EXPECTED_FASTAS = {
    "control_tcr_vs_g12d":    "b72c3e6f7f1d743240e3101235fdf26085b496b96414e81c78d90ee2bea54f7d",
    "control_tcr_vs_wt":      "a4928b133669261e604539a541beff0efc7090ae8c9c6e81c265787415a47894",
    "control_tcr_vs_pepfree": "389fcfe2a0673abd94894e8cf98e6c6086193c48f53805f2d684671449662531",
}


def sha256_str(s):
    return hashlib.sha256(s.encode()).hexdigest()


def sha256_file(p):
    return hashlib.sha256(pathlib.Path(p).read_bytes()).hexdigest()


def validate_fasta(fasta_path):
    """Validate a control FASTA: check sequence hash matches expected."""
    name = fasta_path.stem
    if name not in EXPECTED_FASTAS:
        return False, f"unexpected FASTA: {name}"
    with open(fasta_path) as f:
        lines = f.readlines()
    seq_line = lines[1].strip()
    actual_hash = sha256_str(seq_line)
    if actual_hash != EXPECTED_FASTAS[name]:
        return False, f"hash mismatch: {actual_hash} != {EXPECTED_FASTAS[name]}"
    chains = seq_line.split(":")
    total_aa = sum(len(c) for c in chains)
    return True, f"OK: {len(chains)} chains, {total_aa} aa, hash verified"


def extract_scores(cf_out_dir):
    score_files = sorted(glob.glob(str(cf_out_dir / "*scores*rank*.json")))
    if not score_files:
        score_files = sorted(glob.glob(str(cf_out_dir / "*scores*.json")))
    if not score_files:
        return None

    ranked = []
    for sf in score_files:
        m = re.search(r'rank_(\d+)', pathlib.Path(sf).name)
        rank = int(m.group(1)) if m else 999
        ranked.append((rank, sf))
    ranked.sort()

    all_model_scores = []
    for rank, sf in ranked:
        with open(sf) as f:
            s = json.load(f)
        iptm = s.get("iptm")
        ptm = s.get("ptm")
        plddt = s.get("plddt")
        if isinstance(plddt, list):
            plddt = sum(plddt) / len(plddt) if plddt else None
        m = re.search(r'(model_\d+)', pathlib.Path(sf).name)
        model_name = m.group(1) if m else f"rank_{rank}"
        all_model_scores.append({
            "model": model_name, "rank": rank,
            "iptm": round(iptm, 4) if iptm is not None else None,
            "ptm": round(ptm, 4) if ptm is not None else None,
            "mean_plddt": round(plddt, 2) if plddt is not None else None,
        })

    best_rank, best_file = ranked[0]
    with open(best_file) as f:
        scores = json.load(f)
    iptm = scores.get("iptm")
    ptm = scores.get("ptm")
    plddt = scores.get("plddt")
    if isinstance(plddt, list):
        plddt = sum(plddt) / len(plddt) if plddt else None

    return {
        "iptm": round(iptm, 4) if iptm is not None else None,
        "ptm": round(ptm, 4) if ptm is not None else None,
        "mean_plddt": round(plddt, 2) if plddt is not None else None,
        "best_model": re.search(r'(model_\d+)', pathlib.Path(best_file).name).group(1),
        "n_models": len(ranked),
        "all_models": all_model_scores,
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("fasta_dir", help="Directory containing the 3 control FASTAs")
    ap.add_argument("out_dir", help="Output directory for results")
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    fasta_dir = pathlib.Path(args.fasta_dir)
    out_root = pathlib.Path(args.out_dir)
    out_root.mkdir(parents=True, exist_ok=True)

    fastas = sorted(fasta_dir.glob("control_tcr_vs_*.fasta"))
    if len(fastas) != 3:
        print(f"ERROR: expected 3 control FASTAs, found {len(fastas)}", file=sys.stderr)
        return 1

    print(f"C5 positive-control validation: {len(fastas)} predictions")
    print(f"Settings: {json.dumps(SETTINGS)}")
    print()

    for fa in fastas:
        ok, msg = validate_fasta(fa)
        print(f"  {fa.name}: {msg}")
        if not ok:
            print(f"  ABORT: validation failed", file=sys.stderr)
            return 2

    if args.dry_run:
        print("\nDry run — exiting.")
        return 0

    results = {}
    for fa in fastas:
        name = fa.stem

        design_out = out_root / name
        design_out.mkdir(parents=True, exist_ok=True)
        work = design_out / "work"
        work.mkdir(exist_ok=True)
        cf_out = work / "out"
        cf_out.mkdir(exist_ok=True)

        # Check for existing scores (checkpoint)
        existing_scores = design_out / "scores.json"
        if existing_scores.exists():
            with open(existing_scores) as f:
                s = json.load(f)
            if s.get("iptm") is not None:
                print(f"SKIP {name}: already scored (ipTM={s['iptm']:.4f})")
                results[name] = s
                continue

        work_fa = work / fa.name
        work_fa.write_text(fa.read_text())

        cmd = [
            "colabfold_batch", str(work_fa), str(cf_out),
            "--model-type", SETTINGS["model_type"],
            "--num-recycle", str(SETTINGS["num_recycle"]),
            "--num-models", str(SETTINGS["num_models"]),
            "--rank", SETTINGS["rank"],
            "--msa-mode", SETTINGS["msa_mode"],
        ]

        print(f"\nFOLD {name}...")
        t0 = time.time()
        try:
            p = subprocess.run(cmd, capture_output=True, text=True, timeout=7200)
            rc = p.returncode
            stderr_tail = (p.stderr or "")[-2000:]
        except subprocess.TimeoutExpired:
            rc = -2
            stderr_tail = "TIMEOUT after 7200s"
        except Exception as e:
            rc = -1
            stderr_tail = str(e)
        wall = round(time.time() - t0, 1)

        scores = extract_scores(cf_out) if rc == 0 else None

        if scores and scores.get("iptm") is not None:
            print(f"  OK: ipTM={scores['iptm']:.4f} pTM={scores['ptm']:.4f} "
                  f"pLDDT={scores['mean_plddt']:.1f} wall={wall:.0f}s")
            scores["wall_s"] = wall
            scores["fasta_sha256"] = sha256_file(fa)
            with open(existing_scores, "w") as f:
                json.dump(scores, f, indent=2)
            results[name] = scores
        else:
            print(f"  FAIL: rc={rc} wall={wall:.0f}s")
            if "out of memory" in stderr_tail.lower():
                print(f"  *** OOM — this complex may be too large for this GPU ***")
            with open(design_out / "error.json", "w") as f:
                json.dump({"name": name, "exit_code": rc, "wall_s": wall,
                           "stderr_tail": stderr_tail}, f, indent=2)
            results[name] = {"iptm": None, "error": f"rc={rc}"}

    # Evaluate pass/fail
    print("\n" + "=" * 60)
    print("CONTROL VALIDATION RESULTS")
    print("=" * 60)

    g12d = results.get("control_tcr_vs_g12d", {})
    wt = results.get("control_tcr_vs_wt", {})
    pepfree = results.get("control_tcr_vs_pepfree", {})

    for label, r in [("G12D (bind)", g12d), ("WT (no bind)", wt), ("PepFree (no bind)", pepfree)]:
        iptm = r.get("iptm")
        ptm = r.get("ptm")
        if iptm is not None:
            print(f"  {label}: ipTM={iptm:.4f}, pTM={ptm:.4f}")
        else:
            print(f"  {label}: FAILED — {r.get('error', 'no scores')}")

    g12d_iptm = g12d.get("iptm")
    wt_iptm = wt.get("iptm")
    pf_iptm = pepfree.get("iptm")

    if g12d_iptm is None:
        print("\nVERDICT: INCONCLUSIVE — G12D prediction failed")
        verdict = "INCONCLUSIVE"
    elif wt_iptm is None or pf_iptm is None:
        print("\nVERDICT: INCONCLUSIVE — negative control(s) failed")
        verdict = "INCONCLUSIVE"
    elif g12d_iptm > wt_iptm and g12d_iptm > pf_iptm:
        delta_wt = g12d_iptm - wt_iptm
        delta_pf = g12d_iptm - pf_iptm
        print(f"\nVERDICT: PASS — G12D ipTM is {delta_wt:.4f} above WT, "
              f"{delta_pf:.4f} above peptide-free")
        verdict = "PASS"
    else:
        print(f"\nVERDICT: FAIL — pipeline cannot discriminate G12D from controls")
        print(f"  DO NOT USE THIS PIPELINE FOR DESIGN RANKING")
        verdict = "FAIL"

    summary = {
        "verdict": verdict,
        "timestamp_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "pod_id": os.environ.get("RUNPOD_POD_ID", "unknown"),
        "settings": SETTINGS,
        "results": {k: {"iptm": v.get("iptm"), "ptm": v.get("ptm")} for k, v in results.items()},
    }
    with open(out_root / "control_validation_summary.json", "w") as f:
        json.dump(summary, f, indent=2)

    return 0 if verdict == "PASS" else 1


if __name__ == "__main__":
    sys.exit(main())
