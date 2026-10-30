#!/usr/bin/env python3
"""
C4 FULL RE-SCORE RUNNER — AF2m re-score of designs against P02786-1 homodimer.

Runs ON THE POD. Same settings as the validated canary.
Writes scores.json per design + rescore_run_record.json summary.
Checkpoints after every design so progress survives interruptions.
"""

import argparse
import glob
import hashlib
import json
import os
import pathlib
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

EXPECTED_TARGET_SHA256 = "5345f9f5333fd71cdc0ccfdccab7ef27645c89a437e871af6672470a98b16f93"
EXPECTED_TARGET_LEN = 640


def sha256_str(s):
    return hashlib.sha256(s.encode()).hexdigest()


def sha256_file(p):
    return hashlib.sha256(pathlib.Path(p).read_bytes()).hexdigest()


def validate_fasta(fasta_path):
    with open(fasta_path) as f:
        lines = f.readlines()
    seq_line = lines[1].strip()
    parts = seq_line.split(":")
    if len(parts) != 3:
        return False, f"Expected 3 chains, got {len(parts)}"
    binder, chain_a, chain_b = parts
    if len(chain_a) != EXPECTED_TARGET_LEN:
        return False, f"Chain A length {len(chain_a)} != {EXPECTED_TARGET_LEN}"
    if len(chain_b) != EXPECTED_TARGET_LEN:
        return False, f"Chain B length {len(chain_b)} != {EXPECTED_TARGET_LEN}"
    if chain_a != chain_b:
        return False, "Chain A != Chain B"
    if sha256_str(chain_a) != EXPECTED_TARGET_SHA256:
        return False, "Receptor sha256 mismatch"
    return True, f"OK: binder={len(binder)}aa total={len(binder)+len(chain_a)+len(chain_b)}aa"


def extract_scores(cf_out_dir, name):
    import re
    score_files = sorted(glob.glob(str(cf_out_dir / f"*scores*rank*.json")))
    if not score_files:
        score_files = sorted(glob.glob(str(cf_out_dir / f"*scores*.json")))
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
            "model": model_name,
            "rank": rank,
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

    m = re.search(r'(model_\d+)', pathlib.Path(best_file).name)
    best_model = m.group(1) if m else f"rank_{best_rank}"

    return {
        "iptm": round(iptm, 4) if iptm is not None else None,
        "ptm": round(ptm, 4) if ptm is not None else None,
        "mean_plddt": round(plddt, 2) if plddt is not None else None,
        "best_model": best_model,
        "n_models": len(ranked),
        "all_models": all_model_scores,
    }


def load_checkpoint(out_root):
    record_path = out_root / "rescore_run_record.json"
    if record_path.exists():
        with open(record_path) as f:
            return json.load(f)
    return None


def save_checkpoint(out_root, run_record):
    record_path = out_root / "rescore_run_record.json"
    with open(record_path, "w") as f:
        json.dump(run_record, f, indent=2)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("fasta_dir")
    ap.add_argument("out_dir")
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    fasta_dir = pathlib.Path(args.fasta_dir)
    out_root = pathlib.Path(args.out_dir)
    out_root.mkdir(parents=True, exist_ok=True)

    fastas = sorted(fasta_dir.glob("*_homodimer.fasta"))
    if not fastas:
        print("ERROR: no *_homodimer.fasta files found", file=sys.stderr)
        return 1

    print(f"Re-score: {len(fastas)} designs")
    print(f"Output: {out_root}")
    print(f"Settings: {json.dumps(SETTINGS)}")

    for fa in fastas:
        ok, msg = validate_fasta(fa)
        name = fa.stem.replace("_homodimer", "")
        if not ok:
            print(f"  FAIL {name}: {msg}", file=sys.stderr)
            return 2

    print(f"All {len(fastas)} FASTAs validated.\n")

    if args.dry_run:
        print("Dry run — exiting.")
        return 0

    existing = load_checkpoint(out_root)
    scored_names = set()
    if existing and "designs" in existing:
        for d in existing["designs"]:
            if d.get("scores") and d["scores"].get("iptm") is not None:
                scored_names.add(d["name"])
        print(f"Resuming: {len(scored_names)} already scored\n")

    run_record = existing or {
        "started_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "pod_id": os.environ.get("RUNPOD_POD_ID", "unknown"),
        "gpu": os.environ.get("RUNPOD_GPU_NAME", "unknown"),
        "settings": SETTINGS,
        "target_sha256": EXPECTED_TARGET_SHA256,
        "designs": [],
    }

    existing_design_names = {d["name"] for d in run_record.get("designs", [])}
    total_wall = sum(d.get("wall_s", 0) for d in run_record.get("designs", []))
    n_ok = len(scored_names)
    n_fail = sum(1 for d in run_record.get("designs", []) if d.get("scores") is None and d.get("exit_code") is not None)

    for fa in fastas:
        name = fa.stem.replace("_homodimer", "")

        if name in scored_names:
            continue

        design_out = out_root / name
        design_out.mkdir(parents=True, exist_ok=True)
        work = design_out / "work"
        work.mkdir(exist_ok=True)
        cf_out = work / "out"
        cf_out.mkdir(exist_ok=True)

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

        print(f"FOLD {name} ({n_ok+1}/{len(fastas)})...")
        t0 = time.time()
        try:
            p = subprocess.run(cmd, capture_output=True, text=True, timeout=7200)
            rc = p.returncode
            stderr_tail = (p.stderr or "")[-1000:]
        except subprocess.TimeoutExpired:
            rc = -2
            stderr_tail = "TIMEOUT after 7200s"
        except Exception as e:
            rc = -1
            stderr_tail = str(e)
        t1 = time.time()
        wall = round(t1 - t0, 1)
        total_wall += wall

        scores = extract_scores(cf_out, name) if rc == 0 else None

        rec = {
            "name": name,
            "fasta_sha256": sha256_file(fa),
            "exit_code": rc,
            "wall_s": wall,
            "scores": scores,
        }

        if scores and scores.get("iptm") is not None:
            n_ok += 1
            print(f"  OK {name}: iPTM={scores['iptm']:.4f} pTM={scores['ptm']:.4f} "
                  f"pLDDT={scores['mean_plddt']:.1f} wall={wall:.0f}s [{n_ok}/{len(fastas)}]")
            scores_path = design_out / "scores.json"
            with open(scores_path, "w") as f:
                json.dump({
                    "iptm": scores["iptm"], "ptm": scores["ptm"],
                    "mean_plddt": scores["mean_plddt"],
                    "best_model": scores["best_model"],
                    "all_models": scores["all_models"],
                    "wall_s": wall, "construct": "C_homodimer",
                    "target_sha256": EXPECTED_TARGET_SHA256,
                }, f, indent=2)
        else:
            n_fail += 1
            print(f"  FAIL {name}: rc={rc} wall={wall:.0f}s")
            if "out of memory" in stderr_tail.lower() or "oom" in stderr_tail.lower():
                rec["oom"] = True
                print(f"  *** OOM DETECTED ***")
            with open(design_out / "error.json", "w") as f:
                json.dump({"name": name, "exit_code": rc, "wall_s": wall,
                           "stderr_tail": stderr_tail}, f, indent=2)

        if name not in existing_design_names:
            run_record["designs"].append(rec)
            existing_design_names.add(name)

        run_record["last_updated_utc"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
        run_record["total_wall_s"] = round(total_wall, 1)
        run_record["n_ok"] = n_ok
        run_record["n_fail"] = n_fail
        run_record["n_total"] = len(fastas)
        save_checkpoint(out_root, run_record)

    run_record["finished_utc"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    save_checkpoint(out_root, run_record)

    print(f"\n{'='*60}")
    print(f"RE-SCORE COMPLETE: {n_ok}/{len(fastas)} OK, {n_fail} failed")
    print(f"Total wall: {total_wall:.0f}s ({total_wall/3600:.1f} h)")
    return 0 if n_fail == 0 else 1


if __name__ == "__main__":
    sys.exit(main() or 0)
