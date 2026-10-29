"""ESMFold single-sequence refold for submission candidates.

Uses the public ESM Atlas API for structure prediction from sequence.
Computes pLDDT (predicted confidence) for each design.
"""
import csv
import json
import sys
import time
import urllib.request
import urllib.error
from pathlib import Path

RESULTS = Path(__file__).resolve().parent.parent / "results"
API_URL = "https://api.esmatlas.com/foldSequence/v1/pdb/"
REFOLD_DIR = RESULTS / "e1_refolds"

def fold_sequence(sequence, retries=3, delay=5):
    """Submit sequence to ESMFold API, return PDB text."""
    for attempt in range(retries):
        try:
            req = urllib.request.Request(
                API_URL,
                data=sequence.encode(),
                headers={"Content-Type": "text/plain"},
            )
            with urllib.request.urlopen(req, timeout=120) as resp:
                return resp.read().decode()
        except (urllib.error.URLError, urllib.error.HTTPError, TimeoutError) as e:
            print(f"  Attempt {attempt+1} failed: {e}")
            if attempt < retries - 1:
                time.sleep(delay * (attempt + 1))
    return None

def parse_plddt_from_pdb(pdb_text):
    """Extract per-residue pLDDT from B-factor column of ATOM records."""
    plddts = []
    seen = set()
    for line in pdb_text.splitlines():
        if line.startswith("ATOM") and line[12:16].strip() == "CA":
            res_id = line[22:27].strip()
            if res_id not in seen:
                seen.add(res_id)
                bfactor = float(line[60:66].strip())
                plddts.append(bfactor)
    return plddts

def main():
    REFOLD_DIR.mkdir(parents=True, exist_ok=True)

    # Read reviewed submission sequences (post-Fable/Astra review)
    csv_path = RESULTS / "submission_reviewed_20261023.csv"
    if not csv_path.exists():
        csv_path = RESULTS / "submission_dry_run_20261023.csv"
    designs = []
    with open(csv_path) as f:
        for r in csv.DictReader(f):
            designs.append({"name": r["name"], "sequence": r["sequence"]})

    print(f"Refolding {len(designs)} designs via ESMFold API...")

    results = []
    for i, d in enumerate(designs):
        name = d["name"]
        seq = d["sequence"]
        print(f"[{i+1}/{len(designs)}] {name} ({len(seq)} res)...", end=" ", flush=True)

        pdb_path = REFOLD_DIR / f"{name}.pdb"
        if pdb_path.exists():
            pdb_text = pdb_path.read_text()
            print("(cached)", end=" ")
        else:
            pdb_text = fold_sequence(seq)
            if pdb_text is None:
                print("FAILED")
                results.append({"name": name, "length": len(seq), "mean_plddt": -1,
                                "min_plddt": -1, "plddt_gt70_frac": 0, "status": "API_FAILED"})
                continue
            pdb_path.write_text(pdb_text)

        plddts = parse_plddt_from_pdb(pdb_text)
        if not plddts:
            print("NO PLDDTS")
            results.append({"name": name, "length": len(seq), "mean_plddt": -1,
                            "min_plddt": -1, "plddt_gt70_frac": 0, "status": "PARSE_FAILED"})
            continue

        mean_plddt = sum(plddts) / len(plddts)
        min_plddt = min(plddts)
        gt70_frac = sum(1 for p in plddts if p > 70) / len(plddts)

        status = "PASS" if mean_plddt >= 70 else ("MARGINAL" if mean_plddt >= 50 else "FAIL")
        print(f"pLDDT={mean_plddt:.1f} (min={min_plddt:.1f}, >70: {gt70_frac*100:.0f}%) {status}")

        results.append({
            "name": name,
            "length": len(seq),
            "mean_plddt": f"{mean_plddt:.1f}",
            "min_plddt": f"{min_plddt:.1f}",
            "plddt_gt70_frac": f"{gt70_frac:.3f}",
            "status": status,
        })

        time.sleep(1)  # rate limit courtesy

    # Write results
    out = RESULTS / "e1_refold_results.csv"
    with open(out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["name", "length", "mean_plddt", "min_plddt",
                                           "plddt_gt70_frac", "status"])
        w.writeheader()
        w.writerows(results)

    # Summary
    passed = [r for r in results if r["status"] == "PASS"]
    marginal = [r for r in results if r["status"] == "MARGINAL"]
    failed = [r for r in results if r["status"] in ("FAIL", "API_FAILED", "PARSE_FAILED")]
    print(f"\nSummary: {len(passed)} PASS, {len(marginal)} MARGINAL, {len(failed)} FAIL")
    print(f"Saved to {out}")

if __name__ == "__main__":
    main()
