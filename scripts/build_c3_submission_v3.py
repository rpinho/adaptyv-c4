"""Build final C3 submission v3: 9 BindCraft + 9 variants + 2 RFD4 = 20 designs.

Usage:
    python3 -I scripts/build_c3_submission_v3.py
"""
import csv
from pathlib import Path

RESULTS = Path(__file__).resolve().parent.parent / "results"


def main():
    rows = []

    # 1. BindCraft originals (9 designs from v2 submission)
    with open(RESULTS / "c3_submission_v2_dry_run_20261023.csv") as f:
        for r in csv.DictReader(f):
            if not r["name"].startswith("c3_v_"):
                rows.append(r)
    print(f"BindCraft originals: {len(rows)}")

    # 2. Point-mutation variants (9 designs)
    variant_count = 0
    with open(RESULTS / "c3_submission_v2_dry_run_20261023.csv") as f:
        for r in csv.DictReader(f):
            if r["name"].startswith("c3_v_"):
                rows.append(r)
                variant_count += 1
    print(f"Point-mutation variants: {variant_count}")

    # 3. RFD4 picks (2 designs — ecd_96 PASS + ecd_22 PASS)
    rfd4_picks = {
        "c3_glp1r_ecd_96": None,
        "c3_glp1r_ecd_22": None,
    }
    with open(RESULTS / "c3_rfd4_sequences.csv") as f:
        for r in csv.DictReader(f):
            if r["design"] in rfd4_picks and r["chain"] == "A":
                rfd4_picks[r["design"]] = r["sequence"]

    for name, seq in rfd4_picks.items():
        if seq is None:
            print(f"WARNING: {name} sequence not found!")
            continue
        rows.append({"name": name, "sequence": seq, "molecule_class": "single_chain"})
    print(f"RFD4 designs: {sum(1 for v in rfd4_picks.values() if v)}")

    print(f"\nTotal: {len(rows)} designs")

    # Validate
    for r in rows:
        seq = r["sequence"]
        assert 10 <= len(seq) <= 250, f"{r['name']}: length {len(seq)} out of range"
        assert set(seq) <= set("ACDEFGHIKLMNPQRSTVWY"), f"{r['name']}: invalid chars"

    names = [r["name"] for r in rows]
    assert len(names) == len(set(names)), f"Duplicate names: {[n for n in names if names.count(n) > 1]}"

    out = RESULTS / "c3_submission_v3_20261023.csv"
    with open(out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["name", "sequence", "molecule_class"])
        w.writeheader()
        w.writerows(rows)
    print(f"Saved to {out}")

    print("\nDesign summary:")
    for i, r in enumerate(rows):
        print(f"  {i+1:2d}. {r['name']:<35s} len={len(r['sequence']):3d}  class={r['molecule_class']}")


if __name__ == "__main__":
    main()
