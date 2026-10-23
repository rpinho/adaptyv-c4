"""Build C4 submission v6: replace 2 high-charge designs with charge-optimized 4x variants.

Swaps:
  his015_58  (charge -12, pLDDT 95.6) → c4_v_58_4x  (charge -8, pLDDT 95.7)
  his015_33  (charge -12, pLDDT 94.8) → c4_v_33_4x  (charge -8, pLDDT 94.8)

Usage:
    python3 -I scripts/build_c4_submission_v6.py
"""
import csv
from pathlib import Path

RESULTS = Path(__file__).resolve().parent.parent / "results"

SWAPS = {
    "tfr1_e1_his015_58": "c4_v_58_4x",
    "tfr1_e1_his015_33": "c4_v_33_4x",
}


def main():
    # Load variant sequences
    variant_seqs = {}
    with open(RESULTS / "c4_v5_charge_variants.csv") as f:
        for r in csv.DictReader(f):
            variant_seqs[r["name"]] = r["sequence"]

    # Build v6: start from v5, swap 2 designs
    rows = []
    swapped = 0
    with open(RESULTS / "submission_v5_20261023.csv") as f:
        for r in csv.DictReader(f):
            if r["name"] in SWAPS:
                new_name = SWAPS[r["name"]]
                r["name"] = new_name
                r["sequence"] = variant_seqs[new_name]
                swapped += 1
            rows.append(r)

    print(f"Total: {len(rows)} designs ({swapped} swapped)")

    # Validate
    for r in rows:
        seq = r["sequence"]
        assert 10 <= len(seq) <= 250, f"{r['name']}: length {len(seq)} out of range"
        assert set(seq) <= set("ACDEFGHIKLMNPQRSTVWY"), f"{r['name']}: invalid chars"

    names = [r["name"] for r in rows]
    assert len(names) == len(set(names)), f"Duplicate names"
    assert len(rows) == 20

    def net_charge(seq):
        return seq.count("K") + seq.count("R") - seq.count("D") - seq.count("E")

    out = RESULTS / "submission_v6_20261023.csv"
    with open(out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["name", "sequence", "molecule_class"])
        w.writeheader()
        w.writerows(rows)
    print(f"Saved to {out}")

    print("\nDesign summary:")
    for i, r in enumerate(rows):
        seq = r["sequence"]
        charge = net_charge(seq)
        his_pct = round(seq.count("H") / len(seq) * 100, 1)
        tag = " [VARIANT]" if r["name"].startswith("c4_v_") else ""
        print(f"  {i+1:2d}. {r['name']:<50s} len={len(seq):3d}  charge={charge:+4d}  His={his_pct}%{tag}")

    charges = [net_charge(r["sequence"]) for r in rows]
    print(f"\nCharge range: {min(charges)} to {max(charges)}")


if __name__ == "__main__":
    main()
