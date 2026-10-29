"""Build C4 submission v4: replace 3 high-charge designs with charge-optimized variants.

Swaps:
  his015_62  (charge -16, pLDDT 95.7) → c4_v_62_max  (charge -7, pLDDT 95.6)
  his015_186 (charge -17, pLDDT 92.0) → c4_v_186_max (charge -3, pLDDT 87.9)
  his015_185 (charge -13, pLDDT 73.2) → c4_v_185_4x  (charge -9, pLDDT 77.6)

Usage:
    python3 -I scripts/build_c4_submission_v4.py
"""
import csv
from pathlib import Path

RESULTS = Path(__file__).resolve().parent.parent / "results"

SWAPS = {
    "tfr1_e1_his015_62": "c4_v_62_max",
    "tfr1_e1_his015_186": "c4_v_186_max",
    "tfr1_e1_his015_185": "c4_v_185_4x",
}


def main():
    # Load variant sequences
    variant_seqs = {}
    with open(RESULTS / "c4_charge_variants.csv") as f:
        for r in csv.DictReader(f):
            variant_seqs[r["name"]] = r["sequence"]

    # Build v4: start from v3, swap 3 designs
    rows = []
    swapped = 0
    with open(RESULTS / "submission_v3_20261023.csv") as f:
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

    # Charge summary
    def net_charge(seq):
        return seq.count("K") + seq.count("R") - seq.count("D") - seq.count("E")

    out = RESULTS / "submission_v4_20261023.csv"
    with open(out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["name", "sequence", "molecule_class"])
        w.writeheader()
        w.writerows(rows)
    print(f"Saved to {out}")

    print("\nDesign summary:")
    for i, r in enumerate(rows):
        seq = r["sequence"]
        charge = net_charge(seq)
        tag = " [VARIANT]" if r["name"].startswith("c4_v_") else ""
        print(f"  {i+1:2d}. {r['name']:<40s} len={len(seq):3d}  charge={charge:+4d}{tag}")


if __name__ == "__main__":
    main()
