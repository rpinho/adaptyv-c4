"""Build C4 submission v3: replace 4 weak fill-slots with better candidates.

Keeps the original 16 designs unchanged. Replaces:
  his015_193 (charge -39) → his015_177 (charge 0, pLDDT 95.9)
  his015_117 (tandem repeat) → his015_20 (pLDDT 95.4, charge -2)
  his015_168 (47% hydrophobic) → his015_174 (pLDDT 94.8, His 12.3%)
  his015_145 (generic fill) → plain_98 (pLDDT 94.5, His 10.3%)

Usage:
    python3 -I scripts/build_c4_submission_v3.py
"""
import csv
from pathlib import Path

RESULTS = Path(__file__).resolve().parent.parent / "results"

DROP = {
    "tfr1_e1_his015_193",
    "tfr1_e1_his015_117",
    "tfr1_e1_his015_168",
    "tfr1_e1_his015_145",
}

REPLACEMENTS = [
    "tfr1_e1_his015_177",
    "tfr1_e1_his015_20",
    "tfr1_e1_his015_174",
    "tfr1_e1_plain_98",
]


def main():
    # Load original V2 submission, dropping the 4 fill-slots
    rows = []
    with open(RESULTS / "submission_v2_20261023.csv") as f:
        for r in csv.DictReader(f):
            if r["name"] in DROP:
                continue
            rows.append(r)
    print(f"Kept from V2: {len(rows)}")

    # Load replacement sequences
    rep_seqs = {}
    for csv_path in ["results/e1_plain_sequences.csv", "results/e1_his015_sequences.csv"]:
        with open(RESULTS.parent / csv_path) as f:
            for r in csv.DictReader(f):
                if r["chain"] == "C" and r["design"] in REPLACEMENTS:
                    rep_seqs[r["design"]] = r["sequence"]

    for name in REPLACEMENTS:
        if name not in rep_seqs:
            print(f"WARNING: {name} not found!")
            continue
        rows.append({"name": name, "sequence": rep_seqs[name], "molecule_class": "single_chain"})
    print(f"Added replacements: {len(REPLACEMENTS)}")
    print(f"Total: {len(rows)}")

    # Validate
    for r in rows:
        seq = r["sequence"]
        assert 10 <= len(seq) <= 250, f"{r['name']}: length {len(seq)} out of range"
        assert set(seq) <= set("ACDEFGHIKLMNPQRSTVWY"), f"{r['name']}: invalid chars"

    names = [r["name"] for r in rows]
    assert len(names) == len(set(names)), f"Duplicate names"
    assert len(rows) == 20

    out = RESULTS / "submission_v3_20261023.csv"
    with open(out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["name", "sequence", "molecule_class"])
        w.writeheader()
        w.writerows(rows)
    print(f"Saved to {out}")

    print("\nDesign summary:")
    for i, r in enumerate(rows):
        tag = " [NEW]" if r["name"] in REPLACEMENTS else ""
        print(f"  {i+1:2d}. {r['name']:<40s} len={len(r['sequence']):3d}{tag}")


if __name__ == "__main__":
    main()
