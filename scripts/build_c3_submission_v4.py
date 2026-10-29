"""Build C3 submission v4: drop 2 flagged designs, add 2 better RFD4 picks.

Drops:
  ecd_96  (TM-like, GRAVY +0.72, 55.6% hydrophobic) → c3_glp1r_ecd_172 (pLDDT 89.3, charge 0)
  fff0d0  (MARGINAL refold 68.4, acidic cluster) → c3_glp1r_ecd_171 (pLDDT 88.7, charge -4)

Usage:
    python3 -I scripts/build_c3_submission_v4.py
"""
import csv
from pathlib import Path

RESULTS = Path(__file__).resolve().parent.parent / "results"

DROP = {"c3_glp1r_ecd_96", "c3_bc_l52_fff0d0"}

REPLACEMENTS = {
    "c3_glp1r_ecd_172": None,
    "c3_glp1r_ecd_171": None,
}


def main():
    rows = []
    with open(RESULTS / "c3_submission_v3_20261023.csv") as f:
        for r in csv.DictReader(f):
            if r["name"] in DROP:
                print(f"  DROPPED: {r['name']}")
                continue
            rows.append(r)
    print(f"Kept from v3: {len(rows)}")

    # Load replacement sequences from full pool
    with open(RESULTS / "c3_rfd4_sequences.csv") as f:
        for r in csv.DictReader(f):
            if r["design"] in REPLACEMENTS and r["chain"] == "A":
                REPLACEMENTS[r["design"]] = r["sequence"]

    for name, seq in REPLACEMENTS.items():
        if seq is None:
            print(f"WARNING: {name} not found!")
            continue
        rows.append({"name": name, "sequence": seq, "molecule_class": "single_chain"})
        print(f"  ADDED: {name} ({len(seq)} res)")
    print(f"Total: {len(rows)}")

    # Validate
    for r in rows:
        seq = r["sequence"]
        assert 10 <= len(seq) <= 250, f"{r['name']}: length {len(seq)} out of range"
        assert set(seq) <= set("ACDEFGHIKLMNPQRSTVWY"), f"{r['name']}: invalid chars"

    names = [r["name"] for r in rows]
    assert len(names) == len(set(names)), f"Duplicate names: {[n for n in names if names.count(n) > 1]}"
    assert len(rows) == 20

    def net_charge(seq):
        return seq.count("K") + seq.count("R") - seq.count("D") - seq.count("E")

    out = RESULTS / "c3_submission_v4_20261023.csv"
    with open(out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["name", "sequence", "molecule_class"])
        w.writeheader()
        w.writerows(rows)
    print(f"Saved to {out}")

    print("\nDesign summary:")
    for i, r in enumerate(rows):
        seq = r["sequence"]
        charge = net_charge(seq)
        tag = " [NEW]" if r["name"] in REPLACEMENTS else ""
        print(f"  {i+1:2d}. {r['name']:<35s} len={len(seq):3d}  charge={charge:+4d}{tag}")


if __name__ == "__main__":
    main()
