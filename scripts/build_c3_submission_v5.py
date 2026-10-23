"""Build C3 submission v5: swap 2 redundant l100 single variants for aromatic-rich RFD4.

Drops (per Fable recommendation to reduce l100 variant redundancy):
  c3_v_l100_E79Q     (single variant of c674f8, charge -6)
  c3_v_l100_E93Q     (single variant of c674f8, charge -6)

Keeps: c3_v_l100_E79Q_E93Q (double mutant, charge -5, best optimization)

Adds (WY-rich RFD4 for groove aromatic complementarity):
  c3_glp1r_ecd_29  (pLDDT 93.5, WY 9.4%, 2W+6Y, charge -3)
  c3_glp1r_ecd_12  (pLDDT 84.2, WY 12.0%, 2W+8Y, charge -4)

Usage:
    python3 -I scripts/build_c3_submission_v5.py
"""
import csv
from pathlib import Path

RESULTS = Path(__file__).resolve().parent.parent / "results"

DROP = {"c3_v_l100_E79Q", "c3_v_l100_E93Q"}

REPLACEMENTS = {
    "c3_glp1r_ecd_29": None,
    "c3_glp1r_ecd_12": None,
}


def main():
    rows = []
    with open(RESULTS / "c3_submission_v4_20261023.csv") as f:
        for r in csv.DictReader(f):
            if r["name"] in DROP:
                print(f"  DROPPED: {r['name']}")
                continue
            rows.append(r)
    print(f"Kept from v4: {len(rows)}")

    # Load replacement sequences
    with open(RESULTS / "c3_rfd4_sequences.csv") as f:
        for r in csv.DictReader(f):
            if r["design"] in REPLACEMENTS and r["chain"] == "A":
                REPLACEMENTS[r["design"]] = r["sequence"]

    for name, seq in REPLACEMENTS.items():
        if seq is None:
            print(f"WARNING: {name} not found!")
            continue
        rows.append({"name": name, "sequence": seq, "molecule_class": "single_chain"})
        wy = seq.count("W") + seq.count("Y")
        print(f"  ADDED: {name} ({len(seq)} res, W+Y={wy})")
    print(f"Total: {len(rows)}")

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

    out = RESULTS / "c3_submission_v5_20261023.csv"
    with open(out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["name", "sequence", "molecule_class"])
        w.writeheader()
        w.writerows(rows)
    print(f"Saved to {out}")

    print("\nDesign summary:")
    for i, r in enumerate(rows):
        seq = r["sequence"]
        charge = net_charge(seq)
        wy = seq.count("W") + seq.count("Y")
        wy_pct = round(wy / len(seq) * 100, 1)
        tag = " [NEW]" if r["name"] in REPLACEMENTS else ""
        print(f"  {i+1:2d}. {r['name']:<35s} len={len(seq):3d}  charge={charge:+4d}  WY={wy_pct}%{tag}")

    # Portfolio stats
    all_seqs = [r["sequence"] for r in rows]
    mean_wy = sum(s.count("W") + s.count("Y") for s in all_seqs) / sum(len(s) for s in all_seqs) * 100
    backbones = set()
    for r in rows:
        n = r["name"]
        if n.startswith("c3_v_"):
            backbones.add(n.split("_")[2] + "_" + n.split("_")[3])
        else:
            backbones.add(n)
    print(f"\nPortfolio: {len(backbones)} independent backbones, mean WY={mean_wy:.1f}%")


if __name__ == "__main__":
    main()
