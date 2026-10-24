"""Build C3 submission v6: replace low-WY ecd_22 with aromatic-rich ecd_151.

ecd_22 has only 1.5% WY (1 Tyr in 65 res) — poor aromatic complement for the
W33/W39/Y69-lined GLP1R groove. ecd_151 has 11.5% WY (2W+4Y in 52 res).

Usage:
    python3 -I scripts/build_c3_submission_v6.py
"""
import csv
from pathlib import Path

RESULTS = Path(__file__).resolve().parent.parent / "results"

DROP = {"c3_glp1r_ecd_22"}

REPLACEMENT = {"c3_glp1r_ecd_151": None}


def main():
    rows = []
    with open(RESULTS / "c3_submission_v5_20261023.csv") as f:
        for r in csv.DictReader(f):
            if r["name"] in DROP:
                print(f"  DROPPED: {r['name']} (WY={1.5}%)")
                continue
            rows.append(r)
    print(f"Kept from v5: {len(rows)}")

    with open(RESULTS / "c3_rfd4_sequences.csv") as f:
        for r in csv.DictReader(f):
            if r["design"] == "c3_glp1r_ecd_151" and r["chain"] == "A":
                REPLACEMENT["c3_glp1r_ecd_151"] = r["sequence"]

    for name, seq in REPLACEMENT.items():
        if seq is None:
            print(f"WARNING: {name} not found!")
            continue
        rows.append({"name": name, "sequence": seq, "molecule_class": "single_chain"})
        wy = seq.count("W") + seq.count("Y")
        print(f"  ADDED: {name} ({len(seq)} res, W+Y={wy}, WY={wy/len(seq)*100:.1f}%)")

    assert len(rows) == 20

    for r in rows:
        seq = r["sequence"]
        assert 10 <= len(seq) <= 250
        assert set(seq) <= set("ACDEFGHIKLMNPQRSTVWY")

    names = [r["name"] for r in rows]
    assert len(names) == len(set(names))

    def net_charge(seq):
        return seq.count("K") + seq.count("R") - seq.count("D") - seq.count("E")

    out = RESULTS / "c3_submission_v6_20261023.csv"
    with open(out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["name", "sequence", "molecule_class"])
        w.writeheader()
        w.writerows(rows)
    print(f"Saved to {out}")

    print("\nDesign summary:")
    all_seqs = [r["sequence"] for r in rows]
    for i, r in enumerate(rows):
        seq = r["sequence"]
        charge = net_charge(seq)
        wy = seq.count("W") + seq.count("Y")
        wy_pct = round(wy / len(seq) * 100, 1)
        tag = " [NEW]" if r["name"] in REPLACEMENT else ""
        print(f"  {i+1:2d}. {r['name']:<35s} len={len(seq):3d}  charge={charge:+4d}  WY={wy_pct}%{tag}")

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
