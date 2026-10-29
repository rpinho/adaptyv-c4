"""Build C4 submission v5: drop 4 Fable-flagged designs, add 4 better candidates.

Drops (Fable HIGH severity):
  his015_197 (39-res tandem repeat, domain-swap risk)
  plain_170 (25-res tandem repeat)
  apical_his_70 (VHH, 1 His, can't pH-switch, needs oxidizing expression)
  his015_177 (pI 7.33, precipitation at pH 7.4, only 4 His)

Replacements (all PASS ESMFold, good pH-switching):
  his015_58  (pLDDT 95.6, His 13.5%, Δq 11.8)
  his015_131 (pLDDT 94.9, His 15.8%, Δq 11.7)
  his015_33  (pLDDT 94.8, His 15.3%, Δq 10.3)
  his015_145 (pLDDT 91.4, His 12.3%, Δq 10.3)

Usage:
    python3 -I scripts/build_c4_submission_v5.py
"""
import csv
from pathlib import Path

RESULTS = Path(__file__).resolve().parent.parent / "results"

DROP = {
    "tfr1_e1_his015_197",
    "tfr1_e1_plain_170",
    "tfr1_his030_200m_20261023__tfr1_apical_his_70",
    "tfr1_e1_his015_177",
}

REPLACEMENT_NAMES = [
    "tfr1_e1_his015_58",
    "tfr1_e1_his015_131",
    "tfr1_e1_his015_33",
    "tfr1_e1_his015_145",
]


def main():
    rows = []
    with open(RESULTS / "submission_v4_20261023.csv") as f:
        for r in csv.DictReader(f):
            if r["name"] in DROP:
                print(f"  DROPPED: {r['name']}")
                continue
            rows.append(r)
    print(f"Kept from v4: {len(rows)}")

    # Load replacement sequences from pool CSVs
    rep_seqs = {}
    for csv_name in ["e1_his015_sequences.csv", "e1_plain_sequences.csv"]:
        path = RESULTS / csv_name
        if not path.exists():
            continue
        with open(path) as f:
            for r in csv.DictReader(f):
                if r["chain"] == "C" and r["design"] in REPLACEMENT_NAMES:
                    rep_seqs[r["design"]] = r["sequence"]

    for name in REPLACEMENT_NAMES:
        if name not in rep_seqs:
            print(f"WARNING: {name} not found!")
            continue
        rows.append({"name": name, "sequence": rep_seqs[name], "molecule_class": "single_chain"})
        print(f"  ADDED: {name} ({len(rep_seqs[name])} res)")
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

    out = RESULTS / "submission_v5_20261023.csv"
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
        tag = " [NEW]" if r["name"] in REPLACEMENT_NAMES else ""
        print(f"  {i+1:2d}. {r['name']:<50s} len={len(seq):3d}  charge={charge:+4d}  His={his_pct}%{tag}")


if __name__ == "__main__":
    main()
