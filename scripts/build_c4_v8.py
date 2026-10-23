"""Build C4 v8: replace c4_v_186_max (QQQYQQ prion-like, worst charge variant).

Fable ranked c4_v_186_max worst of 5 charge variants:
  - 15.5% Q+N (far above portfolio mean ~4%)
  - QQQYQQ patch = classic prion-like / self-templating low-complexity signature
  - Net -3 with strong positive swing on acidification → aggregation + non-specific binding
  - Unique value (Δq, His%) covered by 7 other designs

Replacement: his015_165 from pool (pI 5.96, Δq 8.8, His 9.2%, Q+N 4.2%, charge -6).

Usage:
    python3 -I scripts/build_c4_v8.py
"""
import csv
from pathlib import Path

RESULTS = Path(__file__).resolve().parent.parent / "results"

DROP = {"c4_v_186_max"}

REPLACEMENT_NAME = "tfr1_e1_his015_96"


def main():
    replacement_seq = None
    with open(RESULTS / "e1_his015_sequences.csv") as f:
        for r in csv.DictReader(f):
            if r["design"] == REPLACEMENT_NAME and r["chain"] == "C":
                replacement_seq = r["sequence"]
                break
    assert replacement_seq, f"{REPLACEMENT_NAME} not found in pool"

    rows = []
    with open(RESULTS / "submission_v7_20261023.csv") as f:
        for r in csv.DictReader(f):
            if r["name"] in DROP:
                print(f"  DROPPED: {r['name']}")
                rows.append({"name": REPLACEMENT_NAME, "sequence": replacement_seq,
                             "molecule_class": "single_chain"})
                print(f"  ADDED: {REPLACEMENT_NAME} (len={len(replacement_seq)})")
                continue
            rows.append(r)

    assert len(rows) == 20

    for r in rows:
        seq = r["sequence"]
        assert 10 <= len(seq) <= 250
        assert set(seq) <= set("ACDEFGHIKLMNPQRSTVWY")

    names = [r["name"] for r in rows]
    assert len(names) == len(set(names))

    # Save ESMFold candidate
    cand_out = RESULTS / "c4_v8_candidates.csv"
    with open(cand_out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["name", "sequence"])
        w.writeheader()
        w.writerow({"name": REPLACEMENT_NAME, "sequence": replacement_seq})
    print(f"ESMFold candidate: {cand_out}")

    out = RESULTS / "submission_v8_20261023.csv"
    with open(out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["name", "sequence", "molecule_class"])
        w.writeheader()
        w.writerows(rows)
    print(f"Saved v8 to {out}")

    def net_charge(seq):
        return seq.count("K") + seq.count("R") - seq.count("D") - seq.count("E")

    print(f"\nV8 summary ({len(rows)} designs):")
    for i, r in enumerate(rows):
        seq = r["sequence"]
        nc = net_charge(seq)
        hf = seq.count("H") / len(seq) * 100
        qn = (seq.count("Q") + seq.count("N")) / len(seq) * 100
        tag = " [NEW]" if r["name"] == REPLACEMENT_NAME else ""
        print(f"  {i+1:2d}. {r['name'][:45]:45s} len={len(seq):3d} charge={nc:+4d} "
              f"His={hf:4.1f}% Q+N={qn:4.1f}%{tag}")


if __name__ == "__main__":
    main()
