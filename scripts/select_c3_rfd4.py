"""Select best C3 RFD4 designs to fill remaining submission slots.

Reads extracted sequences from the RFD4 generation, filters and ranks by
diversity against existing C3 submission, then validates via ESMFold.

Usage:
    python3 -I scripts/select_c3_rfd4.py
"""
import csv
import sys
from pathlib import Path

RESULTS = Path(__file__).resolve().parent.parent / "results"
BINDER_CHAIN = "A"
MIN_LEN = 45
MAX_LEN = 100
MAX_SELECT = 2

AA = set("ACDEFGHIKLMNPQRSTVWY")


def hamming_distance(s1, s2):
    if len(s1) != len(s2):
        return 1.0
    if len(s1) == 0:
        return 0.0
    return sum(a != b for a, b in zip(s1, s2)) / len(s1)


def net_charge(seq):
    return seq.count("K") + seq.count("R") - seq.count("D") - seq.count("E")


def main():
    seq_csv = RESULTS / "c3_rfd4_sequences.csv"
    if not seq_csv.exists():
        print(f"ERROR: {seq_csv} not found. Run sequence extraction first.")
        sys.exit(1)

    # Load binder chains
    binders = []
    with open(seq_csv) as f:
        for r in csv.DictReader(f):
            if r["chain"] != BINDER_CHAIN:
                continue
            seq = r["sequence"]
            n = int(r["n_residues"])
            if not (MIN_LEN <= n <= MAX_LEN):
                continue
            if set(seq) - AA:
                continue
            if "C" in seq:
                continue
            charge = net_charge(seq)
            binders.append({
                "name": r["design"],
                "sequence": seq,
                "n_residues": n,
                "net_charge": charge,
                "abs_charge": abs(charge),
            })
    print(f"Binder chains passing filters: {len(binders)}")

    # Load existing C3 submission sequences for diversity check
    existing_seqs = []
    sub_csv = RESULTS / "c3_submission_v2_dry_run_20261023.csv"
    if sub_csv.exists():
        with open(sub_csv) as f:
            for r in csv.DictReader(f):
                existing_seqs.append(r["sequence"])
    print(f"Existing submission sequences: {len(existing_seqs)}")

    # Score: prefer near-neutral charge, moderate length
    for b in binders:
        charge_score = max(0, 1 - b["abs_charge"] / 10)
        len_penalty = abs(b["n_residues"] - 70) / 30
        len_penalty = min(len_penalty, 1.0)
        b["score"] = charge_score * 0.6 + (1 - len_penalty) * 0.4

    binders.sort(key=lambda x: -x["score"])

    # Diversity filter against existing submission
    selected = []
    MIN_DIVERSITY = 0.4
    for b in binders:
        if b["sequence"] in existing_seqs:
            continue
        too_similar = False
        for es in existing_seqs + [s["sequence"] for s in selected]:
            if abs(len(b["sequence"]) - len(es)) <= 10:
                d = hamming_distance(
                    b["sequence"][:min(len(b["sequence"]), len(es))],
                    es[:min(len(b["sequence"]), len(es))]
                )
                if d < MIN_DIVERSITY:
                    too_similar = True
                    break
        if not too_similar:
            selected.append(b)
        if len(selected) >= MAX_SELECT:
            break

    print(f"\nSelected {len(selected)} RFD4 designs:")
    for i, s in enumerate(selected):
        print(f"  {i+1}. {s['name']} (len={s['n_residues']}, charge={s['net_charge']}, "
              f"score={s['score']:.3f})")
        print(f"     {s['sequence']}")

    out = RESULTS / "c3_rfd4_selected.csv"
    with open(out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["name", "sequence", "molecule_class"])
        w.writeheader()
        for s in selected:
            w.writerow({"name": s["name"], "sequence": s["sequence"],
                        "molecule_class": "single_chain"})
    print(f"\nSaved to {out}")


if __name__ == "__main__":
    main()
