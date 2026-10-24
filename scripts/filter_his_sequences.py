"""Filter RFD4 extracted sequences: keep binder (chain C) with ≥N His.

Usage:
    python3 -I scripts/filter_his_sequences.py <sequences.csv> [--min-his 2] [--chain C]
"""
import argparse
import csv
import sys


def main():
    p = argparse.ArgumentParser()
    p.add_argument("csv_path")
    p.add_argument("--min-his", type=int, default=2)
    p.add_argument("--chain", default="C")
    p.add_argument("--out", default="")
    args = p.parse_args()

    kept = []
    total = 0
    with open(args.csv_path) as f:
        for r in csv.DictReader(f):
            if r.get("chain") and r["chain"] != args.chain:
                continue
            total += 1
            seq = r["sequence"]
            n_his = seq.count("H")
            if n_his >= args.min_his:
                kept.append({"name": r["design"], "sequence": seq,
                             "n_his": n_his, "length": len(seq)})

    print(f"Total chain {args.chain} sequences: {total}", file=sys.stderr)
    print(f"Kept with ≥{args.min_his} His: {len(kept)}", file=sys.stderr)

    if kept:
        his_counts = [k["n_his"] for k in kept]
        print(f"His range: {min(his_counts)}-{max(his_counts)}, "
              f"mean={sum(his_counts)/len(his_counts):.1f}", file=sys.stderr)

    out = open(args.out, "w") if args.out else sys.stdout
    w = csv.DictWriter(out, fieldnames=["name", "sequence", "n_his", "length"])
    w.writeheader()
    w.writerows(kept)
    if args.out:
        out.close()
        print(f"Saved to {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
