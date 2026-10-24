"""Build C4 submission v8: replace wrong-direction designs with scored acid-site designs.

Reads Wyman-linkage scored results, selects top designs by L (≥ +2 gate),
combines with any keeper designs from v7, outputs submission CSV.

Usage:
    python3 -I scripts/build_submission_v8.py \
        --scored results/c4_acid_scored.json \
        --max 20
"""
import argparse
import csv
import json


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--scored", required=True, help="Wyman-scored JSON from score_c4_acid.py")
    p.add_argument("--max", type=int, default=20, help="Max designs in submission")
    p.add_argument("--gate", type=float, default=2.0, help="Wyman L gate")
    p.add_argument("--out", default="results/submission_v8.csv")
    args = p.parse_args()

    with open(args.scored) as f:
        scored = json.load(f)

    passing = [s for s in scored if s.get("L", 0) >= args.gate]
    passing.sort(key=lambda x: (-x["L"], -x.get("iptm", 0), -x.get("mean_plddt", 0)))

    selected = passing[:args.max]

    print(f"Scored designs: {len(scored)}")
    print(f"Passing L ≥ {args.gate}: {len(passing)}")
    print(f"Selected for submission: {len(selected)}")

    if len(selected) < args.max:
        remaining_slots = args.max - len(selected)
        marginal = [s for s in scored if s.get("L", 0) >= 1.0 and s not in selected]
        marginal.sort(key=lambda x: (-x["L"], -x.get("iptm", 0)))
        backfill = marginal[:remaining_slots]
        if backfill:
            print(f"Backfilling {len(backfill)} designs with L ≥ 1.0")
            selected.extend(backfill)

    if not selected:
        print("WARNING: No designs pass the gate. Using top designs by L regardless.")
        scored.sort(key=lambda x: (-x.get("L", 0), -x.get("iptm", 0)))
        selected = scored[:args.max]

    print(f"\n{'Name':<40s} {'L':>6s} {'ipTM':>6s} {'pLDDT':>6s} {'#H':>3s} {'Len':>4s}")
    print("-" * 70)
    for s in selected:
        print(f"{s['name']:<40s} {s.get('L',0):6.2f} {s.get('iptm',0):6.3f} "
              f"{s.get('mean_plddt',0):6.1f} {s.get('n_his',0):3d} {s.get('length',0):4d}")

    rows = []
    for s in selected:
        rows.append({
            "name": s["name"],
            "sequence": s["sequence"],
            "molecule_class": "single_chain",
        })

    with open(args.out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["name", "sequence", "molecule_class"])
        w.writeheader()
        w.writerows(rows)
    print(f"\nSaved {len(rows)} designs to {args.out}")


if __name__ == "__main__":
    main()
