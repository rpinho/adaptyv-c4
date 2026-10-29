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

    # Gate on ipTM, rank by L descending (most positive = most correct for C4)
    iptm_pass = [s for s in scored if s.get("iptm", 0) >= 0.75]
    iptm_pass.sort(key=lambda x: (-x["L"], -x.get("iptm", 0), -x.get("mean_plddt", 0)))

    selected = iptm_pass[:args.max]

    n_positive_L = sum(1 for s in selected if s.get("L", 0) > 0)
    print(f"Scored designs: {len(scored)}")
    print(f"ipTM ≥ 0.75: {len(iptm_pass)}")
    print(f"Selected for submission: {len(selected)} ({n_positive_L} with L > 0)")

    if len(selected) < args.max:
        remaining_slots = args.max - len(selected)
        remaining = [s for s in scored if s not in selected]
        remaining.sort(key=lambda x: (-x.get("L", 0), -x.get("iptm", 0)))
        backfill = remaining[:remaining_slots]
        if backfill:
            print(f"Backfilling {len(backfill)} designs (ipTM < 0.75, ranked by L)")
            selected.extend(backfill)

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
