"""Rebuild submission after independent reviews (Fable + Astra).

Drops flagged designs, optionally backfills from the candidate pool.
"""
import csv
import json
from pathlib import Path

RESULTS = Path(__file__).resolve().parent.parent / "results"

DROP = {
    "tfr1_e1_his015_193",
    "tfr1_e1_his015_117",
    "tfr1_e1_his015_168",
}

def main():
    with open(RESULTS / "combined_selection_20261023.json") as f:
        data = json.load(f)

    original = data["designs"]
    kept = [d for d in original if d["name"] not in DROP]
    dropped = [d for d in original if d["name"] in DROP]

    print(f"Original: {len(original)} designs")
    print(f"Dropped ({len(dropped)}):")
    for d in dropped:
        print(f"  - {d['name']} ({d['epitope']}, His {d['his_fraction']*100:.1f}%)")
    print(f"Kept: {len(kept)} designs")

    # Check if we can backfill from the candidate pool
    kept_names = {d["name"] for d in kept}
    backfill = []
    pool_path = RESULTS / "e1_selected_candidates.csv"

    # The selected candidates were already the top 17 diverse picks.
    # The 310-candidate pool has more, but we'd need to re-run diversity
    # filtering against the remaining set. For now, just check if any
    # candidates from the pool weren't in the original 17.
    # Actually the original 17 were all from the top of the ranked list,
    # so the next candidates are lower-scored. With 16 designs we're
    # still well within the 1-20 limit. No backfill needed.

    e1 = [d for d in kept if d.get("epitope") == "E1_his707"]
    e2 = [d for d in kept if d.get("epitope") == "E2_apical"]
    print(f"\nFinal split: {len(e2)} E2 + {len(e1)} E1 = {len(kept)} total")

    # Validate
    AA = set("ACDEFGHIKLMNPQRSTVWY")
    for d in kept:
        seq = d["sequence"]
        assert 10 <= len(seq) <= 250, f"{d['name']}: len={len(seq)}"
        bad = set(seq) - AA
        assert not bad, f"{d['name']}: invalid AAs: {bad}"

    out = RESULTS / "combined_selection_reviewed_20261023.json"
    out.write_text(json.dumps({"designs": kept}, indent=2) + "\n")
    print(f"\nSaved to {out}")

    # Also write the reviewed submission CSV
    csv_out = RESULTS / "submission_reviewed_20261023.csv"
    with open(csv_out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["name", "sequence", "molecule_class"])
        w.writeheader()
        for d in kept:
            w.writerow({
                "name": d["name"],
                "sequence": d["sequence"],
                "molecule_class": d["molecule_class"],
            })
    print(f"Submission CSV: {csv_out}")

    print(f"\n{'#':<3} {'Epitope':<12} {'Name':<40} {'Len':<5} {'His%':<7}")
    print("-" * 70)
    for i, d in enumerate(kept):
        print(f"{i+1:<3} {d.get('epitope','?'):<12} {d['name']:<40} "
              f"{len(d['sequence']):<5} {d['his_fraction']*100:5.1f}%")

if __name__ == "__main__":
    main()
