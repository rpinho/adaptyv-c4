"""Build the combined E1+E2 selection JSON for build_submission.py."""
import csv
import json
from pathlib import Path

RESULTS = Path(__file__).resolve().parent.parent / "results"

def main():
    designs = []

    # E2 survivors (apical binders, clash-filtered)
    with open(RESULTS / "clash_filtered_shortlist.csv") as f:
        for r in csv.DictReader(f):
            his_frac = float(r["his_fraction"])
            if his_frac > 0.15:
                print(f"SKIP E2 {r['name']}: His {his_frac*100:.1f}% > 15% cap")
                continue
            designs.append({
                "name": r["name"],
                "sequence": r["sequence"],
                "molecule_class": "single_chain",
                "recommendation_status": "recommended",
                "epitope": "E2_apical",
                "his_fraction": his_frac,
                "clash_verdict": r["clash_04A_verdict"],
            })
            print(f"ADD E2 {r['name']}: {len(r['sequence'])} res, "
                  f"His {his_frac*100:.1f}%, {r['clash_04A_verdict']}")

    # E1 candidates (His707 binders)
    with open(RESULTS / "e1_selected_candidates.csv") as f:
        for r in csv.DictReader(f):
            designs.append({
                "name": r["name"],
                "sequence": r["sequence"],
                "molecule_class": "single_chain",
                "recommendation_status": "recommended",
                "epitope": "E1_his707",
                "his_fraction": float(r["his_fraction"]),
                "source": r["source"],
            })

    print(f"\nTotal designs: {len(designs)} (max 20)")
    e1 = [d for d in designs if d["epitope"] == "E1_his707"]
    e2 = [d for d in designs if d["epitope"] == "E2_apical"]
    print(f"  E1 (His707): {len(e1)}")
    print(f"  E2 (apical): {len(e2)}")

    # Verify all sequences are valid
    AA = set("ACDEFGHIKLMNPQRSTVWY")
    for d in designs:
        seq = d["sequence"]
        assert 10 <= len(seq) <= 250, f"{d['name']}: len={len(seq)}"
        bad = set(seq) - AA
        assert not bad, f"{d['name']}: invalid AAs: {bad}"

    out = RESULTS / "combined_selection_20261023.json"
    out.write_text(json.dumps({"designs": designs}, indent=2) + "\n")
    print(f"\nSaved to {out}")

    # Print summary table
    print(f"\n{'#':<3} {'Epitope':<12} {'Name':<40} {'Len':<5} {'His%':<7}")
    print("-" * 70)
    for i, d in enumerate(designs):
        print(f"{i+1:<3} {d['epitope']:<12} {d['name']:<40} "
              f"{len(d['sequence']):<5} {d['his_fraction']*100:5.1f}%")

if __name__ == "__main__":
    main()
