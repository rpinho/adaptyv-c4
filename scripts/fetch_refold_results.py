"""Fetch ESMFold refold results from Modal volume and apply quality gate."""
import csv
import json
from pathlib import Path

RESULTS = Path(__file__).resolve().parent.parent / "results"

def main():
    csv_path = RESULTS / "e1_refold_results.csv"
    if not csv_path.exists():
        print(f"No refold results at {csv_path}")
        print("Run: modal run modal_esmfold.py::fetch_results")
        return

    results = []
    with open(csv_path) as f:
        for r in csv.DictReader(f):
            results.append(r)

    print(f"{'Name':<45} {'Len':<5} {'pLDDT':<7} {'Min':<7} {'>70%':<7} {'Status'}")
    print("-" * 95)
    for r in results:
        print(f"{r['name']:<45} {r['length']:<5} {r['mean_plddt']:<7} "
              f"{r['min_plddt']:<7} {float(r['plddt_gt70_frac'])*100:5.1f}% {r['status']}")

    passed = [r for r in results if r["status"] == "PASS"]
    marginal = [r for r in results if r["status"] == "MARGINAL"]
    failed = [r for r in results if r["status"] == "FAIL"]
    print(f"\n{len(passed)} PASS, {len(marginal)} MARGINAL, {len(failed)} FAIL")

    if failed:
        print(f"\nFailed designs to drop:")
        for r in failed:
            print(f"  - {r['name']} (pLDDT={r['mean_plddt']})")

    # Build final filtered list
    reviewed = RESULTS / "submission_reviewed_20261023.csv"
    if reviewed.exists():
        designs = []
        fail_names = {r["name"] for r in failed}
        with open(reviewed) as f:
            for r in csv.DictReader(f):
                if r["name"] in fail_names:
                    print(f"  DROPPING {r['name']} (refold FAIL)")
                    continue
                designs.append(r)

        final = RESULTS / "submission_final_20261023.csv"
        with open(final, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=["name", "sequence", "molecule_class"])
            w.writeheader()
            w.writerows(designs)
        print(f"\nFinal submission: {len(designs)} designs → {final}")

if __name__ == "__main__":
    main()
