#!/usr/bin/env python3
"""How much of the v4 score is actually doing work, and the direct seat comparison.

Two questions the headline test cannot answer:

 1. EFFECTIVE WEIGHTS. S = 30P + 25M + 20D + 15C + 10N. If P and M are constant across
    the eligible pool, 55 of the 100 weight points are inert and the ranking is really
    a function of D, C and N only. Measured, not assumed.

 2. DIRECT, IN THE INCUMBENT'S FORMAT. The incumbent reported "better on none, worse on
    two" at raw p<0.05 on the difference of medians. The same tally for the acid-OFF
    seats against the acid-ON seats makes the two directly comparable.

$0, CPU only, read-only.
"""
import json
import pathlib
import subprocess
import sys

import numpy as np

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from c4_v4_falsification_20261026 import (  # noqa: E402
    INDEP_METRICS, HIGHER_IS_BETTER, NPERM, SEED, perm_median_diff, sha_file,
)

OUTDIR = ROOT / "results/c4_ranking_v4_acid_off_20261026"


def main():
    rk = json.loads((OUTDIR / "v4_ranking_acid_off.json").read_text())
    mj = json.loads((ROOT / "results/c4_independent_check_20261026/"
                     "independent_metrics.json").read_text())
    rows = mj["designs"]

    seats = rk["seats"]
    n = len(seats)

    # ---- 1. effective weights ----
    comp = {}
    for k, w in (("P", 30), ("M", 25), ("D", 20), ("C", 15), ("N", 10)):
        vals = np.array([s[k] for s in seats], float)
        comp[k] = {
            "weight": w,
            "seats_min": round(float(vals.min()), 4),
            "seats_max": round(float(vals.max()), 4),
            "seats_median": round(float(np.median(vals)), 4),
            "n_distinct_values": int(len(np.unique(np.round(vals, 6)))),
            "constant_across_seats": bool(np.allclose(vals, vals[0])),
            "weight_points_that_can_discriminate": w if len(np.unique(np.round(vals, 6))) > 1 else 0,
        }
    inert = sum(v["weight"] for v in comp.values() if v["constant_across_seats"])
    live = 100 - inert

    # across the whole eligible pool, is P or M ever non-zero? (they were not in the seats)
    # recompute from the ranking's own notes: P is 0 by rule for every row, M is 0 for
    # every acid-OFF row because the mouse audit covers the acid-ON seats only.
    comp["effective_weights"] = {
        "total_weight": 100,
        "inert_weight_points": inert,
        "live_weight_points": live,
        "inert_fraction": round(inert / 100, 4),
        "reading": ("S is effectively 20D + 15C + 10N on this pool: P cannot exceed 0 "
                    "(no controlled second pH estimator) and M is 0 for every acid-OFF "
                    "design (mouse audit covers the incumbent acid-ON 20 only). The "
                    "rule's two heaviest components are silent here."),
    }

    # ---- 2. direct comparison, incumbent format ----
    seat_shas = {s["sequence_sha256"] for s in seats}
    by_sha = {r["sequence_sha256"]: r for r in rows}
    A = [by_sha[s["sequence_sha256"]] for s in seats if s["sequence_sha256"] in by_sha]
    B = [r for r in rows if r["in_v4"]]
    rng = np.random.default_rng(SEED)
    direct = {}
    for m in INDEP_METRICS:
        a = [r[m] for r in A if r.get(m) is not None]
        b = [r[m] for r in B if r.get(m) is not None]
        diff, med_a, p = perm_median_diff(a, b, rng)
        direct[m] = {
            "acidOFF_median": round(med_a, 4),
            "acidON_median": round(float(np.median(b)), 4),
            "median_diff": round(diff, 4), "p": round(p, 5),
            "higher_is_better": HIGHER_IS_BETTER[m],
            "acidOFF_better": bool((diff > 0) == HIGHER_IS_BETTER[m]),
            "p_bonferroni_13": round(min(1.0, p * 13), 5),
        }
    sig = [m for m in INDEP_METRICS if direct[m]["p"] < 0.05]
    n_better_raw = sum(1 for m in sig if direct[m]["acidOFF_better"])
    n_worse_raw = len(sig) - n_better_raw
    n_better_any = sum(1 for m in INDEP_METRICS if direct[m]["acidOFF_better"])
    sigb = [m for m in INDEP_METRICS if direct[m]["p_bonferroni_13"] < 0.05]

    payload = {
        "created_at_utc": subprocess.check_output(
            ["date", "-u", "+%Y-%m-%dT%H:%M:%SZ"]).decode().strip(),
        "purpose": "Effective component weights and the direct seat-set comparison",
        "component_audit": comp,
        "direct_acidOFF_vs_acidON": direct,
        "direct_tally": {
            "n_metrics": 13,
            "n_significant_raw_p05": len(sig), "significant_metrics": sig,
            "n_significant_acidOFF_better": n_better_raw,
            "n_significant_acidOFF_worse": n_worse_raw,
            "n_metrics_acidOFF_directionally_better": n_better_any,
            "n_metrics_acidOFF_directionally_worse": 13 - n_better_any,
            "n_significant_bonferroni13": len(sigb), "bonferroni_metrics": sigb,
            "incumbent_format": (f"acid-OFF seats vs acid-ON seats: better on {n_better_raw}, "
                                 f"worse on {n_worse_raw} at raw p<0.05, "
                                 f"{len(sigb)} after Bonferroni"),
        },
        "inputs": {"ranking_sha256": sha_file(OUTDIR / "v4_ranking_acid_off.json"),
                   "metrics_sha256": sha_file(
                       ROOT / "results/c4_independent_check_20261026/independent_metrics.json")},
    }
    o = OUTDIR / "component_audit_v4.json"
    if o.exists():
        sys.exit(f"REFUSING TO OVERWRITE {o}")
    o.write_text(json.dumps(payload, indent=1))

    print("EFFECTIVE WEIGHTS")
    for k in ("P", "M", "D", "C", "N"):
        v = comp[k]
        print(f"  {k} weight {v['weight']:2d}  seats {v['seats_min']}..{v['seats_max']} "
              f"median {v['seats_median']}  distinct {v['n_distinct_values']} "
              f"{'CONSTANT (inert)' if v['constant_across_seats'] else 'discriminating'}")
    print(f"  -> inert {inert}/100 weight points = {inert}% ; live {live}")
    print("\nDIRECT: acid-OFF seats vs incumbent acid-ON seats (incumbent format)")
    for m in sorted(INDEP_METRICS, key=lambda m: direct[m]["p"])[:6]:
        v = direct[m]
        print(f"  {m:24s} {v['acidOFF_median']:>9} vs {v['acidON_median']:>9} "
              f"p={v['p']:<8} bonf={v['p_bonferroni_13']:<8} acidOFF_better={v['acidOFF_better']}")
    print("  ", payload["direct_tally"]["incumbent_format"])
    print("\nwrote", o.relative_to(ROOT), sha_file(o))


if __name__ == "__main__":
    main()
