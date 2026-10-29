#!/usr/bin/env python3
"""Positive control: reproduce the INCUMBENT acid-ON result with this lane's own code.

The falsification test above is only interpretable if it agrees with the run it is
being compared to. The incumbent (results/c4_independent_check_20261026/analysis.json)
found: 20 acid-ON seats vs 322 unselected, better on NONE, worse on TWO
(core_charged_frac p=0.0329, frac_acidic_DE p=0.0497), nothing Bonferroni-significant.

If this script reproduces those numbers, the v4 acid-OFF null is a real replication of
a known result rather than a broken test. If it does not, that is a bug in this lane
and must be reported as such.

$0, CPU only, read-only.
"""
import json
import pathlib
import sys

import numpy as np

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from c4_v4_falsification_20261026 import (  # noqa: E402
    INDEP_METRICS, HIGHER_IS_BETTER, NPERM, SEED, perm_median_diff, sha_file,
)

OUTDIR = ROOT / "results/c4_ranking_v4_acid_off_20261026"


def main():
    mj = json.loads((ROOT / "results/c4_independent_check_20261026/"
                     "independent_metrics.json").read_text())
    rows = mj["designs"]
    v4 = [r for r in rows if r["in_v4"]]
    rest = [r for r in rows if not r["in_v4"]]

    rng = np.random.default_rng(SEED)
    out = {}
    for m in INDEP_METRICS:
        a = [r[m] for r in v4 if r.get(m) is not None]
        b = [r[m] for r in rest if r.get(m) is not None]
        diff, med_a, p = perm_median_diff(a, b, rng)
        out[m] = {"v4_median": round(med_a, 4), "rest_median": round(float(np.median(b)), 4),
                  "median_diff": round(diff, 4), "p": round(p, 5),
                  "higher_is_better": HIGHER_IS_BETTER[m],
                  "v4_better": bool((diff > 0) == HIGHER_IS_BETTER[m])}

    ref = json.loads((ROOT / "results/c4_independent_check_20261026/analysis.json").read_text())
    refq = ref["Q2_v4_vs_unselected"]

    agree, disagree = [], []
    for m in INDEP_METRICS:
        mine, theirs = out[m], refq[m]
        ok_p = abs(mine["p"] - theirs["perm_p_two_sided"]) < 0.02
        ok_b = mine["v4_better"] == theirs["v4_better"]
        (agree if (ok_p and ok_b) else disagree).append(
            {"metric": m, "my_p": mine["p"], "ref_p": theirs["perm_p_two_sided"],
             "my_better": mine["v4_better"], "ref_better": theirs["v4_better"],
             "verdict_matches": ok_p and ok_b})

    # The reference claim "better on none, worse on two" refers to the metrics that reach
    # p<0.05 raw, not to the direction of every median. Compare like with like.
    n_better = sum(1 for m in INDEP_METRICS if out[m]["v4_better"])
    n_worse = len(INDEP_METRICS) - n_better
    sig = [m for m in INDEP_METRICS if out[m]["p"] < 0.05]
    n_sig_better = sum(1 for m in sig if out[m]["v4_better"])
    n_sig_worse = len(sig) - n_sig_better
    payload = {
        "created_at_utc": __import__("subprocess").check_output(
            ["date", "-u", "+%Y-%m-%dT%H:%M:%SZ"]).decode().strip(),
        "purpose": "Positive control: my implementation must reproduce the incumbent result",
        "reference": "results/c4_independent_check_20261026/analysis.json Q2_v4_vs_unselected",
        "reference_sha256": sha_file(
            ROOT / "results/c4_independent_check_20261026/analysis.json"),
        "n_v4": len(v4), "n_rest": len(rest),
        "my_result": out,
        "reproduced": {
            "n_metrics_matching": len(agree),
            "n_metrics_disagreeing": len(disagree),
            "disagreements": disagree,
            "my_n_metrics_directionally_better": n_better,
            "my_n_metrics_directionally_worse": n_worse,
            "my_n_significant_raw_p05": len(sig),
            "my_significant_metrics": sig,
            "my_n_significant_better": n_sig_better,
            "my_n_significant_worse": n_sig_worse,
            "reference_claim": "better on none, worse on two (significant at raw p<0.05)",
            "matches_reference_claim": n_sig_better == 0 and n_sig_worse == 2,
        },
    }
    o = OUTDIR / "incumbent_positive_control.json"
    if o.exists():
        sys.exit(f"REFUSING TO OVERWRITE {o}")
    o.write_text(json.dumps(payload, indent=1))
    print(f"metrics matching reference: {len(agree)}/13")
    print(f"direction: {n_better} better / {n_worse} worse (context only)")
    print(f"significant raw p<0.05: {len(sig)} -> {n_sig_better} better / {n_sig_worse} worse")
    print(f"  {sig}")
    print("reference claim: none better, two worse ->",
          "REPRODUCED" if (n_sig_better == 0 and n_sig_worse == 2) else "DIFFERS")
    if disagree:
        print("DISAGREEMENTS:", json.dumps(disagree, indent=1))
    print("wrote", o.relative_to(ROOT), sha_file(o))


if __name__ == "__main__":
    main()
