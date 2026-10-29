#!/usr/bin/env python3
"""Robustness of the v4 null: alternative statistics, and acid-OFF vs acid-ON seats.

The headline test uses the difference of MEDIANS (matching the incumbent). A null on
one statistic can be a median artifact, so this re-tests the same two comparisons with
three more:

  MEAN          difference of means (sensitive to tails)
  CLIFF         Cliff's delta (rank-based effect size, robust, no distributional
                assumption); |delta| < 0.147 is "negligible" (Romano et al. thresholds)
  AUC           probability a random seat exceeds a random rejected design

and then asks the question the rule change actually poses:

  DIRECT        the 20 acid-OFF seats vs the 20 incumbent acid-ON seats. The sign fix
                is only worth something if the corrected set is measurably different
                from the set it replaces.

All p-values are permutation, 20,000 draws, two-sided, same seed policy as the main test.
Bonferroni over 13 metrics per family.

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
    INDEP_METRICS, HIGHER_IS_BETTER, NPERM, SEED, sha_file,
)

OUTDIR = ROOT / "results/c4_ranking_v4_acid_off_20261026"


def perm_stat(a, b, fn, center, rng, nperm=NPERM):
    """Two-sided permutation test of fn against its NULL CENTER.

    v1 of this function compared |fn_perm - obs| to |obs|, which centers the null on the
    observed value. For a statistic whose null spread is small relative to the effect
    (Cliff's delta on n=20 vs 293 has spread ~0.03) that drives EVERY p to the floor
    1/(nperm+1) and reports all 13 metrics as significant. The correct comparison is
    |fn_perm - center| >= |obs - center|, where center is the value fn takes under no
    effect: 0 for mean difference and Cliff's delta, 0.5 for AUC.
    The buggy run was discarded before being reported; this is the corrected version.
    """
    a, b = np.asarray(a, float), np.asarray(b, float)
    a, b = a[~np.isnan(a)], b[~np.isnan(b)]
    if len(a) < 2 or len(b) < 2:
        return None, None
    obs = fn(a, b) - center
    pooled = np.concatenate([a, b])
    na = len(a)
    cnt = 0
    for _ in range(nperm):
        rng.shuffle(pooled)
        if abs(fn(pooled[:na], pooled[na:]) - center) >= abs(obs) - 1e-12:
            cnt += 1
    return obs + center, (cnt + 1) / (nperm + 1)


def d_mean(a, b):
    return float(np.mean(a) - np.mean(b))


def d_cliff(a, b):
    """Cliff's delta = P(a>b) - P(a<b)."""
    gt = sum(float((x > b).sum()) for x in a)
    lt = sum(float((x < b).sum()) for x in a)
    return (gt - lt) / (len(a) * len(b))


def d_auc(a, b):
    gt = sum(float((x > b).sum()) for x in a)
    eq = sum(float((x == b).sum()) for x in a)
    return (gt + 0.5 * eq) / (len(a) * len(b))


def main():
    mj = json.loads((ROOT / "results/c4_independent_check_20261026/"
                     "independent_metrics.json").read_text())
    rows = mj["designs"]
    rk = json.loads((OUTDIR / "v4_ranking_acid_off.json").read_text())
    seat_shas = {s["sequence_sha256"] for s in rk["seats"]}
    inelig = {x["name"] for x in rk["ineligible"]}

    seats = [r for r in rows if r["sequence_sha256"] in seat_shas]
    rest_all = [r for r in rows if r["sequence_sha256"] not in seat_shas]
    rest_inelig = [r for r in rest_all if r["name"] in inelig]
    unseated = [r for r in rest_all if r["name"] not in inelig]
    incumbent = [r for r in rows if r["in_v4"]]

    rng = np.random.default_rng(SEED + 77)
    stats = {"mean": d_mean, "cliffs_delta": d_cliff, "auc": d_auc}
    CENTERS = {"mean": 0.0, "cliffs_delta": 0.0, "auc": 0.5}

    # ---- controls on THIS implementation's null centering ----
    rng_n = np.random.default_rng(SEED + 99)
    noise = rng_n.normal(size=len(rows))
    nseats = len(seats)
    ctrl = {}
    for name, fn in stats.items():
        _, p_noise = perm_stat(noise[:nseats], noise[nseats:], fn, CENTERS[name], rng_n)
        planted_a = np.ones(nseats) + rng_n.normal(scale=0.01, size=nseats)
        planted_b = np.zeros(len(rest_all)) + rng_n.normal(scale=0.01, size=len(rest_all))
        _, p_planted = perm_stat(planted_a, planted_b, fn, CENTERS[name], rng_n)
        ctrl[name] = {"noise_p": round(p_noise, 5), "noise_silent": p_noise > 0.05,
                      "planted_p": p_planted, "planted_fires": p_planted <= 2 / (NPERM + 1)}
    all_ctrl_ok = all(v["noise_silent"] and v["planted_fires"] for v in ctrl.values())

    families = {}
    for label, B in (("A_vs_all_v4_ineligible", rest_inelig),
                     ("B_vs_eligible_unseated", unseated),
                     ("DIRECT_acidOFF_vs_acidoN_seats", incumbent)):
        blk = {}
        for m in INDEP_METRICS:
            a = [r[m] for r in seats if r.get(m) is not None]
            b = [r[m] for r in B if r.get(m) is not None]
            row = {}
            for name, fn in stats.items():
                obs, p = perm_stat(a, b, fn, CENTERS[name], rng)
                if obs is None:
                    row[name] = None
                    continue
                row[name] = {"value": round(obs, 5), "p": round(p, 5)}
            # direction of benefit for each statistic (Cliff/AUC > 0.5 = seats higher)
            hi = HIGHER_IS_BETTER[m]
            better = {}
            for name in ("mean", "cliffs_delta"):
                if row.get(name):
                    better[name] = bool((row[name]["value"] > 0) == hi)
            if row.get("auc"):
                better["auc"] = bool((row["auc"]["value"] > 0.5) == hi)
            blk[m] = {"higher_is_better": hi, "stats": row, "seats_better": better,
                      "min_p_across_stats": round(min(v["p"] for v in row.values()
                                                      if v), 5)}
        min_p = min(v["min_p_across_stats"] for v in blk.values())
        n_any_sig_raw = sum(1 for v in blk.values() if v["min_p_across_stats"] < 0.05)
        n_any_sig_bonf = sum(1 for v in blk.values() if v["min_p_across_stats"] * 13 < 0.05)
        families[label] = {
            "n_seats": len(seats), "n_comparison": len(B),
            "per_metric": blk,
            "min_p_over_all_metric_stat_combinations": min_p,
            "n_metrics_sig_raw_p05_any_stat": n_any_sig_raw,
            "n_metrics_sig_bonferroni13_any_stat": n_any_sig_bonf,
            "bonferroni_also_over_3_stats_39_tests": round(min(1.0, min_p * 39), 5),
            "max_abs_cliffs_delta": round(max(
                abs(v["stats"]["cliffs_delta"]["value"]) for v in blk.values()
                if v["stats"].get("cliffs_delta")), 4),
        }

    payload = {
        "created_at_utc": subprocess.check_output(
            ["date", "-u", "+%Y-%m-%dT%H:%M:%SZ"]).decode().strip(),
        "purpose": "Is the v4 null a median artifact? And is the corrected set different "
                   "from the set it replaces?",
        "inputs": {"metrics_sha256": sha_file(
                       ROOT / "results/c4_independent_check_20261026/independent_metrics.json"),
                   "ranking_sha256": sha_file(OUTDIR / "v4_ranking_acid_off.json")},
        "statistics": {"median": "see falsification_v4.json", "mean": "difference of means",
                       "cliffs_delta": "P(a>b)-P(a<b); |d|<0.147 negligible, <0.33 small",
                       "auc": "P(random seat > random rejected)"},
        "n_permutations": NPERM,
        "controls": ctrl,
        "all_controls_fired": all_ctrl_ok,
        "families": families,
    }
    o = OUTDIR / "robustness_v4.json"
    if o.exists():
        sys.exit(f"REFUSING TO OVERWRITE {o}")
    o.write_text(json.dumps(payload, indent=1))

    print("controls:", json.dumps({k: [v["noise_silent"], v["planted_fires"]]
                                   for k, v in ctrl.items()}))
    print("all controls fired:", all_ctrl_ok)
    for label, f in families.items():
        print(f"\n== {label}  (seats {f['n_seats']} vs {f['n_comparison']}) ==")
        print(f"   min p over 13 metrics x 3 stats = {f['min_p_over_all_metric_stat_combinations']}"
              f"  | Bonferroni over 39 = {f['bonferroni_also_over_3_stats_39_tests']}")
        print(f"   metrics significant raw p<0.05 (any stat): {f['n_metrics_sig_raw_p05_any_stat']}/13")
        print(f"   metrics significant Bonferroni-13 (any stat): "
              f"{f['n_metrics_sig_bonferroni13_any_stat']}/13")
        print(f"   max |Cliff's delta| = {f['max_abs_cliffs_delta']} "
              f"({'negligible' if f['max_abs_cliffs_delta'] < 0.147 else 'non-negligible'})")
        worst = sorted(f['per_metric'].items(), key=lambda kv: kv[1]['min_p_across_stats'])[:3]
        for m, v in worst:
            print(f"     {m:24s} min p {v['min_p_across_stats']}  "
                  f"delta={v['stats']['cliffs_delta']['value']:+.4f}  "
                  f"auc={v['stats']['auc']['value']:.3f}")
    print("\nwrote", o.relative_to(ROOT), sha_file(o))


if __name__ == "__main__":
    main()
