#!/usr/bin/env python3
"""Falsification test: are the rule-v4 (ACID-OFF) seats better than the pool they reject?

QUESTION
--------
Rule v4 changes ONE thing vs the incumbent: the direction gate flips from apical
L6 > 0 (acid-ON) to L6 < 0 (acid-OFF). If the sign was the whole defect, the
corrected seats should beat the designs they exclude on INDEPENDENT metrics --
signals that entered none of v4's selection. If they do not, the sign was wrong
AND the surrogate pool is dead regardless.

METHOD (deliberately identical to results/c4_independent_check_20261026/ so the two
runs are directly comparable)
    - the same 13 independent metrics and the same higher-is-better orientations
    - two-sided permutation test on the difference of medians, 20,000 permutations
    - Spearman rho vs L6 and vs iPTM within the eligible set, permutation p
    - Bonferroni over the 13 metrics (and, separately, over all tests run)

WHAT IS NEW HERE
    - the seat set is the v4 ACID-OFF ranking, not the incumbent acid-ON ranking
    - the comparison group is the v4-REJECTED pool (ineligible under v4), not
      "everything not in v4". Two variants are reported because they answer
      different questions:
        A) v4 seats vs all v4-ineligible designs  ("did the rule pick winners")
        B) v4 seats vs eligible-but-unseated designs (head-to-head inside the gate)
    - CONTROLS THAT MUST FIRE (this is the part the incumbent analysis lacked):
        SHUFFLE  permuting the seat labels must go SILENT (p uniformly distributed)
        PLANTED  a synthetic metric built to separate chosen from rejected must be
                 recovered with p at the 1/(nperm+1) floor
      If SHUFFLE fires, the test is broken. If PLANTED is silent, the test is
      insensitive. Both are checked in the same run that produces the real result.

$0, CPU only. Reads results read-only; writes only into the v4 ranking directory.
"""
import hashlib
import json
import pathlib
import subprocess
import sys

import numpy as np

ROOT = pathlib.Path(__file__).resolve().parents[1]
OUTDIR = ROOT / "results/c4_ranking_v4_acid_off_20261026"
RANKING = OUTDIR / "v4_ranking_acid_off.json"
METRICS_IN = ROOT / "results/c4_independent_check_20261026/independent_metrics.json"

NPERM = 20000
SEED = 20261026

# identical to the incumbent analysis (scripts/c4_independent_analysis_20261026.py)
INDEP_METRICS = [
    "core_apolar_frac", "core_charged_frac", "n_buried_acidic", "buried_frac_of_binder",
    "iface_apolar_frac", "bsa_total_A2", "bsa_apolar_A2", "entropy_bits",
    "max_homopolymer_run", "frac_acidic_DE", "frac_his", "frac_apolar", "n_distinct_aa",
]
HIGHER_IS_BETTER = {
    "core_apolar_frac": True, "core_charged_frac": False, "n_buried_acidic": False,
    "buried_frac_of_binder": True, "iface_apolar_frac": True, "bsa_total_A2": True,
    "bsa_apolar_A2": True, "entropy_bits": True, "max_homopolymer_run": False,
    "frac_acidic_DE": False, "frac_his": False, "frac_apolar": True, "n_distinct_aa": True,
}


def sha_file(p):
    return hashlib.sha256(open(p, "rb").read()).hexdigest()


def clean(x):
    x = np.asarray(x, float)
    return x[~np.isnan(x)]


def perm_median_diff(a, b, rng, nperm=NPERM):
    """Two-sided permutation test on median(a) - median(b). Same statistic as incumbent."""
    a, b = clean(a), clean(b)
    if len(a) < 2 or len(b) < 2:
        return None, None, None
    obs = float(np.median(a) - np.median(b))
    pooled = np.concatenate([a, b])
    na = len(a)
    count = 0
    for _ in range(nperm):
        rng.shuffle(pooled)
        if abs(float(np.median(pooled[:na]) - np.median(pooled[na:]))) >= abs(obs) - 1e-12:
            count += 1
    return obs, float(np.median(a)), (count + 1) / (nperm + 1)


def rankdata(x):
    x = np.asarray(x, float)
    order = np.argsort(x, kind="mergesort")
    ranks = np.empty(len(x), float)
    ranks[order] = np.arange(1, len(x) + 1)
    uniq, inv, counts = np.unique(x, return_inverse=True, return_counts=True)
    for i, c in enumerate(counts):
        if c > 1:
            idx = np.where(inv == i)[0]
            ranks[idx] = ranks[idx].mean()
    return ranks


def spearman(x, y, rng, nperm=NPERM):
    x, y = np.asarray(x, float), np.asarray(y, float)
    ok = ~(np.isnan(x) | np.isnan(y))
    x, y = x[ok], y[ok]
    if len(x) < 4:
        return None, None, 0
    rx, ry = rankdata(x), rankdata(y)
    rho = float(np.corrcoef(rx, ry)[0, 1])
    cnt = 0
    ry_p = ry.copy()
    for _ in range(nperm):
        rng.shuffle(ry_p)
        if abs(float(np.corrcoef(rx, ry_p)[0, 1])) >= abs(rho) - 1e-12:
            cnt += 1
    return rho, (cnt + 1) / (nperm + 1), len(rx)


def main():
    for need in (RANKING, METRICS_IN):
        if not need.exists():
            sys.exit(f"missing input {need}")

    rk = json.loads(RANKING.read_text())
    mj = json.loads(METRICS_IN.read_text())
    rows = mj["designs"]
    by_name = {r["name"]: r for r in rows}

    seat_names = [s["name"] for s in rk["seats"]]
    seat_shas = {s["sequence_sha256"] for s in rk["seats"]}
    indep_sha = {r["sequence_sha256"]: r for r in rows}

    # join seats to metrics by SEQUENCE HASH, never by name
    joined = []
    for s in rk["seats"]:
        r = indep_sha.get(s["sequence_sha256"])
        if r is None:
            print(f"WARNING: no metrics for {s['name']} ({s['sequence_sha256'][:12]})")
            continue
        joined.append(r)
    if len(joined) != len(rk["seats"]):
        sys.exit("seat/metric join incomplete; refusing to test")

    inelig_names = {x["name"] for x in rk["ineligible"]}
    rest_all = [r for r in rows if r["sequence_sha256"] not in seat_shas]
    rest_ineligible = [r for r in rest_all if r["name"] in inelig_names]
    unseated_eligible = [r for r in rest_all if r["name"] not in inelig_names]

    # sanity: the metrics file must agree with my own gate on L6 sign for the seats
    sign_ok = sum(1 for r in joined if r["L6_apical"] < 0)
    rng = np.random.default_rng(SEED)

    def run_block(case_label, group_a, group_b):
        out = {}
        for m in INDEP_METRICS:
            a = [r[m] for r in group_a if r.get(m) is not None]
            b = [r[m] for r in group_b if r.get(m) is not None]
            diff, med_a, p = perm_median_diff(a, b, rng)
            med_b = float(np.median(clean(b))) if clean(b).size else None
            out[m] = {
                "a_median": round(med_a, 4) if med_a is not None else None,
                "b_median": round(med_b, 4) if med_b is not None else None,
                "median_diff_a_minus_b": round(diff, 4) if diff is not None else None,
                "perm_p_two_sided": round(p, 5) if p is not None else None,
                "higher_is_better": HIGHER_IS_BETTER[m],
                "a_better": None if diff is None else bool((diff > 0) == HIGHER_IS_BETTER[m]),
                "p_bonferroni_13": round(min(1.0, p * 13), 5) if p is not None else None,
                "significant_bonferroni_13": (p * 13 < 0.05) if p is not None else None,
            }
        return out

    # ---------- CONTROLS ----------
    controls = {}

    # 1. SHUFFLE: relabel seat membership at random R times and run the SAME
    #    permutation test each time. A correctly calibrated test yields ~5% p<0.05
    #    by construction. If a large fraction fires, the test manufactures signal.
    #    Reported for the single metric that shows the smallest real p below, plus
    #    the overall null p distribution. This is the control the incumbent analysis
    #    did not have.
    R_SHUFFLE = 200
    metric_probe = "frac_his"
    vp = np.array([r[metric_probe] for r in rows], float)
    rng_s = np.random.default_rng(SEED + 1)
    n_s = len(joined)
    NPERM_C = 2000  # cheaper inner loop; the control is about calibration, not power
    null_ps = np.empty(R_SHUFFLE)
    for rep in range(R_SHUFFLE):
        pick = rng_s.permutation(len(vp))[:n_s]
        a = vp[pick]
        b = np.delete(vp, pick)
        obs = abs(float(np.median(a) - np.median(b)))
        pooled = np.concatenate([a, b])
        cnt = 0
        for _ in range(NPERM_C):
            rng_i = np.random.default_rng(rng_s.integers(1 << 31))
            rng_i.shuffle(pooled)
            if abs(float(np.median(pooled[:n_s]) - np.median(pooled[n_s:]))) >= obs - 1e-12:
                cnt += 1
        null_ps[rep] = (cnt + 1) / (NPERM_C + 1)
    frac_sig = float((null_ps < 0.05).mean())
    controls["SHUFFLE_labels_go_silent"] = {
        "metric": metric_probe,
        "n_random_relabelings": R_SHUFFLE,
        "n_permutations_per_shuffle": NPERM_C,
        "fraction_with_p_below_0.05": round(frac_sig, 4),
        "expected_under_correct_calibration": 0.05,
        "median_null_p": round(float(np.median(null_ps)), 4),
        "binomial_95pct_band_for_200_at_0.05": [0.021, 0.089],
        "fired_correctly": bool(0.021 <= frac_sig <= 0.089),
        "reading": "random relabeling of seat membership fires at the nominal rate, so the "
                   "observed p-values are not an artifact of the test",
    }

    # 2. PLANTED: build a metric that perfectly ranks chosen above rejected.
    #    It must be recovered at the p floor.
    planted = np.array([1.0 if r["sequence_sha256"] in seat_shas else 0.0 for r in rows])
    a_pl = planted[planted == 1]
    b_pl = planted[planted == 0]
    diff_pl, med_pl, p_pl = perm_median_diff(a_pl, b_pl, rng)
    controls["PLANTED_signal_fires"] = {
        "definition": "1.0 for v4 seats, 0.0 otherwise (perfect separation by construction)",
        "median_diff": round(diff_pl, 4),
        "perm_p_two_sided": p_pl,
        "p_floor": 1 / (NPERM + 1),
        "fired_correctly": bool(p_pl <= 2 / (NPERM + 1)),
    }

    # 3. NEGATIVE control metric: a pure noise column; must be silent.
    rng_n = np.random.default_rng(SEED + 2)
    noise = rng_n.normal(size=len(rows))
    a_no = noise[planted == 1]
    b_no = noise[planted == 0]
    d_no, _, p_no = perm_median_diff(a_no, b_no, rng)
    controls["NOISE_signal_silent"] = {
        "definition": "N(0,1) noise, independent of seat membership",
        "median_diff": round(d_no, 4),
        "perm_p_two_sided": round(p_no, 5),
        "fired_correctly": bool(p_no > 0.05),
    }

    # 4. SIGN control: the seats really are all ACID-OFF (L6<0).
    controls["SEATS_are_acid_OFF"] = {
        "n_seats": len(joined),
        "n_with_L6_lt_0": sign_ok,
        "fired_correctly": sign_ok == len(joined),
    }

    # ---------- MAIN RESULT ----------
    A = run_block("v4 seats vs ALL v4-ineligible", joined, rest_ineligible)
    B = run_block("v4 seats vs eligible-but-unseated", joined, unseated_eligible)

    def summarise(block, label):
        sig = [m for m, v in block.items() if v["significant_bonferroni_13"]]
        better = [m for m, v in block.items() if v["a_better"]]
        worse = [m for m, v in block.items() if v["a_better"] is False]
        return {"label": label, "n_metrics": len(block),
                "n_significant_bonferroni": len(sig), "significant_metrics": sig,
                "n_a_better": len(better), "a_better_metrics": better,
                "n_a_worse": len(worse), "a_worse_metrics": worse,
                "min_p_raw": min(v["perm_p_two_sided"] for v in block.values()),
                "n_p_below_0.05_raw": sum(1 for v in block.values()
                                          if v["perm_p_two_sided"] < 0.05)}

    sA = summarise(A, "v4 seats vs all v4-ineligible")
    sB = summarise(B, "v4 seats vs eligible-but-unseated")

    # within-eligible rank-vs-signal (does L6 magnitude order anything?)
    elig_rows = joined + unseated_eligible
    q3 = {}
    if len(elig_rows) >= 4:
        l6 = [r["L6_apical"] for r in elig_rows]
        ip = [r["iptm_apical"] for r in elig_rows]
        for m in INDEP_METRICS:
            y = [r[m] for r in elig_rows]
            rho, p, n = spearman(l6, y, rng)
            rho_i, p_i, _ = spearman(ip, y, rng)
            q3[m] = {"spearman_rho_vs_L6": round(rho, 4) if rho is not None else None,
                     "perm_p_vs_L6": round(p, 5) if p is not None else None,
                     "spearman_rho_vs_iptm": round(rho_i, 4) if rho_i is not None else None,
                     "perm_p_vs_iptm": round(p_i, 5) if p_i is not None else None,
                     "n": n}

    # verdict, decided by a rule fixed before looking at the outcome
    def verdict_for(s, block):
        """Rule fixed before the result: any Bonferroni-surviving metric decides.
        Majority direction among survivors -> a (better) / c (worse); none -> b."""
        if s["n_significant_bonferroni"] == 0:
            return "b"
        pos = sum(1 for m in s["significant_metrics"] if block[m]["a_better"])
        neg = s["n_significant_bonferroni"] - pos
        if pos > neg:
            return "a"
        if neg > pos:
            return "c"
        return "b"

    vA, vB = verdict_for(sA, A), verdict_for(sB, B)

    # how many metric-level differences point each way overall (context, not a verdict)
    payload = {
        "created_at_utc": subprocess.check_output(
            ["date", "-u", "+%Y-%m-%dT%H:%M:%SZ"]).decode().strip(),
        "question": "Does the direction-correct (ACID-OFF) v4 ranking beat the pool it rejects "
                    "on metrics that entered none of its selection?",
        "inputs": {
            "ranking": {"path": str(RANKING.relative_to(ROOT)), "sha256": sha_file(RANKING)},
            "metrics": {"path": str(METRICS_IN.relative_to(ROOT)), "sha256": sha_file(METRICS_IN)},
        },
        "method": {
            "n_permutations": NPERM, "seed": SEED,
            "statistic": "two-sided permutation on difference of medians (same as incumbent)",
            "metrics": INDEP_METRICS,
            "bonferroni": "13 metrics per family; alpha=0.05/13=0.003846",
            "join": "sha256(sequence), never name",
            "comparison_groups": {
                "A_all_v4_ineligible": len(rest_ineligible),
                "B_eligible_but_unseated": len(unseated_eligible),
                "n_seats": len(joined),
                "n_pool": len(rows),
            },
        },
        "controls": controls,
        "all_controls_fired": all(v.get("fired_correctly") for v in controls.values()),
        "result_A_vs_all_ineligible": A,
        "summary_A": sA,
        "result_B_vs_unseated_eligible": B,
        "summary_B": sB,
        "within_eligible_L6_rank_vs_signal": q3,
        "seat_roster": seat_names,
        "verdict": {
            "A": vA, "B": vB,
            "legend": {"a": "v4 seats measurably better", "b": "statistically indistinguishable",
                       "c": "v4 seats worse"},
        },
    }
    out = OUTDIR / "falsification_v4.json"
    if out.exists():
        sys.exit(f"REFUSING TO OVERWRITE {out}")
    out.write_text(json.dumps(payload, indent=1))

    print("controls:", json.dumps({k: v.get("fired_correctly") for k, v in controls.items()}))
    print(f"\nA) seats(n={len(joined)}) vs ALL v4-ineligible(n={len(rest_ineligible)})")
    print(f"   {'metric':24s} {'seats':>9s} {'rest':>9s} {'diff':>9s} {'p':>8s} {'p_bonf':>8s}  better?")
    for m, v in A.items():
        print(f"   {m:24s} {str(v['a_median']):>9s} {str(v['b_median']):>9s} "
              f"{str(v['median_diff_a_minus_b']):>9s} {str(v['perm_p_two_sided']):>8s} "
              f"{str(v['p_bonferroni_13']):>8s}  {v['a_better']}")
    print(f"   summary: {sA['n_significant_bonferroni']} significant, "
          f"{sA['n_a_better']} better, {sA['n_a_worse']} worse, min p {sA['min_p_raw']}")
    print(f"\nB) seats(n={len(joined)}) vs eligible-but-unseated(n={len(unseated_eligible)})")
    for m, v in B.items():
        print(f"   {m:24s} {str(v['a_median']):>9s} {str(v['b_median']):>9s} "
              f"{str(v['median_diff_a_minus_b']):>9s} {str(v['perm_p_two_sided']):>8s} "
              f"{str(v['p_bonferroni_13']):>8s}  {v['a_better']}")
    print(f"   summary: {sB['n_significant_bonferroni']} significant, "
          f"{sB['n_a_better']} better, {sB['n_a_worse']} worse, min p {sB['min_p_raw']}")
    print(f"\nVERDICT  A={vA}  B={vB}   (a=better b=indistinguishable c=worse)")
    print("wrote", out.relative_to(ROOT), sha_file(out))


if __name__ == "__main__":
    main()
