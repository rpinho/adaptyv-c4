#!/usr/bin/env python3
"""
Analyse the independent foldability/interface metrics against the v4 ranking.

Reads results/c4_independent_check_20261026/independent_metrics.json and asks:

 Q1  Does the whole 342-design pool look like a foldable protein, measured against
     the natural human TfR1 reference computed in the same run?
 Q2  Do the 20 v4-selected rows differ from the 322 unselected ones on signals
     that entered none of v4's selection?
 Q3  Within the 40 eligible designs (iPTM>=0.60 and apical L6>0), does the L6 rank
     that assigned the seats predict any independent quality signal?

Statistics are permutation-based and implemented here (no scipy dependency):
 - two-sided permutation test on the difference in medians
 - Spearman rho via Pearson correlation on ranks, with a permutation p-value

CPU only. Read-only with respect to all measurement inputs.
"""

import argparse
import hashlib
import json
import os
import sys

import numpy as np

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RNG = np.random.default_rng(20261026)
NPERM = 20000

# Signals that entered NONE of v4's selection (iPTM, apical L6, 2.0 A clash gate).
INDEP_METRICS = [
    "core_apolar_frac",
    "core_charged_frac",
    "n_buried_acidic",
    "buried_frac_of_binder",
    "iface_apolar_frac",
    "bsa_total_A2",
    "bsa_apolar_A2",
    "entropy_bits",
    "max_homopolymer_run",
    "frac_acidic_DE",
    "frac_his",
    "frac_apolar",
    "n_distinct_aa",
]

HIGHER_IS_BETTER = {
    "core_apolar_frac": True,
    "core_charged_frac": False,
    "n_buried_acidic": False,
    "buried_frac_of_binder": True,
    "iface_apolar_frac": True,
    "bsa_total_A2": True,
    "bsa_apolar_A2": True,
    "entropy_bits": True,
    "max_homopolymer_run": False,
    "frac_acidic_DE": False,
    "frac_his": False,
    "frac_apolar": True,
    "n_distinct_aa": True,
}


def perm_median_diff(a, b, nperm=NPERM):
    """Two-sided permutation test on the difference in medians (a - b)."""
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    a = a[~np.isnan(a)]
    b = b[~np.isnan(b)]
    if len(a) < 2 or len(b) < 2:
        return None, None, None
    obs = float(np.median(a) - np.median(b))
    pooled = np.concatenate([a, b])
    na = len(a)
    count = 0
    for _ in range(nperm):
        RNG.shuffle(pooled)
        d = np.median(pooled[:na]) - np.median(pooled[na:])
        if abs(d) >= abs(obs) - 1e-12:
            count += 1
    p = (count + 1) / (nperm + 1)
    return obs, float(np.median(a)), p


def rankdata(x):
    """Average ranks, ties shared."""
    x = np.asarray(x, float)
    order = np.argsort(x, kind="mergesort")
    ranks = np.empty(len(x), float)
    ranks[order] = np.arange(1, len(x) + 1)
    # average ties
    uniq, inv, counts = np.unique(x, return_inverse=True, return_counts=True)
    for i, c in enumerate(counts):
        if c > 1:
            idx = np.where(inv == i)[0]
            ranks[idx] = ranks[idx].mean()
    return ranks


def spearman(x, y, nperm=NPERM):
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    ok = ~(np.isnan(x) | np.isnan(y))
    x, y = x[ok], y[ok]
    if len(x) < 4:
        return None, None, 0
    rx, ry = rankdata(x), rankdata(y)
    rho = float(np.corrcoef(rx, ry)[0, 1])
    cnt = 0
    ry_perm = ry.copy()
    for _ in range(nperm):
        RNG.shuffle(ry_perm)
        r = np.corrcoef(rx, ry_perm)[0, 1]
        if abs(r) >= abs(rho) - 1e-12:
            cnt += 1
    return rho, (cnt + 1) / (nperm + 1), len(x)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--inp", default="results/c4_independent_check_20261026/independent_metrics.json")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    path = os.path.join(REPO, args.inp)
    payload = json.load(open(path))
    rows = payload["designs"]
    controls = payload["controls"]

    out_path = os.path.join(REPO, args.out)
    if os.path.exists(out_path):
        sys.exit(f"REFUSING TO OVERWRITE: {out_path}")

    nat_core = controls["natural_core_reference_11qc_A"]
    nat_seq = controls["NAT_SEQ_human_apical_189_383"]
    nat_iface_hetero = controls["POS_HETERO_1suv_C_vs_A"]
    nat_iface_homo = controls["POS_HOMO_11qc_A_vs_B"]

    v4 = [r for r in rows if r["in_v4"]]
    rest = [r for r in rows if not r["in_v4"]]
    eligible = [r for r in rows if r["iptm_apical"] >= 0.60 and r["L6_apical"] > 0]

    report = {
        "created_at_utc": os.popen("date -u +%Y-%m-%dT%H:%M:%SZ").read().strip(),
        "input": args.inp,
        "input_sha256": hashlib.sha256(open(path, "rb").read()).hexdigest(),
        "n_pool": len(rows),
        "n_v4": len(v4),
        "n_rest": len(rest),
        "n_eligible_iptm_and_L6": len(eligible),
        "natural_reference": {
            "core_apolar_frac_11qc_A": nat_core["core_apolar_frac"],
            "core_charged_frac_11qc_A": nat_core["core_charged_frac"],
            "buried_frac_11qc_A": nat_core["buried_frac_of_chain"],
            "entropy_bits_apical": nat_seq["entropy_bits"],
            "max_homopolymer_run_apical": nat_seq["max_homopolymer_run"],
            "frac_acidic_DE_apical": nat_seq["frac_acidic_DE"],
            "frac_his_apical": nat_seq["frac_his"],
            "frac_apolar_apical": nat_seq["frac_apolar"],
            "iface_apolar_frac_real_hetero_1suv": nat_iface_hetero["iface_apolar_frac"],
            "iface_apolar_frac_real_homodimer_11qc": nat_iface_homo["iface_apolar_frac"],
        },
    }

    # ---- Q1: pool vs natural reference ----
    q1 = {}
    for m in INDEP_METRICS:
        vals = np.array([r[m] for r in rows if r[m] is not None], float)
        q1[m] = {
            "pool_median": round(float(np.median(vals)), 4),
            "pool_p10": round(float(np.percentile(vals, 10)), 4),
            "pool_p90": round(float(np.percentile(vals, 90)), 4),
            "pool_min": round(float(vals.min()), 4),
            "pool_max": round(float(vals.max()), 4),
            "n": int(len(vals)),
        }
    report["Q1_pool_distribution"] = q1

    # fraction of pool worse than the natural reference on key axes
    def frac_worse(metric, ref, higher_better):
        vals = np.array([r[metric] for r in rows if r[metric] is not None], float)
        if higher_better:
            return round(float((vals < ref).mean()), 4), int((vals < ref).sum()), len(vals)
        return round(float((vals > ref).mean()), 4), int((vals > ref).sum()), len(vals)

    report["Q1_fraction_of_pool_worse_than_natural"] = {
        "core_apolar_frac_below_natural_0.738": frac_worse("core_apolar_frac", nat_core["core_apolar_frac"], True),
        "core_charged_frac_above_natural_0.056": frac_worse("core_charged_frac", nat_core["core_charged_frac"], False),
        "entropy_below_natural_4.095": frac_worse("entropy_bits", nat_seq["entropy_bits"], True),
        "frac_acidic_DE_above_natural_0.092": frac_worse("frac_acidic_DE", nat_seq["frac_acidic_DE"], False),
        "frac_his_above_natural_0.021": frac_worse("frac_his", nat_seq["frac_his"], False),
        "frac_apolar_below_natural_0.508": frac_worse("frac_apolar", nat_seq["frac_apolar"], True),
        "iface_apolar_frac_below_real_hetero_0.546": frac_worse(
            "iface_apolar_frac", nat_iface_hetero["iface_apolar_frac"], True),
    }

    # ---- Q2: v4 20 vs unselected 322 ----
    q2 = {}
    for m in INDEP_METRICS:
        a = [r[m] for r in v4 if r[m] is not None]
        b = [r[m] for r in rest if r[m] is not None]
        diff, med_a, p = perm_median_diff(a, b)
        med_b = float(np.median(np.asarray(b, float))) if b else None
        q2[m] = {
            "v4_median": round(med_a, 4) if med_a is not None else None,
            "rest_median": round(med_b, 4) if med_b is not None else None,
            "median_diff_v4_minus_rest": round(diff, 4) if diff is not None else None,
            "perm_p_two_sided": round(p, 5) if p is not None else None,
            "higher_is_better": HIGHER_IS_BETTER[m],
            "v4_better": (
                None if diff is None
                else bool((diff > 0) == HIGHER_IS_BETTER[m])
            ),
        }
    report["Q2_v4_vs_unselected"] = q2

    # ---- Q3: within the 40 eligible, does L6 rank predict independent quality? ----
    q3 = {}
    if len(eligible) >= 4:
        l6 = [r["L6_apical"] for r in eligible]
        ip = [r["iptm_apical"] for r in eligible]
        for m in INDEP_METRICS:
            y = [r[m] for r in eligible]
            rho, p, n = spearman(l6, y)
            rho_i, p_i, _ = spearman(ip, y)
            q3[m] = {
                "spearman_rho_vs_L6": round(rho, 4) if rho is not None else None,
                "perm_p_vs_L6": round(p, 5) if p is not None else None,
                "spearman_rho_vs_iptm": round(rho_i, 4) if rho_i is not None else None,
                "perm_p_vs_iptm": round(p_i, 5) if p_i is not None else None,
                "n": n,
            }
    report["Q3_within_eligible_L6_rank_vs_independent_signal"] = q3

    # ---- per-row table for the 20 v4 seats ----
    keep = ["name", "v4_rank", "L6_apical", "iptm_apical", "core_apolar_frac",
            "core_charged_frac", "n_buried_acidic", "iface_apolar_frac",
            "bsa_total_A2", "entropy_bits", "max_homopolymer_run",
            "frac_acidic_DE", "frac_his", "frac_apolar", "sequence_sha256"]
    report["v4_rows"] = [
        {k: (r[k][:12] if k == "sequence_sha256" else r[k]) for k in keep}
        for r in sorted(v4, key=lambda r: r["v4_rank"])
    ]

    with open(out_path, "w") as fh:
        json.dump(report, fh, indent=1)
    print(f"wrote {out_path}")

    # ---- console summary ----
    print("\n== Q1  pool (n=342) vs natural human TfR1 reference ==")
    for k, v in report["Q1_fraction_of_pool_worse_than_natural"].items():
        print(f"  {k:52s} {v[1]:4d}/{v[2]}  ({v[0]*100:.1f}%)")

    print("\n== Q2  v4 (n=20) vs unselected (n=322), permutation test ==")
    print(f"  {'metric':24s} {'v4':>9s} {'rest':>9s} {'diff':>9s} {'p':>8s}  v4 better?")
    for m, v in q2.items():
        print(f"  {m:24s} {str(v['v4_median']):>9s} {str(v['rest_median']):>9s} "
              f"{str(v['median_diff_v4_minus_rest']):>9s} {str(v['perm_p_two_sided']):>8s}  {v['v4_better']}")

    print(f"\n== Q3  within {len(eligible)} eligible: does L6 rank predict independent quality? ==")
    print(f"  {'metric':24s} {'rho(L6)':>9s} {'p':>8s} {'rho(iPTM)':>10s} {'p':>8s}")
    for m, v in q3.items():
        print(f"  {m:24s} {str(v['spearman_rho_vs_L6']):>9s} {str(v['perm_p_vs_L6']):>8s} "
              f"{str(v['spearman_rho_vs_iptm']):>10s} {str(v['perm_p_vs_iptm']):>8s}")


if __name__ == "__main__":
    main()
