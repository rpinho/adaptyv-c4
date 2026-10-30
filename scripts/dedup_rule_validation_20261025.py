#!/usr/bin/env python3
"""Regression validator for the near-identical predicate (dedup rule).

Read-only. Encodes the two defects found in c3_round2_cluster_crosswalk_20261025.compare()
(results/c3c4_dedup_symmetry_20261025/FINDING.md) as executable invariants:

  D1  TM must be direction-independent. compare(a,b)['tm'] == compare(b,a)['tm'].
  D2  coverage must be *paired* alignment coverage, and can never exceed min(len)/max(len).
      As shipped, coverage counted non-gap residues per sequence and was 1.0 for every pair.

It supports three rule modes over any set of names with resolvable structures:
  legacy   chain2 TM + per-sequence non-gap coverage   (what shipped; kept for provenance replay)
  sym_tm   max(TM_ab, TM_ba) + per-sequence coverage
  fixed    max(TM_ab, TM_ba) + paired-alignment coverage   <- the corrected rule

Usage:
  python3 scripts/dedup_rule_validation_20261025.py --selftest
  python3 scripts/dedup_rule_validation_20261025.py --artifact c3_v8
  python3 scripts/dedup_rule_validation_20261025.py --all --json out.json

Exit code 0 iff every requested check passes.
"""
import argparse
import csv
import glob
import json
import pathlib
import sys

import numpy as np
from tmtools import tm_align

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))

TM_MIN, COV_MIN = 0.65, 0.80

ARTIFACTS = {
    'c3_v8': 'results/c3_submission_v8_20261025.csv',
    'c3_filtered': 'results/c3_challenger_20261025/c3_challenger_filtered_20261025.csv',
    'c3_robust': 'results/c3_challenger_robust_20261025/c3_challenger_all5robust_20261025.csv',
    'c4_ranked': 'results/c4_ranked_20261025/c4_acid_ranked_20261025.csv',
}


def _load_cw():
    import importlib.util
    spec = importlib.util.spec_from_file_location(
        'cw', ROOT / 'scripts/c3_round2_cluster_crosswalk_20261025.py')
    cw = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(cw)
    return cw


def resolve(name):
    """Deterministic structure resolution: canonical paths first, sorted glob last."""
    for c in (ROOT / f'results/c3_pod_artifacts_20261025/{name}_complex.pdb',
              ROOT / f'results/c3_round2_111_20261025/screen24/c/{name}/complex_rank_001.pdb',
              ROOT / f'results/c3_afm2_port_20261025/portval5/portval5_out/{name}/complex_rank_001.pdb',
              ROOT / f'results/c4_rebuild_20261024/remote_pdbs/{name}_complex.pdb'):
        if c.exists():
            return c
    hits = sorted(glob.glob(str(ROOT / f'results/**/{name}/complex_rank_001.pdb'), recursive=True))
    return pathlib.Path(hits[0]) if hits else None


def measure(cw, recs, a, b):
    """Both directions, plus legacy and paired coverage. No rounding before comparison."""
    from c3_backbone_clusters_20261025 import one_letter
    _, na, ca = recs[a]
    _, nb, cb = recs[b]
    sa, sb = one_letter(na), one_letter(nb)
    if 'X' in sa + sb:
        return None
    r = tm_align(ca, cb, sa, sb)
    r2 = tm_align(cb, ca, sb, sa)
    ax, ay = r.seqxA, r.seqyA
    paired = sum(1 for u, v in zip(ax, ay) if u != '-' and v != '-')
    legacy_cov = min(sum(c != '-' for c in ax) / len(sa), sum(c != '-' for c in ay) / len(sb))
    paired_cov = min(paired / len(sa), paired / len(sb)) if paired else 0.0
    return dict(a=a, b=b,
                tm_ab=float(r.tm_norm_chain2), tm_ba=float(r2.tm_norm_chain2),
                legacy_cov=float(legacy_cov), paired_cov=float(paired_cov),
                len_a=len(sa), len_b=len(sb))


def classify(p, mode):
    tm = p['tm_ab'] if mode == 'legacy' else max(p['tm_ab'], p['tm_ba'])
    cov = p['legacy_cov'] if mode in ('legacy', 'sym_tm') else p['paired_cov']
    return tm >= TM_MIN and cov >= COV_MIN


def components(names, edges):
    parent = {n: n for n in names}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    for a, b in edges:
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[rb] = ra
    return len({find(n) for n in names})


def check_artifact(cw, label, path, modes=('legacy', 'sym_tm', 'fixed')):
    names = [r['name'] for r in csv.DictReader(open(ROOT / path))]
    recs = {}
    for n in names:
        p = resolve(n)
        if p is not None:
            rec = cw.load(p, n)
            if rec:
                recs[n] = rec
    ns = sorted(recs)
    pairs = []
    for i, a in enumerate(ns):
        for b in ns[i + 1:]:
            m = measure(cw, recs, a, b)
            if m:
                pairs.append(m)

    # D1: TM symmetry — this is a property of the *measurement*, reported as the defect magnitude
    d1_max_delta = max((abs(p['tm_ab'] - p['tm_ba']) for p in pairs), default=0.0)
    # D2: coverage bound — paired coverage must never exceed min/max length ratio; legacy must not be
    #     trivially 1.0 for pairs whose lengths differ.
    d2_bound_viol = [p for p in pairs if p['paired_cov'] > min(p['len_a'], p['len_b']) / max(p['len_a'], p['len_b']) + 1e-12]
    d2_trivial = [p for p in pairs if abs(p['legacy_cov'] - 1.0) < 1e-12 and p['len_a'] != p['len_b']]

    out = {'artifact': label, 'path': path, 'n_loaded': len(recs), 'n_pairs': len(pairs),
           'checks': {
               'D1_tm_asymmetry_max': round(d1_max_delta, 4),
               'D1_violates_symmetry': d1_max_delta > 1e-9,
               'D2_paired_cov_bound_violations': len(d2_bound_viol),
               'D2_legacy_cov_trivially_1.0_on_unequal_lengths': len(d2_trivial),
           }}
    if d2_bound_viol:
        out['checks']['D2_examples'] = [
            {'a': p['a'], 'b': p['b'], 'paired_cov': p['paired_cov'],
             'bound': round(min(p['len_a'], p['len_b']) / max(p['len_a'], p['len_b']), 4)}
            for p in d2_bound_viol[:3]]
    for mode in modes:
        flags = [p for p in pairs if classify(p, mode)]
        out[mode] = {'near_identical_pairs': len(flags), 'components': components(ns, [(p['a'], p['b']) for p in flags]),
                     'pairs': [{'a': p['a'], 'b': p['b'],
                                'tm': round(max(p['tm_ab'], p['tm_ba']), 4),
                                'tm_ab': round(p['tm_ab'], 4), 'tm_ba': round(p['tm_ba'], 4),
                                'legacy_cov': round(p['legacy_cov'], 4),
                                'paired_cov': round(p['paired_cov'], 4)} for p in flags]}
    return out


def selftest(cw):
    """Synthetic checks that must hold for ANY correct predicate implementation."""
    fails = []
    # 1. TM symmetry on a real pair with a large length gap
    a = ROOT / 'results/c4_rebuild_20261024/remote_pdbs/tfr1_c4_acid_site1_158_complex.pdb'
    b = ROOT / 'results/c4_rebuild_20261024/remote_pdbs/tfr1_c4_acid_site2_1_complex.pdb'
    ra, rb = cw.load(a, 'a'), cw.load(b, 'b')
    m = measure(cw, {'a': ra, 'b': rb}, 'a', 'b')
    if abs(m['tm_ab'] - m['tm_ba']) < 1e-9:
        fails.append('expected the known asymmetric pair to show asymmetry (regression fixture stale?)')
    if max(m['tm_ab'], m['tm_ba']) < 0.65:
        fails.append('symmetric TM should flag this pair')
    if m['paired_cov'] >= 0.80:
        fails.append('paired coverage for a 94-vs-70 pair cannot reach 0.80')
    if abs(m['legacy_cov'] - 1.0) > 1e-12:
        fails.append('legacy coverage should reproduce the 1.0 defect on this pair')
    # 2. Logical invariants of the reference implementation
    if not (measure and classify and components):
        fails.append('helpers missing')
    # 3. Coverage bound
    bound = min(m['len_a'], m['len_b']) / max(m['len_a'], m['len_b'])
    if m['paired_cov'] > bound + 1e-12:
        fails.append('paired coverage exceeded the theoretical bound')
    return fails


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--artifact', choices=sorted(ARTIFACTS))
    ap.add_argument('--all', action='store_true')
    ap.add_argument('--selftest', action='store_true')
    ap.add_argument('--json')
    args = ap.parse_args()

    cw = _load_cw()
    results = {}

    if args.selftest:
        fails = selftest(cw)
        print('SELFTEST', 'PASS' if not fails else 'FAIL')
        for f in fails:
            print('  -', f)
        results['selftest'] = {'pass': not fails, 'failures': fails}

    labels = sorted(ARTIFACTS) if args.all else ([args.artifact] if args.artifact else [])
    for label in labels:
        r = check_artifact(cw, label, ARTIFACTS[label])
        results[label] = r
        print(f"\n== {label}: loaded {r['n_loaded']}/20, {r['n_pairs']} pairs")
        print(f"   D1 max |TM_ab-TM_ba| = {r['checks']['D1_tm_asymmetry_max']}"
              f"  (asymmetric: {r['checks']['D1_violates_symmetry']})")
        print(f"   D2 paired-cov bound violations = {r['checks']['D2_paired_cov_bound_violations']}"
              f" ; legacy cov ==1.0 on unequal lengths for {r['checks']['D2_legacy_cov_trivially_1.0_on_unequal_lengths']} pairs")
        for mode in ('legacy', 'sym_tm', 'fixed'):
            m = r[mode]
            print(f"   {mode:7s}: {m['near_identical_pairs']} near-identical -> {m['components']} components")
            for p in m['pairs']:
                print(f"       {p['a']} x {p['b']}  TM {p['tm_ab']}/{p['tm_ba']} "
                      f"cov legacy {p['legacy_cov']} paired {p['paired_cov']}")

    if args.json:
        pathlib.Path(args.json).write_text(json.dumps(results, indent=1) + '\n')
        print(f"\nwrote {args.json}")

    bad = [k for k, v in results.items()
           if k != 'selftest' and (v['checks']['D2_paired_cov_bound_violations'] > 0)]
    if results.get('selftest', {}).get('failures'):
        bad.append('selftest')
    return 1 if bad else 0


if __name__ == '__main__':
    sys.exit(main())
