#!/usr/bin/env python3
"""Is the FULL corrected predicate symmetric, or only its TM term? ($0, read-only)

Astra's open objection: the corrected code takes TM from both directions but coverage from the forward
alignment only, so the predicate may still be directional. This measures it on EVERY pair of the 200-key
C3 universe rather than asserting it.

Tests, per pair:
  T1  tm_ab reversed == tm_ba            (TM term symmetry)
  T2  paired_columns(ab) == paired_columns(ba)   (coverage term invariance)
  T3  coverage_ab == coverage_ba
  T4  full predicate verdict equal in both directions

Writes results/c3c4_dedup_symmetry_20261025/predicate_symmetry.json
"""
import glob
import importlib.util
import json
import pathlib
import sys

from tmtools import tm_align

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
from c3_backbone_clusters_20261025 import one_letter  # noqa: E402

CROSS = ROOT / 'results/c3_r3_20261025/cross_tm_20261025.json'
OUT = ROOT / 'results/c3c4_dedup_symmetry_20261025/predicate_symmetry.json'
TM_MIN, COV_MIN = 0.65, 0.80


def _load_cw():
    spec = importlib.util.spec_from_file_location(
        'cw', ROOT / 'scripts/c3_round2_cluster_crosswalk_20261025.py')
    cw = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(cw)
    return cw


def find_pdb(n):
    for c in (ROOT / f'results/c3_pod_artifacts_20261025/{n}_complex.pdb',
              ROOT / f'results/c3_round2_111_20261025/backbone_pdb/{n}.pdb',
              ROOT / f'results/c3_round2_111_20261025/screen24/c/{n}/complex_rank_001.pdb'):
        if c.exists():
            return c
    hits = [x for x in sorted(glob.glob(str(ROOT / f'results/**/{n}*.pdb'), recursive=True))
            if 'backbone_pdb' in x or 'complex' in x]
    return pathlib.Path(hits[0]) if hits else None


def main():
    cw = _load_cw()
    cross = json.loads(CROSS.read_text())['pairs']
    nodes = sorted({a for a in cross} | {b for a in cross for b in cross[a]})
    recs = {}
    for n in nodes:
        f = find_pdb(n)
        if f:
            r = cw.load(f, n)
            if r:
                recs[n] = r
    ns = sorted(recs)
    print(f'nodes {len(nodes)}, loaded {len(recs)}')

    t1 = t2 = t3 = t4 = 0
    n = 0
    examples = {'t1': [], 't2': [], 't3': [], 't4': []}
    for i, a in enumerate(ns):
        for b in ns[i + 1:]:
            _, na, ca = recs[a]
            _, nb, cb = recs[b]
            sa, sb = one_letter(na), one_letter(nb)
            if 'X' in sa + sb:
                continue
            r = tm_align(ca, cb, sa, sb)
            r2 = tm_align(cb, ca, sb, sa)
            p1 = sum(1 for u, v in zip(r.seqxA, r.seqyA) if u != '-' and v != '-')
            p2 = sum(1 for u, v in zip(r2.seqxA, r2.seqyA) if u != '-' and v != '-')
            c1 = min(p1 / len(sa), p1 / len(sb)) if p1 else 0.0
            c2 = min(p2 / len(sa), p2 / len(sb)) if p2 else 0.0
            v1 = (r.tm_norm_chain2 >= TM_MIN) and (c1 >= COV_MIN)
            v2 = (r2.tm_norm_chain2 >= TM_MIN) and (c2 >= COV_MIN)
            n += 1
            if abs(r.tm_norm_chain2 - r2.tm_norm_chain2) > 1e-9:
                t1 += 1
                if len(examples['t1']) < 3:
                    examples['t1'].append([a, b, r.tm_norm_chain2, r2.tm_norm_chain2])
            if p1 != p2:
                t2 += 1
                if len(examples['t2']) < 3:
                    examples['t2'].append([a, b, p1, p2])
            if abs(c1 - c2) > 1e-9:
                t3 += 1
                if len(examples['t3']) < 3:
                    examples['t3'].append([a, b, c1, c2])
            if v1 != v2:
                t4 += 1
                if len(examples['t4']) < 5:
                    examples['t4'].append([a, b, {'tm': r.tm_norm_chain2, 'cov': c1},
                                           {'tm': r2.tm_norm_chain2, 'cov': c2}])

    payload = {'pairs_tested': n,
               'T1_directional_tm': t1, 'T2_differing_paired_columns': t2,
               'T3_directional_coverage': t3, 'T4_directional_final_verdict': t4,
               'examples': examples,
               'verdict': ('predicate symmetric' if t4 == 0 else 'predicate NOT symmetric')}
    OUT.write_text(json.dumps(payload, indent=1) + '\n')
    print(json.dumps({k: v for k, v in payload.items() if k != 'examples'}, indent=1))
    return 0


if __name__ == '__main__':
    sys.exit(main())
