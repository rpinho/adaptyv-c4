#!/usr/bin/env python3
"""Fully symmetric version of the near-identical predicate, and its effect on the C3 edge set. ($0)

v1 of the "corrected" rule took `max(TM_ab, TM_ba)` but coverage from the FORWARD alignment only.
Astra flagged that this leaves the predicate directional. Measured on all 19,900 pairs: 17 pairs have
differing paired-column counts and 67 pairs get different final verdicts by direction.

v2 defines both terms symmetrically:
    TM       = max(tm_norm_chain2(ab), tm_norm_chain2(ba))
    coverage = max(paired_columns(ab), paired_columns(ba)) / max(len_a, len_b)
This is invariant under argument swap by construction (both terms are maxima over the same two
alignments, and lengths are order-free).

Reports the edge set under v1 and v2 and the component counts, so the effect is measured.
Writes results/c3c4_dedup_symmetry_20261025/rule_v2.json
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
OUT = ROOT / 'results/c3c4_dedup_symmetry_20261025/rule_v2.json'
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


def comps(nodes, edges):
    parent = {n: n for n in nodes}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    for a, b in edges:
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[rb] = ra
    groups = {}
    for n in nodes:
        groups.setdefault(find(n), []).append(n)
    return groups


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

    v1, v2, flips = [], [], []
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
            sym_tm = max(r.tm_norm_chain2, r2.tm_norm_chain2)
            cov1 = min(p1 / len(sa), p1 / len(sb)) if p1 else 0.0
            cov2 = max(p1, p2) / max(len(sa), len(sb))
            e1 = sym_tm >= TM_MIN and cov1 >= COV_MIN
            e2 = sym_tm >= TM_MIN and cov2 >= COV_MIN
            if e1:
                v1.append((a, b, round(sym_tm, 4), round(cov1, 4)))
            if e2:
                v2.append((a, b, round(sym_tm, 4), round(cov2, 4)))
            if e1 != e2:
                flips.append([a, b, round(sym_tm, 4), round(cov1, 4), round(cov2, 4)])

    s1 = {(a, b) for a, b, _, _ in v1}
    s2 = {(a, b) for a, b, _, _ in v2}
    g2 = comps(ns, s2)
    g1 = comps(ns, s1)

    payload = {
        'rule_v1': 'max TM + FORWARD-alignment coverage (directional in 67/19900 pairs)',
        'rule_v2': 'max TM + max(paired cols)/max(len) -> symmetric by construction',
        'v1_edges': len(s1), 'v2_edges': len(s2),
        'edges_changed': len(s1 ^ s2),
        'flip_pairs': flips,
        'v1_components': len(g1), 'v2_components': len(g2),
        'v1_largest': max(len(g) for g in g1.values()),
        'v2_largest': max(len(g) for g in g2.values()),
        'v2_edges_list': v2,
    }
    OUT.write_text(json.dumps(payload, indent=1) + '\n')
    print(json.dumps({k: v for k, v in payload.items() if k not in ('flip_pairs', 'v2_edges_list')}, indent=1))
    for f in flips:
        print('  flip:', f)
    return 0


if __name__ == '__main__':
    sys.exit(main())
