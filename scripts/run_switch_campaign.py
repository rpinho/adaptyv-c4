#!/usr/bin/env python
"""CPU switch installation on one best-L scored candidate per ranked C4 backbone.

A source must belong to the 25-backbone ranked set. This intentionally uses the
best measured candidate for that backbone, which may differ from its ranked seq0.
Historical linkage excludes chain B; new linkage uses the complete AB target.
All provenance, exclusions, plans, failures, linkage terms and costs go to files.
"""
import argparse
import csv
import json
import re
from pathlib import Path

from install_switch import make_plan, run_plan, save_json
from ph_linkage import score_one

ROOT = Path(__file__).resolve().parents[1]


def substrates(ranked_csv, linkage_csv, pdb_dir):
    with open(ranked_csv) as fh:
        ranked = {re.sub(r'_seq\d+$', '', r['design']): r for r in csv.DictReader(fh)}
    with open(linkage_csv) as fh:
        measured = sorted(csv.DictReader(fh), key=lambda r: (float(r['L']), r['design']))
    selected = {}
    for row in measured:
        backbone = re.sub(r'_candidate\d+\.pdb$', '', row['design'])
        if backbone not in ranked or backbone in selected:
            continue
        pdb = Path(pdb_dir) / row['design']
        if not pdb.is_file():
            raise ValueError('missing source PDB: %s' % pdb)
        selected[backbone] = dict(backbone=backbone, source_pdb=str(pdb.resolve()),
                                   historical_L=float(row['L']), ranked_design=ranked[backbone]['design'])
    if len(selected) != len(ranked):
        raise ValueError('scored representatives missing for %s' % sorted(set(ranked) - set(selected)))
    return list(selected.values())


def score_model(result, directory):
    """Score each retained twin and its matched reference, with the complete target."""
    cache = {}
    for item in result['pair_references'] + [r for r in result['variants'] if r['kept']]:
        name = item['pdb']
        scored = score_one(str(directory / name), result['binder_chain'], result['target_chain'],
                           6.0, 7.4, 0.0, verbose=False)
        if scored is None:
            raise RuntimeError('PROPKA failed: %s' % name)
        # Preserve full residue terms, including the adverse acid tail.
        item['linkage'] = scored
        cache[name] = scored
    for row in result['variants']:
        if row['kept']:
            if row['kind'] == 'switch':
                expected = {(s['chain'], s['resnum']) for s in row['sites']}
                observed = {(r['chain'], r['resnum']) for r in row['linkage']['rows']
                            if r['resname'] == 'HIS'}
                missing = sorted(expected - observed)
                if missing:
                    # PROPKA can suppress a ring covalently coupled to an N terminus.
                    # Keep the partial diagnostic, but never rank it as complete L.
                    row['partial_linkage'] = row.pop('linkage')
                    row['linkage_status'] = 'incomplete_selected_histidines'
                    row['missing_histidines'] = [dict(chain=c, resnum=n) for c, n in missing]
                    continue
            row['linkage_status'] = 'complete'
            row['delta_L_vs_reference'] = round(
                row['linkage']['L'] - cache[row['reference']]['L'], 3)
    return result


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--ranked-csv', type=Path, default=ROOT / 'results/pods/c4_E1_human/3_Ranked/!_Ranked.csv')
    ap.add_argument('--linkage-csv', type=Path, default=ROOT / 'results/v3_all_linkage.csv')
    ap.add_argument('--pdb-dir', type=Path, default=ROOT / 'results/v3_all_pdb')
    ap.add_argument('--out-dir', type=Path, required=True)
    ap.add_argument('--limit', type=int)
    ap.add_argument('--max-positions', type=int, default=8)
    ap.add_argument('--max-pairs', type=int, default=1)
    ap.add_argument('--max-ddg', type=float, default=3.0)
    ap.add_argument('--plan-only', action='store_true')
    a = ap.parse_args()
    selected = substrates(a.ranked_csv, a.linkage_csv, a.pdb_dir)
    if a.limit is not None:
        if a.limit < 1:
            ap.error('--limit must be positive')
        selected = selected[:a.limit]
    out = a.out_dir
    out.mkdir(parents=True, exist_ok=True)
    campaign = dict(status='ranking_only' if a.plan_only else 'running',
                    selection='one best historical-L candidate per ranked backbone; ascending L',
                    historical_target_chains='A (chain B omitted)', target_chains='AB', binder_chain='C',
                    linkage_comparison='Use relaxed AB reference; historical L is ordering only.',
                    min_shift=0.0,
                    seed=20261021, max_ddg=a.max_ddg,
                    max_pairs=a.max_pairs, max_positions=a.max_positions,
                    runpod_session_gpu_spend_usd=0, runpod_last_recorded_total_usd=33.87,
                    runpod_pods='both left stopped; no pod API actions', designs=[])
    save_json(out / 'campaign.json', campaign)
    for rank, source in enumerate(selected, 1):
        directory = out / ('%02d_%s' % (rank, source['backbone']))
        directory.mkdir(parents=True, exist_ok=True)
        plan = make_plan(source['source_pdb'], 'C', 'AB', a.max_positions, 2, a.max_pairs)
        plan.update(source)
        save_json(directory / 'plan.json', dict(plan, status='ranking_only'))
        entry = dict(rank=rank, **source, directory=str(directory.resolve()))
        if a.plan_only:
            entry.update(status='ranking_only', pairs=plan['pairs'])
        else:
            manifest = directory / (Path(source['source_pdb']).stem + '_variants.json')
            try:
                if manifest.exists():
                    result = json.loads(manifest.read_text())
                    # Resume only matching, completed modelling; never reuse stale parameters.
                    for key in ('source_sha256', 'code_sha256', 'target_chain', 'binder_chain', 'max_positions',
                                'sites_per_variant', 'max_pairs', 'pairs'):
                        if result[key] != plan[key]:
                            raise ValueError('resume mismatch: %s' % key)
                    if result['status'] != 'completed' or result['max_ddg'] != a.max_ddg:
                        raise ValueError('incomplete or incompatible manifest; inspect before resuming')
                else:
                    result = run_plan(plan, directory, a.max_ddg)
                if result['status'] == 'completed':
                    if not result.get('linkage_completed'):
                        result = score_model(result, directory)
                        result['linkage_completed'] = True
                        save_json(manifest, result)
                    entry.update(status='completed', variants=result['variants'],
                                 pair_references=result['pair_references'])
                else:
                    entry.update(status=result['status'], errors=result['errors'])
            except Exception as exc:
                entry.update(status='failed', error=repr(exc))
        campaign['designs'].append(entry)
        save_json(out / 'campaign.json', campaign)
        print('%02d/%02d %s %s' % (rank, len(selected), source['backbone'], entry['status']), flush=True)
    if not a.plan_only:
        campaign['status'] = ('completed' if all(r['status'] == 'completed' for r in campaign['designs'])
                              else 'completed_with_errors')
    save_json(out / 'campaign.json', campaign)
    if campaign['status'] == 'completed_with_errors':
        raise SystemExit(1)


if __name__ == '__main__':
    main()
