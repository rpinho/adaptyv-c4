#!/usr/bin/env python
"""Reproduce the C4 parent queue and arithmetic sensitivity analysis; no modelling.

Counterfactual pKa changes hold every other contribution fixed. They identify
leverage to investigate, not mutations or experimental/predictive affinity ratios.
"""
import csv
import hashlib
import json
from pathlib import Path

from switch_sites import contribution

ROOT = Path(__file__).resolve().parents[1]


def write_table(path, rows):
    path.with_suffix('.json').write_text(json.dumps(rows, indent=2) + '\n')
    with path.with_suffix('.csv').open('w', newline='') as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main():
    source = ROOT / 'results/switch_install/campaign/campaign.json'
    campaign = json.loads(source.read_text())
    assert campaign['status'] == 'completed'
    rows, parents = [], []
    for design in campaign['designs']:
        candidates = []
        for variant in design['variants']:
            if variant['kind'] != 'switch' or not variant.get('linkage'):
                continue
            control = next(v for v in design['variants']
                           if v['pair'] == variant['pair'] and v['kind'] == 'control')
            sites = {(s['chain'], s['resnum']) for s in variant['sites']}
            for term in variant['linkage']['rows']:
                if term['resname'] != 'HIS' or (term['chain'], term['resnum']) not in sites:
                    continue
                old = contribution(term['pka_free'], term['pka_bound'])
                scenarios = [('raise_free_only', 7.4, term['pka_bound']),
                             ('deepen_bound_only', term['pka_free'], 3.5),
                             ('both', 7.4, 3.5)]
                for label, pf, pb in scenarios:
                    rows.append(dict(
                        backbone=design['backbone'], pair=variant['pair'], source_pdb=variant['pdb'],
                        site=term['resnum'], base_L=variant['linkage']['L'], control_kept=control['kept'],
                        scenario=label, original_free=term['pka_free'], original_bound=term['pka_bound'],
                        hypothetical_free=pf, hypothetical_bound=pb,
                        hypothetical_L=round(variant['linkage']['L'] - old + contribution(pf, pb), 3),
                        caveat='Arithmetic sensitivity only: all other contributions held fixed; '
                               'not a designed mutation or prediction.'))
            if control['kept'] and control.get('linkage'):
                candidates.append((variant['linkage']['L'], variant, control))
        if not candidates:
            continue
        _, variant, control = min(candidates, key=lambda v: v[0])
        parents.append(dict(
            backbone=design['backbone'], pair=variant['pair'], L=variant['linkage']['L'],
            ddG=variant['ddG'], matched_QQ_L=control['linkage']['L'],
            HH_minus_QQ=round(variant['linkage']['L'] - control['linkage']['L'], 3),
            variant_path=str(Path(design['directory']) / variant['pdb']),
            source_path=design['source_pdb'],
            status='priority for context/ensemble validation; not submission selected'))
    parents.sort(key=lambda r: r['L'])
    for rank, parent in enumerate(parents, 1):
        parent['priority_rank'] = rank
    out = ROOT / 'results/next_campaign'
    out.mkdir(exist_ok=True)
    write_table(out / 'parent_priorities', parents)
    write_table(out / 'sensitivity', rows)
    metadata = dict(source=str(source), source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
                    script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                    parent_ids=len(parents), sensitivity_scenarios=len(rows),
                    modelling_performed=False, runpod_new_gpu_spend_usd=0)
    (out / 'provenance.json').write_text(json.dumps(metadata, indent=2) + '\n')
    print('%d parent IDs; %d sensitivity scenarios; no models generated' % (len(parents), len(rows)))


if __name__ == '__main__':
    main()
