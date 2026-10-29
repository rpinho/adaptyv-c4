#!/usr/bin/env python
"""Write reviewable CSVs and a factual report from completed CPU switch manifests."""
import argparse
import csv
import json
from pathlib import Path


def write_csv(path, rows):
    if not rows:
        return
    keys = list(dict.fromkeys(k for row in rows for k in row))
    with path.open('w', newline='') as fh:
        writer = csv.DictWriter(fh, fieldnames=keys)
        writer.writeheader()
        writer.writerows(rows)


def summarize(out):
    campaign = json.loads((out / 'campaign.json').read_text())
    rows, residues = [], []
    for design in campaign['designs']:
        controls = {v['pair']: v for v in design.get('variants', []) if v['kind'] == 'control'}
        switches_by_pair = {v['pair']: v for v in design.get('variants', []) if v['kind'] == 'switch'}
        for v in design.get('variants', []):
            linkage = v.get('linkage', {})
            selected = {(s['chain'], s['resnum']) for s in v['sites']}
            his = [r for r in linkage.get('rows', []) if r['resname'] == 'HIS'
                   and (r['chain'], r['resnum']) in selected]
            ctrl = controls.get(v['pair'], {})
            row = dict(rank=design['rank'], backbone=design['backbone'], pair=v['pair'],
                       sites=';'.join('%s%s%d' % (s['chain'], s['wt'], s['resnum']) for s in v['sites']),
                       kind=v['kind'], mutation=v['mut'], status=v['status'],
                       linkage_status=v.get('linkage_status', 'not_scored'),
                       historical_L=design['historical_L'], reference_dG=v['base_dG'],
                       dG=v.get('dG'), ddG=v.get('ddG'), kept=v['kept'],
                       L=linkage.get('L'), delta_L_vs_reference=v.get('delta_L_vs_reference'),
                       control_kept=ctrl.get('kept', False),
                       matched_pair_retained=bool(switches_by_pair.get(v['pair'], {}).get('kept', False)
                                                 and ctrl.get('kept', False)),
                       switch_gate=bool(v['kind'] == 'switch' and v['kept'] and linkage.get('switch_ok')),
                       selected_His_free_pKa=';'.join('%s%d:%.2f' % (r['chain'], r['resnum'], r['pka_free']) for r in his),
                       selected_His_bound_pKa=';'.join('%s%d:%.2f' % (r['chain'], r['resnum'], r['pka_bound']) for r in his),
                       selected_His_free_ge_7=sum(r['pka_free'] >= 7 for r in his),
                       positive_residue_sum=round(sum(r['contribution'] for r in linkage.get('rows', []) if r['contribution'] > 0), 3),
                       pdb=v['pdb'])
            if v['kind'] == 'switch' and linkage and ctrl.get('linkage'):
                row['L_HH_minus_QQ'] = round(linkage['L'] - ctrl['linkage']['L'], 3)
            rows.append(row)
            for r in linkage.get('rows', []):
                residues.append(dict(backbone=design['backbone'], pair=v['pair'], kind=v['kind'], **r))
        for reference in design.get('pair_references', []):
            for r in reference.get('linkage', {}).get('rows', []):
                residues.append(dict(backbone=design['backbone'], pair=reference['pair'], kind='reference', **r))
    write_csv(out / 'variants.csv', rows)
    write_csv(out / 'residues.csv', residues)
    switches = [r for r in rows if r['kind'] == 'switch']
    retained = [r for r in switches if r['kept'] and r['L'] is not None]
    retained.sort(key=lambda r: r['L'])
    gates = [r for r in retained if r['switch_gate']]
    passed = [r for r in gates if r['matched_pair_retained']]
    complete = sum(d['status'] == 'completed' for d in campaign['designs'])
    text = [
        '# CPU switch-installation results', '',
        'Campaign status: **%s**. Completed backbones: **%d/%d**.' % (campaign['status'], complete, len(campaign['designs'])), '',
        '- Configured cap: %d ranked pairs per backbone from %d shortlisted sites.' % (
            campaign.get('max_pairs', 1), campaign.get('max_positions', 8)),
        '- HH pairs attempted: %d; passing neutral-energy filters: %d; complete linkage scores: %d.' % (
            len(switches), sum(r['kept'] for r in switches), len(retained)),
        '- HH passing energy filters but unscorable for complete installed-site linkage: %d.' % sum(
            r['kept'] and r['L'] is None for r in switches),
        '- Complete-linkage HH with a retained QQ twin: %d.' % sum(r['matched_pair_retained'] for r in retained),
        '- Retained HH with L <= -2: %d; with a retained QQ twin: %d.' % (len(gates), len(passed)),
        '- HH pairs with both selected His free-state pKa >= 7: %d.' % sum(r['selected_His_free_ge_7'] >= 2 for r in retained),
        '- Runpod GPU spend this session: **$0.00**. Last recorded historical pod spend: **$33.87**. '
        'Both pods left stopped. Persistent-volume storage not freshly reconciled.', '',
        'References and variants use target AB and binder C. Historical ordering used A only. '
        'The whole complex is first repacked and side chains minimised with the backbone fixed; '
        'all variants then use the same local protocol as their reference. '
        f'Retention requires negative reference and variant dG and ddG <= {campaign.get("max_ddg", 3.0):g} REU. '
        'ddG is Rosetta REU; L is PROPKA Wyman linkage at pH 6.0 versus 7.4. '
        'Neither establishes experimental affinity or the assayed pH ratio.', '',
        '## Best retained HH variants', '',
        '| Backbone | Sites | L | ΔL vs relaxed AB reference | ddG (REU) | QQ retained | Selected free pKa |',
        '|---|---|---:|---:|---:|---|---|',
    ]
    for r in retained[:10]:
        text.append('| %s | %s | %.3f | %+.3f | %+.2f | %s | %s |' % (
            r['backbone'].split('_forced_')[-1], r['sites'], r['L'], r['delta_L_vs_reference'],
            r['ddG'], r['control_kept'], r['selected_His_free_pKa']))
    if not retained:
        text.append('')
        text.append('No retained HH result available; inspect manifest failures before drawing conclusions.')
    text += ['', '## Limits and stop point', '',
             'The configured ranked-pair cap is a screen, not an exhaustive combinatorial search. '
             'The coordinate ranking is a heuristic; selected acidic residues can themselves be removed by mutation '
             'and neighbouring sites can interact. Measured mutant pKa and the full acid budget, not pair-score sums, '
             'must decide whether the mechanism improved.', '',
             'Ranking used source coordinates before preparation. Rotamer changes can materially change PROPKA pKa; '
             'single-conformer sensitivity limits absolute pKa and energy claims.', '',
             'No APBS estimator, orthogonal refolding, mouse validation, account creation or submission was run. '
             'This CPU campaign stops after its configured backbone/pair count. No rented pod is needed.', '',
             'Files: `campaign.json` (full state), per-design manifests and PDBs, '
             '`variants.csv` (including rejected twins), `residues.csv` (both tails of every scored linkage sum).', '']
    (out / 'RESULTS.md').write_text('\n'.join(text))


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('out_dir', type=Path)
    summarize(ap.parse_args().out_dir)
