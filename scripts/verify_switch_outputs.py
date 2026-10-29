#!/usr/bin/env python
"""Audit completed switch outputs against their source sequence and backbone."""
import argparse
import hashlib
import json
import math
from pathlib import Path

from install_switch import AA1, save_json


def coordinates(path):
    residues, bb = {}, {}
    for line in path.read_text().splitlines():
        if not line.startswith('ATOM'):
            continue
        key = (line[21], int(line[22:26]))
        residues[key] = AA1[line[17:20].strip()]
        atom = line[12:16].strip()
        if atom in {'N', 'CA', 'C', 'O'}:
            bb[key + (atom,)] = tuple(float(line[k:k+8]) for k in (30, 38, 46))
    return residues, bb


def verify(out):
    campaign = json.loads((out / 'campaign.json').read_text())
    assert campaign['status'] == 'completed', 'campaign must complete first'
    checks = []
    for design in campaign['designs']:
        directory = Path(design['directory'])
        manifest = directory / (Path(design['source_pdb']).stem + '_variants.json')
        result = json.loads(manifest.read_text())
        assert result['status'] == 'completed' and result['linkage_completed']
        source = Path(result['source_pdb'])
        assert hashlib.sha256(source.read_bytes()).hexdigest() == result['source_sha256']
        for script, digest in result['code_sha256'].items():
            assert hashlib.sha256(Path(__file__).with_name(script).read_bytes()).hexdigest() == digest
        residues, backbone = coordinates(source)
        files = [result['preparation']] + result['pair_references']
        files += [r for r in result['variants'] if r['kept']]
        max_distance = 0.0
        for item in files:
            pdb = directory / item['pdb']
            assert hashlib.sha256(pdb.read_bytes()).hexdigest() == item['sha256']
            expected = dict(residues)
            if item.get('kind') in ('switch', 'control'):
                for site, aa in zip(item['sites'], item['mut']):
                    assert site['chain'] == result['binder_chain']
                    expected[(site['chain'], site['resnum'])] = aa
                assert item['dG'] < 0 and item['base_dG'] < 0
                assert math.isfinite(item['ddG']) and item['ddG'] <= result['max_ddg']
                assert abs(item['ddG'] - (item['dG'] - item['base_dG'])) < 1e-8
            actual, other_backbone = coordinates(pdb)
            assert expected == actual, 'unexpected sequence/target change in %s' % pdb
            assert set(backbone) <= set(other_backbone)
            displacement = max(math.dist(backbone[k], other_backbone[k]) for k in backbone)
            assert displacement <= 0.002, 'backbone moved in %s' % pdb
            max_distance = max(max_distance, displacement)
            if 'linkage' in item:
                linkage = item['linkage']
                terms = linkage['rows']
                assert math.isfinite(linkage['L'])
                assert abs(sum(r['contribution'] for r in terms) - linkage['L']) <= .0005 * (len(terms) + 1) + 1e-8
                if item.get('kind') == 'switch':
                    sites = {(s['chain'], s['resnum']) for s in item['sites']}
                    his = [r for r in terms if r['resname'] == 'HIS' and (r['chain'], r['resnum']) in sites]
                    assert len(his) == len(sites), 'selected His missing from linkage terms'
            elif item.get('kept'):
                assert item.get('linkage_status') == 'incomplete_selected_histidines'
                assert item.get('missing_histidines') and 'partial_linkage' in item
                assert 'delta_L_vs_reference' not in item
        checks.append(dict(backbone=design['backbone'], files_verified=len(files),
                           max_backbone_displacement_A=max_distance,
                           unchanged_target_sequence=True, mutations_verified=True,
                           source_output_and_code_sha256_verified=True,
                           energy_gates_and_linkage_sums_verified=True))
    save_json(out / 'integrity.json', dict(status='passed', designs=checks))
    print('Integrity passed for %d backbones, %d PDBs' % (
        len(checks), sum(r['files_verified'] for r in checks)))


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('out_dir', type=Path)
    verify(ap.parse_args().out_dir)
