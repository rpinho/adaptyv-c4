"""Inventory existing monomers; create sequence-verified comparison inputs."""
import csv
import json
import sys
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / 'scripts'))
import novelty_level_predictor as novelty
from Bio.PDB.MMCIF2Dict import MMCIF2Dict

HERE = Path(__file__).resolve().parent


def main():
    if (HERE / 'STRUCTURES.json').exists():
        raise FileExistsError('Do not overwrite registered comparison inputs')
    inputs = ROOT / 'results/novelty_predictor/analysis_inputs.csv'
    rows = list(csv.DictReader(inputs.open()))
    by_sequence = defaultdict(list)
    for row in rows:
        by_sequence[row['sequence']].append(row['name'])
    matches = defaultdict(list)
    inventory = []
    files = sorted(p for p in (ROOT / 'results').rglob('*.cif')
                   if 'monomer' in str(p.relative_to(ROOT)).lower()
                   and 'novelty_predictor' not in str(p.relative_to(ROOT)))
    for path in files:
        entry = dict(path=str(path.relative_to(ROOT)), sha256=novelty.digest(path), chains=[])
        try:
            data = MMCIF2Dict(str(path))
            chains = sorted(set(data.get('_atom_site.auth_asym_id',
                                         data.get('_atom_site.label_asym_id', []))))
            for chain in chains:
                residues = novelty.ca_residues(path, chain)
                sequence = ''.join(r[1] for r in residues)
                names = by_sequence.get(sequence, [])
                entry['chains'].append(dict(chain=chain, length=len(residues),
                                            sequence=sequence, matched_names=names))
                for name in names:
                    matches[name].append(dict(name=name, status='available',
                        path=str(path), sha256=entry['sha256'], chain=chain,
                        domains=[dict(start=1, end=len(residues),
                            source='Whole-chain single-domain assumption for short minibinder; NOT three-predictor consensus')]))
        except Exception as error:
            entry['error'] = str(error)
        inventory.append(entry)
    if any(len(values) != 1 for values in matches.values()):
        raise ValueError('Multiple matching monomers: selection not preregistered')
    baseline_path = ROOT / 'results/novelty_predictor/exact_score/run/predictions.json'
    baseline = json.loads(baseline_path.read_text())
    # The positive control is the sole nonretired analysis input, identified by
    # input scope, not by any measured score in this run.
    retired = {r['name'] for r in csv.DictReader(
        (ROOT / 'results/diversity_review/retired_eight_metric_pool.csv').open())}
    controls = [r['name'] for r in rows if r['name'] not in retired]
    if len(controls) != 1:
        raise ValueError('Expected exactly one native control')
    existing = {r['name']: r for r in baseline['designs']}
    structures = []
    for row in rows:
        name = row['name']
        if name in matches:
            structures.append(matches[name][0])
        elif name == controls[0]:
            structures.append(existing[name]['structure_evidence']['input'])
        else:
            structures.append(dict(name=name, status='unavailable',
                reason='No exact full-chain sequence match in saved monomer CIF inventory; not a new folding failure'))
    novelty.write_json(HERE / 'INVENTORY.json', dict(
        input_sha256=novelty.digest(inputs), preregistration_sha256=novelty.digest(HERE/'PREREGISTRATION.md'),
        frozen_baseline_sha256=novelty.digest(baseline_path),
        files_examined=len(inventory), parse_error_count=sum('error' in x for x in inventory),
        matched_designs=sorted(matches), native_control=controls[0], files=inventory))
    novelty.write_json(HERE / 'STRUCTURES.json', dict(structures=structures))
    print(json.dumps(dict(files_examined=len(inventory), matched_designs=len(matches),
                          native_controls=len(controls))))


if __name__ == '__main__':
    main()
