"""Verify the paired CPU run and report coordinate sensitivity without selection."""
import csv
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / 'scripts'))
import novelty_level_predictor as novelty
from run_context_ensemble import fit

HERE = Path(__file__).resolve().parent


def main():
    inventory = json.loads((HERE / 'INVENTORY.json').read_text())
    baseline_path = ROOT / 'results/novelty_predictor/exact_score/run/predictions.json'
    run_path = HERE / 'run/predictions.json'
    baseline = json.loads(baseline_path.read_text())
    run = json.loads(run_path.read_text())
    assert novelty.digest(baseline_path) == inventory['frozen_baseline_sha256']
    assert novelty.digest(HERE / 'PREREGISTRATION.md') == inventory['preregistration_sha256']
    assert baseline['input_csv_sha256'] == run['input_csv_sha256'] == inventory['input_sha256']
    assert run['structure_manifest_sha256'] == novelty.digest(HERE / 'STRUCTURES.json')
    assert run['predictor_sha256'] == baseline['predictor_sha256'] == novelty.digest(ROOT/'scripts/novelty_level_predictor.py')
    assert run['structure_database'] == baseline['structure_database']
    for key in ('prefilter_modes', 'hit_cap_per_discovery_mode', 'exact_tmscore_export', 'evalue'):
        assert run['structural_search'][key] == baseline['structural_search'][key]
    assert run['structural_search']['exact_tmscore_export'] == 1
    assert novelty.digest(HERE/'run/structure/hits.tsv') == run['structural_search']['raw_hits_sha256']
    old = {r['name']: r for r in baseline['designs']}
    new = {r['name']: r for r in run['designs']}
    assert set(old) == set(new)
    for name in old:
        for key in ('sequence_sha256', 'sequence_evidence_class', 'best_detected_identity',
                    'sequence_database_gaps', 'best_hit', 'hit_count', 'top_identity_hits'):
            assert old[name][key] == new[name][key], (name, key)
    control = inventory['native_control']
    assert new[control]['local_predicted_level'] == old[control]['local_predicted_level'] == 1
    assert new[control]['structure_evidence'] == old[control]['structure_evidence']
    assert run['structural_search']['query_domains'] == len(inventory['matched_designs']) + 1
    pairs = []
    for name in inventory['matched_designs']:
        before, after = old[name], new[name]
        a, b = before['structure_evidence'], after['structure_evidence']
        rr = []
        for evidence in (a, b):
            spec = evidence['input']
            assert novelty.digest(spec['path']) == spec['sha256']
            rr.append(novelty.ca_residues(Path(spec['path']), spec['chain']))
        assert ''.join(r[1] for r in rr[0]) == ''.join(r[1] for r in rr[1])
        assert a['input']['domains'] == b['input']['domains']
        _, _, rmsd = fit([r[2] for r in rr[1]], [r[2] for r in rr[0]])
        ah, bh = a['domains'][0]['best_hit'], b['domains'][0]['best_hit']
        pairs.append(dict(name=name, length=before['length'], sequence_sha256=before['sequence_sha256'],
            complex_source=a['input'], monomer_source=b['input'], paired_ca_rmsd_angstrom=rmsd,
            complex_best_query_tm=ah['qtmscore'], monomer_best_query_tm=bh['qtmscore'],
            complex_best_target=ah['target'], monomer_best_target=bh['target'],
            complex_local_level=before['local_predicted_level'], monomer_local_level=after['local_predicted_level'],
            level_changed=before['local_predicted_level'] != after['local_predicted_level'],
            organizer_level_established=False))
    summary = dict(status='completed_cpu_only', pairs=pairs,
        changed_level_count=sum(p['level_changed'] for p in pairs), paired_design_count=len(pairs),
        positive_control=dict(name=control, local_level=new[control]['local_predicted_level'],
            query_tm=new[control]['structure_evidence']['domains'][0]['best_hit']['qtmscore']),
        files_examined=inventory['files_examined'], parse_error_count=inventory['parse_error_count'],
        unmatched_design_count=len(old)-len(pairs)-1,
        full_pool_missing_complex_coordinates_recovered=sum(
            old[p['name']]['local_predicted_level'] is None for p in pairs),
        structural_search=run['structural_search'],
        input_fingerprints=dict(preregistration=novelty.digest(HERE/'PREREGISTRATION.md'),
            inventory=novelty.digest(HERE/'INVENTORY.json'), structures=novelty.digest(HERE/'STRUCTURES.json'),
            baseline=novelty.digest(baseline_path), run=novelty.digest(run_path)),
        checks=dict(sequence_evidence_unchanged=True, paired_full_sequences_equal=True,
            unchanged_database_and_search_settings=True, unchanged_positive_control=True,
            predictor_unchanged=True, coordinate_fingerprints_verified=True,
            all_inputs_accounted=True),
        limitations=['Four related retired designs from one parent; not a representative sample',
            'Saved monomer prediction is not organizer ESMFold/ESMFold2-Fast prediction',
            'Whole-chain domain assumption, historical structural snapshot, finite discovery, missing patents',
            'Residual standalone TM-align numeric disagreement is not resolved by this experiment',
            'No official novelty level or binding/physical qualification established'],
        spend=dict(runpod_recorded_usd=56.5899880444631, our_modal_estimate_usd=10.718947158590089,
            combined_recorded_usd=67.30893520305318, new_gpu_or_cloud_cpu_usd=0.0))
    novelty.write_json(HERE/'COMPARISON.json', summary)
    keys = ['name', 'length', 'paired_ca_rmsd_angstrom', 'complex_best_query_tm',
            'monomer_best_query_tm', 'complex_local_level', 'monomer_local_level', 'level_changed']
    with (HERE/'COMPARISON.csv').open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=keys)
        writer.writeheader()
        writer.writerows({key: p[key] for key in keys} for p in pairs)
    print(json.dumps(dict(paired_designs=len(pairs), changed_levels=summary['changed_level_count'],
                         elapsed_seconds=run['structural_search']['elapsed_seconds'],
                         recovered_missing_coordinates=summary['full_pool_missing_complex_coordinates_recovered'])))


if __name__ == '__main__':
    main()
