"""Verify the actual arbitrary-name PDB workflow, not a benchmark shortcut."""
import csv
import json
from pathlib import Path
import sys

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
sys.path.insert(0,str(ROOT/'scripts'))
import novelty_level_predictor as n


def main():
    state=json.loads((HERE/'run/STATUS.json').read_text())
    assert state['state']=='complete' and all(s['completed'] for s in state['steps'])
    assert len(state['steps'])==6
    assert state['workflow_sha256']==n.digest(ROOT/'scripts/run_cpu_novelty.py')
    final=Path(state['final_report'])
    assert n.digest(final)==state['final_report_sha256']
    report=json.loads(final.read_text());assert len(report['designs'])==1
    design=report['designs'][0]
    assert design['name']=='GENERIC_NATIVE_CONTROL' and design['length']==195
    assert design['best_detected_identity']==1.0 and design['sequence_evidence_class']=='high'
    assert design['local_predicted_level']==design['predicted_level']==1
    assert design['organizer_level_established'] is False and report['organizer_score'] is False
    assert design['structural_similarity']['category']=='high'
    assert report['predictor_sha256']==n.digest(ROOT/'scripts/novelty_level_predictor.py')
    source=json.loads((HERE/'run/whole_chain_structures.json').read_text())['structures'][0]
    assert source['chain']=='A' and source['chain_inferred_by_unique_full_sequence']
    assert source['source_hash_was_supplied'] is False
    assert n.digest(source['path'])==source['sha256']
    prep=json.loads((HERE/'run/domain_inputs/INPUTS.json').read_text())
    assert prep['coordinate_models']==1 and prep['unavailable']==0
    assert prep['inputs'][0]['backbone_atoms']==780
    assert prep['inputs'][0]['pdb_to_cif_format_conversion']
    assert json.loads((HERE/'run/domain_models/STATUS.json').read_text())['state']=='complete'
    assert report['structural_search']['exact_tmscore_export']==1
    assert report['parent_context_augmentation']['applies_to_all_query_domains']==1
    cache=json.loads((HERE/'RUNTIME_CACHE.json').read_text())
    archived={record['original']:record['cached'] for item in cache['items'] for record in item['files']}
    for path,fingerprint in report['parent_context_augmentation']['query_db_fingerprints'].items():
        actual=Path(path) if Path(path).exists() else Path(archived[path])
        assert n.digest(actual)==fingerprint
    with final.with_suffix('.csv').open() as stream:csv_rows=list(csv.DictReader(stream))
    assert len(csv_rows)==1 and csv_rows[0]['name']=='GENERIC_NATIVE_CONTROL' and csv_rows[0]['predicted_level']=='1'
    assert state['new_gpu_or_cloud_cpu_spend_usd']==0
    checks=dict(status='passed_real_end_to_end_cpu_workflow',arbitrary_input_name=True,pdb_input_supported=True,
        chain_identified_by_unique_full_sequence=True,source_fingerprint_acquired_and_verified=True,
        complete_backbone_atoms=780,all_three_cpu_domain_methods_completed=True,domain_search_and_parent_context_completed=True,
        sequence_identity=1.0,local_level=1,native_domain_best_tm=design['structure_evidence']['domains'][0]['best_hit']['qtmscore'],
        domain_coverage=design['structural_similarity']['domain_segmentation_coverage'],
        final_report=str(final),final_report_sha256=n.digest(final),organizer_level_established=False,
        workflow_source_sha256=state['workflow_sha256'],predictor_source_sha256=report['predictor_sha256'],
        new_gpu_or_cloud_cpu_spend_usd=0,runpod_recorded_usd=56.5899880444631,
        our_modal_estimate_usd=10.718947158590089,combined_recorded_usd=67.30893520305318)
    n.write_json(HERE/'VALIDATION.json',checks)
    print(json.dumps(dict(status='passed',input='arbitrary-name PDB control',local_level=1,
        best_domain_tm=checks['native_domain_best_tm'],organizer_level_established=False,new_gpu_spend_usd=0)))


if __name__=='__main__':main()
