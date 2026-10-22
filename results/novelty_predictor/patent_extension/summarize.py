"""Verify real patent evidence and compare unchanged retired-pool local calls."""
from collections import Counter
import csv
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT/'scripts'))
import novelty_level_predictor as n


def main():
    scope=json.loads((HERE/'INPUT_SCOPE.json').read_text())
    acquisition=json.loads((HERE/'ACQUISITION.json').read_text())
    before=json.loads((ROOT/'results/novelty_predictor/exact_score/run/predictions.json').read_text())
    after=json.loads((HERE/'run/predictions.json').read_text())
    records=n.load_designs(HERE/'analysis_inputs.csv')
    assert n.digest(HERE/'analysis_inputs.csv')==after['input_csv_sha256']==scope['analysis_inputs_sha256']
    assert n.digest(HERE/'DATABASES.json')==after['database_manifest_sha256']==scope['database_manifest_sha256']
    assert n.digest(HERE/'STRUCTURES.json')==after['structure_manifest_sha256']==scope['structure_manifest_sha256']
    assert n.digest(HERE/'PROTOCOL.md')==scope['protocol_sha256']
    assert n.digest(ROOT/'scripts/novelty_level_predictor.py')==before['predictor_sha256']==after['predictor_sha256']
    assert len(records)==len(after['designs'])==75
    old={r['name']:r for r in before['designs']}
    new={r['name']:r for r in after['designs']}
    assert set(new)-set(old)=={scope['added_control']}
    searches={d['id']:d for d in after['database_searches']}
    assert set(searches)==set(n.REQUIRED_DATABASES)
    assert all(d['search_status']=='completed' for d in searches.values())
    assert searches['patents']['scope_complete'] is False
    assert searches['patents']['fasta_sha256']==acquisition['fasta_sha256']
    for d in before['database_searches']:
        if d['id']=='patents': continue
        assert searches[d['id']]['fasta_sha256']==d['fasta_sha256']
    sequence_root=HERE/'sequence_and_partial_structure_run'
    assert n.digest(sequence_root/'predictions.json')==after['sequence_report_sha256']
    assert n.digest(ROOT/'results/novelty_predictor/exact_score/run/predictions.json')==scope['inherited_structure_report_sha256']
    assert after['structural_search']['query_domains']==scope['expected_coordinate_domains']==68
    expected_coordinates={r['name'] for r in before['designs'] if r['structure_evidence']['status']=='coordinates_verified'}
    actual_coordinates={r['name'] for r in after['designs'] if r['structure_evidence']['status']=='coordinates_verified'}
    assert expected_coordinates==actual_coordinates, 'Full coordinate scope must match frozen baseline'
    for database,d in searches.items():
        assert n.digest(sequence_root/(database+'.hits.tsv'))==d['raw_hits_sha256']
    queries={f'q{i:06d}':r for i,r in enumerate(records,1)}
    hits=n.parse_hits(sequence_root/'patents.hits.tsv','patents',queries)
    control_query=next(q for q,r in queries.items() if r['name']==scope['added_control'])
    source_id=acquisition['positive_control']['source_header'].split()[0]
    control_hits=[h for h in hits[control_query] if h['fident']==1 and h['qstart']==1
                  and h['qend']==h['qlen']==len(queries[control_query]['sequence'])
                  and h['nident']==h['alnlen']==h['qlen'] and h['target']==source_id]
    assert control_hits, 'Registered patent source sequence must match itself fully at 100% identity'
    assert before['structure_database']==after['structure_database']
    assert before['structural_search']['exact_tmscore_export']==after['structural_search']['exact_tmscore_export']==1
    changed_structure=[]; changed_level=[]; patent_matches=[]
    retired={r['name'] for r in csv.DictReader((ROOT/'results/diversity_review/retired_eight_metric_pool.csv').open())}
    for name,r in old.items():
        assert r['sequence_sha256']==new[name]['sequence_sha256']
        if r.get('structure_evidence')!=new[name].get('structure_evidence'):
            changed_structure.append(name)
        if r.get('local_predicted_level')!=new[name].get('local_predicted_level'):
            changed_level.append(dict(name=name, old_local_level=r.get('local_predicted_level'),
                                     new_local_level=new[name].get('local_predicted_level')))
    for q,hs in hits.items():
        name=queries[q]['name']
        if name in retired and hs:
            patent_matches.append(dict(name=name, hit_count=len(hs), best_hit=n.identity_ranked_hits(hs)[0]))
    counts=Counter(str(new[name].get('local_predicted_level')) for name in retired)
    positive=next(name for name in old if name not in retired)
    assert new[positive]['local_predicted_level']==1
    assert all(new[name]['organizer_level_established'] is False for name in new)
    result=dict(status='completed_cpu_only_public_patent_extension', retired_pool_local_level_counts=dict(counts),
        retired_designs=len(retired), patent_matched_retired_designs=len(patent_matches),
        patent_matches=patent_matches, changed_local_levels=changed_level,
        changed_structural_evidence_names=changed_structure,
        native_positive_control=dict(name=positive, local_level=new[positive]['local_predicted_level'],
            query_tm=new[positive]['structure_evidence']['domains'][0]['best_hit']['qtmscore']),
        patent_positive_control=dict(name=scope['added_control'], source_id=source_id, full_identity_hit=control_hits[0]),
        sequence_searches=[dict(id=k, seconds=v['elapsed_seconds'], raw_hit_rows=v['hit_rows']) for k,v in searches.items()],
        structural_search=after['structural_search'], patent_scope_complete=False,
        organizer_levels_established=0, production_source_unchanged=True,
        source_records=acquisition['exported_records'], source_residues=acquisition['exported_residues'],
        report_sha256=n.digest(HERE/'run/predictions.json'),
        spend=dict(runpod_recorded_usd=56.5899880444631, our_modal_estimate_usd=10.718947158590089,
            combined_recorded_usd=67.30893520305318, new_gpu_or_cloud_cpu_usd=0.0))
    n.write_json(HERE/'VALIDATION.json',result)
    print(json.dumps(dict(retired_local_levels=dict(counts), patent_matched_retired_designs=len(patent_matches),
        changed_local_levels=len(changed_level), changed_structural_evidence=len(changed_structure),
        patent_control_full_identity_verified=True, organizer_levels_established=0)))


if __name__=='__main__': main()
