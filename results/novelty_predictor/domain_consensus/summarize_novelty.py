"""Check real domain-search evidence, then compare both fixed local policies."""
from collections import Counter
import argparse
import csv
import json
from pathlib import Path
import sys

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
sys.path.insert(0,str(ROOT/'scripts'))
import novelty_level_predictor as n


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--parent-context',action='store_true')
    args=ap.parse_args()
    old_path=ROOT/'results/novelty_predictor/patent_extension/run/predictions.json'
    run_folder=HERE/('parent_context_run' if args.parent_context else 'novelty_run')
    new_path=run_folder/'predictions.json'
    old=json.loads(old_path.read_text());new=json.loads(new_path.read_text())
    assert old['input_csv_sha256']==new['input_csv_sha256']
    assert new['sequence_report_sha256']==n.digest(old_path)
    assert new['predictor_sha256']==n.digest(ROOT/'scripts/novelty_level_predictor.py')
    assert new['structure_manifest_sha256']==n.digest(HERE/'STRUCTURES.json')
    assert old['structure_database']==new['structure_database']
    assert old['constants']==new['constants']
    for key in ('prefilter_modes','hit_cap_per_discovery_mode','exact_tmscore_export','evalue'):
        assert old['structural_search'][key]==new['structural_search'][key]
    assert new['structural_search']['exact_tmscore_export']==1
    assert old['structural_search']['query_domains']==68 and new['structural_search']['query_domains']==70
    assert n.digest(run_folder/'structure/hits.tsv')==new['structural_search']['raw_hits_sha256']
    raw_domains=json.loads((run_folder/'structure/query_domains.json').read_text())
    hits=n.parse_structure_hits(run_folder/'structure/hits.tsv',raw_domains)
    old_by_name={r['name']:r for r in old['designs']};new_by_name={r['name']:r for r in new['designs']}
    assert set(old_by_name)==set(new_by_name) and len(new_by_name)==75
    retired={r['name'] for r in csv.DictReader((ROOT/'results/diversity_review/retired_eight_metric_pool.csv').open())}
    comparisons=[];changed=[];verified_domains=0
    for before,after in zip(old['designs'],new['designs']):
        assert before['name']==after['name']
        for field in ('sequence_sha256','best_hit','top_identity_hits','best_detected_identity',
                      'hit_count','sequence_evidence_class','sequence_database_gaps'):
            assert before[field]==after[field],(after['name'],field)
        assert after['organizer_level_established'] is False
        if before['structure_evidence']['status']!='coordinates_verified':
            assert after['structure_evidence']['status']=='coordinates_unavailable'
            continue
        a,b=before['structure_evidence'],after['structure_evidence']
        assert {k:v for k,v in a['input'].items() if k!='domains'}=={k:v for k,v in b['input'].items() if k!='domains'}
        assert n.digest(b['input']['path'])==b['input']['sha256']
        matches=[];details=[]
        for d in b['domains']:
            hs=hits[d['query_id']]
            assert d['best_hit']==(hs[0] if hs else None)
            positions=n.domain_positions(d,after['length'])
            assert positions==d['query_residue_positions']
            best=hs[0] if hs else None
            matches.append(dict(**n.domain_coordinates(d),tm_score=best['qtmscore'] if best else 0.,normalization='query_domain_length'))
            details.append(dict(positions=positions,query_length=len(positions),consensus_quality=d['consensus_quality'],
                best_query_tm=best['qtmscore'] if best else None,best_target=best['target'] if best else None))
            verified_domains+=1
        measured=n.structural_class(after['length'],matches)
        assert measured['category']==after['structural_similarity']['category']
        assert measured['high_coverage']==after['structural_similarity']['high_coverage']
        assert measured['moderate_coverage']==after['structural_similarity']['moderate_coverage']
        seq=after['sequence_evidence_class'];seq='low' if seq=='no_significant_hit' else seq
        assert n.general_level(seq,measured['category'])==after['local_predicted_level']
        item=dict(name=after['name'],length=after['length'],retired_design=after['name'] in retired,
            whole_chain_local_level=before['local_predicted_level'],consensus_domain_local_level=after['local_predicted_level'],
            changed=before['local_predicted_level']!=after['local_predicted_level'],
            domain_count=len(details),domain_coverage=after['structural_similarity']['domain_segmentation_coverage'],
            high_similarity_coverage=measured['high_coverage'],moderate_similarity_coverage=measured['moderate_coverage'],
            whole_chain_best_query_tm=a['domains'][0]['best_hit']['qtmscore'],domain_matches=details)
        comparisons.append(item)
        if item['changed']: changed.append(item)
    assert verified_domains==70 and len(comparisons)==68
    native=new_by_name['KNOWN_TFR1_POSITIVE_CONTROL'];assert native['local_predicted_level']==1
    assert new_by_name['KNOWN_PATENT_SOURCE_POSITIVE_CONTROL']['local_predicted_level'] is None
    domain_counts=Counter(str(new_by_name[name].get('local_predicted_level')) for name in retired)
    baseline_counts=Counter(str(old_by_name[name].get('local_predicted_level')) for name in retired)
    result=dict(status='completed_local_three_method_domain_comparison',
        retired_design_count=73,whole_chain_local_counts=dict(baseline_counts),consensus_domain_local_counts=dict(domain_counts),
        changed_retired_calls=sum(item['retired_design'] for item in changed),changed_calls=changed,comparisons=comparisons,
        native_control=dict(local_level=1,length=native['length'],structure=native['structural_similarity'],
            domains=native['structure_evidence']['domains']),
        actual_foldseek_sequences=json.loads((HERE/'QUERY_DATABASE_CHECK.json').read_text())['actual_foldseek_entries_verified'],
        verified_query_domains=verified_domains,search=new['structural_search'],
        organizer_levels_established=0,local_policy='Checked consensus domains for every input; no favorable per-design policy selection',
        invariants=dict(sequence_evidence_unchanged=True,full_coordinate_sources_unchanged=True,
            reference_database_unchanged=True,named_thresholds_unchanged=True,search_settings_unchanged=True,
            domain_query_sequences_and_exact_score_denominators_verified=True),
        source_sha256=dict(old_report=n.digest(old_path),new_report=n.digest(new_path),
            comparison_policy=n.digest(HERE/'COMPARISON_POLICY.md'),indexing_checks=n.digest(HERE/'INDEXING_CHECKS.json')),
        spend=dict(runpod_recorded_usd=56.5899880444631,our_modal_estimate_usd=10.718947158590089,
            combined_recorded_usd=67.30893520305318,new_gpu_or_cloud_cpu_usd=0))
    if args.parent_context:
        previous_path=HERE/'novelty_run/predictions.json'
        previous=json.loads(previous_path.read_text())
        provenance=new['parent_context_augmentation']
        assert provenance['base_report_sha256']==n.digest(previous_path)
        assert provenance['parent_report_sha256']==n.digest(old_path)
        assert provenance['script_sha256']==n.digest(ROOT/'scripts/augment_novelty_parent_context.py')
        assert provenance['applies_to_all_query_domains']==70
        assert n.digest(run_folder/'structure/parent_candidates.tsv')==provenance['raw_seed_tsv_sha256']
        assert n.digest(run_folder/'structure/parent_hits.tsv')==provenance['raw_parent_hits_sha256']
        for path,fingerprint in provenance['query_db_fingerprints'].items():assert n.digest(path)==fingerprint
        prior_domains=json.loads((HERE/'novelty_run/structure/query_domains.json').read_text())
        assert raw_domains==prior_domains
        prior_hits=n.parse_structure_hits(HERE/'novelty_run/structure/hits.tsv',prior_domains)
        no_hit_before=sum(not prior_hits[q] for q in prior_domains)
        no_hit_after=sum(not hits[q] for q in prior_domains)
        for q in prior_domains:
            if prior_hits[q]:assert hits[q][0]['qtmscore']>=prior_hits[q][0]['qtmscore']
        previous_names={r['name']:r for r in previous['designs']}
        changed_from_domain_only=[]
        for design in new['designs']:
            earlier=previous_names[design['name']]
            assert design['sequence_sha256']==earlier['sequence_sha256']
            if 'input' in design['structure_evidence']:
                assert design['structure_evidence']['input']==earlier['structure_evidence']['input']
            else:
                assert design['structure_evidence']==earlier['structure_evidence']
            before_level=earlier.get('local_predicted_level');after_level=design.get('local_predicted_level')
            if before_level is not None:assert after_level<=before_level
            if before_level!=after_level:
                changed_from_domain_only.append(dict(name=design['name'],before=before_level,after=after_level))
        mechanisms=Counter()
        for item in comparisons:
            if item['consensus_domain_local_level']!=4 or not item['retired_design']:continue
            if any((d['best_query_tm'] or 0)>=.5 for d in item['domain_matches']):
                mechanisms['matching_domains_below_whole_query_coverage_rule']+=1
            elif any(d['best_query_tm'] is None for d in item['domain_matches']):
                mechanisms['no_candidate_match_finite_search']+=1
            else:mechanisms['returned_domain_scores_below_moderate_TM']+=1
        result.update(parent_context_augmentation=provenance,domain_no_hit_before=no_hit_before,domain_no_hit_after=no_hit_after,
            changed_from_domain_only=changed_from_domain_only,local4_evidence_mechanisms=dict(mechanisms),
            prior_best_scores_preserved=True,local_level_cannot_increase_with_more_detected_similarity_verified=True)
    output=HERE/('PARENT_CONTEXT_VALIDATION.json' if args.parent_context else 'VALIDATION.json')
    n.write_json(output,result)
    columns=['name','length','whole_chain_local_level','consensus_domain_local_level','changed','domain_count',
             'domain_coverage','high_similarity_coverage','moderate_similarity_coverage','whole_chain_best_query_tm','domain_matches']
    with (HERE/('PARENT_CONTEXT_COMPARISON.csv' if args.parent_context else 'COMPARISON.csv')).open('w',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=columns);writer.writeheader()
        for item in comparisons:
            row={k:item[k] for k in columns};row['domain_matches']=json.dumps(row['domain_matches'],separators=(',',':'))
            writer.writerow(row)
    print(json.dumps(dict(retired_whole_chain=dict(baseline_counts),retired_consensus_domains=dict(domain_counts),
        changed_retired_calls=result['changed_retired_calls'],verified_query_domains=verified_domains,
        native_control_local_level=1,organizer_levels_established=0,
        domain_no_hit_after=result.get('domain_no_hit_after'),parent_context=args.parent_context)))


if __name__=='__main__':main()
