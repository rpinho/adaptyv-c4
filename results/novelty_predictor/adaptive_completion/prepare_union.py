"""Deduplicate actual domain+whole-parent candidates for unresolved local levels."""
import collections
from datetime import datetime,timezone
import json
from pathlib import Path
import time

from native import HERE,ROOT,CACHE,INCOMING,BASE,DOMAIN,PARENT,BINARY,n,CandidateDB,lookup,files,source_audit


def main():
    if (HERE/'CANDIDATE_UNION.json').exists():raise FileExistsError('Candidate union already frozen')
    started=time.monotonic()
    check=json.loads((HERE/'SEED_CHECK.json').read_text());assert check['state']=='passed'
    source=source_audit();assert source['binary_sha256']==check['source']['binary_sha256']
    base=json.loads(BASE.read_text());validation=json.loads((INCOMING/'VALIDATION.json').read_text())
    assert validation['state']=='complete' and n.digest(BASE)==validation['final_report_sha256']
    assert base['structural_search']['positive_evidence_only']
    assert base['input_csv_sha256']==n.digest(INCOMING/'analysis_inputs.csv')
    by_name={r['name']:r for r in base['designs']}
    known={name for name,r in by_name.items() if r.get('local_predicted_level') is not None}
    assert all(by_name[name]['local_predicted_level'] in (1,2) for name in known)
    domains=json.loads((DOMAIN/'query_domains.json').read_text())
    keys=lookup(DOMAIN/'querydb');assert set(keys)==set(domains)
    expected=base['parent_context_augmentation']['query_db_fingerprints']
    assert files(DOMAIN/'querydb')==expected
    parent_checkpoint=PARENT/'DISCOVERY_CHECKPOINT.json'
    parent_report=json.loads((INCOMING/'run_bounded/whole_chain/predictions.json').read_text())
    assert parent_report['structural_search']['reused_discovery_checkpoint_sha256']==n.digest(parent_checkpoint)
    saved=json.loads(parent_checkpoint.read_text())
    for relative,fingerprint in saved['files'].items():assert n.digest(PARENT/relative)==fingerprint
    parent_domains=json.loads((PARENT/'query_domains.json').read_text())
    parent_keys=lookup(PARENT/'querydb')
    parent_by_name={}
    for q,domain in parent_domains.items():
        name=domain['name'];assert name not in parent_by_name
        assert n.domain_positions(domain,by_name[name]['length'])==list(range(1,by_name[name]['length']+1))
        parent_by_name[name]=q
    assert set(parent_by_name)==set(by_name)
    reference=base['structure_database'];target=Path(reference['prefix'])
    assert reference==json.loads((ROOT/'results/novelty_predictor/structural/DATABASE.json').read_text())
    for suffix,fingerprint in reference['files'].items():assert n.digest(str(target)+suffix)==fingerprint
    selected={q:d for q,d in domains.items() if d['name'] not in known}
    old_raw=BASE.parent/'structure/hits.tsv'
    assert n.digest(old_raw)==base['structural_search']['raw_hits_sha256']
    positive=n.parse_structure_hits(old_raw,domains)
    wanted={h['target'] for q in selected for h in positive[q]}
    positive_target_keys={}
    with Path(str(target)+'.lookup').open() as stream:
        for line in stream:
            key,name,_=line.rstrip('\n').split('\t')
            if name in wanted:
                if name in positive_target_keys:raise ValueError('Ambiguous reference header')
                positive_target_keys[name]=int(key)
    assert set(positive_target_keys)==wanted
    domain_candidates=CandidateDB(DOMAIN/'discovered')
    parent_candidates=CandidateDB(PARENT/'discovered')
    seeds=CACHE/'union_candidates.tsv';per_domain={};total_raw=total_unique=0
    with seeds.open('w') as stream:
        for q,domain in sorted(selected.items(),key=lambda item:keys[item[0]]):
            own_rows=domain_candidates.records(keys[q])
            parent_q=parent_by_name[domain['name']]
            parent_rows=parent_candidates.records(parent_keys[parent_q])
            own={int(row.split(b'\t',1)[0]) for row in own_rows}
            inherited={int(row.split(b'\t',1)[0]) for row in parent_rows}
            refs=own|inherited
            assert {positive_target_keys[h['target']] for h in positive[q]}<=refs,'Existing positive evidence absent from actual candidate union'
            for ref in sorted(refs):stream.write(f'{keys[q]}\t{ref}\t0\t0\n')
            total_raw+=len(own_rows)+len(parent_rows);total_unique+=len(refs)
            per_domain[q]=dict(name=domain['name'],native_query_key=keys[q],own_discovery_rows=len(own_rows),
                parent_discovery_rows=len(parent_rows),unique_domain_targets=len(own),
                unique_parent_targets=len(inherited),unique_union_targets=len(refs))
    domain_candidates.close();parent_candidates.close()
    result=dict(state='prepared',created_at=datetime.now(timezone.utc).isoformat(),
        base_report=str(BASE),base_report_sha256=n.digest(BASE),source=source,
        input_csv_sha256=base['input_csv_sha256'],reference_database=reference,
        query_db=str(DOMAIN/'querydb'),query_db_files=files(DOMAIN/'querydb'),
        domain_discovery_files={suffix:n.digest(DOMAIN/('discovered'+suffix)) for suffix in ('','.index','.dbtype')},
        parent_discovery_checkpoint_sha256=n.digest(parent_checkpoint),
        query_domains=str(DOMAIN/'query_domains.json'),query_domains_sha256=n.digest(DOMAIN/'query_domains.json'),
        preserved_proven_point_inputs=len(known),unresolved_inputs=len(by_name)-len(known),
        inputs_with_unresolved_accepted_domains=len({d['name'] for d in selected.values()}),
        unresolved_query_domains=len(selected),no_accepted_domain_inputs=sum(not r.get('structure_evidence',{}).get('domains') for r in by_name.values()),
        original_selected_candidate_rows=total_raw,unique_candidate_pairs=total_unique,
        seed_tsv=str(seeds),seed_tsv_sha256=n.digest(seeds),per_domain=per_domain,
        elapsed_seconds=time.monotonic()-started,script_sha256=n.digest(__file__),
        new_gpu_or_remote_cpu_spend_usd=0,
        policy='All independent-domain plus complete whole-parent discovery targets for unresolved inputs, unique per query/target. Retain proven points unchanged.')
    n.write_json(HERE/'CANDIDATE_UNION.json',result)
    print(json.dumps({k:result[k] for k in ('preserved_proven_point_inputs','inputs_with_unresolved_accepted_domains','unresolved_query_domains','original_selected_candidate_rows','unique_candidate_pairs')}))


if __name__=='__main__':main()
