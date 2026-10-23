#!/usr/bin/env python3
"""CPU-only additional domain TM-align candidates from frozen whole-chain hits.

Preserves original alignments and imports the production novelty rule evaluator.
All output remains provisional; this is not an official novelty score.
"""
import argparse
import copy
from datetime import datetime,timezone
import json
import os
from pathlib import Path
import shutil
import signal
import subprocess
import time

import novelty_level_predictor as n


def report_and_hits(path):
    report=json.loads(path.read_text())
    folder=path.parent/'structure'
    raw=folder/'hits.tsv'; domains=json.loads((folder/'query_domains.json').read_text())
    if n.digest(raw)!=report['structural_search']['raw_hits_sha256']:
        raise ValueError('Raw structure evidence fingerprint differs')
    return report,domains,n.parse_structure_hits(raw,domains),raw


def main(argv=None):
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--base-report',type=Path,required=True)
    ap.add_argument('--parent-report',type=Path,required=True)
    ap.add_argument('--query-db',type=Path)
    ap.add_argument('--foldseek',default='foldseek')
    ap.add_argument('--threads',type=int,default=4)
    ap.add_argument('--timeout-stage',type=int,default=300)
    ap.add_argument('--out',type=Path,required=True)
    args=ap.parse_args(argv)
    if not 1<=args.threads<=4 or args.timeout_stage<=0:ap.error('CPU budget:1–4 threads and positive stage timeout')
    if args.out.exists():raise FileExistsError('Fresh output directory required')
    binary=shutil.which(args.foldseek)
    if not binary:raise FileNotFoundError('Foldseek required; no substitute aligner')
    base,domains,base_hits,base_raw=report_and_hits(args.base_report)
    parent,parent_domains,parent_hits,parent_raw=report_and_hits(args.parent_report)
    if base.get('parent_context_augmentation'):raise ValueError('No repeat augmentation of the same evidence')
    if base['input_csv_sha256']!=parent['input_csv_sha256'] or base['structure_database']!=parent['structure_database']:
        raise ValueError('Parent must use identical full sequences and reference snapshot')
    if base['predictor_sha256']!=n.digest(n.__file__):raise ValueError('Base predictor source differs; inspect provenance first')
    if base['structural_search']['exact_tmscore_export']!=1 or parent['structural_search']['exact_tmscore_export']!=1:
        raise ValueError('Exact exported TM evidence required')
    old_by_name={r['name']:r for r in parent['designs']}
    base_by_name={r['name']:r for r in base['designs']}
    if set(old_by_name)!=set(base_by_name):raise ValueError('Input names differ')
    target=Path(base['structure_database']['prefix']).expanduser()
    for suffix,fingerprint in base['structure_database']['files'].items():
        if n.digest(str(target)+suffix)!=fingerprint:raise ValueError('Reference DB fingerprint differs')
    query=args.query_db or args.base_report.parent/'structure/querydb'
    query=query.resolve()
    query_keys={}
    with Path(str(query)+'.lookup').open() as stream:
        for line in stream:
            key,name,_=line.rstrip('\n').split('\t');query_keys[name]=int(key)
    if set(query_keys)!=set(domains):raise ValueError('Query DB names differ from domain evidence')
    query_files={str(p):n.digest(p) for p in sorted(query.parent.glob(query.name+'*')) if p.is_file()}
    candidates={};wanted=set()
    for q,d in domains.items():
        current=base_by_name[d['name']]; previous=old_by_name[d['name']]
        if current['sequence_sha256']!=previous['sequence_sha256'] or current['length']!=previous['length']:
            raise ValueError('Protein sequence differs')
        a=current['structure_evidence']['input'];b=previous['structure_evidence']['input']
        if {k:v for k,v in a.items() if k!='domains'}!={k:v for k,v in b.items() if k!='domains'}:
            raise ValueError('Parent coordinate source differs')
        if n.digest(a['path'])!=a['sha256']:raise ValueError('Coordinate fingerprint differs')
        previous_domains=previous['structure_evidence']['domains']
        if len(previous_domains)!=1 or n.domain_positions(previous_domains[0],current['length'])!=list(range(1,current['length']+1)):
            raise ValueError('Parent evidence must cover the complete same protein')
        refs={h['target'] for h in parent_hits[previous_domains[0]['query_id']]}
        candidates[q]=refs;wanted|=refs
    target_keys={}
    with Path(str(target)+'.lookup').open() as stream:
        for line in stream:
            key,name,_=line.rstrip('\n').split('\t')
            if name in wanted:
                if name in target_keys:raise ValueError('Ambiguous target identifier')
                target_keys[name]=int(key)
    if set(target_keys)!=wanted:raise ValueError('Parent targets absent from verified reference lookup')
    args.out.mkdir(parents=True); folder=args.out/'structure';folder.mkdir()
    n.write_json(folder/'query_domains.json',domains)
    seeds=folder/'parent_candidates.tsv'
    with seeds.open('w') as stream:
        for q in sorted(candidates,key=query_keys.get):
            for ref in sorted(candidates[q],key=target_keys.get):
                stream.write(f'{query_keys[q]}\t{target_keys[ref]}\t0\t0\n')
    seed_db,aligned=folder/'parent_prefilter',folder/'parent_tmaligned'
    extra=folder/'parent_hits.tsv'
    commands=[
        [binary,'tsv2db',str(seeds),str(seed_db),'--output-dbtype','7'],
        [binary,'tmalign',str(query),str(target),str(seed_db),str(aligned),
         '--threads',str(args.threads),'--tmalign-fast','0','--tmscore-threshold','0',
         '--tmscore-threshold-mode','1','--tmalign-hit-order','1','-a','1'],
        [binary,'convertalis',str(query),str(target),str(aligned),str(extra),
         '--threads',str(args.threads),'--exact-tmscore','1','--format-output',','.join(n.FOLDSEEK_FIELDS)]]
    start=time.monotonic();steps=[]
    for i,cmd in enumerate(commands):
        step=dict(command=cmd,started_at=datetime.now(timezone.utc).isoformat());steps.append(step)
        n.write_json(args.out/'STATUS.json',dict(state='running',steps=steps,new_gpu_or_cloud_cpu_spend_usd=0))
        with (folder/f'stage_{i}.log').open('w') as log:
            process=subprocess.Popen(cmd,stdout=log,stderr=subprocess.STDOUT,start_new_session=True,
                env={**os.environ,'CUDA_VISIBLE_DEVICES':'','OMP_NUM_THREADS':str(args.threads)})
            step['pid']=process.pid
            try:code=process.wait(timeout=args.timeout_stage)
            except subprocess.TimeoutExpired:
                os.killpg(process.pid,signal.SIGTERM)
                try:process.wait(timeout=5)
                except subprocess.TimeoutExpired:os.killpg(process.pid,signal.SIGKILL);process.wait()
                raise
            if code:raise subprocess.CalledProcessError(code,cmd)
        step['completed']=True
    combined=folder/'hits.tsv'
    with combined.open('wb') as dst:
        for source in (base_raw,extra):
            with source.open('rb') as src:shutil.copyfileobj(src,dst)
    hits=n.parse_structure_hits(combined,domains)
    result=copy.deepcopy(base)
    for design in result['designs']:
        if design['structure_evidence']['status']!='coordinates_verified':continue
        matches=[]
        for d in design['structure_evidence']['domains']:
            hs=hits[d['query_id']]
            d.update(best_hit=hs[0] if hs else None,top_hits=hs[:10],hit_count=len(hs))
            matches.append(dict(**n.domain_coordinates(d),tm_score=hs[0]['qtmscore'] if hs else 0.,normalization='query_domain_length'))
        measured=n.structural_class(design['length'],matches)
        measured.update(domain_segmentation_coverage=sum(len(n.domain_positions(d,design['length'])) for d in matches)/design['length'],
            evidence_status='matches_measured' if any(d['best_hit'] for d in design['structure_evidence']['domains']) else 'no_candidate_match_detected_not_zero_TM',
            tm_boundary_precision_caveat=any((d['best_hit'] or {}).get('qtmscore') in (.5,.8) for d in design['structure_evidence']['domains']))
        observed=design['sequence_evidence_class'];observed='low' if observed=='no_significant_hit' else observed
        level=n.general_level(observed,measured['category'])
        design.update(structural_similarity=measured,local_predicted_level=level,predicted_level=level,
            possible_levels_if_structure_policy_accepted=n.possible_levels(observed,design['sequence_database_gaps'],measured['category']),
            local_gate_assessment='provisional_local_level_'+str(level),organizer_level_established=False)
    for q,old_hits in base_hits.items():
        if old_hits and hits[q][0]['qtmscore']<old_hits[0]['qtmscore']:raise ValueError('Original best evidence was lost')
    provenance=dict(base_report_sha256=n.digest(args.base_report),parent_report_sha256=n.digest(args.parent_report),
        parent_raw_hits_sha256=n.digest(parent_raw),base_raw_hits_sha256=n.digest(base_raw),query_db_prefix=str(query),
        query_db_fingerprints=query_files,target_lookup_sha256=n.digest(str(target)+'.lookup'),
        additional_target_pairs=sum(map(len,candidates.values())),unique_parent_reference_targets=len(wanted),
        applies_to_all_query_domains=len(candidates),commands=commands,elapsed_seconds=time.monotonic()-start,
        raw_parent_hits_sha256=n.digest(extra),raw_seed_tsv_sha256=n.digest(seeds),
        script_sha256=n.digest(__file__),policy='Add every same-parent whole-chain hit to every domain; retain original domain alignments; deduplicate by best query TM')
    result.update(parent_context_augmentation=provenance,generated_at=datetime.now(timezone.utc).isoformat())
    result['structural_search'].update(raw_hits_sha256=n.digest(combined),
        raw_alignment_rows=sum(1 for _ in combined.open()),hit_rows=sum(map(len,hits.values())),
        parent_context_candidate_augmentation=True)
    n.save_report(result,args.out)
    n.write_json(args.out/'STATUS.json',dict(state='complete',steps=steps,new_gpu_or_cloud_cpu_spend_usd=0))


if __name__=='__main__':main()
