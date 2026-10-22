"""Rescore the exact saved alignment and test its effect as an optimizer seed."""
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import time

HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('audit_standalone',HERE/'run_check.py')
audit=importlib.util.module_from_spec(spec);spec.loader.exec_module(audit)


def main():
    from Bio.Data.PDBData import protein_letters_3to1 as aa
    def sequence(path):
        lines=[line for line in Path(path).read_text().splitlines()
               if line.startswith('ATOM') and line[12:16].strip()=='CA']
        return ''.join(aa.get(line[17:20].strip(),'X') for line in lines)
    installation=json.loads((HERE/'INSTALLATION.json').read_text())
    binary=installation['binary']
    assert audit.sha(binary)==installation['binary_sha256']
    pairs=json.loads((HERE/'PAIRS.json').read_text())['pairs']
    references={r['database_key']:r for r in json.loads((HERE/'REFERENCE_MANIFEST.json').read_text())['records']}
    folder=HERE/'saved_alignment';folder.mkdir(exist_ok=True)
    planned=dict(registered_at=datetime.now(timezone.utc).isoformat(),trigger='Default standalone best0.4286 vs bundled0.476; same broad class but >0.01 difference',
        procedure='All27 same pairs: validate full query/target sequences against reported aligned substrings; pad only unaligned termini with gaps; run standalone-I fixed correspondence and-i starting-seed options; retain all outcomes and original free runs.',
        interpretation='Fixed recovery diagnoses scoring/coordinate consistency for a known alignment. Seeded recovery with lower free default implicates optimizer initialization. Any category crossing is reported, never threshold-adjusted.',new_gpu_spend_usd=0)
    (HERE/'SAVED_ALIGNMENT_PLAN.json').write_text(json.dumps(planned,indent=2)+'\n')
    jobs=[]
    for i,pair in enumerate(pairs,1):
        hit=pair['hit'];q=sequence(pair['query_path']);ref=references[pair['database_key']];t=sequence(ref['path'])
        assert len(q)==hit['qlen'] and len(t)==hit['tlen']
        qs,qe,ts,te=hit['qstart']-1,hit['qend'],hit['tstart']-1,hit['tend']
        assert hit['qaln'].replace('-','')==q[qs:qe],(i,'query sequence mismatch')
        assert hit['taln'].replace('-','')==t[ts:te],(i,'target sequence mismatch')
        qa=q[:qs]+'-'*ts+hit['qaln']+q[qe:]+'-'*len(t[te:])
        ta='-'*qs+t[:ts]+hit['taln']+'-'*len(q[qe:])+t[te:]
        assert len(qa)==len(ta) and qa.replace('-','')==q and ta.replace('-','')==t
        aln=folder/f'{i:03d}.fasta';aln.write_text('>query\n'+qa+'\n>reference\n'+ta+'\n')
        for flag,mode in [('-I','fixed'),('-i','seeded')]:
            jobs.append((i,pair,ref,aln,flag,mode))
    def run(job):
        i,pair,ref,aln,flag,mode=job
        cmd=[binary,pair['query_path'],ref['path'],flag,str(aln)]
        out=subprocess.run(cmd,capture_output=True,text=True,check=True,timeout=30,
            env={**os.environ,'OMP_NUM_THREADS':'1','CUDA_VISIBLE_DEVICES':''})
        path=folder/f'{i:03d}_{mode}.txt';path.write_text(out.stdout)
        scores=audit.parse_output(out.stdout)
        assert scores['structure_1']['length']==pair['hit']['qlen']
        return dict(pair_index=i,mode=mode,query_name=pair['query_name'],role=pair['role'],target=pair['hit']['target'],
            foldseek_query_TM=pair['hit']['qtmscore'],query_TM=scores['structure_1']['tm_score'],
            target_TM=scores['structure_2']['tm_score'],command=cmd,
            alignment_sha256=audit.sha(aln),output_sha256=audit.sha(path))
    start=time.monotonic()
    with ThreadPoolExecutor(max_workers=4) as pool:results=list(pool.map(run,jobs))
    retired=[r for r in results if r['role']=='remaining_local4']
    summary=dict(status='completed',comparisons=len(results),query_and_reference_sequences_verified=True,
        fixed_best_query_TM=max(r['query_TM'] for r in retired if r['mode']=='fixed'),
        seeded_best_query_TM=max(r['query_TM'] for r in retired if r['mode']=='seeded'),
        fixed_maximum_absolute_difference=max(abs(r['query_TM']-r['foldseek_query_TM']) for r in retired if r['mode']=='fixed'),
        seeded_pairs_at_least_0p5=sum(r['query_TM']>=.5 for r in retired if r['mode']=='seeded'),
        elapsed_seconds=time.monotonic()-start,results=results,new_gpu_spend_usd=0)
    (HERE/'SAVED_ALIGNMENT_RESULTS.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps({k:v for k,v in summary.items() if k!='results'}))


if __name__=='__main__':main()
