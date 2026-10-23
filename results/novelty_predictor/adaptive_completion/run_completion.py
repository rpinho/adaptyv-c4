"""Finite detached CPU completion of unresolved structural comparisons."""
import ast
import argparse
import copy
from datetime import datetime,timezone
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import time

from native import HERE,ROOT,CACHE,BASE,BINARY,n,files,run,source_audit

FULL_ALIGNMENT_TIMEOUT=8*60*60
EXPORT_TIMEOUT=2*60*60


def stage_timeout(stage):
    return {0:1200,1:FULL_ALIGNMENT_TIMEOUT,2:EXPORT_TIMEOUT}[stage]


def export_resume(previous,commands,union_sha256,aligned):
    """Resume only a failed export backed by a recorded completed alignment."""
    if previous.get('state')!='failed' or previous.get('current_stage')!=2:
        raise ValueError('Export resume requires a failed export stage')
    if previous.get('candidate_union_sha256')!=union_sha256:
        raise ValueError('Candidate union changed since alignment')
    steps=previous.get('steps',[])
    if len(steps)!=3 or any(not steps[i].get('completed') for i in (0,1)):
        raise ValueError('Completed alignment receipt missing')
    for i in (0,1):
        if steps[i]['command']!=list(map(str,commands[i])):
            raise ValueError('Cached alignment command/input paths differ')
    # The failed export may use a different destination, but no scoring option
    # or source database may change during recovery.
    old=list(steps[2]['command']);new=list(map(str,commands[2]))
    old[5]=new[5]
    if old!=new:raise ValueError('Export scoring settings changed')
    from native import CandidateDB
    cached=CandidateDB(aligned)
    cached.close()
    fingerprints=files(aligned)
    return copy.deepcopy(steps[:2]),dict(files=fingerprints,completed_alignment_receipt=steps[1],
        scope='Reused completed native alignment; only exact-score export is repeated')


def unchanged_rules():
    old=subprocess.check_output(['git','show','6353dc8:scripts/novelty_level_predictor.py'],cwd=ROOT,text=True)
    current=Path(n.__file__).read_text()
    names={'sequence_class','general_level','domain_intervals','domain_coordinates','domain_positions',
        'structural_class','possible_levels','apply_positive_evidence_bounds'}
    def relevant(text):
        found={}
        for node in ast.parse(text).body:
            if isinstance(node,ast.FunctionDef) and node.name in names:
                found[node.name]=ast.dump(node,include_attributes=False)
            elif isinstance(node,ast.Assign):
                for target in node.targets:
                    if isinstance(target,ast.Name) and target.id.startswith(('SEQUENCE_','STRUCTURE_','MMSEQS_','FOLDSEEK_')):
                        found[target.id]=ast.dump(node.value,include_attributes=False)
        return found
    before,after=relevant(old),relevant(current)
    if before!=after:raise ValueError('Classification rules or search constants changed')
    return dict(reference_revision='6353dc8',semantic_rule_hash=hashlib.sha256(json.dumps(after,sort_keys=True).encode()).hexdigest(),
        checked_functions=sorted(names),current_predictor_sha256=n.digest(n.__file__))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--resume-export',action='store_true')
    args=parser.parse_args()
    status=HERE/'STATUS.json'
    previous=json.loads(status.read_text()) if status.exists() else None
    if previous and not args.resume_export:raise FileExistsError('This finite completion has already started; use --resume-export only for a failed export')
    if args.resume_export and not previous:raise FileNotFoundError('No failed export receipt to resume')
    state=dict(state='verifying_inputs',pid=os.getpid(),started_at=datetime.now(timezone.utc).isoformat(),
        new_gpu_or_remote_cpu_spend_usd=0,generator=False,numerical_threads=4)
    start=time.monotonic();steps=[];status_published=False
    try:
        union=json.loads((HERE/'CANDIDATE_UNION.json').read_text())
        assert union['state']=='prepared' and n.digest(BASE)==union['base_report_sha256']
        assert source_audit()['binary_sha256']==union['source']['binary_sha256']
        assert files(Path(union['query_db']))==union['query_db_files']
        assert n.digest(union['seed_tsv'])==union['seed_tsv_sha256']
        rules=unchanged_rules()
        target=Path(union['reference_database']['prefix'])
        for suffix,expected in union['reference_database']['files'].items():assert n.digest(str(target)+suffix)==expected
        output=CACHE/'full_remaining'
        if not args.resume_export:output.mkdir()
        prefilter=output/'union_prefilter';aligned=output/'full_aligned'
        hits=output/('hits_reexport.tsv' if args.resume_export else 'hits.tsv')
        if hits.exists() or list(output.glob(hits.name+'.*')):
            raise FileExistsError('Export destination already exists; preserve prior artifacts')
        commands=[
            [BINARY,'tsv2db',union['seed_tsv'],prefilter,'--output-dbtype','7'],
            [BINARY,'tmalign',union['query_db'],target,prefilter,aligned,'--threads','4','--tmalign-fast','0',
             '--tmscore-threshold','0','--tmscore-threshold-mode','1','--tmalign-hit-order','1',
             '--max-accept','2147483647','--max-rejected','2147483647','-a','1'],
            [BINARY,'convertalis',union['query_db'],target,aligned,hits,'--threads','4','--exact-tmscore','1',
             '--format-output',','.join(n.FOLDSEEK_FIELDS)]]
        first_stage=0
        if args.resume_export:
            steps,receipt=export_resume(previous,commands,n.digest(HERE/'CANDIDATE_UNION.json'),aligned)
            archive=HERE/'attempt_01_export_timeout';archive.mkdir()
            for name in ('STATUS.json','PROCESS_STATUS.json','stage_2.log'):
                path=HERE/name
                if path.exists():(archive/name).write_bytes(path.read_bytes())
            n.write_json(HERE/'EXPORT_RESUME.json',dict(recorded_at=datetime.now(timezone.utc).isoformat(),
                failed_status_sha256=n.digest(archive/'STATUS.json'),alignment_cache=receipt,
                export_timeout_seconds=EXPORT_TIMEOUT,script_sha256=n.digest(__file__),
                new_gpu_or_remote_cpu_spend_usd=0))
            state.update(resumed_export_only=True,prior_elapsed_seconds=previous['elapsed_seconds'],
                reused_alignment_files=receipt['files'])
            first_stage=2
        for i,command in enumerate(commands):
            if i<first_stage:continue
            state.update(state='running',current_stage=i,steps=steps,rule_compatibility=rules,
                candidate_union_sha256=n.digest(HERE/'CANDIDATE_UNION.json'),unique_candidate_pairs=union['unique_candidate_pairs'])
            n.write_json(status,state)
            status_published=True
            run(command,HERE/('stage_'+str(i)+'.log'),stage_timeout(i),steps)
        state.update(state='alignments_complete',steps=steps,raw_hits=str(hits),raw_hits_sha256=n.digest(hits),
            elapsed_seconds=time.monotonic()-start);n.write_json(status,state)
        finalizer=HERE/'finalize.py'
        if finalizer.exists():
            with (HERE/'finalize.log').open('w') as log:
                subprocess.run([sys.executable,str(finalizer)],stdout=log,stderr=subprocess.STDOUT,check=True,timeout=1800)
            state['state']='complete'
        else:state['state']='alignments_complete_report_pending'
        state.update(finished_at=datetime.now(timezone.utc).isoformat(),elapsed_seconds=time.monotonic()-start)
        n.write_json(status,state);print(json.dumps(state))
    except Exception as error:
        state.update(state='failed',error=repr(error),elapsed_seconds=time.monotonic()-start,steps=steps)
        n.write_json(status if status_published or previous is None else HERE/'EXPORT_RESUME_REJECTED.json',state);raise


if __name__=='__main__':main()
