"""Finite detached CPU completion of unresolved structural comparisons."""
import ast
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
    status=HERE/'STATUS.json'
    if status.exists():raise FileExistsError('This finite completion has already started')
    state=dict(state='verifying_inputs',pid=os.getpid(),started_at=datetime.now(timezone.utc).isoformat(),
        new_gpu_or_remote_cpu_spend_usd=0,generator=False,numerical_threads=4)
    n.write_json(status,state);start=time.monotonic();steps=[]
    try:
        union=json.loads((HERE/'CANDIDATE_UNION.json').read_text())
        assert union['state']=='prepared' and n.digest(BASE)==union['base_report_sha256']
        assert source_audit()['binary_sha256']==union['source']['binary_sha256']
        assert files(Path(union['query_db']))==union['query_db_files']
        assert n.digest(union['seed_tsv'])==union['seed_tsv_sha256']
        rules=unchanged_rules()
        target=Path(union['reference_database']['prefix'])
        for suffix,expected in union['reference_database']['files'].items():assert n.digest(str(target)+suffix)==expected
        output=CACHE/'full_remaining';output.mkdir()
        prefilter=output/'union_prefilter';aligned=output/'full_aligned';hits=output/'hits.tsv'
        commands=[
            [BINARY,'tsv2db',union['seed_tsv'],prefilter,'--output-dbtype','7'],
            [BINARY,'tmalign',union['query_db'],target,prefilter,aligned,'--threads','4','--tmalign-fast','0',
             '--tmscore-threshold','0','--tmscore-threshold-mode','1','--tmalign-hit-order','1',
             '--max-accept','2147483647','--max-rejected','2147483647','-a','1'],
            [BINARY,'convertalis',union['query_db'],target,aligned,hits,'--threads','4','--exact-tmscore','1',
             '--format-output',','.join(n.FOLDSEEK_FIELDS)]]
        for i,command in enumerate(commands):
            state.update(state='running',current_stage=i,steps=steps,rule_compatibility=rules,
                candidate_union_sha256=n.digest(HERE/'CANDIDATE_UNION.json'),unique_candidate_pairs=union['unique_candidate_pairs'])
            n.write_json(status,state)
            run(command,HERE/('stage_'+str(i)+'.log'),FULL_ALIGNMENT_TIMEOUT if i==1 else 1200,steps)
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
        n.write_json(status,state);raise


if __name__=='__main__':main()
