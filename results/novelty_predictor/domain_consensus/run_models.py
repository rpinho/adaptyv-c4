"""Finite sequential local CPU control or complete supplied-coordinate census."""
from datetime import datetime,timezone
import hashlib
import json
import os
from pathlib import Path
import shutil
import signal
import subprocess
import sys
import time

HERE=Path(__file__).resolve().parent
CACHE=Path.home()/'.cache/adaptyv-c4/novelty_domain_source_audit'
PY=CACHE/'venv/bin/python'
SOURCE=CACHE/'ted-pinned'


def sha(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def verify_rows(path,expected):
    rows=[line.split('\t') for line in path.read_text().splitlines() if line]
    assert len(rows)==len(expected), 'Missing/extra predictor rows; exit zero is insufficient'
    assert {r[0] for r in rows}==set(expected)
    for row in rows:
        query,md5,nres,ndom,chopping,score=row
        spec=expected[query]
        assert int(nres)==spec['length']
        assert md5==spec['sequence_md5']
        assert 0<=float(score)<=1
        if chopping in ('NULL','NO_SS'):
            assert int(ndom)==0
            continue
        domains=chopping.split(','); assert len(domains)==int(ndom)
        occupied=set()
        for domain in domains:
            for fragment in domain.split('_'):
                lo,hi=map(int,fragment.split('-'))
                assert 1<=lo<=hi<=int(nres), 'Out-of-range domain coordinates'
                span=set(range(lo,hi+1))
                assert not span&occupied,'Overlapping output domains/fragments'
                occupied|=span
    return len(rows)


def main():
    phase=sys.argv[1]
    if phase not in ('control','all'): raise ValueError('Expected control or all')
    if phase=='all':
        preflight=json.loads((HERE/'control/STATUS.json').read_text())
        assert preflight['state']=='complete'
        assert json.loads((HERE/'INDEXING_CHECKS.json').read_text())['all_passed']
    output=HERE/phase
    output.mkdir()
    inputs=output/'inputs'; inputs.mkdir()
    manifest=json.loads((HERE/'INPUTS.json').read_text())
    available=[r for r in manifest['inputs'] if r['status']=='full_backbone_verified'
               and (phase=='all' or r['name']=='KNOWN_TFR1_POSITIVE_CONTROL')]
    expected={}
    for item in available:
        assert sha(item['path'])==item['sha256']
        target=inputs/(item['query']+'.pdb'); shutil.copy2(item['path'],target)
        expected[item['query']]=dict(**item)
    import csv
    sequences={r['name']:r['sequence'] for r in csv.DictReader((HERE.parents[2]/'results/novelty_predictor/analysis_inputs.csv').open())}
    for value in expected.values(): value['sequence_md5']=hashlib.md5(sequences[value['name']].encode()).hexdigest()
    targets=output/'targets.txt'; targets.write_text(''.join(str(inputs/(q+'.pdb'))+'\n' for q in expected))
    state=dict(state='prepared',phase=phase,controller_pid=os.getpid(),inputs=len(expected),steps=[],new_gpu_or_cloud_cpu_spend_usd=0)
    def save(): (output/'STATUS.json').write_text(json.dumps(state,indent=2)+'\n')
    def run(name,command,cwd=None):
        step=dict(name=name,command=list(map(str,command)),started_at=datetime.now(timezone.utc).isoformat())
        state['steps'].append(step); state['state']='running'; save(); start=time.monotonic()
        with (output/(name+'.log')).open('w') as log:
            process=subprocess.Popen(step['command'],cwd=cwd,env={**os.environ,'PATH':str(PY.parent)+os.pathsep+os.environ['PATH'],
                'CUDA_VISIBLE_DEVICES':'','OMP_NUM_THREADS':'4','OPENBLAS_NUM_THREADS':'4','MKL_NUM_THREADS':'4','VECLIB_MAXIMUM_THREADS':'4'},
                stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
            step['pid']=process.pid; save()
            try:
                code=process.wait(timeout=1200)
            except subprocess.TimeoutExpired:
                os.killpg(process.pid,signal.SIGTERM)
                try: process.wait(timeout=5)
                except subprocess.TimeoutExpired:
                    os.killpg(process.pid,signal.SIGKILL); process.wait()
                raise
            if code: raise subprocess.CalledProcessError(code,step['command'])
        step.update(completed=True,elapsed_seconds=time.monotonic()-start); save()
    try:
        # Chainsaw takes a directory; run it before UniDoc creates chopped PDBs.
        # Merizo still precedes UniDoc as required by the published dependency.
        for method in ('chainsaw','merizo','unidoc'):
            raw=output/('raw_'+method+'.txt')
            filtered=output/('chopping_'+method+'.txt')
            if method=='chainsaw':
                command=[PY,HERE/'cpu_entry.py',method,'--structure_directory',inputs,'--output',raw,
                         '--save_dir',output/'chainsaw_cache']
            elif method=='merizo':
                command=[PY,HERE/'cpu_entry.py',method,'-l',targets,'--out',raw,'-d','cpu']
            else:
                command=[PY,CACHE/'native_build/unidoc/Run_UniDoc_from_scratch_structure_afdb.py',
                    '-l',targets,'--out',raw,'--inherit_chopping',output/'chopping_merizo.txt']
            run(method,command)
            run('filter_'+method,[PY,SOURCE/'scripts/filter_domains.py',raw,'-o',filtered,
                '--offset_resi','1' if method=='chainsaw' else '0'])
            count=verify_rows(filtered,expected)
            state['steps'][-1]['verified_rows']=count; save()
        assert json.loads((HERE/'INDEXING_CHECKS.json').read_text())['all_passed']
        run('consensus',[PY,CACHE/'consensus_indexing_compat/get_consensus.py','-c',
            *[output/('chopping_'+method+'.txt') for method in ('merizo','chainsaw','unidoc')],
            '-o',output/'raw_consensus.tsv'])
        run('filter_consensus',[PY,CACHE/'consensus_indexing_compat/filter_domains_consensus.py',output/'raw_consensus.tsv',
            '-o',output/'consensus.tsv'])
        rows=[line.split('\t') for line in (output/'consensus.tsv').read_text().splitlines() if line]
        assert len(rows)==len(expected) and {r[0] for r in rows}==set(expected)
        for item in available: assert sha(inputs/(item['query']+'.pdb'))==item['sha256'],'Original model input was mutated'
        state.update(state='complete',original_inputs_unchanged=True,
            outputs_sha256={p.name:sha(p) for p in output.glob('*.txt')},consensus_sha256=sha(output/'consensus.tsv'))
        save(); print(json.dumps(dict(phase=phase,state='complete',inputs=len(expected))))
    except Exception as error:
        state.update(state='failed',error=repr(error)); save(); raise


if __name__=='__main__': main()
