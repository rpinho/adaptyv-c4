"""Read completed uncompressed native candidate keys; execute finite CPU tools."""
from datetime import datetime,timezone
import collections
import hashlib
import json
import os
from pathlib import Path
import signal
import subprocess
import sys
import time

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
CACHE=Path.home()/'.cache/adaptyv-c4/novelty_adaptive_completion'
INCOMING=ROOT/'results/novelty_predictor/incoming_rfd4_20261022'
BASE=INCOMING/'run_bounded/parent_context/predictions.json'
DOMAIN=INCOMING/'run_bounded/domain_search/structure'
PARENT=INCOMING/'run/whole_chain/structure'
BINARY=Path.home()/'.local/opt/foldseek/foldseek/bin/foldseek'
sys.path.insert(0,str(ROOT/'scripts'))
import novelty_level_predictor as n


def lookup(prefix):
    result={}
    with Path(str(prefix)+'.lookup').open() as stream:
        for line in stream:
            key,name,_=line.rstrip('\n').split('\t')
            if name in result:raise ValueError('Duplicate native query name')
            result[name]=int(key)
    return result


class CandidateDB:
    def __init__(self,prefix):
        self.prefix=Path(prefix)
        if int.from_bytes(Path(str(prefix)+'.dbtype').read_bytes(),'little')!=5:
            raise ValueError('Only normal, nonextended, uncompressed alignment result DBs supported')
        if not self.prefix.is_file():raise ValueError('Expected completed merged single-file candidate DB')
        self.stream=self.prefix.open('rb');size=self.prefix.stat().st_size;self.index={}
        with Path(str(prefix)+'.index').open() as stream:
            for line in stream:
                key,offset,length=map(int,line.split())
                if key in self.index or not (offset>=0 and length>=1 and offset+length<=size):
                    raise ValueError('Invalid/duplicate candidate index')
                self.index[key]=(offset,length)

    def records(self,key):
        if key not in self.index:return []
        offset,length=self.index[key];self.stream.seek(offset);block=self.stream.read(length)
        if len(block)!=length or not block.endswith(b'\0'):raise ValueError('Incomplete native candidate entry')
        rows=block[:-1].splitlines()
        for row in rows:
            if not row.split(b'\t',1)[0].isdigit():raise ValueError('Non-numeric native target key')
        return rows

    def targets(self,key):return {int(row.split(b'\t',1)[0]) for row in self.records(key)}
    def close(self):self.stream.close()


def files(prefix):
    return {str(p):n.digest(p) for p in sorted(prefix.parent.glob(prefix.name+'*')) if p.is_file()}


def stream_best(paths,domains,keep=10):
    """Exact top observed unique targets; memory bounded independently of rows.

    Once a target leaves the top set its old rank cannot beat a later, stronger
    top-set cutoff. A later row that qualifies therefore contains its own best
    relevant score. Equal-score rows keep the first payload, as in the full parser.
    """
    retained={q:{} for q in domains};worst={};counts={}
    rank=lambda r:(-r['qtmscore'],r['target'])
    for path in paths:
        count=collections.Counter();counts[str(path)]=count
        for row in n.iter_structure_hits(path,domains):
            q=row['query'];count[q]+=1;top=retained[q];old=top.get(row['target'])
            if old is not None:
                if row['qtmscore']<=old['qtmscore']:continue
                top[row['target']]=row
            elif len(top)<keep:top[row['target']]=row
            elif rank(row)<worst[q]:
                top[row['target']]=row
                del top[max(top,key=lambda target:rank(top[target]))]
            else:continue
            if len(top)==keep:worst[q]=max(map(rank,top.values()))
    return {q:sorted(top.values(),key=rank) for q,top in retained.items()},counts


def run(command,log,timeout,steps):
    step=dict(command=list(map(str,command)),started_at=datetime.now(timezone.utc).isoformat());steps.append(step)
    begin=time.monotonic()
    with log.open('w') as stream:
        process=subprocess.Popen(step['command'],stdout=stream,stderr=subprocess.STDOUT,
            env={**os.environ,'CUDA_VISIBLE_DEVICES':'','OMP_NUM_THREADS':'4','OPENBLAS_NUM_THREADS':'4'},start_new_session=True)
        step['pid']=process.pid
        n.write_json(log.parent/'PROCESS_STATUS.json',dict(state='running',steps=steps,new_gpu_or_remote_cpu_spend_usd=0))
        try:code=process.wait(timeout=timeout)
        except subprocess.TimeoutExpired:
            os.killpg(process.pid,signal.SIGTERM)
            try:process.wait(timeout=5)
            except subprocess.TimeoutExpired:os.killpg(process.pid,signal.SIGKILL);process.wait()
            raise
        if code:raise subprocess.CalledProcessError(code,step['command'])
    step.update(completed=True,elapsed_seconds=time.monotonic()-begin)
    n.write_json(log.parent/'PROCESS_STATUS.json',dict(state='complete',steps=steps,new_gpu_or_remote_cpu_spend_usd=0))


def source_audit():
    revision='463739e0014a1549a527de589102cde98f802f37'
    version=subprocess.check_output([str(BINARY),'version'],text=True).strip()
    if version!=revision:raise ValueError('Installed Foldseek commit differs from audited source')
    tree=json.loads((CACHE/'foldseek_tree.json').read_text())
    assert tree['sha']==revision
    names={'src/strucclustutils/tmalign.cpp':'tmalign.cpp','src/commons/TMaligner.cpp':'TMaligner.cpp',
        'src/commons/TMaligner.h':'TMaligner.h','lib/tmalign/TMalign.cpp':'TMalign_library.cpp'}
    verified={}
    for original,local in names.items():
        item=next(x for x in tree['tree'] if x['path']==original);data=(CACHE/local).read_bytes()
        assert hashlib.sha1(b'blob '+str(len(data)).encode()+b'\0'+data).hexdigest()==item['sha']
        verified[original]=dict(git_blob=item['sha'],sha256=hashlib.sha256(data).hexdigest(),
            url='https://github.com/steineggerlab/foldseek/blob/'+revision+'/'+original)
    return dict(commit=revision,binary_sha256=n.digest(BINARY),source_files=verified,
        rationale='Native TM-align candidate parser reads target keys only; full pairwise align receives coordinates/sequences, not discovery seeds; no RNG found in the audited TMalign library.')
