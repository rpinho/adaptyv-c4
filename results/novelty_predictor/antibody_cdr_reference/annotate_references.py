"""Resume finite local annotation of existing reference sequences, never generation."""
import collections
from datetime import datetime,timezone
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import random
import sqlite3
import subprocess
import sys
import time
import zlib

from prepare_pool import HERE,ROOT,CACHE,DB,AUDIT,sha,table_hash

os.environ.update(CUDA_VISIBLE_DEVICES='',OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',
    MKL_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1')
sys.path.insert(0,str(CACHE/'runtime'))
import numpy as np
import torch
from anarcii import Anarcii

CHUNK_SIZE=128
REFERENCE_ANNOTATION_MAX_LENGTH=2000
MAX_RUN_SECONDS=8*60*60


def save(path,value):
    temporary=path.with_suffix('.tmp')
    temporary.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n');temporary.replace(path)


def main():
    os.nice(10)
    pool=json.loads((HERE/'POOL.json').read_text())
    calibration=json.loads((AUDIT/'RUNTIME.json').read_text())
    state=dict(state='verifying',pid=os.getpid(),started_at=datetime.now(timezone.utc).isoformat(),
        cpu=True,numerical_threads=1,seed=0,mode='accuracy',seq_type='antibody',chunk_size=CHUNK_SIZE,
        max_run_seconds=MAX_RUN_SECONDS,reference_length_limit=REFERENCE_ANNOTATION_MAX_LENGTH,
        pool_sha256=sha(HERE/'POOL.json'),script_sha256=sha(__file__),
        new_gpu_or_remote_cpu_spend_usd=0,generates_or_mutates_sequences=False)
    status=HERE/'ANNOTATION_STATUS.json';save(status,state)
    connection=sqlite3.connect(DB)
    start=time.monotonic()
    try:
        assert table_hash(connection,'SELECT id,sequence_sha256,canonical FROM heavy ORDER BY id')==pool['immutable_heavy_rows_sha256']
        assert table_hash(connection,'SELECT q,heavy FROM query_heavy ORDER BY q,heavy')==pool['immutable_query_heavy_rows_sha256']
        for package,version in calibration['versions'].items():
            assert importlib.metadata.version(package)==version,package+' version changed'
        for field in ('source_files','weights'):
            for relative,expected in calibration[field].items():
                assert sha(CACHE/'runtime/anarcii'/relative)==expected,relative+' changed'
        random.seed(0);np.random.seed(0);torch.manual_seed(0)
        torch.set_num_threads(1);torch.set_num_interop_threads(1)
        model=Anarcii(seq_type='antibody',mode='accuracy',batch_size=32,cpu=True,ncpu=1)
        assert str(model.device)=='cpu' and torch.get_num_threads()==1
        connection.execute('CREATE INDEX IF NOT EXISTS heavy_annotation_state ON heavy(annotation_state,id)')
        connection.execute("UPDATE heavy SET annotation_state='excluded_noncanonical',error='Original ambiguity retained; no sequence repair or invented CDR' WHERE canonical=0 AND annotation_state='pending'")
        connection.execute("UPDATE heavy SET annotation_state='excluded_length_budget',error='Beyond declared finite model input budget' WHERE length(sequence)>? AND annotation_state='pending'",(REFERENCE_ANNOTATION_MAX_LENGTH,))
        connection.commit()
        state['state']='running';save(status,state)
        while True:
            batch=connection.execute("SELECT id,sequence,sequence_sha256 FROM heavy WHERE annotation_state='pending' ORDER BY id LIMIT ?",(CHUNK_SIZE,)).fetchall()
            if not batch:break
            if time.monotonic()-start>MAX_RUN_SECONDS:
                state.update(state='finite_cpu_budget_reached',remaining=connection.execute("SELECT count(*) FROM heavy WHERE annotation_state='pending'").fetchone()[0]);save(status,state);return
            for key,sequence,expected in batch:
                assert hashlib.sha256(sequence.encode()).hexdigest()==expected,'Reference sequence changed'
            inputs={f'r{key:09d}':sequence for key,sequence,_ in batch}
            begin=time.monotonic();raw=model.number(inputs)
            assert set(raw)==set(inputs),'Missing/extra reference output'
            updates=[]
            for key,sequence,_ in batch:
                result=raw[f'r{key:09d}'];kind=result['chain_type']
                text=''.join(aa for _,aa in (result.get('numbering') or []) if aa!='-')
                positions=[tuple(position) for position,aa in (result.get('numbering') or [])]
                assigned=kind=='H' and not result.get('error')
                span_ok=bool(assigned and text and text==sequence[result['query_start']:result['query_end']+1]
                    and len(positions)==len(set(positions)) and all(1<=p[0]<=128 for p in positions))
                mapped={p[0]:aa for p,aa in (result.get('numbering') or []) if not p[1].strip()}
                flanks=bool(span_ok and mapped.get(104,'-')!='-' and mapped.get(118,'-')!='-')
                loop=''.join(aa for (position,insertion),aa in (result.get('numbering') or [])
                    if 105<=position<=117 and aa!='-') if flanks else None
                raw_bytes=json.dumps(result,separators=(',',':'),allow_nan=False).encode()
                error=result.get('error') or (None if span_ok else 'No verified source-matched H annotation')
                updates.append(('done',kind,result.get('score'),error,zlib.compress(raw_bytes),
                    hashlib.sha256(raw_bytes).hexdigest(),loop,int(flanks),int(span_ok),key))
            connection.executemany('UPDATE heavy SET annotation_state=?,chain_type=?,score=?,error=?,raw=?,raw_sha256=?,cdrh3=?,flanks_complete=?,span_verified=? WHERE id=?',updates)
            connection.commit()
            counts=dict(connection.execute('SELECT annotation_state,count(*) FROM heavy GROUP BY annotation_state'))
            state.update(states=counts,last_chunk_seconds=time.monotonic()-begin,
                last_reference_id=batch[-1][0],elapsed_seconds=time.monotonic()-start,
                updated_at=datetime.now(timezone.utc).isoformat())
            save(status,state)
        connection.execute('PRAGMA wal_checkpoint(TRUNCATE)')
        state.update(state='complete',states=dict(connection.execute('SELECT annotation_state,count(*) FROM heavy GROUP BY annotation_state')),
            valid_complete_h_loops=connection.execute("SELECT count(*) FROM heavy WHERE annotation_state='done' AND flanks_complete=1 AND span_verified=1 AND cdrh3 IS NOT NULL").fetchone()[0],
            elapsed_seconds=time.monotonic()-start,finished_at=datetime.now(timezone.utc).isoformat())
        save(status,state)
        scorer=HERE/'score_evidence.py'
        if scorer.exists():
            with (HERE/'scoring.log').open('w') as stream:
                subprocess.run([sys.executable,str(scorer)],stdout=stream,stderr=subprocess.STDOUT,check=True,timeout=1800)
            state['conditional_report_complete']=True;save(status,state)
        print(json.dumps(state))
    except Exception as error:
        state.update(state='failed',error=repr(error),elapsed_seconds=time.monotonic()-start);save(status,state);raise
    finally:
        connection.close()


if __name__=='__main__':main()
