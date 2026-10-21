#!/usr/bin/env python3
"""Detached local queue: after verified pod stop, run the CPU follow-up."""
import json
import argparse
import os
from pathlib import Path
import subprocess
import time

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results/gpu_campaign'


def main(out=OUT,followup='results/gpu_followup'):
    state=dict(status='waiting_for_verified_pod_stop',pid=os.getpid(),started=time.time())
    path=out/'continuation.json'
    while time.time()-state['started']<8000:
        path.write_text(json.dumps(state,indent=2)+'\n')
        status=json.loads((out/'supervisor.json').read_text())
        if status.get('stop_verified'):
            state['status']='cpu_followup_running';path.write_text(json.dumps(state,indent=2)+'\n')
            with (out/'cpu_followup.log').open('ab') as log:
                process=subprocess.run([str(ROOT/'.venv-rosetta/bin/python'),'-u',
                                        str(ROOT/'scripts/postprocess_gpu_campaign.py'),
                                        '--source',str(out/'mirrored'),'--out',followup],cwd=ROOT,
                                       stdin=subprocess.DEVNULL,stdout=log,stderr=subprocess.STDOUT,timeout=7350)
            state.update(status='cpu_followup_finished' if process.returncode==0 else 'cpu_followup_failed',
                         returncode=process.returncode,finished=time.time())
            path.write_text(json.dumps(state,indent=2)+'\n');return
        if status['status']=='STOP_FAILED_REQUIRES_ATTENTION':
            state['status']='blocked_on_verified_pod_stop';break
        time.sleep(30)
    else:state['status']='supervisor_wait_timeout'
    path.write_text(json.dumps(state,indent=2)+'\n')


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--state-dir',default=str(OUT))
    parser.add_argument('--followup',default='results/gpu_followup');args=parser.parse_args()
    main(Path(args.state_dir).resolve(),args.followup)
