#!/usr/bin/env python3
"""Mirror a bounded campaign and stop its POD on exit/failure/budget expiry.

Runs locally, detached from the interactive turn. Uses existing OAuth only for
the named migrated pod. It never terminates a pod, starts another, or changes
account settings. API failures remain explicit and are retried, never called a stop.
"""
import datetime
import argparse
import json
import os
from pathlib import Path
import subprocess
import time
import traceback

from runpod_control import request

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results/gpu_campaign'
POD='2p9x2qea1ietbw'
SSH=['ssh','-i',str(Path.home()/'.ssh/runpod_ed25519'),'-p','13524',
     '-o','BatchMode=yes','-o','ConnectTimeout=10','-o','ServerAliveInterval=10',
     '-o','ServerAliveCountMax=2','root@157.157.221.29']
WORKER=1349
REMOTE_PROJECT='/workspace/runs/c4_full_context_apical_v1'
REMOTE_LOG='/workspace/logs/c4_full_context_apical_v1.log'
PRIOR_SPEND=0.0
TOTAL_CAP=3.58


def save(path,data):
    tmp=path.with_suffix(path.suffix+'.tmp');tmp.write_text(json.dumps(data,indent=2)+'\n');tmp.replace(path)


def decode(result):
    if result.get('isError'):raise RuntimeError(str(result))
    text='\n'.join(c['text'] for c in result.get('content',[]) if c.get('type')=='text')
    return json.loads(text)


def inspect():
    data=decode(request('tools/call',dict(name='get-pod',arguments={'id':POD})))
    return {k:v for k,v in data.items() if k in ['id','name','status','cost','gpu','startedAt','createdAt']}


def mirror():
    import shlex
    dest=OUT/'mirrored';dest.mkdir(exist_ok=True)
    remote=shlex.join(SSH[:-1])
    cp=subprocess.run(['rsync','-az','--timeout=25','-e',remote,
                       '--exclude=compile_cache/','--exclude=*.npz','--exclude=*.pkl',
                       SSH[-1]+':'+REMOTE_PROJECT+'/',str(dest)+'/'],
                      capture_output=True,text=True,timeout=50)
    log=subprocess.run(SSH+['tail -100 '+REMOTE_LOG],
                       capture_output=True,text=True,timeout=25)
    (OUT/'worker_tail.log').write_text(log.stdout+log.stderr)
    return dict(rsync_exit=cp.returncode,stderr=cp.stderr[-1000:],tail_exit=log.returncode)


def main():
    state=dict(status='monitoring',pid=os.getpid(),pod=POD,worker_pid=WORKER,started=time.time(),
               checks=0,ssh_failures=0,mirrors=[],errors=[],rate_source='fresh get-pod cost',
               max_compute_usd=TOTAL_CAP,prior_attempt_estimated_compute_usd=PRIOR_SPEND,stop_verified=False)
    status=OUT/'supervisor.json';save(status,state)
    pod={'cost':1.80}
    billing_start=state['started']-600
    deadline=state['started']
    reason='unknown'
    try:
        pod=inspect();save(OUT/'live_pod.json',pod)
        billing_start=datetime.datetime.fromisoformat(pod['startedAt'].replace('Z','+00:00')).timestamp()
        # Reserve three minutes for failed calls/retries inside the $3.58 bound.
        deadline=billing_start+((TOTAL_CAP-PRIOR_SPEND)/float(pod['cost']))*3600-180
        state.update(billing_start=billing_start,hourly_compute_usd=pod['cost'],stop_deadline=deadline)
        while time.time()<deadline:
            check=subprocess.run(SSH+[
                f'ps -p {WORKER} -o stat=; nvidia-smi --query-gpu=utilization.gpu,memory.used --format=csv,noheader'],
                capture_output=True,text=True,timeout=30)
            state['checks']+=1;state['last_check']=time.time()
            state['estimated_compute_usd']=round((time.time()-billing_start)*float(pod['cost'])/3600,4)
            state['total_estimated_compute_usd']=round(PRIOR_SPEND+state['estimated_compute_usd'],4)
            state['last_gpu_worker_state']=check.stdout.strip()
            if check.returncode:
                state['ssh_failures']+=1
            else:
                state['ssh_failures']=0
                lines=check.stdout.splitlines()
                # ps output is a process state; without it first line is GPU metrics.
                if not lines or not lines[0].strip() or ',' in lines[0] or lines[0].startswith('Z'):
                    reason='worker_exited';break
            if state['ssh_failures']>=3:
                reason='three_consecutive_ssh_failures';break
            try:state['mirrors'].append(dict(time=time.time(),**mirror()))
            except Exception as exc:state['errors'].append(dict(stage='mirror',error=repr(exc)))
            save(status,state)
            time.sleep(min(30,max(0,deadline-time.time())))
        else:reason='compute_cap'
    except BaseException as exc:
        reason='supervisor_exception';state['errors'].append(dict(error=repr(exc),traceback=traceback.format_exc()))
    finally:
        state.update(status='stopping_pod',stop_reason=reason);save(status,state)
        if time.time()<deadline:
            try:state['mirrors'].append(dict(time=time.time(),**mirror()))
            except Exception as exc:state['errors'].append(dict(stage='final_mirror',error=repr(exc)))
        # Retry actual pod stop and verify through independent get-pod read.
        for attempt in range(6):
            try:
                try:current=inspect()
                except Exception as exc:
                    state['errors'].append(dict(stage='pre_stop_read',error=repr(exc)))
                    current={'status':'unknown'}
                if current['status']!='EXITED':
                    result=request('tools/call',dict(name='pod-action',arguments={'id':POD,'body':{'action':'stop'}}))
                    if result.get('isError'):state['errors'].append(dict(stage='stop_reply',result=result))
                current=inspect();save(OUT/'live_pod.json',current)
                if current['status']=='EXITED':
                    state.update(status='pod_stopped',stop_verified=True,stopped_at=time.time())
                    state['estimated_compute_usd']=round((time.time()-billing_start)*float(pod['cost'])/3600,4)
                    state['total_estimated_compute_usd']=round(PRIOR_SPEND+state['estimated_compute_usd'],4)
                    save(status,state);break
            except Exception as exc:state['errors'].append(dict(stage='stop',attempt=attempt,error=repr(exc)))
            save(status,state);time.sleep(10)
        if not state['stop_verified']:
            state['status']='STOP_FAILED_REQUIRES_ATTENTION';save(status,state)


def configure(path):
    global OUT,POD,SSH,WORKER,REMOTE_PROJECT,REMOTE_LOG,PRIOR_SPEND,TOTAL_CAP
    import re
    cfg=json.loads(Path(path).read_text());OUT=(ROOT/cfg['state_dir']).resolve();OUT.mkdir(exist_ok=True,parents=True)
    POD=cfg['pod'];WORKER=int(cfg['worker_pid']);REMOTE_PROJECT=cfg['remote_project'];REMOTE_LOG=cfg['remote_log']
    if not re.fullmatch(r'/workspace/runs/[A-Za-z0-9_.-]+',REMOTE_PROJECT):raise ValueError('unsafe project path')
    if not re.fullmatch(r'/workspace/logs/[A-Za-z0-9_.-]+',REMOTE_LOG):raise ValueError('unsafe log path')
    SSH=['ssh','-i',str(Path.home()/'.ssh/runpod_ed25519'),'-p',str(int(cfg['ssh_port'])),
         '-o','BatchMode=yes','-o','ConnectTimeout=10','-o','ServerAliveInterval=10',
         '-o','ServerAliveCountMax=2','root@'+cfg['ssh_host']]
    PRIOR_SPEND=float(cfg.get('prior_estimated_compute_usd',0));TOTAL_CAP=float(cfg.get('total_compute_cap_usd',3.58))


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--config');args=parser.parse_args()
    if args.config:configure(args.config)
    main()
