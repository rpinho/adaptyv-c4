#!/usr/bin/env python3
"""Wait without renting compute, then run one predeclared native-target control.

Only resumes this campaign's existing pod. Never creates accounts/resources,
changes GPU class, exceeds the quoted rate, or launches duplicate workers.
The existing supervisor owns pod shutdown once its first status is verified.
"""
import argparse
import json
import os
from pathlib import Path
import subprocess
import sys
import time

from runpod_control import request

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results/gpu_campaign_v3'
POD='2p9x2qea1ietbw'
RATE_LIMIT=1.79
TOTAL_CAP=3.58
PRIOR=1.7952


class CapacityUnavailable(Exception):pass


def save(data):
    temporary=OUT/'queue.json.tmp';temporary.write_text(json.dumps(data,indent=2)+'\n')
    temporary.replace(OUT/'queue.json')


def call(name,arguments):
    result=request('tools/call',dict(name=name,arguments=arguments))
    text='\n'.join(c['text'] for c in result.get('content',[]) if c.get('type')=='text')
    if result.get('isError'):
        if 'not enough free GPUs' in text:raise CapacityUnavailable('host GPU occupied')
        raise RuntimeError(name+' failed: '+text[:600])
    return json.loads(text)


def stop_verified():
    for _ in range(4):
        try:
            call('pod-action',{'id':POD,'body':{'action':'stop'}})
            if call('get-pod',{'id':POD})['status']=='EXITED':return True
        except Exception:pass
        time.sleep(5)
    return False


def launch(state):
    live=call('get-pod',{'id':POD})
    if live['status']!='EXITED':
        state['reason']='pod already active; leave it alone';return False
    price=call('get-gpu-type',{'id':'NVIDIA A100-SXM4-80GB'})['price']['secure']
    if price>RATE_LIMIT:raise RuntimeError('catalog rate exceeds prequoted limit')
    requested=False;delegated=False
    try:
        requested=True
        try:call('pod-action',{'id':POD,'body':{'action':'start'}})
        except CapacityUnavailable:
            requested=False;raise
        state.update(status='booting',hourly_usd=price);save(state)
        for _ in range(24):
            live=call('get-pod',{'id':POD})
            if live['status']=='RUNNING' and live.get('ssh',{}).get('direct'):
                if live.get('gpu',{}).get('count',0)!=1:raise RuntimeError('GPU allocation missing')
                if live['cost']>RATE_LIMIT:raise RuntimeError('pod rate exceeds prequoted limit')
                break
            time.sleep(5)
        else:raise RuntimeError('direct SSH did not become ready within boot bound')
        direct=live['ssh']['direct'];host=direct['host'];port=direct['port']
        ssh=['ssh','-i',str(Path.home()/'.ssh/runpod_ed25519'),'-p',str(port),
             '-o','BatchMode=yes','-o','ConnectTimeout=10','-o','StrictHostKeyChecking=accept-new','root@'+host]
        for local,remote in [(OUT/'campaign.json','/workspace/runs/c4_native_apical_ablation_v3.json'),
                             (ROOT/'results/gpu_campaign_v2/human_11qc_full_monomer.pdb',
                              '/workspace/runs/human_11qc_full_monomer.pdb')]:
            command=['scp','-i',str(Path.home()/'.ssh/runpod_ed25519'),'-P',str(port),
                     '-o','BatchMode=yes','-o','ConnectTimeout=10','-o','StrictHostKeyChecking=accept-new',
                     str(local),'root@'+host+':'+remote]
            subprocess.run(command,check=True,capture_output=True,text=True,timeout=45)
        payload='''import subprocess,os,json
from pathlib import Path
name=subprocess.check_output(['nvidia-smi','--query-gpu=name','--format=csv,noheader'],text=True).strip()
assert name=='NVIDIA A100-SXM4-80GB',name
env=dict(os.environ,BINDCRAFT_AF2_PARAMS='/workspace/af2_params',BINDCRAFT_WORKERS_PER_GPU='1',BINDCRAFT_DESIGN_WORKERS='1')
log=open('/workspace/logs/c4_native_apical_ablation_v3.log','ab')
p=subprocess.Popen(['/workspace/src/BindCraft2/.venv/bin/bindcraft','design','/workspace/runs/c4_native_apical_ablation_v3.json'],cwd='/workspace/src/BindCraft2',env=env,stdin=subprocess.DEVNULL,stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
print(json.dumps(dict(worker_pid=p.pid)))
'''
        result=subprocess.run(ssh+['python -'],input=payload,capture_output=True,text=True,check=True,timeout=30)
        worker=json.loads(result.stdout)['worker_pid']
        config=dict(state_dir=str(OUT.relative_to(ROOT)),pod=POD,worker_pid=worker,
                    ssh_host=host,ssh_port=port,remote_project='/workspace/runs/c4_native_apical_ablation_v3',
                    remote_log='/workspace/logs/c4_native_apical_ablation_v3.log',
                    prior_estimated_compute_usd=PRIOR,total_compute_cap_usd=TOTAL_CAP)
        config_path=OUT/'supervisor_config.json';config_path.write_text(json.dumps(config,indent=2)+'\n')
        with (OUT/'supervisor.log').open('ab') as log:
            monitor=subprocess.Popen([sys.executable,'-u',str(ROOT/'scripts/supervise_gpu_campaign.py'),
                                      '--config',str(config_path)],cwd=ROOT,stdin=subprocess.DEVNULL,
                                     stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
        for _ in range(30):
            path=OUT/'supervisor.json'
            if path.exists():
                status=json.loads(path.read_text())
                if status.get('pid')==monitor.pid and status.get('worker_pid')==worker:
                    delegated=True;break
            if monitor.poll() is not None:raise RuntimeError('supervisor exited at startup')
            time.sleep(1)
        if not delegated:raise RuntimeError('supervisor did not acknowledge ownership')
        with (OUT/'continuation.log').open('ab') as log:
            queue=subprocess.Popen([sys.executable,'-u',str(ROOT/'scripts/continue_after_gpu.py'),
                                    '--state-dir',str(OUT),'--followup','results/gpu_followup_v3'],
                                   cwd=ROOT,stdin=subprocess.DEVNULL,stdout=log,stderr=subprocess.STDOUT,
                                   start_new_session=True)
        state.update(status='modeling_delegated_to_supervisor',worker_pid=worker,monitor_pid=monitor.pid,
                     continuation_pid=queue.pid,delegated_at=time.time());save(state)
        return True
    finally:
        if requested and not delegated:
            verified=stop_verified();state.update(status='launch_failed_pod_stopped' if verified else
                                                  'STOP_FAILED_REQUIRES_ATTENTION',stop_verified=verified)
            save(state)


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--wait-hours',type=float,default=6)
    args=parser.parse_args();OUT.mkdir(exist_ok=True,parents=True)
    state=dict(status='waiting_for_capacity',pid=os.getpid(),started=time.time(),pod=POD,attempts=0,
               max_hourly_usd=RATE_LIMIT,prior_compute_estimate=PRIOR,total_compute_cap_usd=TOTAL_CAP,
               wait_deadline=time.time()+args.wait_hours*3600)
    while time.time()<state['wait_deadline']:
        if (OUT/'CANCEL').exists():state['status']='cancelled';break
        state.update(status='checking_capacity',attempts=state['attempts']+1,last_attempt=time.time());save(state)
        try:
            if launch(state):return
        except CapacityUnavailable:state['reason']='host GPU occupied; no compute started'
        except Exception as exc:
            state['error']=repr(exc)
            if state.get('status')!='STOP_FAILED_REQUIRES_ATTENTION':state['status']='failed'
            save(state);return
        state.update(status='waiting_for_capacity',next_attempt=time.time()+600);save(state)
        while time.time()<min(state['next_attempt'],state['wait_deadline']):
            if (OUT/'CANCEL').exists():state['status']='cancelled';save(state);return
            time.sleep(5)
    else:state['status']='capacity_wait_expired'
    save(state)


if __name__=='__main__':main()
