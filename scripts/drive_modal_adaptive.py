"""Run declared experiments, select the next from actual failures, release every GPU."""
import csv
import fcntl
import json
from pathlib import Path
import time

import drive_modal_campaign as driver
import modal_search_policy as policy

ROOT=Path(__file__).resolve().parents[1]
STATE=ROOT/'results/modal_campaign_v5'
BUDGET=27.90
MAX_BATCHES=6
PRIOR_MODAL_RUNTIME=2.87600756690237+0.2559073828754234


def outcome(folder):
    campaign=folder/'mirrored/campaign'
    table=campaign/'3_Ranked/!_Ranked.csv'
    accepted=len(list(csv.DictReader(table.open()))) if table.exists() else 0
    trajectories=campaign/'1_Trajectories/!_Trajectories.csv'
    stages=[r['terminated'] for r in csv.DictReader(trajectories.open())] if trajectories.exists() else []
    return accepted,stages


def abort_after_driver_error(folder):
    record=folder/'supervisor.json'
    snapshot=json.loads(record.read_text()) if record.exists() else {}
    if not snapshot.get('call_id'):
        return dict(gpu_released_verified=False,error='No confirmed call ID; inspect app before any next allocation')
    driver.modal.FunctionCall.from_id(snapshot['call_id']).cancel(terminate_containers=True)
    function=driver.modal.Function.from_name('adaptyv-c4-native-control','run_campaign')
    for _ in range(30):
        stats=function.get_current_stats()
        if stats.num_total_runners==stats.num_running_inputs==stats.backlog==0:
            return dict(gpu_released_verified=True,call_id=snapshot['call_id'])
        time.sleep(5)
    return dict(gpu_released_verified=False,call_id=snapshot['call_id'])


def main():
    STATE.mkdir(parents=True,exist_ok=True)
    driver.CANCEL_ROOT=STATE
    with (STATE/'search.lock').open('w') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        path=STATE/'search.json'
        state=json.loads(path.read_text()) if path.exists() else dict(
            started=time.time(),profile=policy.PROFILES[0],completed=[],accepted=0,
            runtime_estimate_upper_usd=PRIOR_MODAL_RUNTIME,budget_usd=BUDGET,max_batches=MAX_BATCHES)
        if state.get('status') and state['status']!='active':
            print(json.dumps(state),flush=True)
            return
        state['status']='active'
        driver.atomic(path,state)
        for index in range(len(state['completed']),MAX_BATCHES):
            if (STATE/'CANCEL').exists():
                state.update(status='cancelled');break
            # Reserve the complete next two-hour envelope before allocation.
            if state['runtime_estimate_upper_usd']+7800/3600*driver.HOURLY_MAX>BUDGET:
                state.update(status='runtime_budget_reserve_reached');break
            profile=state['profile']
            run_name=profile if index==0 else profile+f'_replicate_{index+1}'
            folder=STATE/f'batch_{index+1}'
            driver.STATE=folder
            state.update(current_batch=index+1,current_run=run_name)
            driver.atomic(path,state)
            try:
                result=driver.main(run_name,{'campaign_seed':20261021+index})
            except Exception as exc:
                state.update(status='stopped_on_driver_exception',error=repr(exc))
                try:
                    state['emergency_release']=abort_after_driver_error(folder)
                except Exception as shutdown_exc:
                    state['emergency_release']=dict(gpu_released_verified=False,error=repr(shutdown_exc))
                driver.atomic(path,state)
                break
            # Submission-to-release wall time conservatively includes queue/startup.
            elapsed=result.get('released',time.time())-result['submitted']
            spend=elapsed/3600*driver.HOURLY_MAX
            state['runtime_estimate_upper_usd']+=spend
            if not result.get('gpu_released_verified') or result.get('status')!='returned' or result.get('result',{}).get('status')!='completed':
                state.update(status='stopped_on_execution_or_release_failure',last_result=result);break
            if result.get('cpu_followup_returncode')!=0:
                state.update(status='stopped_on_cpu_followup_failure',last_result=result);break
            accepted,stages=outcome(folder)
            state['accepted']+=accepted
            next_choice,reason=policy.next_profile(profile,accepted,stages)
            state['completed'].append(dict(run_name=run_name,folder=str(folder),accepted=accepted,
                                           terminated_stages=stages,runtime_estimate_upper_usd=spend,next_reason=reason))
            state.update(profile=next_choice,updated=time.time())
            driver.atomic(path,state)
            if state['accepted']>=3:
                state.update(status='scaffold_set_ready_for_review');break
            if next_choice is None:
                state.update(status=reason);break
        else:
            state['status']='finite_search_limit'
        state['updated']=time.time()
        driver.atomic(path,state)
        print(json.dumps(state),flush=True)


if __name__=='__main__':
    main()
