"""Five declared native-input seed replicates; finite, no heuristic switching."""
import fcntl
import json
import os
from pathlib import Path
import time

import drive_modal_campaign as driver
from drive_modal_adaptive import abort_after_driver_error
from campaign_admission import read_outcome,followup_path
from modal_search_policy import REPAIRED

ROOT=Path(__file__).resolve().parents[1]
STATE=ROOT/'results/modal_native_replicates_v14'
PRIOR_MODAL_USD=6.384988167371454
LIMIT_SECONDS=2100
RESERVE_SECONDS=2700
SEEDS=tuple(range(20261022,20261027))


def require_idle():
    snapshot={}
    for name in ('run_campaign','frame_diagnostic'):
        stats=driver.modal.Function.from_name('adaptyv-c4-native-control',name).get_current_stats()
        snapshot[name]=dict(runners=stats.num_total_runners,inputs=stats.num_running_inputs,backlog=stats.backlog)
        if stats.num_total_runners or stats.num_running_inputs or stats.backlog:
            raise RuntimeError('Another input/container exists; no new spawn')
    return dict(checked=time.time(),functions=snapshot)


def main():
    STATE.mkdir(parents=True,exist_ok=True)
    driver.CANCEL_ROOT=STATE
    with (STATE/'series.lock').open('w') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        path=STATE/'series.json'
        state=json.loads(path.read_text()) if path.exists() else dict(
            started=time.time(),finished=[],batches={},seeds=list(SEEDS),
            max_batches=len(SEEDS),runtime_planning_budget_usd=len(SEEDS)*RESERVE_SECONDS/3600*driver.HOURLY_MAX,
            prior_modal_runtime_estimate_usd=PRIOR_MODAL_USD,
            accounting_scope='Completed batches; add active supervisor submission elapsed at HOURLY_MAX for live estimate')
        if state.get('status') in {'completed','cancelled','stopped_for_review','stopped_on_failure'}:
            return state
        if state['seeds']!=list(SEEDS) or state['max_batches']!=5:raise ValueError('Queue declaration changed')
        state.update(status='active',pid=os.getpid());driver.atomic(path,state)
        for index,seed in enumerate(SEEDS,2):
            key=str(index)
            if index in state['finished']:continue
            if (STATE/'CANCEL').exists():state['status']='cancelled';break
            spent=sum(b['runtime_upper_usd'] for b in state['batches'].values())
            reserve=RESERVE_SECONDS/3600*driver.HOURLY_MAX
            if spent+reserve>state['runtime_planning_budget_usd']+1e-6:
                state['status']='stopped_for_review';state['reason']='runtime_budget_reserve';break
            run=REPAIRED+f'_replicate_{index}'
            driver.STATE=STATE/f'replicate_{index}'
            state.update(current_run=run,current_seed=seed,updated=time.time());driver.atomic(path,state)
            prior=driver.STATE/'supervisor.json'
            result=json.loads(prior.read_text()) if prior.exists() else {}
            try:
                if not (result.get('gpu_released_verified') and result.get('cpu_followup_returncode')==0):
                    if not result.get('call_id'):
                        state['pre_spawn_idle_snapshot']=require_idle();driver.atomic(path,state)
                    result=driver.main(run,{'campaign_seed':seed},max_runtime_seconds=LIMIT_SECONDS)
                if not result.get('gpu_released_verified') or result.get('status')!='returned' or result.get('result',{}).get('status')!='completed':
                    raise RuntimeError('Model execution or actual GPU release failed')
                if result.get('cpu_followup_returncode')!=0:raise RuntimeError('CPU follow-up failed')
                outcome=read_outcome(driver.STATE/'mirrored/campaign',followup_path(ROOT,run))
                elapsed=result['released']-result['submitted']
                state['batches'][key]=dict(seed=seed,run=run,runtime_upper_usd=elapsed/3600*driver.HOURLY_MAX,outcome=outcome)
                state['runtime_upper_usd']=sum(b['runtime_upper_usd'] for b in state['batches'].values())
                state['cumulative_modal_runtime_upper_usd']=PRIOR_MODAL_USD+state['runtime_upper_usd']
                if not outcome['followup_completed_without_errors']:raise RuntimeError('Scientific evidence requires review')
                state['finished'].append(index)
                if outcome['native_context_admitted_trajectory_ids']:
                    state['status']='stopped_for_review';state['reason']='native_context_outputs_ready';break
                driver.atomic(path,state)
            except Exception as exc:
                state.update(status='stopped_on_failure',error=repr(exc))
                # Driver normally handles shutdown; also handle exceptions outside its loop.
                current=json.loads(prior.read_text()) if prior.exists() else {}
                if current.get('call_id') and not current.get('gpu_released_verified'):
                    try:
                        state['emergency_shutdown']=abort_after_driver_error(driver.STATE)
                        if state['emergency_shutdown'].get('gpu_released_verified'):
                            current.update(gpu_released_verified=True,released=time.time())
                            driver.atomic(prior,current)
                    except Exception as stop_error:state['shutdown_error']=repr(stop_error)
                if current.get('submitted'):
                    state['batches'][key]=dict(seed=seed,run=run,failed=True,
                        runtime_upper_usd=((current.get('released') or time.time())-current['submitted'])/3600*driver.HOURLY_MAX,
                        gpu_released_verified=current.get('gpu_released_verified',False))
                    state['runtime_upper_usd']=sum(b['runtime_upper_usd'] for b in state['batches'].values())
                    state['cumulative_modal_runtime_upper_usd']=PRIOR_MODAL_USD+state['runtime_upper_usd']
                break
        else:state['status']='completed'
        state['updated']=time.time();driver.atomic(path,state)
        return state


if __name__=='__main__':print(json.dumps(main()),flush=True)
