"""Bounded invocation of Claude's existing app; contains no generation model.

--execute is used only after Ricardo approves the quoted GPU rate. The finished
novelty report must already be committed. Default mode performs local checks.
"""
import argparse
import csv
from datetime import datetime, timezone
import hashlib
import importlib.util
import json
from pathlib import Path
import subprocess
import re
import sys
import time

ROOT=Path(__file__).resolve().parents[2]
SOURCE=Path('/Users/ricardo/.worktrees/adaptyv-c4/claude-submission/modal_rfd4.py')
EXPECTED_SOURCE='3554f595c0c20a98aff02356d81800828d22619e569acda99777a4f1b348e609'
MODAL=ROOT/'.venv-modal/bin/modal'
ARMS=[('tfr1_c200m_20261023','tfr1_apical_binder'),('tfr1_his200m_20261023','tfr1_apical_his')]
FUNCTION_TIMEOUT=9000
GPU_USD_PER_SECOND=.000542


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def apps_released(snapshot,app_ids):
    if not app_ids:return False
    return all(row.get('state')=='stopped' and str(row.get('tasks'))=='0'
               for row in snapshot if row.get('app_id') in app_ids)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--execute',action='store_true')
    parser.add_argument('--sweep',choices=('initial','higher-his'),default='initial')
    args=parser.parse_args()
    arms=ARMS;expected_source=EXPECTED_SOURCE;output_name='rfd4_geom_20261023'
    if args.sweep=='higher-his':
        prepared=json.loads((ROOT/'results/generation_pivot/HIGHER_PRIOR_PLAN.json').read_text())
        expected_source=prepared['launcher_sha256'];arms=prepared['arms'];output_name='rfd4_his_sweep_20261023'
        assert len(arms)==2 and all(len(arm)==2 and all(re.fullmatch(r'[A-Za-z0-9_]+',s) for s in arm) for arm in arms)
    assert sha(SOURCE)==expected_source,'Canonical launcher changed; recheck it before spending'
    plan=dict(source=str(SOURCE),source_sha256=expected_source,arms=arms,replicates_each=200,
              batch_size=8,seed=42,metrics='geom',function_timeout_seconds=FUNCTION_TIMEOUT,
              max_gpu_containers=2,gpu_usd_per_hour=GPU_USD_PER_SECOND*3600,
              cpu_ram_charges_additional=True,contains_generator=False)
    if not args.execute:
        print(json.dumps(plan,indent=2));return
    if args.sweep=='higher-his':
        approval=json.loads((ROOT/'results/generation_pivot/HIGHER_PRIOR_GPU_APPROVAL.json').read_text())
        assert approval['approved_by']=='Ricardo' and approval['approved_arms']==[a[0] for a in arms]
        assert approval['model_invoker_sha256']==sha(__file__)
        checked_inputs=ROOT/'results/generation_pivot/HIGHER_PRIOR_VALIDATION.json'
        assert sha(checked_inputs)==prepared['cpu_validation_sha256']
        assert json.loads(checked_inputs.read_text())['state']=='passed_cpu_materialization'
        first=ROOT/'results/rfd4_geom_20261023'
        completed=json.loads((first/'STATUS.json').read_text())
        assert completed['state']=='outputs_ready_for_local_validation' and completed['app_stop_confirmed']
        snapshot=json.loads((first/'apps_after_stop.json').read_text())
        assert apps_released(snapshot,completed['apps_started']),'Initial GPU allocations must be released first'
        with (first/'ranking/shortlist.csv').open(newline='') as stream:
            assert list(csv.DictReader(stream)),'Rank the first measured arms before this sweep'
    report=ROOT/'results/novelty_predictor/adaptive_completion/predictions.json'
    validation=report.parent/'VALIDATION.json'
    checked=json.loads(validation.read_text())
    assert checked['state']=='complete' and checked['predictions_sha256']==sha(report)
    committed=subprocess.check_output(['git','show','HEAD:'+str(report.relative_to(ROOT))],cwd=ROOT)
    assert hashlib.sha256(committed).hexdigest()==sha(report),'Commit resolved novelty levels before generation'
    output=ROOT/'results'/output_name;output.mkdir()
    state=dict(state='initializing',started_at=datetime.now(timezone.utc).isoformat(),plan=plan,
               generator_source_unchanged=True,apps_started=[],calls=[],new_runpod_spend_usd=0)
    def save():(output/'STATUS.json').write_text(json.dumps(state,indent=2)+'\n')
    save();app_id=None;started=time.monotonic();error=None
    try:
        spec=importlib.util.spec_from_file_location('modal_rfd4',SOURCE)
        module=importlib.util.module_from_spec(spec)
        sys.path.insert(0,str(SOURCE.parent));sys.modules['modal_rfd4']=module
        spec.loader.exec_module(module)
        # Reuse the canonical App and Function, with bounded invocation resources.
        with module.modal.enable_output(),module.app.run(detach=False):
            app_id=module.app.app_id;state['apps_started']=[app_id];save()
            function=module.generate.with_options(timeout=FUNCTION_TIMEOUT,max_containers=2,
                                                  scaledown_window=2,retries=0)
            calls=[]
            for name,input_spec in arms:
                call=function.spawn(num_replicates=200,batch_size=8,seed=42,
                                    out_name=name,metrics='geom',spec=input_spec)
                calls.append(call);state['calls'].append(dict(out_name=name,spec=input_spec,call_id=call.object_id));save()
            state['state']='running';save()
            for call,record in zip(calls,state['calls']):
                call.get(timeout=max(1,FUNCTION_TIMEOUT+600-(time.monotonic()-started)))
                record['completed_at']=datetime.now(timezone.utc).isoformat();save()
            state['state']='generation_finished_exporting_sequences';save()
            for name,_ in arms:
                module.export_sequences.with_options(timeout=300,max_containers=1,scaledown_window=2,retries=0).remote(run=name)
            state['state']='outputs_ready_for_local_validation';save()
    except BaseException as caught:
        error=caught;state.update(state='failed',error=repr(caught));save()
    finally:
        if app_id:
            stopped=False
            for attempt in range(3):
                try:
                    result=subprocess.run([str(MODAL),'app','stop',app_id],capture_output=True,text=True,timeout=60)
                    (output/('stop_'+str(attempt)+'.log')).write_text(result.stdout+result.stderr)
                    if result.returncode==0:stopped=True;break
                except Exception as stop_error:
                    (output/('stop_'+str(attempt)+'.log')).write_text(repr(stop_error))
                time.sleep(5)
            state['app_stop_confirmed']=stopped
            try:
                snapshot=subprocess.run([str(MODAL),'app','list','--json'],capture_output=True,text=True,timeout=60)
                (output/'apps_after_stop.json').write_text(snapshot.stdout)
                if snapshot.returncode==0:
                    found=[row for row in json.loads(snapshot.stdout) if row.get('app_id')==app_id]
                    state['app_absent_from_active_listing']=not found
                    state['app_allocation_release_verified']=stopped and apps_released(json.loads(snapshot.stdout),[app_id])
            except Exception as snapshot_error:
                state['release_snapshot_error']=repr(snapshot_error)
            if not stopped and error is None:error=RuntimeError('App stop was not confirmed')
        state.update(finished_at=datetime.now(timezone.utc).isoformat(),elapsed_wall_seconds=time.monotonic()-started,
                     account_wide_bill_verified=False,recorded_prior_runpod_usd=56.5899880444631,
                     recorded_prior_own_modal_estimate_usd=10.718947158590089)
        state['new_modal_gpu_time_cost_upper_estimate_usd']=2*state['elapsed_wall_seconds']*GPU_USD_PER_SECOND if app_id else 0
        state['estimate_scope']='GPU only, using two maximum containers over whole controller wall time; CPU/RAM/build/storage excluded'
        save()
    if error:raise error


if __name__=='__main__':main()
