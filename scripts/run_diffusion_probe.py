"""Persistent finite generation driver; CPU context audit after actual release."""
from pathlib import Path
import dataclasses
import subprocess
import json
import hashlib
import drive_modal_campaign as driver
from drive_modal_adaptive import abort_after_driver_error

ROOT=Path(__file__).resolve().parents[1]
STATE=ROOT/'results/diffusion_v18'

def stop_project_app():
    app=driver.modal.App.lookup('adaptyv-c4-native-control')
    command=[str(ROOT/'.venv-modal/bin/modal'),'app','stop',app.app_id,'--yes']
    result=subprocess.run(command,capture_output=True,text=True,timeout=60)
    (STATE/'force_app_stop.log').write_text(result.stdout+result.stderr);result.check_returncode()
    rows=json.loads(subprocess.check_output([str(ROOT/'.venv-modal/bin/modal'),'app','list','--json'],text=True,timeout=30))
    match=[r for r in rows if r.get('app_id')==app.app_id]
    if len(match)!=1 or match[0]['state']!='stopped' or match[0]['tasks']!='0':raise RuntimeError('Actual project-app shutdown unverified')
    return dict(gpu_released_verified=True,app_stopped_verified=True,app_id=app.app_id,app_state=match[0])

if __name__=='__main__':
    STATE.mkdir(exist_ok=True)
    if not (STATE/'CPU_PREFLIGHT.json').exists():raise RuntimeError('Fresh CPU image proof required')
    proof=json.loads((STATE/'CPU_RUNTIME_PROOF.json').read_text())
    if proof.get('status')!='passed' or proof['modal_runner_python']!=[3,11]:raise RuntimeError('Actual Modal CPU bootstrap must pass before GPU allocation')
    if proof['generation_worker_sha256']!=hashlib.sha256((ROOT/'scripts/diffusion_worker.py').read_bytes()).hexdigest():raise RuntimeError('Actual CPU bootstrap does not match current generation worker')
    attempt=json.loads((STATE/'CURRENT_RESOURCE_ATTEMPT.json').read_text())
    run=attempt.get('run_name','c4_diffusion_v18')
    if run not in ('c4_diffusion_v18','c4_diffusion_v18_retry4'):raise ValueError('Unprepared retry identity')
    driver.STATE=STATE;driver.CANCEL_ROOT=STATE;driver.HOURLY_MAX=3.198528
    if not (STATE/'supervisor.json').exists():
        stats={n:dataclasses.asdict(driver.modal.Function.from_name('adaptyv-c4-native-control',n).get_current_stats())
               for n in ('run_campaign','frame_diagnostic','sequence_refold','diffusion_generate')}
        if any(s['num_total_runners'] or s['num_running_inputs'] or s['backlog'] for s in stats.values()):raise RuntimeError('Existing compute prevents allocation')
        driver.atomic(STATE/'pre_spawn_idle.json',stats)
    try:
        result=driver.main(run,{},function_name='diffusion_generate',perform_cpu_followup=False,max_runtime_seconds=2400,max_startup_seconds=300)
        if not result.get('gpu_released_verified'):
            driver.atomic(STATE/'emergency_release.json',stop_project_app())
            raise RuntimeError('Per-call release failed; stopped owned project app')
        if result.get('gpu_released_verified') and result.get('result',{}).get('status')=='completed':
            subprocess.run([str(ROOT/'.venv-rosetta/bin/python'),str(ROOT/'scripts/audit_diffusion_backbone.py')],check=True)
            if json.loads((STATE/'BACKBONE_AUDIT.json').read_text()).get('sequence_design_allowed'):
                subprocess.run([str(ROOT/'.venv-rosetta/bin/python'),str(ROOT/'scripts/prepare_diffusion_sequences.py')],check=True)
    except BaseException:
        try:release=abort_after_driver_error(STATE,function_name='diffusion_generate')
        except Exception as exc:release=dict(gpu_released_verified=False,error=repr(exc))
        if not release.get('gpu_released_verified'):release=stop_project_app()
        driver.atomic(STATE/'emergency_release.json',release)
        raise
