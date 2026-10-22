"""Persistent finite generation driver; CPU context audit after actual release."""
from pathlib import Path
import dataclasses
import subprocess
import drive_modal_campaign as driver
from drive_modal_adaptive import abort_after_driver_error

ROOT=Path(__file__).resolve().parents[1]
STATE=ROOT/'results/diffusion_v18'

if __name__=='__main__':
    STATE.mkdir(exist_ok=True)
    if not (STATE/'CPU_PREFLIGHT.json').exists():raise RuntimeError('Fresh CPU image proof required')
    driver.STATE=STATE;driver.CANCEL_ROOT=STATE;driver.HOURLY_MAX=3.198528
    if not (STATE/'supervisor.json').exists():
        stats={n:dataclasses.asdict(driver.modal.Function.from_name('adaptyv-c4-native-control',n).get_current_stats())
               for n in ('run_campaign','frame_diagnostic','sequence_refold','diffusion_generate')}
        if any(s['num_total_runners'] or s['num_running_inputs'] or s['backlog'] for s in stats.values()):raise RuntimeError('Existing compute prevents allocation')
        driver.atomic(STATE/'pre_spawn_idle.json',stats)
    try:
        result=driver.main('c4_diffusion_v18',{},function_name='diffusion_generate',perform_cpu_followup=False,max_runtime_seconds=2400)
        if result.get('gpu_released_verified') and result.get('result',{}).get('status')=='completed':
            subprocess.run([str(ROOT/'.venv-rosetta/bin/python'),str(ROOT/'scripts/audit_diffusion_backbone.py')],check=True)
    except BaseException:
        driver.atomic(STATE/'emergency_release.json',abort_after_driver_error(STATE,function_name='diffusion_generate'))
        raise
