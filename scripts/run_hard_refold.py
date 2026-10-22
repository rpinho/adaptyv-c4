"""One supervised allocation; no CPU pH follow-up or automatic redesign."""
from pathlib import Path
import dataclasses
import json
import subprocess
import drive_modal_campaign as driver
from drive_modal_adaptive import abort_after_driver_error

ROOT=Path(__file__).resolve().parents[1]
STATE=ROOT/'results/hard_sequence_refold_v15'

if __name__=='__main__':
    driver.STATE=STATE;driver.CANCEL_ROOT=STATE
    if not (STATE/'supervisor.json').exists():
        snapshot={name:dataclasses.asdict(driver.modal.Function.from_name('adaptyv-c4-native-control',name).get_current_stats())
                  for name in ('run_campaign','frame_diagnostic','sequence_refold')}
        if any(s['num_total_runners'] or s['num_running_inputs'] or s['backlog'] for s in snapshot.values()):raise RuntimeError('Compute not idle')
        driver.atomic(STATE/'pre_spawn_idle.json',snapshot)
    try:
        result=driver.main('c4_hard_sequence_refold_v15',{},function_name='sequence_refold',perform_cpu_followup=False,max_runtime_seconds=2400)
        if result.get('gpu_released_verified') and result.get('result',{}).get('status')=='completed':
            subprocess.run([str(ROOT/'.venv-rosetta/bin/python'),str(ROOT/'scripts/audit_hard_refolds.py')],check=True)
    except BaseException:
        release=abort_after_driver_error(STATE,function_name='sequence_refold')
        driver.atomic(STATE/'emergency_release.json',release)
        raise
