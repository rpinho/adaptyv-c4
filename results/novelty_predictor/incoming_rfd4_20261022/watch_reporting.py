"""Finite detached observer: report this existing CPU run; never launch compute."""
from datetime import datetime, timedelta, timezone
import argparse
import json
import os
from pathlib import Path
import signal
import subprocess
import sys
import time

HERE = Path(__file__).resolve().parent


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--workflow-dir',type=Path,default=HERE/'run')
    parser.add_argument('--state-file',type=Path,default=HERE/'REPORTING_STATUS.json')
    args=parser.parse_args()
    signal.signal(signal.SIGHUP, signal.SIG_IGN)
    deadline = datetime.now(timezone.utc) + timedelta(hours=12)
    state = dict(state='watching_existing_cpu_run', pid=os.getpid(),
        started_at=datetime.now(timezone.utc).isoformat(), deadline=deadline.isoformat(),
        launches_gpu_or_remote_compute=False, launches_or_repeats_analysis=False)
    status = args.state_file
    def save():
        state['updated_at'] = datetime.now(timezone.utc).isoformat()
        temp = status.with_suffix('.tmp'); temp.write_text(json.dumps(state, indent=2) + '\n'); temp.replace(status)
    save()
    try:
        while datetime.now(timezone.utc) < deadline:
            try:
                workflow = json.loads((args.workflow_dir / 'STATUS.json').read_text())
            except (FileNotFoundError, json.JSONDecodeError):
                time.sleep(15); continue
            state['workflow_state'] = workflow['state']
            state['workflow_controller_pid'] = workflow['controller_pid']
            state['current_step'] = (workflow.get('steps') or [{}])[-1].get('name')
            save()
            if workflow['state'] == 'complete':
                with (HERE / 'finish_reporting.log').open('w') as log:
                    subprocess.run([sys.executable, str(HERE / 'finish_reporting.py'),'--workflow-dir',str(args.workflow_dir)],
                        stdout=log, stderr=subprocess.STDOUT, check=True, timeout=300)
                state['state'] = 'report_complete'; save(); return
            if workflow['state'] == 'failed':
                raise RuntimeError('Existing workflow failed: ' + str(workflow.get('error')))
            try:
                os.kill(workflow['controller_pid'], 0)
            except ProcessLookupError:
                raise RuntimeError('Controller exited before a completed/failed checkpoint; no blind retry')
            time.sleep(15)
        raise TimeoutError('Finite 12-hour reporting window ended; no new compute was launched')
    except Exception as error:
        state.update(state='reporting_failed', error=repr(error)); save(); raise


if __name__ == '__main__':
    main()
