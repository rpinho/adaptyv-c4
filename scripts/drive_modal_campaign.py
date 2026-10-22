"""Resume-safe local driver: one remote call, mirrored evidence, scale-to-zero check."""
import dataclasses
import fcntl
import json
import os
from pathlib import Path
import subprocess
import time

import modal
from campaign_admission import followup_path

ROOT = Path(__file__).resolve().parents[1]
STATE = ROOT / "results/modal_campaign_v3"
CANCEL_ROOT = STATE
RUN = "c4_native_apical_ablation_v3"
HOURLY_MAX = 3.9492 + 4 * 0.04716 + 64 * 0.007992


def atomic(path, value):
    temp = path.with_suffix(path.suffix + ".tmp")
    temp.write_text(json.dumps(value, indent=2) + "\n")
    temp.replace(path)


def mirror(vol, run_name=RUN):
    """Download small durable scientific outputs; heavy trajectories stay in the volume."""
    for item in vol.listdir(run_name, recursive=True):
        if item.type != modal.types.FileEntryType.FILE or item.size > 25 * 1024**2:
            continue
        relative = Path(item.path).relative_to(run_name)
        if relative.suffix not in {".json", ".jsonl", ".log", ".txt", ".csv", ".cif", ".pdb", ".trb"}:
            continue
        target = STATE / "mirrored" / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        if target.exists() and target.stat().st_size == item.size and relative.name not in {"status.json", "worker.log"}:
            continue
        temp = target.with_suffix(target.suffix + ".download")
        with temp.open("wb") as out:
            for chunk in vol.read_file(item.path):
                out.write(chunk)
        temp.replace(target)


def main(run_name=RUN, overrides=None, function_name='run_campaign', perform_cpu_followup=True,max_runtime_seconds=7800):
    STATE.mkdir(parents=True, exist_ok=True)
    with (STATE / "driver.lock").open("w") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        status_path = STATE / "supervisor.json"
        status = json.loads(status_path.read_text()) if status_path.exists() else {}
        function = modal.Function.from_name("adaptyv-c4-native-control", function_name)
        vol = modal.Volume.from_name("adaptyv-c4-results")
        if not status.get("call_id"):
            if status.get("spawn_intent"):
                raise RuntimeError("Prior spawn has uncertain outcome; inspect Modal before retrying")
            status.update(spawn_intent=time.time(), pid=os.getpid(), status="spawning")
            atomic(status_path, status)
            call = function.spawn(run_name, overrides or {})
            status.update(call_id=call.object_id, submitted=time.time(), status="submitted")
            atomic(status_path, status)
        else:
            call = modal.FunctionCall.from_id(status["call_id"])
        status["pid"] = os.getpid()
        print(json.dumps(status), flush=True)
        deadline = status["submitted"] + max_runtime_seconds
        completed = False
        while not completed:
            if (STATE / "CANCEL").exists() or (CANCEL_ROOT / "CANCEL").exists() or time.time() > deadline:
                call.cancel(terminate_containers=True)
                status.update(status="cancelled", cancelled=time.time())
                completed = True
            else:
                try:
                    status["result"] = call.get(timeout=0)
                    status.update(status="returned", returned=time.time())
                    completed = True
                except TimeoutError:
                    pass
                except Exception as exc:
                    status.update(status="remote_failed", error=repr(exc))
                    # A failed input should not leave a paid container behind.
                    call.cancel(terminate_containers=True)
                    completed = True
            try:
                mirror(vol, run_name)
                remote_path = STATE / "mirrored/status.json"
                if remote_path.exists():
                    remote = json.loads(remote_path.read_text())
                    status["remote"] = remote
                    status["compute_estimate_upper_usd"] = remote.get("elapsed_seconds", 0) / 3600 * HOURLY_MAX
                    status["estimate_excludes"] = "image build, container startup/teardown, storage, other apps"
            except Exception as exc:
                status["last_mirror_error"] = repr(exc)
            try:
                status["function_stats"] = dataclasses.asdict(function.get_current_stats())
            except Exception as exc:
                status["last_stats_error"] = repr(exc)
            status["checked"] = time.time()
            atomic(status_path, status)
            if not completed:
                time.sleep(30)

        # Actual runner termination, not just worker exit, is the completion gate.
        for _ in range(30):
            stats = function.get_current_stats()
            status["function_stats"] = dataclasses.asdict(stats)
            if stats.num_total_runners == 0 and stats.num_running_inputs == 0 and stats.backlog == 0:
                status.update(gpu_released_verified=True, released=time.time())
                break
            time.sleep(5)
        else:
            status.update(gpu_released_verified=False, status="REQUIRES_CONTAINER_SHUTDOWN")
        atomic(status_path, status)
        mirror(vol, run_name)
        if status.get("gpu_released_verified") and perform_cpu_followup:
            followup = followup_path(ROOT, run_name)
            followup.mkdir(exist_ok=True)
            with (followup / "postprocess.log").open("a") as log:
                result = subprocess.run([
                    str(ROOT / ".venv-rosetta/bin/python"),
                    str(ROOT / "scripts/postprocess_gpu_campaign.py"),
                    "--source", str(STATE / "mirrored/campaign"), "--out", str(followup)
                ], cwd=ROOT, stdout=log, stderr=subprocess.STDOUT)
            status["cpu_followup_returncode"] = result.returncode
            status["cpu_followup_directory"] = str(followup)
            atomic(status_path, status)
        print(json.dumps(status), flush=True)
        return status


if __name__ == "__main__":
    main()
