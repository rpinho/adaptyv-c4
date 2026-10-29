"""Three bounded seed replicates with actual GPU release and CPU checks between them."""
import fcntl
import argparse
import json
import os
from pathlib import Path
import time

import drive_modal_campaign as driver

ROOT = Path(__file__).resolve().parents[1]
STATE = ROOT / "results/modal_campaign_v3"


def main():
    STATE.mkdir(parents=True, exist_ok=True)
    driver.CANCEL_ROOT = STATE
    with (STATE / "series.lock").open("w") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        path = STATE / "series.json"
        state = json.loads(path.read_text()) if path.exists() else {
            "started": time.time(), "finished_batches": [], "max_batches": 3,
            "runtime_upper_hourly_usd": driver.HOURLY_MAX,
            "runtime_envelope_usd": 6 * driver.HOURLY_MAX,
        }
        state.update(pid=os.getpid(), status="active")
        driver.atomic(path, state)
        for index in range(3):
            if index in state["finished_batches"]:
                continue
            if (STATE / "CANCEL").exists():
                state["status"] = "cancelled"
                break
            run_name = driver.RUN if index == 0 else driver.RUN + f"_replicate_{index + 1}"
            driver.STATE = STATE if index == 0 else STATE / f"replicate_{index + 1}"
            state.update(current_batch=index, current_run=run_name)
            driver.atomic(path, state)
            result = driver.main(run_name, {"campaign_seed": 20261021 + index})
            state.setdefault("batch_runtime_estimates_upper_usd", {})[str(index)] = result.get("compute_estimate_upper_usd", 0)
            state["runtime_estimate_upper_usd"] = sum(state["batch_runtime_estimates_upper_usd"].values())
            if not result.get("gpu_released_verified") or result.get("status") != "returned":
                state.update(status="stopped_on_execution_failure", last_result=result)
                break
            if result.get("result", {}).get("status") != "completed" or result.get("cpu_followup_returncode") != 0:
                state.update(status="stopped_on_model_or_followup_failure", last_result=result)
                break
            state["finished_batches"].append(index)
            driver.atomic(path, state)
        else:
            state["status"] = "series_completed"
        state["updated"] = time.time()
        driver.atomic(path, state)
        print(json.dumps(state), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--state-dir", default="results/modal_campaign_v3")
    parser.add_argument("--run-name", default=driver.RUN)
    args = parser.parse_args()
    STATE = (ROOT / args.state_dir).resolve()
    driver.RUN = args.run_name
    main()
