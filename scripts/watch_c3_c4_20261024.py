#!/usr/bin/env python3
"""Record bounded Modal progress and stop when both authorized runs resolve."""
import csv
import json
import pathlib
import subprocess
import time

ROOT = pathlib.Path(".")
MODAL = "modal"
LOG = ROOT / "results/runner_monitor_20261024.jsonl"
APPS = {"c3": "ap-yoYiuB1uFgZY7n807wWB8g", "c4": "ap-x8qUMEomla0Dhcf7x8Zzb8"}


def modal(*args):
    result = subprocess.run([MODAL, *args], capture_output=True, text=True, timeout=45)
    if result.returncode:
        if "No such file or directory" in result.stderr and args[:2] == ("volume", "ls"):
            return "[]"
        raise RuntimeError(result.stderr.strip()[-1000:])
    return result.stdout


def snapshot():
    apps = {app["App ID"]: app for app in json.loads(modal("app", "list", "--json"))}
    state = {"time": time.strftime("%Y-%m-%dT%H:%M:%S%z"), "apps": {}}
    for key, app_id in APPS.items():
        app = apps.get(app_id)
        state["apps"][key] = app or {"App ID": app_id, "State": "not_listed"}
    for key, volume, prefix in (
        ("c3", "c3-af2m-runs", "c3_ecd_complex_20261024"),
        ("c4", "rfd4-runs", "boltz2_c4_apical"),
    ):
        rows = json.loads(modal("volume", "ls", "--json", volume, prefix + "/"))
        suffix = "scores.json" if key == "c3" else "_result.json"
        files = [row["Filename"] for row in rows if row["Type"] == "file"]
        state[key + "_entries"] = len(rows)
        if key == "c3":
            with (ROOT / "results/c3_screen_input_20261024.csv").open() as source:
                names = [row["name"] for row in csv.DictReader(source)]
            score_files = []
            for name in names:
                listing = json.loads(modal("volume", "ls", "--json", volume,
                                           f"{prefix}/{name}/scores.json"))
                score_files.extend(row["Filename"] for row in listing if row["Type"] == "file")
            state[key + "_count"] = len(score_files)
            state["c3_consolidated"] = any(name.endswith("consolidated_scores.json") for name in files)
        else:
            state[key + "_count"] = sum(name.endswith(suffix) for name in files)
            state["c4_final_result"] = any(name.endswith("tfr1_c4_acid_site1_41_result.json") for name in files)
    return state


with LOG.open("a") as log:
    while True:
        try:
            state = snapshot()
            log.write(json.dumps(state) + "\n")
            log.flush()
            print(json.dumps(state), flush=True)
            c3 = state["apps"]["c3"]
            c4 = state["apps"]["c4"]
            active = lambda app: app.get("Stopped at") is None and app.get("State") != "stopped"
            if state.get("c4_final_result") or not active(c4):
                if not active(c3):
                    break
        except Exception as error:
            event = {"time": time.strftime("%Y-%m-%dT%H:%M:%S%z"), "error": str(error)}
            log.write(json.dumps(event) + "\n")
            log.flush()
            print(json.dumps(event), flush=True)
        time.sleep(300)
