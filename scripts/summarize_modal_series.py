"""Summarize finished native-target controls without treating rejected folds as binders."""
import csv
import datetime
import json
from pathlib import Path
import re

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "results/modal_campaign_v3"


def main():
    batches = []
    for folder in [BASE, BASE / "replicate_2", BASE / "replicate_3"]:
        supervisor = json.loads((folder / "supervisor.json").read_text())
        rows = list(csv.DictReader((folder / "mirrored/campaign/1_Trajectories/!_Trajectories.csv").open()))
        log = (folder / "mirrored/worker.log").read_text()
        accepted = re.search(r"campaign done: (\d+) accepted design", log)
        readings = re.findall(r"(passed|rejected) at (\w+) design stage\s+i_pTM=([\d.]+)\s+pLDDT=([\d.]+)", log)
        batches.append({
            "folder": str(folder.relative_to(ROOT)), "trajectories": len(rows),
            "accepted": int(accepted[1]) if accepted else None,
            "stopped_at": [r["terminated"] for r in rows],
            "autotuned": [r["autotuned"] for r in rows],
            "stage_readings": [dict(outcome=a, stage=b, iptm=float(c), plddt=float(d)) for a,b,c,d in readings],
            "gpu_released_verified": supervisor["gpu_released_verified"],
            "runtime_estimate_upper_usd": supervisor["compute_estimate_upper_usd"],
            "finished": supervisor["remote"]["finished"],
        })
    report = dict(batches=batches, trajectories=sum(b["trajectories"] for b in batches),
                  accepted=sum(b["accepted"] for b in batches),
                  modal_runtime_estimate_upper_usd=sum(b["runtime_estimate_upper_usd"] for b in batches),
                  runpod_supervised_runtime_estimate_usd=1.7952,
                  billing_excludes="CPU image build, startup/teardown, storage, older campaigns",
                  gpu_released_verified=all(b["gpu_released_verified"] for b in batches))
    (BASE / "summary.json").write_text(json.dumps(report, indent=2) + "\n")
    finished = datetime.datetime.fromtimestamp(max(b["finished"] for b in batches), datetime.timezone.utc)
    lines = ["# Native-target control: completed, no accepted designs", "",
             f"Three independent seeds produced {report['accepted']} accepted designs from {report['trajectories']} trajectories. Five stopped at screen; one passed screen (pLDDT0.63) and stopped at refine (pLDDT0.55). All autotuned fields were empty. No ProteinMPNN redesign or pKa-scoring candidate was produced.", "",
             f"The last worker finished {finished.isoformat()}. Each batch independently verified zero GPU runners before its CPU follow-up. A fresh Modal stats read at2026-10-21 23:58 UTC also reported zero runners, inputs and backlog. Nothing was modeling at that checkpoint.", "",
             f"Estimated Modal runtime upper total **${report['modal_runtime_estimate_upper_usd']:.4f}**; prior newly supervised Runpod runtime **$1.7952**. Setup, startup/teardown, storage and older campaigns are separate; these are runtime estimates, not reconciled bills.", "",
             "| Batch | Trajectories | Accepted | Stop stages | Runtime upper estimate |", "|---|---:|---:|---|---:|"]
    for index,b in enumerate(batches,1):
        lines.append(f"| {index} | {b['trajectories']} | {b['accepted']} | {', '.join(b['stopped_at'])} | ${b['runtime_estimate_upper_usd']:.4f} |")
    lines += ["", "The control does not support a claim of improved scaffold quality. Six attempts are too few to establish that the epitope is impossible, but they do not justify repeatedly renting the same unchanged configuration. Saved trajectory structures make diagnosis possible without another GPU. Full native target, opposite-subunit restoration, clash checks and final quality gates remain required.", "",
              "The next bounded hypothesis tests upstream's large-complex initialization on the same full native monomer, patch, lengths and seeds. It injects the target coordinates at the start of the structure module; the old run used zero initialization. This is a computational setup hypothesis, not a demonstrated improvement. See ../modal_campaign_v4/DECISION.md."]
    (BASE / "RESULTS.md").write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
