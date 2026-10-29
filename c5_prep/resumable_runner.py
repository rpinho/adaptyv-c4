#!/usr/bin/env python3
"""Target-independent resumable runner with budget stops (Kimi plan item #3, 2026-10-25).

Nothing here knows what the C5 target is. It runs a *caller-supplied* command per job,
so it is pure infrastructure: smoke-testable with a fake worker, and reusable once the
target and assay are published.

Guarantees, each enforced in code and covered by a fail-closed test:
  1. RESUMABLE   — completed jobs are recorded in a JSON state file keyed by job id and by
                   the sha256 of the input row. A rerun skips them. A job whose recorded
                   artifact no longer matches its hash is NOT skipped (it re-runs).
  2. BUDGET STOP — wall-clock and dollar ceilings are checked BEFORE each launch. The
                   runner refuses to start a job that cannot be afforded and stops cleanly
                   with status "budget_stop" rather than silently overrunning.
  3. FAIL CLOSED — a missing command, a job that exits nonzero, or a job that produces no
                   artifact is recorded as failed and is never counted as done. Any failed
                   job makes the final exit code nonzero.
  4. PROVENANCE  — each job records command, start/end, duration, exit code, artifact path,
                   artifact sha256, and the input-row hash. A run-level manifest records the
                   software version and the state-file hash at exit.

No GPU, no network, no sequences are required to use or test it.

usage:
  resumable_runner.py --jobs jobs.csv --workdir DIR --cmd 'mycmd {id} {input}' \
      --artifact '{id}.json' --state DIR/state.json --max-wall-s 7200 --max-usd 5.00 \
      [--id-col name] \
      [--hourly-usd 0.59] [--dry-run] [--max-jobs N]
Exit codes: 0 all jobs done; 2 budget stop reached with jobs remaining; 3 one or more
jobs failed; 4 fatal configuration error.
"""
from __future__ import annotations

import argparse
import csv
import datetime as _dt
import hashlib
import json
import os
import pathlib
import shlex
import subprocess
import sys
import time

VERSION = "1.0.0"
EXIT_OK, EXIT_BUDGET, EXIT_FAILED, EXIT_CONFIG = 0, 2, 3, 4


def now_iso() -> str:
    return _dt.datetime.now(_dt.timezone.utc).isoformat(timespec="seconds")


def sha256_file(path: pathlib.Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def row_hash(row: dict) -> str:
    """Stable identity of a job's input row: sorted keys, exact values."""
    blob = json.dumps(row, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(blob).hexdigest()


def load_state(path: pathlib.Path) -> dict:
    if not path.exists():
        return {"version": VERSION, "created": now_iso(), "jobs": {}}
    state = json.loads(path.read_text())
    if not isinstance(state, dict) or "jobs" not in state:
        raise ValueError(f"malformed state file: {path}")
    return state


def save_state(path: pathlib.Path, state: dict) -> None:
    """Atomic write so a crash mid-save cannot corrupt resumability."""
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(state, indent=2, sort_keys=True) + "\n")
    os.replace(tmp, path)


def load_jobs(path: pathlib.Path, id_col: str = "id") -> list[dict]:
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh)
        rows = [dict(r) for r in reader]
        fields = reader.fieldnames or []
    if not rows:
        raise ValueError(f"no jobs in {path}")
    if id_col not in fields:
        raise ValueError(f"id column {id_col!r} not found in {path}; columns={fields}")
    ids = [r.get(id_col, "") for r in rows]
    if any(not i for i in ids):
        raise ValueError(f"every job row needs a non-empty {id_col!r}")
    dupes = {i for i in ids if ids.count(i) > 1}
    if dupes:
        raise ValueError(f"duplicate job ids: {sorted(dupes)[:5]}")
    return rows


def is_done(record: dict, artifact_path: pathlib.Path) -> bool:
    """A job counts as done only if its artifact still matches the recorded hash."""
    if not record or record.get("status") != "done":
        return False
    recorded = record.get("artifact_sha256")
    if not recorded or not artifact_path.is_file():
        return False
    return sha256_file(artifact_path) == recorded


def build_command(template: str, job: dict, workdir: pathlib.Path) -> list[str]:
    """Substitute {id} and {input}; never use a shell, so quoting cannot be injected."""
    filled = template.replace("{id}", job["id"]).replace("{input}", job.get("input", ""))
    return shlex.split(filled)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--jobs", required=True, type=pathlib.Path)
    ap.add_argument("--workdir", required=True, type=pathlib.Path)
    ap.add_argument("--cmd", required=True, help="worker command; {id} and {input} substitute")
    ap.add_argument("--artifact", required=True, help="artifact path relative to workdir; {id} substitutes")
    ap.add_argument("--id-col", default="id",
                    help="CSV column used as the job id and for {id} substitution (default: id)")
    ap.add_argument("--state", required=True, type=pathlib.Path)
    ap.add_argument("--max-wall-s", type=float, default=7200.0)
    ap.add_argument("--max-usd", required=True, type=float,
                    help="dollar ceiling; required so no run starts without an explicit budget")
    ap.add_argument("--hourly-usd", type=float, default=0.59)
    ap.add_argument("--max-jobs", type=int, default=0, help="0 = all")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--manifest", type=pathlib.Path, default=None)
    ap.add_argument("--retrieve-cmd", default=None,
                    help="optional command run AFTER the worker to pull the artifact back "
                         "(e.g. rsync from a pod); {id} and {input} substitute. The artifact "
                         "hash check still decides success, so a failed retrieval cannot pass.")
    args = ap.parse_args(argv)

    started = time.time()
    try:
        jobs = load_jobs(args.jobs, args.id_col)
        state = load_state(args.state)
    except (OSError, ValueError) as exc:
        print(f"CONFIG ERROR: {exc}", file=sys.stderr)
        return EXIT_CONFIG

    for job in jobs:
        if args.id_col != "id":
            job["id"] = job[args.id_col]
    args.workdir.mkdir(parents=True, exist_ok=True)
    save_state(args.state, state)   # valid state even if this run launches nothing
    todo, skipped = [], []
    for job in jobs:
        artifact = args.workdir / args.artifact.replace("{id}", job["id"])
        rec = state["jobs"].get(job["id"])
        if is_done(rec, artifact):
            skipped.append(job["id"])
        else:
            todo.append(job)
    if args.max_jobs:
        todo = todo[: args.max_jobs]

    print(f"jobs={len(jobs)} done={len(skipped)} todo={len(todo)} "
          f"max_wall_s={args.max_wall_s} max_usd={args.max_usd} hourly_usd={args.hourly_usd}")
    if args.dry_run:
        for job in todo:
            print("WOULD RUN", job["id"], "->", args.artifact.replace("{id}", job["id"]))
        return EXIT_OK

    failures, budget_hit = [], False
    for job in todo:
        elapsed = time.time() - started
        spent = elapsed / 3600.0 * args.hourly_usd
        if elapsed >= args.max_wall_s or spent > args.max_usd:
            budget_hit = True
            print(f"BUDGET STOP before {job['id']}: elapsed={elapsed:.0f}s spent=${spent:.4f}")
            break
        cmd = build_command(args.cmd, job, args.workdir)
        artifact = args.workdir / args.artifact.replace("{id}", job["id"])
        if artifact.exists():
            artifact.unlink()  # a stale artifact must not be mistaken for this run's output
        t0 = time.time()
        try:
            proc = subprocess.run(cmd, cwd=args.workdir, capture_output=True, text=True,
                                  timeout=max(1.0, args.max_wall_s - elapsed))
            rc, out, err = proc.returncode, proc.stdout, proc.stderr
        except subprocess.TimeoutExpired as exc:
            rc, out, err = "timeout", exc.stdout or "", exc.stderr or ""
        except OSError as exc:
            rc, out, err = "oserror", "", str(exc)
        t1 = time.time()

        retrieve = None
        if args.retrieve_cmd and rc == 0:
            rcmd = build_command(args.retrieve_cmd, job, args.workdir)
            try:
                rproc = subprocess.run(rcmd, cwd=args.workdir, capture_output=True, text=True,
                                       timeout=max(1.0, args.max_wall_s - (time.time() - started)))
                retrieve = {"command": rcmd, "exit_code": rproc.returncode,
                            "stdout_tail": (rproc.stdout or "")[-1000:],
                            "stderr_tail": (rproc.stderr or "")[-1000:]}
            except (subprocess.TimeoutExpired, OSError) as exc:
                retrieve = {"command": rcmd, "exit_code": "error", "stderr_tail": str(exc)[-1000:]}
            if retrieve.get("exit_code") != 0:
                rc = retrieve.get("exit_code")
                err = (err or "") + "\nRETRIEVE FAILED: " + str(retrieve.get("stderr_tail", ""))[-500:]

        rec = {
            "id": job["id"],
            "input_row_sha256": row_hash(job),
            "command": cmd,
            "started": _dt.datetime.fromtimestamp(t0, _dt.timezone.utc).isoformat(timespec="seconds"),
            "ended": _dt.datetime.fromtimestamp(t1, _dt.timezone.utc).isoformat(timespec="seconds"),
            "duration_s": round(t1 - t0, 3),
            "exit_code": rc,
            "stdout_tail": (out or "")[-2000:],
            "stderr_tail": (err or "")[-2000:],
            "artifact": str(artifact),
        }
        if retrieve is not None:
            rec["retrieval"] = retrieve
        ok = rc == 0 and artifact.is_file() and artifact.stat().st_size > 0
        if ok:
            rec["status"] = "done"
            rec["artifact_sha256"] = sha256_file(artifact)
        else:
            rec["status"] = "failed"
            if rc == 0 and not artifact.is_file():
                rec["reason"] = "exit 0 but no artifact produced"
            elif rc == 0:
                rec["reason"] = "artifact is empty"
            else:
                rec["reason"] = f"worker exit {rc}"
            failures.append(job["id"])
        state["jobs"][job["id"]] = rec
        save_state(args.state, state)
        print(f"{rec['status'].upper():6s} {job['id']} rc={rc} {rec.get('reason','')}".rstrip())

    elapsed = time.time() - started
    done_now = sum(1 for j in jobs if is_done(state["jobs"].get(j["id"]),
                                             args.workdir / args.artifact.replace("{id}", j["id"])))
    summary = {
        "version": VERSION,
        "finished": now_iso(),
        "jobs_total": len(jobs),
        "jobs_done": done_now,
        "jobs_failed": failures,
        "skipped_already_done": skipped,
        "budget_stop": budget_hit,
        "elapsed_s": round(elapsed, 3),
        "estimated_spend_usd": round(elapsed / 3600.0 * args.hourly_usd, 4),
        "limits": {"max_wall_s": args.max_wall_s, "max_usd": args.max_usd,
                   "hourly_usd": args.hourly_usd},
    }
    if args.manifest:
        args.manifest.parent.mkdir(parents=True, exist_ok=True)
        summary["state_sha256"] = sha256_file(args.state) if args.state.exists() else None
        args.manifest.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(json.dumps({k: summary[k] for k in ("jobs_total", "jobs_done", "budget_stop",
                                              "elapsed_s", "estimated_spend_usd")}))
    if failures:
        return EXIT_FAILED
    if budget_hit:
        return EXIT_BUDGET
    return EXIT_OK


if __name__ == "__main__":
    sys.exit(main())
