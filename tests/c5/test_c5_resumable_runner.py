#!/usr/bin/env python3
"""Fail-closed tests for c5_prep/resumable_runner.py. No GPU, no network, no sequences.

A fake worker is generated per test so the runner's contract is tested without touching
any real pipeline. Every guarantee claimed in the module docstring gets a test.
"""
import json, pathlib, subprocess, sys, tempfile, textwrap

ROOT = pathlib.Path(__file__).resolve().parents[2]
RUNNER = ROOT / "c5_prep" / "resumable_runner.py"
res = []


def check(name, cond, detail=""):
    res.append(bool(cond))
    print(("PASS " if cond else "FAIL ") + name + ("" if cond else f"  {detail}"))


def run(args, **kw):
    return subprocess.run([sys.executable, str(RUNNER), *args], capture_output=True, text=True, **kw)


def make_fake(tmp, body, name="fake_worker.py"):
    """Write a fake worker. `name` must be unique when two fakes coexist in one run."""
    p = tmp / name
    p.write_text("#!/usr/bin/env python3\nimport sys, pathlib\n" + textwrap.dedent(body))
    return p


JOBS = "id,input\nj1,a\nj2,b\nj3,c\n"

with tempfile.TemporaryDirectory() as td:
    tmp = pathlib.Path(td)
    jobs = tmp / "jobs.csv"; jobs.write_text(JOBS)

    # --- 1. happy path: all jobs run, artifacts hashed, manifest written
    worker = make_fake(tmp, """
        p = pathlib.Path(sys.argv[1] + '.out')
        p.write_text('result for ' + sys.argv[1])
        print('made', p.name)
    """)
    st = tmp / "state1.json"; man = tmp / "manifest1.json"
    r = run(["--jobs", str(jobs), "--workdir", str(tmp / "w1"), "--cmd",
             f"{sys.executable} {worker} {{id}} {{input}}", "--artifact", "{id}.out",
             "--state", str(st), "--manifest", str(man), "--max-usd", "5.00"])
    check("happy path exit 0", r.returncode == 0, r.stdout + r.stderr)
    state = json.loads(st.read_text())
    check("happy path records 3 done", sum(1 for v in state["jobs"].values() if v["status"] == "done") == 3)
    check("happy path records artifact hash",
          all(len(v.get("artifact_sha256", "")) == 64 for v in state["jobs"].values()))
    check("happy path records input row hash",
          all(len(v.get("input_row_sha256", "")) == 64 for v in state["jobs"].values()))
    m = json.loads(man.read_text())
    check("manifest records state hash", len(m.get("state_sha256", "")) == 64)
    check("manifest reports jobs_done=3", m.get("jobs_done") == 3)

    # --- 2. resume: a second run launches nothing and still exits 0
    r2 = run(["--jobs", str(jobs), "--workdir", str(tmp / "w1"), "--cmd",
              f"{sys.executable} {worker} {{id}} {{input}}", "--artifact", "{id}.out",
              "--state", str(st), "--max-usd", "5.00"])
    check("resume exit 0", r2.returncode == 0, r2.stdout + r2.stderr)
    check("resume launches nothing", "todo=0" in r2.stdout, r2.stdout)
    check("resume reports all done", '\"jobs_done\": 3' in r2.stdout, r2.stdout)

    # --- 3. stale artifact: artifact modified after completion -> NOT treated as done
    (tmp / "w1" / "j2.out").write_text("tampered-content-that-differs")
    prior_hash = json.loads(st.read_text())["jobs"]["j2"]["artifact_sha256"]
    r3 = run(["--jobs", str(jobs), "--workdir", str(tmp / "w1"), "--cmd",
              f"{sys.executable} {worker} {{id}} {{input}}", "--artifact", "{id}.out",
              "--state", str(st), "--max-usd", "5.00"])
    check("tampered artifact re-runs that job", "todo=1" in r3.stdout, r3.stdout)
    check("tampered artifact run exit 0", r3.returncode == 0)
    state3 = json.loads(st.read_text())
    import hashlib as _h
    now_hash = _h.sha256((tmp / "w1" / "j2.out").read_bytes()).hexdigest()
    check("tampered artifact was regenerated", "tampered" not in (tmp / "w1" / "j2.out").read_text())
    check("state hash matches the regenerated artifact",
          state3["jobs"]["j2"]["artifact_sha256"] == now_hash)
    check("regenerated job is done again", state3["jobs"]["j2"]["status"] == "done")

with tempfile.TemporaryDirectory() as td:
    tmp = pathlib.Path(td)
    jobs = tmp / "jobs.csv"; jobs.write_text(JOBS)

    # --- 4. worker failure: nonzero exit -> failed, nonzero runner exit, not "done"
    bad = make_fake(tmp, "sys.exit(7)")
    st = tmp / "state_bad.json"
    r = run(["--jobs", str(jobs), "--workdir", str(tmp / "wb"), "--cmd",
             f"{sys.executable} {bad} {{id}} {{input}}", "--artifact", "{id}.out",
             "--state", str(st), "--max-usd", "5.00"])
    check("worker failure -> runner exit 3", r.returncode == 3, r.stdout + r.stderr)
    state = json.loads(st.read_text())
    check("failed job not marked done", all(v["status"] == "failed" for v in state["jobs"].values()))

    # --- 5. exit 0 but no artifact -> still a failure (fail closed)
    silent = make_fake(tmp, "pass")
    st2 = tmp / "state_silent.json"
    r = run(["--jobs", str(jobs), "--workdir", str(tmp / "ws"), "--cmd",
             f"{sys.executable} {silent} {{id}} {{input}}", "--artifact", "{id}.out",
             "--state", str(st2), "--max-usd", "5.00"])
    check("no artifact -> exit 3", r.returncode == 3, r.stdout + r.stderr)
    check("no artifact records a reason", "no artifact" in r.stdout, r.stdout)

    # --- 6. budget stop: max_jobs=1 with 3 jobs -> exit 2, only 1 launched
    good = make_fake(tmp, "pathlib.Path(sys.argv[1] + '.out').write_text('ok')")
    st3 = tmp / "state_budget.json"
    r = run(["--jobs", str(jobs), "--workdir", str(tmp / "wbu"), "--cmd",
             f"{sys.executable} {good} {{id}} {{input}}", "--artifact", "{id}.out",
             "--state", str(st3), "--max-jobs", "1", "--max-usd", "5.00"])
    check("max_jobs caps launches", "todo=1" in r.stdout, r.stdout)
    check("max_jobs leaves remaining jobs pending",
          sum(1 for v in json.loads(st3.read_text())["jobs"].values() if v["status"] == "done") == 1)
    r = run(["--jobs", str(jobs), "--workdir", str(tmp / "wbu"), "--cmd",
             f"{sys.executable} {good} {{id}} {{input}}", "--artifact", "{id}.out",
             "--state", str(st3), "--max-usd", "5.00"])
    check("second capped run resumes the remaining 2", "todo=2" in r.stdout and r.returncode == 0,
          r.stdout)

    # --- 7. wall ceiling already exceeded -> refuses to launch anything
    st4 = tmp / "state_wall.json"
    r = run(["--jobs", str(jobs), "--workdir", str(tmp / "ww"), "--cmd",
             f"{sys.executable} {good} {{id}} {{input}}", "--artifact", "{id}.out",
             "--state", str(st4), "--max-wall-s", "0", "--max-usd", "5.00"])
    check("zero wall budget -> exit 2", r.returncode == 2, r.stdout + r.stderr)
    check("zero wall budget launched nothing", "BUDGET STOP" in r.stdout, r.stdout)
    check("zero wall budget leaves valid state", json.loads(st4.read_text())["jobs"] == {})

    # --- 8. duplicate ids -> config error, no jobs launched
    dup = tmp / "dup.csv"; dup.write_text("id,input\nj1,a\nj1,b\n")
    r = run(["--jobs", str(dup), "--workdir", str(tmp / "wd"), "--cmd", "true",
             "--artifact", "{id}.out", "--state", str(tmp / "sd.json"), "--max-usd", "5.00"])
    check("duplicate ids -> exit 4", r.returncode == 4, r.stdout + r.stderr)

    # --- 9. malformed state file -> config error (never silently discards history)
    bad_state = tmp / "bad_state.json"; bad_state.write_text('{"nope": 1}')
    r = run(["--jobs", str(jobs), "--workdir", str(tmp / "wz"), "--cmd", "true",
             "--artifact", "{id}.out", "--state", str(bad_state), "--max-usd", "5.00"])
    check("malformed state -> exit 4", r.returncode == 4, r.stdout + r.stderr)

    # --- 10. dry run launches nothing and writes no state
    st5 = tmp / "state_dry.json"
    r = run(["--jobs", str(jobs), "--workdir", str(tmp / "wdr"), "--cmd", "true",
             "--artifact", "{id}.out", "--state", str(st5), "--max-usd", "5.00", "--dry-run"])
    check("dry run exit 0", r.returncode == 0, r.stdout + r.stderr)
    check("dry run launches nothing", "BUDGET STOP" not in r.stdout and "WOULD RUN" in r.stdout)
    check("dry run records no completed jobs",
          json.loads(st5.read_text())["jobs"] == {})
    check("dry run lists all 3", r.stdout.count("WOULD RUN") == 3, r.stdout)

    # --- 11. command is executed without a shell (no injection, no globbing)
    inj = make_fake(tmp, "pathlib.Path(sys.argv[1] + '.out').write_text('ok')")
    st6 = tmp / "state_inj.json"
    r = run(["--jobs", str(jobs), "--workdir", str(tmp / "wi"), "--cmd",
             f"{sys.executable} {inj} {{id}} '$(touch pwned)'", "--artifact", "{id}.out",
             "--state", str(st6), "--max-usd", "5.00"])
    check("shell metacharacters are not evaluated", not (tmp / "wi" / "pwned").exists())
    check("metachar run still completes jobs", r.returncode == 0, r.stdout + r.stderr)

failed = res.count(False)
print(f"\n{len(res) - failed}/{len(res)} passed")
# (original summary moved to the end so retrieve-cmd tests execute)

# --- 12. --retrieve-cmd: runs after a successful worker and its failure must block success
#     The artifact hash check remains the final gate. These tests exist because the
#     retrieval path was added without coverage and its whole purpose is fail-closed.
with tempfile.TemporaryDirectory() as td:
    tmp = pathlib.Path(td)
    jobs = tmp / "jobs.csv"; jobs.write_text(JOBS)

    # 12a. worker writes artifact; retriever runs and succeeds -> done, retrieval recorded
    w = make_fake(tmp, "pathlib.Path(sys.argv[1] + '.out').write_text('result')")
    ret_ok = make_fake(tmp, """
        # retriever proves {id} substitution reached it by writing a sidecar
        pathlib.Path(sys.argv[1] + '.retrieved').write_text('pulled ' + sys.argv[1])
    """, name="fake_retriever_ok.py")
    st = tmp / "st_ret_ok.json"
    r = run(["--jobs", str(jobs), "--workdir", str(tmp / "wr"), "--cmd",
             f"{sys.executable} {w} {{id}} {{input}}", "--artifact", "{id}.out",
             "--state", str(st), "--max-usd", "5.00",
             "--retrieve-cmd", f"{sys.executable} {ret_ok} {{id}}"])
    check("retrieve ok: exit 0", r.returncode == 0, r.stdout + r.stderr)
    state = json.loads(st.read_text())
    check("retrieve ok: all 3 done",
          sum(1 for v in state["jobs"].values() if v["status"] == "done") == 3)
    check("retrieve ok: retrieval recorded with exit_code 0",
          all(v.get("retrieval", {}).get("exit_code") == 0 for v in state["jobs"].values()))
    check("retrieve ok: retriever actually ran with {id} substituted",
          all((tmp / "wr" / f"{j}.retrieved").exists() for j in ("j1", "j2", "j3")))

    # 12b. worker succeeds, retriever FAILS -> job must NOT be done; runner exits nonzero
    w2 = make_fake(tmp, "pathlib.Path(sys.argv[1] + '.out').write_text('result')")
    ret_bad = make_fake(tmp, "sys.exit(9)", name="fake_retriever_bad.py")
    st2 = tmp / "st_ret_bad.json"
    r = run(["--jobs", str(jobs), "--workdir", str(tmp / "wrb"), "--cmd",
             f"{sys.executable} {w2} {{id}} {{input}}", "--artifact", "{id}.out",
             "--state", str(st2), "--max-usd", "5.00",
             "--retrieve-cmd", f"{sys.executable} {ret_bad} {{id}}"])
    check("retrieve fail: runner exit nonzero", r.returncode != 0, r.stdout + r.stderr)
    state2 = json.loads(st2.read_text())
    check("retrieve fail: no job marked done",
          sum(1 for v in state2["jobs"].values() if v["status"] == "done") == 0)
    check("retrieve fail: reason names the retrieval failure",
          all("RETRIEVE FAILED" in v.get("stderr_tail", "") for v in state2["jobs"].values()))
    check("retrieve fail: retrieval exit_code preserved as 9",
          all(v.get("retrieval", {}).get("exit_code") == 9 for v in state2["jobs"].values()))

    # 12c. retriever runs but leaves no artifact -> still failure (hash gate is final)
    w3 = make_fake(tmp, "pass")
    ret_noop = make_fake(tmp, "pass", name="fake_retriever_noop.py")
    st3 = tmp / "st_ret_noop.json"
    r = run(["--jobs", str(jobs), "--workdir", str(tmp / "wrn"), "--cmd",
             f"{sys.executable} {w3} {{id}} {{input}}", "--artifact", "{id}.out",
             "--state", str(st3), "--max-usd", "5.00",
             "--retrieve-cmd", f"{sys.executable} {ret_noop} {{id}}"])
    check("retrieve no-artifact: runner exit 3", r.returncode == 3, r.stdout + r.stderr)
    check("retrieve no-artifact: nothing marked done",
          sum(1 for v in json.loads(st3.read_text())["jobs"].values() if v["status"] == "done") == 0)

    # 12d. worker fails -> retriever must NOT run (no retrieval key recorded)
    w4 = make_fake(tmp, "sys.exit(5)")
    ret_marker = make_fake(tmp, "pathlib.Path('RETRIEVER_RAN').write_text('x')",
                           name="fake_retriever_marker.py")
    st4 = tmp / "st_ret_skip.json"
    r = run(["--jobs", str(jobs), "--workdir", str(tmp / "wrs"), "--cmd",
             f"{sys.executable} {w4} {{id}} {{input}}", "--artifact", "{id}.out",
             "--state", str(st4), "--max-usd", "5.00",
             "--retrieve-cmd", f"{sys.executable} {ret_marker} {{id}}"])
    check("retriever skipped when worker fails: exit nonzero", r.returncode != 0)
    check("retriever skipped when worker fails: no RETRIEVER_RAN marker",
          not (tmp / "wrs" / "RETRIEVER_RAN").exists())
    check("retriever skipped when worker fails: no retrieval record",
          all("retrieval" not in v for v in json.loads(st4.read_text())["jobs"].values()))

failed2 = res.count(False)
print(f"\nTOTAL {len(res) - failed2}/{len(res)} passed")
sys.exit(1 if failed2 else 0)
