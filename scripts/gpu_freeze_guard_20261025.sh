#!/bin/zsh
# C3 GPU freeze guard. Closes any C4 GPU spend for the remainder of the C3 sprint.
# Idempotent. Stops Modal c4-apical apps and RunPod pods named ^adaptyv-c4 only, and kills local
# launchers. Non-C4 pods (e.g. a peer session's live C3 job) are LOGGED and LEFT RUNNING.
# Logs every action. Does not touch C3 artifacts, frozen files, the app-server, or anything else.
set -u
cd .
LOG=results/orchestrator_20261025/gpu_freeze_guard.log
mkdir -p results/orchestrator_20261025
ts(){ date '+%Y-%m-%d %H:%M:%S'; }

# 1) kill local C4 launchers before they can re-create anything
for p in $(pgrep -f 'modal_c4_apical_canary' 2>/dev/null); do
  kill -TERM "$p" 2>/dev/null && echo "$(ts) killed launcher $p" >> "$LOG"
done

# 2) stop every live c4-apical Modal app
live=$(timeout 60 .venv-modal/bin/modal app list 2>/dev/null | python3 -c '
import re,sys
t=sys.stdin.read()
for line in t.splitlines():
    if "c4-apica" in line or "c4-apical" in line:
        m=re.match(r"\s*(ap-[A-Za-z0-9]+)", line)
        if m and "stopped" not in line and "stopping" not in line:
            print(m.group(1))
' 2>/dev/null | sort -u)
for app in ${(f)live}; do
  [ -z "$app" ] && continue
  timeout 90 .venv-modal/bin/modal app stop -y "$app" >/dev/null 2>&1 \
    && echo "$(ts) modal stop $app" >> "$LOG"
done

# 3) stop any *C4* RunPod pod that is not EXITED.
#    SCOPED ON PURPOSE (2026-10-25 20:21, orchestrator cycle): the original version stopped ANY
#    non-EXITED pod. That would have killed another session's live C3 measurement, which is exactly
#    the destructive overshoot to avoid. Only pods whose name matches the C4 baseline prefix are
#    stopped; everything else is logged and left running so the owner can see it.
python3 scripts/runpod_control.py --tool list-pods 2>/dev/null | python3 -c '
import json,sys,re
try:
    d=json.load(sys.stdin); t=json.loads(d["content"][0]["text"])
except Exception:
    sys.exit(0)
for p in t.get("pods",[]):
    if p.get("status") == "EXITED":
        continue
    name = p.get("name") or ""
    if re.match(r"^adaptyv-c4", name):
        print("STOP", p["id"], name)
    else:
        print("SKIP", p["id"], name)
' 2>/dev/null | while read -r action pod name; do
  [ -z "$pod" ] && continue
  if [ "$action" = "STOP" ]; then
    python3 scripts/runpod_control.py --tool pod-action --arguments "{\"id\":\"$pod\",\"body\":{\"action\":\"stop\"}}" >/dev/null 2>&1 \
      && echo "$(ts) runpod stop $pod ($name)" >> "$LOG"
  else
    echo "$(ts) runpod leave-running $name ($pod) [not C4-scoped]" >> "$LOG"
  fi
done

# 4) heartbeat: prove the guard ran, even when it had nothing to stop.
#    Without this the log is silent on a healthy run, which is indistinguishable from a dead guard.
echo "$(ts) heartbeat ok" >> "$LOG"
