# Active bounded modeling campaign

**Historical attempt:** this full-dimer pilot failed with OOM. The supervisor
stopped the pod and independently verified EXITED; estimated compute $0.3855.
Its CPU follow-up correctly found no accepted scaffold. **Current activity is
in `../gpu_campaign_v2/README.md`, under the same combined $3.58 compute cap.**

Authorization: Ricardo's explicit instruction to keep modeling and use GPU funds
without waiting for another prompt. Accounts and entries remain Ricardo-only.
Everything is on `codex/install-two-state-switch`; no KB/main writes.

Pod **2p9x2qea1ietbw**, migrated from the stopped C3 pod after both original hosts
could not reacquire GPUs. A100-SXM4-80GB confirmed by nvidia-smi. MCP get-pod cost
**$1.79/hour**; console displayed$1.80/hour. Both original pods remain stopped.
Fresh balance $132.71, pre-start idle/storage rate $0.019/hour. New compute cap
**$3.58**, three minutes reserved for stop retries. Runtime estimates are not
settled billing; `billing_read.json` can lag.

The first pilot uses the entire resolved 11QC human dimer (1278 residues) and a
compact accessible apical patch A205,A206,A208,A209,A210. One worker, at most
3 trajectories or 2 accepted designs. Default backbone/confidence filters remain;
accessible-termini objective added for the binder's capture geometry. Only
accepted native-target refolds can enter CPU pKa scoring; lysinated trajectory
structures are forbidden there. Glycans, assay scar and capture tag are not
modeled; mouse cross-reactivity and independent folding validation are pending.

Live files:

- `supervisor.json`: current worker/GPU evidence, estimated compute and stop state.
- `worker_tail.log`: mirrored GPU log, including compilation/memory failures.
- `mirrored/`: campaign artifacts copied off the pod, excluding compile caches.
- `continuation.json`: automatic local CPU follow-up after verified pod shutdown.
- `../gpu_followup/status.json`: appears when that CPU queue begins.
- `../full_context_audit/RESULTS.md`: measured reason for rebuilding generation.
- `../free_networks_v3/RESULTS.md`: independently prepared free-module scout.

Supervisor PID and worker PID are in the live JSON; inspect before acting, since
PIDs can change. The supervisor stops the pod on worker exit, three consecutive
SSH failures, a supervisor exception, or the cap, and checks EXITED with a fresh
API read. `STOP_FAILED_REQUIRES_ATTENTION` means actual pod shutdown was not
verified. Do not equate a stopped Python worker with stopped billing.

Emergency pod stop using the existing authenticated Runpod MCP session:

```sh
python3 scripts/runpod_control.py --tool pod-action --arguments '{"id":"2p9x2qea1ietbw","body":{"action":"stop"}}'
python3 scripts/runpod_control.py --tool get-pod --arguments '{"id":"2p9x2qea1ietbw"}'
```

No credentials are copied into this repo; the control script reads the existing
local OAuth entry and sends it only to its matching Runpod endpoint. Never commit
raw get-pod records containing environment values. `live_pod.json` is filtered.

Shutdown regression checks: four mocked tests cover exit, cap, mirror failure,
and initial API-read failure. Native-context admission tests reject cropped and
lysinated inputs. Run the full suite with `.venv-rosetta/bin/python -m unittest
discover -s tests -v`.
