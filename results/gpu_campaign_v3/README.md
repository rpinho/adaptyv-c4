# Next native-target control: automatic capacity queue

The previous GPU pilot completed three trajectories with no accepted design.
Its pod was stopped and independently verified EXITED. The CPU continuation's
startup race is repaired and replayed; there is no accepted scaffold to score.
See `../gpu_campaign_v2/RESULTS.md`.

Both attempted existing host restarts were refused for lack of free GPUs. The
global catalog has A100-SXM stock in other US regions, but the cached environment
and weights are on the Iceland network volume. The Iceland inventory has no
A100 availability in this read. No new GPU resource or account was created.

**Live state is `queue.json`.** The local queue checks the existing migrated
pod at ten-minute intervals for at most six hours. Waiting rents no compute.
It launches only the predeclared control in `DECISION.md`, on the same A100 at
no more than the quoted **$1.79/hour**, and only if the pod was EXITED. It leaves
an already-running pod alone. GPU/SSH boot failure causes actual pod shutdown.
There is no need for another user prompt or approval.

On success it verifies supervisor ownership and starts the repaired CPU
continuation. Live GPU progress then moves to `supervisor.json`; the supervisor
stops the POD on exit/failure/cap and verifies EXITED. Combined compute cap remains
**$3.58**, deducting the prior **$1.7952** runtime estimate. No more than two new
trajectories. Billing may lag; network storage is additional.

`capacity_wait_expired`, `failed`, or `STOP_FAILED_REQUIRES_ATTENTION` are explicit
end states, not claims of continued modeling. Waiting for capacity is not a GPU
model run. All 21 regressions pass, including rate/ownership/allocation guards
and the follow-up startup race.

To cancel an unlaunched queue, create `results/gpu_campaign_v3/CANCEL`. If the
queue already delegated a GPU run, use the existing Runpod stop/verification
commands in the v2 README; do not merely kill its Python worker. Account creation
and challenge submission remain Ricardo-only. No main or KB writes.
