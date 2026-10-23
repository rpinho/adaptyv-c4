# Generation pivot readiness — 23 October 2026

The latest Ricardo instruction transfers work to the existing generation lane
after the cached novelty export is recovered, validated and committed. This
record prepares that transition; no generator or GPU has been started.

Use the canonical launcher on `claude/submission-and-control`, commit `e8fdfa9`,
in `/Users/ricardo/.worktrees/adaptyv-c4/claude-submission`. Its tracked files are
clean. Launcher SHA256:
`3554f595c0c20a98aff02356d81800828d22619e569acda99777a4f1b348e609`.
Do not build another generator or exporter. Ricardo retains accounts/submission.

First paid stage, approved by Ricardo after pricing in GPU_APPROVAL.json: the two existing 200-input
arms, batch8, seed42, inline `geom`, unchanged input specs. Fresh destinations
`tfr1_c200m_20261023` and `tfr1_his200m_20261023` are absent from the volume
listing. The read-only Modal app listing showed no running generation app.

Use the existing function with invocation-specific timeout9000s,
max_containers1, retries0 and short idle scale-down. At campaign completion or
stop, stop its app and verify no running inputs/allocations; stopping local
workers is insufficient. No Runpod pod is needed. Preserve app/call IDs,
start/finish times, output receipts and compute-cost scope. No detached
unbounded app or automatic retry/campaign successor.

Live public [Modal pricing](https://modal.com/pricing), checked23 October:
L40S $0.000542/s = **$1.9512 per GPU-hour**; two concurrent GPUs $3.9024/hour.
CPU and host memory are additional metered charges. Prior handoff estimated
about$4 per200 designs, or about$8 for these two arms; this is an estimate from
previous runs, not a guaranteed cost or measured geom-campaign runtime.
9000s per arm bounds GPU execution time, not the full account bill.
[Invocation options](https://modal.com/docs/reference/modal.Function#with_options)
are supported by both the installed SDK signature and primary documentation.

Read-only metric preflight downloaded the existing four-row inline metric
control, `geom_control_metrics.csv`. It contains actual hotspot obedience,
clash fractions, unreasonable-bond fractions and conditioning adherence.
It contains no iPTM, pLDDT, pAE, interface ddG or developability scores.
The existing ranker weights old confidence column names; feeding these raw
metrics into that policy would not establish a meaningful ranking. Join
verified sequences to actual metrics and explicitly adapt/report the ranking
policy before using it. Missing scores must stay unknown. No ranking or
submission-readiness claim is made by this preflight.

After these measured arms and their ranking, the requested higher-His arms
and restored-context audit remain next. No changed conditioning or metric
implementation is claimed tested here. Novelty analysis is sufficient to
handoff after the current completion; do not start another novelty campaign.

Recorded own-lane spend: Runpod$56.589988 + Modal estimate$10.718947 = about
$67.31. New GPU/cloud compute$0.00. Claude's approximately$15 and the complete
account bill remain separately unverified.
