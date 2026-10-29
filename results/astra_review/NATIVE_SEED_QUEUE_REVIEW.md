# Native seed queue review — 2026-10-22

## Verdict

**No remaining launch blocker found in the reviewed finite controller after the corrections below.** This supports at most five declared native-input seed replicates, conditional on the installed-image CPU preflights and existing resource checks. It is not evidence of a geometry cure, a favorable hit rate, or authorization for an open-ended search. Any admitted output stops the queue for review.

Reviewed source hashes:

- `scripts/run_native_replicates.py`: `8ad0a7398965bc1122df2bcb185b6806acd76a23956aee6b34f0d267b32d6f35`
- `scripts/monitored_design.py`: `0474a7635610de77cd0514bf0598c4dbe230575fb51882a252a7f2bb09e84f3e`
- `scripts/modal_campaign.py`: `47e56b7e86c4b491dc01c37645ccc585a50d560ff9ac6cec52d645cefc712921`

I also inspected the shared driver, admission evaluator, policy builder and queue/monitor tests. This review allocated no GPU and changed no implementation. Parent reports V12 completed/rejected and actual GPU release verified; prior cumulative Modal runtime estimate$6.3849882, Runpod$35.6652 including historical usage, consultations$1.1923, with extras unreconciled.

## Scientific isolation and seed routing

I executed the actual local `make_config()` for seeds20261022–20261026 and compared each resulting dictionary with archived **actual v9** `results/modal_repaired_native_v9/mirrored/campaign.json`. Every difference was confined to `campaign_seed`, `campaign_name` and `project_folder`. Native target, original65–75 length distribution, frame-initialization setting, cutoff20Å, weights, filters and `max_trajectories=1` remain unchanged. Forced targeting is false.

The first-trajectory CPU key construction matches the pinned execution path: campaign folds trajectory number1 into the campaign PRNG key; `run_trajectory()` splits that key and passes its first child to initialization. The new replicate-name/seed checks map indices2–6 to the declared five seeds and allow the original length distribution instead of incorrectly demanding69 residues. Observation does not draw another random key or change the gradient returned to the optimizer.

This is a study of **stochastic variation under the fixed sampling distribution**. Length and model/sequence sampling may vary as consequences of the seed. It does not isolate sequence initialization while holding every realized trajectory at69 residues. Different PRNG seeds are the operational independence mechanism; do not infer biological diversity or independent successful scaffolds from seed count alone.

The image build invokes each actual declared seed's CPU initialization before any GPU function allocation. Parent reports all five passed, with lengths69,67,68,74,75, complete639-residue target, zero sequence differences and valid frames. I verified the wiring/source but did not independently rerun the installed-image stack. Preserve these actual reports and first-step sequence hashes with each run. Parent reports73 local tests passing; I inspected the new tests rather than rerunning the full suite.

## Controller and corrected operational boundaries

- Finite seed tuple, one trajectory per call, one series lock and per-call driver locks prevent ordinary duplicate queue launches. Completed/cancelled/review/failure states do not automatically restart. Resume validates the seed declaration and maximum batch count; confirmed existing call IDs are resumed rather than respawned. Uncertain spawn intent remains a stop requiring inspection.
- The new `require_idle()` checks **both** deployed GPU functions for zero runners, running inputs and backlog before each genuinely new spawn, including after resumed/cached prior batches. It fails closed on a busy cloud state. A historical `gpu_released_verified` flag alone is no longer the pre-spawn evidence. This is an observation immediately before launch, not a distributed lock against unrelated concurrent launchers; no other launchers should run during this study.
- Normal continuation requires successful remote completion, actual scale-to-zero verification, CPU follow-up exit success and error-free scientific evidence. Native-context output stops for review regardless of any later engineering-advance flag. Scientific errors, worker failures, cancellation and release failures stop successors. The strict native admission evaluator remains in use.
- Remote monitored worker budget is30 minutes (termination begins with its existing60-second margin); local monitoring cancels after35 minutes, with polling and shutdown overhead. Neither is an exact invoice ceiling under an API outage. The separately declared45-minute **planning reserve** is$3.486996 per call/$17.43498 for five. Completed submission-to-verified-release runtime is charged to the planning ledger, and another full reserve must fit before the next call.
- The original35-minute reserve omitted polling/teardown overhead. The revised reserve and honest planning label address this. Failure paths now account for submitted runtime and persist emergency release evidence. If actual release cannot be confirmed, the estimate is only through the observation time and the queue stops; resource status still needs resolution. Active runtime is explicitly additional to the completed-batch ledger.

No claim is made that image build, storage or consultation charges are covered by this runtime planning figure. Historical Runpod usage must not be relabeled as new queue cost.

## Interpretation and remaining limits

Five seeds can reveal whether the observed off-patch/malformed outcome repeats under this declared configuration; they cannot establish a reliable success probability. Near-patch overlaps remain failures. Lower scalar losses and higher model confidence do not override physical geometry, native chemistry, full receptor context or unchanged final gates.

The monitor still aborts exact/near collapse and nonfinite output, while less extreme malformed geometry is logged rather than silently treated as valid. That behavior is unchanged from the reference and does not lower admission requirements. A successful native-input search is still computational evidence, not measured affinity, pH selectivity or mouse cross-reactivity.

The domain adapter remains a separate CPU-only alternative. Nothing in this queue routes into it, changes profiles in response to failure, admits a reduced-domain score or extends the predeclared five-seed study automatically.
