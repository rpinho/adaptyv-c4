# Overnight gap, CPU repair and resumed optimization — 2026-10-22

## Accountability / timeline

The last Oct 21 GPU diagnostic started at 20:26:51 EDT and finished at
20:28:34 EDT. Its GPU release was independently verified at **20:28:54 EDT**.
At Oct 22 **07:04:26 EDT**, both Modal functions still had zero runners, inputs
and backlog. There was no overnight design-model run or new candidate.

The agent completed evening consultation / application work, then left no
continuing research job running. This did not satisfy Ricardo's instruction
to keep improving models. Rosalind identity verification was a separate task
and did not require suspending the authorized research. No claim of overnight
modeling, training, CPU research or autonomous progress should be made.

## CPU work performed after the morning question

Implemented scientific outcome accounting in `scripts/campaign_admission.py`:
generator rows, native-context admitted trajectories and pH-engineering proxies
are separate. Admission requires completed error-free scientific follow-up,
verified native sequence, the existing backbone/NAG clearances, and a matching
canonical structure SHA. Exit code zero alone cannot establish admission.
Unknown/conflicting trajectory identities fail closed. Sequence siblings count
once, both before bounded CPU parent selection and across controller batches.
These identities do not establish structural diversity or binding.

Legacy raw-row search states cannot resume into new allocations. Stage labels
alone no longer trigger automatic length or confidence-weight changes. An
admitted structure still requires acquisition review before more sampling.
The old v5 controller remains cancelled; it was not restarted.

All **55 local regressions passed**; see `regression_tests.log`. Replaying the
saved native batches with the new accounting still yields zero admitted
trajectories, as expected (`historical_outcomes.json`). No acceptance thresholds
were relaxed and no old result was promoted.

Astra's new CPU/source review is in `../astra_review/ACQUISITION_REVIEW.md`.
It found no hotspot-mapping or interface-loss-routing bug in the pinned source.
The next missing experiment is actual optimization with the repaired frames;
the previous v8 result was just one initial forward/gradient.

## Single repaired-native optimization test

New run `c4_repaired_native_v9` preserves v5's seed20261021, first length69,
complete native target in the audited local frame, loss/optimizer settings,
forced_targeting=false and acceptance gates. It permits **one trajectory**,
with no automatic successor. The scientific change from collapsed v5 is the
initial/recycled frame repair. It is not a completed matched endpoint comparison
until results exist. No claim of exact byte-for-byte first-forward equivalence
is made; the first real input sequence-array SHA/shape/dtype is now recorded.

CPU image preflight prepares the real first target/binder state and checks its
repaired frames using the actual upstream converter before allocation. The
GPU worker wraps the installed CLI only to observe gradient calls, record
physical requested-patch distance / CA geometry, and abort on nonfinite loss,
gradients or collapse. Compile-only calls pass through. Actual outputs are
returned unchanged to the original optimizer. Runtime module/monitor hashes
are recorded; physical observations are not candidate acceptance gates.

Worker runtime cap30min (termination begins one minute before), local monitor
cap35min, underlying function max2h as a final infrastructure fallback. The
local monitor cancels the exact call on its deadline and independently verifies
actual zero runners/inputs/backlog before CPU follow-up. Single-use containers,
zero retries, no Runpod start and no automatic successor.

Current launch/result state belongs in `../modal_repaired_native_v9/supervisor.json`;
this decision file alone is not evidence that a GPU is running or a run passed.

## Spend

Before this new run: newly supervised Runpod runtime estimate **$1.7952**;
cumulative Modal runtime estimate **$4.0531838**; external consultation reported
usage **$1.1922999**. Historical campaigns, CPU image build and setup/storage
remain separate/unreconciled.

[Modal pricing](https://modal.com/pricing), checked Oct22: H100GPU
$3.9492/hour; A10080GB fallback $2.4984/hour. Requested4physicalCPUcores and64GiB
add at most$0.700128/hour. H100 runtime upper **$4.649328/hour**; the35min local
envelope is approximately **$2.7121**, plus setup/storage. Existing continued
modeling/spending authorization applies; no repeated permission wait.
