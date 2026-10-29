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

## Completed v9 / conditional v10 follow-up

V9 completed50 actual gradient rounds on H100 with finite gradients and no exact
collapse. It yielded **0 accepted designs**, rejected at screen for pLDDT0.54
(i_pTM0.17). No round approached the requested patch: closest CA distance38.512A
was the initial prediction; final54.747A. The lowest-loss round was24, not the
final round. Minimum adjacent-CA distance over the run fell to0.376A, so
noncollapse must not be described as a physically valid fold. Raw round evidence
and summary: `../modal_repaired_native_v9/mirrored/physical_steps.jsonl`
and `../modal_repaired_native_v9/RESULTS.json`.

Actual GPU release was verified at07:24:09EDT with zero runners/inputs/backlog;
CPU follow-up returned0 and explicitly no accepted scaffolds. Conservative
submission-through-release runtime estimate **$0.5836719**, Modal cumulative
**$4.6368557**. This morning run is separate from the overnight gap.

Astra supported a conditional **single-variable objective test**, not a promised
fix or proven superior hit-rate strategy: narrow interface_contacts cutoff20→8A.
Its exact CPU loss check found this removes the low-cost19A loophole, but does
not establish a new far-distance gradient or solve acquisition. Therefore the
motivation is stricter proximity, not rescuing absent/saturated guidance.

Prepared `c4_interface8_v10` only after the completed off-patch v9 result. Same
native target/local frame, first seed/length, repaired frames, all other loss
weights and acceptance gates; one trajectory and identical runtime/release caps.
Monitoring now also records already-weighted scalar loss components. These
describe numerical loss reduction, not gradient dominance. Raw loss values
cannot be compared across different objective cutoffs as if unchanged.
All **56** local regressions pass including the exact configuration delta.
New image deployment builds on CPU before any new allocation; live state for
this prospective test will be `../modal_interface8_v10/supervisor.json`.

V10 completion:50 gradients,0 accepted, screen rejectionpLDDT.31/i_pTM.20.
Closest requested-patch distance38.806A was again the initial round; final48.779A.
Actual GPU release and CPU no-acceptance follow-up verified. Runtime upper
$0.5815800, Modal cumulative$5.2184357. Next single targeting-heuristic rationale
and explicit surrogate/native-restoration limits are in
`FORCED_TARGETING_DECISION.md`; all57 local regressions pass. CPU image build
must finish before any new allocation. No early surrogate result may qualify.
# Latest completed checkpoint: 2026-10-22 07:45 EDT

All three morning tests completed: v9/v10/v11 each50 gradients, each0 accepted,
all pLDDT screen rejection. V11 .36pLDDT/.13iPTM; forced surrogate sequence375
changes in every observed call. Its1.255A patch CA minimum is overlap, not
binding, and no native restoration occurred. See the run's RESULTS.json.
Actual resource release verified07:45:29EDT (zero runners/inputs/backlog).
Modal cumulative runtime upper estimate$5.8023551; newly supervised Runpod
$1.7952 plus historical$33.87; consultations$1.1922999. Other billing remains
unreconciled. No overnight execution is claimed. No automatic next GPU job.
