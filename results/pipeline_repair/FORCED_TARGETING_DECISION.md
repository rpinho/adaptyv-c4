# Conditional targeting-heuristic test after v9/v10 acquisition failures

V9 and v10 each completed50 gradient calls with no exact collapse and no
accepted design. Closest requested-hotspot CA distance was the initial round:
38.512A for v9,38.806A for v10. V10 ended48.779A away, rejected at screen for
pLDDT0.31(i_pTM0.20). Narrowing the distogram contact cutoff did not solve
physical acquisition in this one-seed test. This does not establish a universal
hit-rate comparison or a reason to lower acceptance gates.

Actual v10 resource release and CPU follow-up were verified before preparing
another allocation. Conservative v10 submission-through-release runtime
estimate$0.5815800; cumulative Modal runtime estimate **$5.2184357**.
Newly supervised Runpod estimate$1.7952 and consultation reported usage$1.1922999
remain separate; setup/storage/older campaigns unreconciled.

## One different, declared hypothesis

Test `forced_targeting=true` as the **only scientific configuration difference
from repaired-native v9**, with v9's original contact cutoff20A, original losses,
weights, seed/first length, complete local-frame native input structure and all
acceptance gates. One trajectory,30min worker cap /35min monitor cap, actual
resource release before CPU follow-up, no automatic successor. H100 runtime
upper$4.649328/hour including requested CPU/RAM, same existing authorization.

The hypothesis is that upstream's targeting heuristic can focus optimization
on the intended patch where native-sequence search has not acquired it. This is
not a proven fix or evidence of superior yield. No new COM or coordinate loss
is introduced in this test. The existing COM objective targets the whole
receptor centroid and is unsuitable as a local-patch substitute.

## Explicit representation and admission limits

This **is not native receptor chemistry throughout search**. The pinned
implementation replaces exposed target residues outside a10A hotspot shell
with lysines during early optimization while retaining the complete target
backbone. The prior exact source audit found371 sequence changes on this input.
Early screen/refine/anneal results and failures are surrogate-target predictions
and must never become native-context candidates or enter pKa engineering.

Astra traced native restoration to `harden.prepare`, before harden gradients and
filters. Target-schedule cached entries are replaced and predictions cleared;
later mutation polish, refolds and final filters use restored target states.
This source flow does not replace independent verification. Final output must
pass existing CPU native sequence/numbering, restored full-dimer backbone/NAG
clearance, canonical-file hash and scientific admission checks. Generator rows
alone are never qualifying candidates. No assay-target mutations are proposed.

CPU preflight now separately prepares native and forced states and records
their actual sequence differences from the native input. Per-gradient monitor
records the same difference count, making any eventual restoration observable.
The input PDB, native reference and full-backbone context are preserved. The
source/monitor hashes and initial sequence-array hash are recorded.

Preflight correction before allocation: initialization alone precedes the
upstream forcing operation, so the first preflight incorrectly reported0
changes in both modes. It was corrected to execute the actual forcing callback
before auditing frames/identity and to require the declared representation.
The corrected CPU build passed: native0 changes; the current local-frame forced
input **375** changes, complete639-residue backbone observed in both. The older
371 count is historical for a different input representation, not this measured
prepared state. No GPU was allocated until this corrected check passed.

All57 local regressions pass, including exact forced-vs-native configuration
delta and rejection of lysinated inputs by context admission. Details and source
lines: `../astra_review/ACQUISITION_REVIEW.md`.

Current launch/result evidence, if allocated after CPU build succeeds, belongs
in `../modal_forced_native_restore_v11/supervisor.json`. This decision alone is
not evidence of a running job, native-target restoration or success.
