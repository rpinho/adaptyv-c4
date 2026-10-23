# Actual CPU annotation: 400 existing designs

Ran run_annotation.py with ANARCII2.0.8, antibody/accuracy, CPU, one numerical
thread, seed0 and shipped fingerprinted weights. No sequence generation,
mutation, folding, GPU or cloud worker. Inputs are all400 existing RFD4 rows;
the separate general-scale structural novelty run is still in progress.

| Existing arm | Inputs | H annotation | K annotation | L annotation | Unassigned |
| --- | ---: | ---: | ---: | ---: | ---: |
| tfr1_c200 | 200 | 34 | 10 | 5 | 151 |
| tfr1_his200 | 200 | 25 | 6 | 5 | 164 |
| Total | 400 | 59 | 16 | 10 | 315 |

All85 assigned annotations map exactly to their reported contiguous source
sequence span. All IMGT labels are unique and within positions1–128, retaining
insertion labels.50 of59 H annotations have numbered residues at both CDR3
flanks104/118;9 do not. A CDR-region substring is not automatically a complete
CDRH3. The unassigned315 are not established non-antibodies.

Actual model-call times, excluding environment setup/imports: controls0.224747s;
400-input batch3.772071s. Native input/output hashes, model/source versions and
weights are in RUNTIME.json and VALIDATION.json. Raw outputs remain in
control_raw.json and batch_raw.json; ANNOTATIONS.json holds per-query evidence.

Controls: publisher12E8 H→H (score30.764790), L→K (29.945995); native195-residue
receptor→unassigned (13.203511). The first adapter expected null instead of the
publisher's `F` failure sentinel. That failed assertion and raw controls are
preserved; SCHEMA_CORRECTION.md documents the representation-only fix. No model
threshold or scientific control condition was adjusted.

## What this does and does not establish

These are local antigen-receptor sequence annotations, not the organizer's
structural type calls. H does not distinguish VHH from VH, and a light-domain
annotation does not determine the submission molecule class. Exact organizer
types established:0. Antibody-specific novelty levels measured:0; no reference
CDRH3 comparison has yet run. No binding/assay success or qualified design is
established.

The organizer uses a different novelty rule for antibody fragments. Therefore
general-scale scores alone cannot determine official acceptance/rejection for
these annotated queries. [Published scoring description](https://www.adaptyvbio.com/blog/novelty).
Both the fixed general-scale evidence and this separate scope audit are retained;
neither policy is selected because it produces a more favorable score.

## Spend

Recorded Runpod $56.589988 + own Modal estimate $10.718947 =~$67.31, unchanged.
New GPU/cloud compute:$0.00. Separate generation-lane/account billing excluded.
