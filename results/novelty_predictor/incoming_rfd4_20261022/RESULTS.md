# Existing RFD4 outputs: CPU novelty evidence

Source: HANDOFF_FROM_CLAUDE.md, Claude branch e8fdfa9. Both previously generated
200-row arm exports were retrieved through read-only Modal Volume operations.
No Modal app, generator, cloud worker or GPU was started. Source listings and
byte-for-byte exported CSVs are retained here; coordinates live in the local
cache recorded by DATA_PLAN.json and the supplied coordinate manifest.

## Fixed analysis policy, recorded before score inspection

Analyze all 400 rows, with no functional, histidine, metric or favorable-score
preselection. Match each row to its exact generated CIF filename, require exactly
one chain matching the complete input sequence, fingerprint the source, and
require complete backbone coordinates in the existing workflow. Do not assume
the designed-chain letter. Extra conditioning CIFs are not query structures.

Use scripts/run_cpu_novelty.py with the same pinned five public sequence sources,
historical PDB structural reference, three CPU domain methods, corrected consensus,
exact TM-score export and all-parent-context follow-up as the delivered workflow.
The structural stage execution budget is increased to 1200 seconds for this larger
batch. This changes a process time limit, not any scientific threshold or search
setting. Fixed named sequence/TM/coverage equality rules remain unchanged.

These are supplied generator coordinates, not independently predicted monomer
structures. Novelty evidence does not establish folding, receptor binding,
pH selectivity, assay success or submission eligibility. Every numerical level
is provisional; exact organizer levels established remain zero. Database scope,
historical snapshots, finite retrieval and organizer domain-policy gaps remain.
This lane does not generate, mutate, rank, emit entries or submit designs.

## Results

| Existing arm | Rows | Local 1 | Local 2 | Local 3 | Local 4 | Unestablished |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| tfr1_c200 | 200 | 40 | 40 | 0 | 0 | 120 |
| tfr1_his200 | 200 | 33 | 30 | 0 | 0 | 137 |

Actual command and source hashes: SOURCE.json. Final report: run_bounded/parent_context/predictions.json.
Per-design identities, alignment lengths, coverage and domain scores: EVIDENCE.csv.
Input order/sequence/coordinate fingerprints and all five raw sequence-search hashes were verified; JSON/CSV levels agree and the production local rules recompute exactly.
These are general-scale provisional predictions only. Antibody/CDR typing is unimplemented; actual antibody-database match evidence and its consequence are in ANTIBODY_SCOPE.md and ANTIBODY_SCOPE_EXAMPLE.json. No official acceptance/rejection, binding, pH-selectivity or assay result is established.
BOUNDED_POLICY.md applies: this is positive evidence with unfinished alignments, not exhaustive best structural similarity. Unresolved levels are withheld; no local3/4 is certified by missing matches.

## Spend

Recorded tracked lane: Runpod $56.589988 + own Modal estimate $10.718947 =
$67.308935 (~$67.31), unchanged. Claude's handoff separately estimates about $15
for its lane; that is not independently verified account billing. New GPU/cloud
compute spend here: $0.00. The Volume operations only read/download existing files.
Storage and any account-wide or other-lane activity are not measured by this audit.
