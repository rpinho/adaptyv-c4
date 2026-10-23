# Conditional antibody-scale evidence: CPU continuation plan

Registered before any reference CDR identity result. The sequence-annotation
audit found59 H-domain annotations;50 have both IMGT104/118 flanks. These50 can
support a complete annotated-CDRH3 comparison. Keep the other350 inputs explicit
as unestablished for this particular analysis, including9 incomplete H domains,
26 light annotations and315 unassigned. This is an annotation-quality condition,
not selection by favorable novelty score. No official structural typing is inferred.

Use the already completed MMseqs PLAbDab hits for these same full-sequence
queries, with their existing input/hash/search settings. Resolve every returned
reference identifier to its exact sequence in the fingerprinted public FASTA.
Then recover **every** associated heavy sequence from the original paired CSV,
not just the representative pair retained by sequence deduplication. Include
matching unpaired H sequences from the original unpaired CSV. Keep counts of
unpaired light hits without a heavy partner; do not invent partners.

Annotation references will use the same pinned ANARCII2.0.8 CPU/accuracy model,
one numerical thread and seed0. Stream unique heavy sequences in finite chunks,
save raw annotations/checkpoints, verify numbered source spans and IMGT labels,
and require both CDR3 flanks for a complete reference loop. Noncanonical,
unassigned, partial or ambiguous references are counted as gaps, never repaired
or silently discarded from a claim of complete reference coverage. The source
unpaired `numbered` field was inspected but is not automatically assumed to be
the same coordinate/numbering convention as the current model.

Compute exact edit distances on the existing annotated loops. Primary local
identity convention:1-Levenshtein/max(query_length,reference_length). Report
sensitivity to query-length and reference-length denominators (clamped at0),
and normalized indel/LCS similarity2*LCS/(query_length+reference_length).
Each convention independently takes its best reference. The organizer's exact
edit-distance normalization is unspecified, so none is labeled the official one.

Use named antibody-scale boundaries: CDR0.70 is familiar, CDR0.95 is known;
global0.70/0.95 are on the familiar/known side. Apply the published antibody
truth table only **conditional on that scoring regime**, with the already
measured global sequence identity. Never choose whichever general-vs-antibody
policy gives a better score. If reference annotations are incomplete, allow
unseen references to increase similarity; withhold a numeric point unless every
completion gives the same level. Empty reference sets have null identity.

Controls: an independently identified publisher heavy-chain loop compared to
itself must give identity1 under every convention; the existing input must remain
unchanged. Unequal-length toy-string and exact-boundary checks will verify
normalization/threshold behavior. These checks are software controls, not designs.

Deliver per-query measured identities, closest reference evidence, reference
coverage/gaps, conditional levels/intervals, and unresolved typing/normalization
scope. No generation, sequence mutation, contact optimization, ranking, export
of submissions, upload, account action or GPU. No binding or assay conclusion.
New GPU/cloud compute:$0.00. Recorded Runpod$56.589988 +own Modal$10.718947
=~$67.31; separate generation/account billing excluded.
