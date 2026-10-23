# Measured indexing correction and PDB formatting repair

The original pinned TED consensus code fails a basic identity invariant: three
identical one-based `1-100` domains produce `2-100`, and three identical `25-74`
domains produce `26-75`. A discontinuous `1-40_61-100` also shifts and loses its
terminal boundary. These are synthetic software fixtures, not candidate results.

The cause is the call to `domstr_to_assignment(c, nres, offset=1)` even though
that function already converts one-based input using the slice `start-1:end`.
The local compatibility copy changes exactly that call's offset from 1 to 0.
The pinned original stays unchanged. INDEXING_CHECKS.json records both source
hashes, the exact one-line difference, and negative/positive fixture outputs.
No IoU, domain-size, fragment-size or novelty threshold changes.

Six consensus fixtures now pass: full chain, interior interval, two-method
agreement, discontinuous domain, disjoint disagreement, and no assigned domains.
Five checks exercise the actual published domain filter: 24 versus 25 residues,
four versus five residue fragments, and Chainsaw's zero-to-one-based conversion.
This is an explicitly corrected local TED-style implementation; the organizer's
versions and treatment of this issue remain unknown. It is not evidence that
the organizer shares the bug or that any measured candidate level changed.

Separately, our first PDB writer emitted a bare `TER` footer. Bio.PDB and the two
neural predictors accepted it, but pdb-tools' residue selector reads residue
metadata from TER and failed before the UniDoc classifier ran. The writer now
emits a complete TER record for the actual last residue. All 68 original ATOM
blocks are byte-for-byte preserved (TER_ATOM_PRESERVATION.json); no sequence,
coordinate, model confidence or physical geometry was repaired. Earlier inputs,
failed control and logs remain available for diagnosis.

The coordinate-rounding check also now separates PDB text rounding from
Bio.PDB's float32 parser quantization rather than mixing both under one bound.
Both quantities are measured in INPUTS.json. This changes a representation
check, not a scientific threshold or coordinates. The CPU wrapper fixes Torch
and NumPy seeds to zero for local reproducibility and records CPU/threads/version
in each model log. These settings are local choices, not identified organizer
settings.

CPU preflight and complete-census results are separate from these software
fixtures; do not claim a working three-method census solely because fixtures
pass. New GPU/cloud CPU cost is $0.00.
