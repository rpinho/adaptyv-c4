# Domain adapter review — 2026-10-22

Scope: CPU/source review of `scripts/restore_apical_context.py` SHA256 `d1838b93b2312402e36e952c3eec77e4c96551b0eeb38df885c4ce7282518279` and `tests/test_apical_restoration.py` SHA256 `c9f26da6e8199491e449214044e86ea0bf5530ff019c5357cbb0624a03cb32c3`. No root code, GPU state or data file changed. Additional adversarial checks used in-memory copies of Bio.PDB models, an in-memory PDB writer and mocked file lookup/hash boundaries; they are logic tests, not a claim to have rerun the parent's complete on-disk test suite.

## Verdict

The **rigid-restoration calculation and chain layout are correct for the present native template and ordinary standard-residue inputs**. Its staging-only labels are appropriate. It is not a native complex refold, admission or binding result. Two input-integrity gaps should be fixed before connecting this adapter to a future search/refold workflow.

## Findings

### 1. Native reference is not independently verified — moderate priority

`native_model()` (lines25–31) verifies residue numbers122–760 only, despite its error message referring to identity. It does not check either chain against the expected native sequence, full backbone completeness, finite coordinates or a pinned reference hash. `audit_domain()` then treats this same file as the sequence authority.

Reproduction: with an in-memory native template carrying A142 GLY→ALA, an otherwise valid domain/binder fixture restores successfully and writes A142 ALA. The audit still declares native sequence verification. The present real file is unchanged and native; this exposes a fail-open boundary, not evidence that today's results are mutated.

Missing finite-coordinate validation also matters: a NaN in omitted target backbone or observed NAG coordinates can yield NaN clearance; `NaN < 2` is false, so the current comparison does not reliably reject it. Require finite reference coordinates and finite computed clearances, not merely nonempty arrays.

Fix before integration: verify A/B sequence independently against the canonical UniProt122–760 reference; assert complete native backbone, coverage/insertion-code policy and finite coordinates. Pin or verify the approved11QC template digest as provenance. Retain downstream native admission independently.

### 2. Input filtering can silently change the binder — moderate priority

`protein_residues()` and `ProteinOnly` (lines21–22,60–62) exclude unrecognized residues. The binder's size/geometry checks operate only on that filtered list, and the writer silently drops the rest.

Reproduction: append one UNK residue to an otherwise accepted30-residue distant binder. Restoration accepts31 input residues and writes30 residues. The source/output hash difference records bytes but does not explain that the molecule changed. This violates an atom/sequence-preserving staging contract.

Reject unsupported protein residue types and ambiguous alternate conformers explicitly, or implement a declared validated normalization policy. Assert that every input binder residue/atom intended for preservation survives writing, with the same identities/order and only the intended rigid coordinate transform. Use strict PDB parsing and reject duplicate/ambiguous identifiers rather than suppressing warnings from a permissive parser. Generated standard-residue BindCraft outputs are not known to exercise this defect, but the boundary should fail closed.

### 3. Binder geometry and glycans need precise staging language — limitation

The adapter checks binder backbone completeness and finiteness, but not consecutive peptide geometry, local Cα spacing or self-overlap. An in-memory fixture with a binder N atom displaced20Å passes restoration when distant from the receptor. That is not an admission bypass because `qualified_candidate=False` and native refolding remain required. Given the observed malformed binders, record a binder-geometry audit and preferably reject broken peptide chains before spending on refolding; do not describe a successful staging call as physical validation.

All six observed NAG residues in current11QC are correctly included in the clearance check in the same native frame as the transformed binder. **They are then excluded from the written PDB by `ProteinOnly`.** State “observed NAG clearance checked; glycans omitted from output” explicitly. Current wording mentions only unresolved glycans. Future refold/final admission must repeat full-context/glycan checks on the new pose; this pre-refold clearance is not transferable evidence.

## Verified strengths and remaining tests

- Domain author IDs189–383, sequence, insertion codes, N/CA/C/O, finite coordinates and peptide connectivity are checked. Existing `cif2pdb.py` preserves `auth_seq_id`/`auth_asym_id` when present; actual search output still needs that round-trip test before integration. Chain layout is deliberately targetA/binderB in and nativeA/B/binderC out.
- The row-vector Kabsch fit uses all195 Cα atoms, correct handedness correction and translation. The same transform is applied to the binder while native receptor A/B remain fixed. Both overall-domain and patch RMSD plus maximum displacement are checked.
- With the actual template, an in-memory PDB write/read restored all10,114 native A/B protein heavy-atom coordinates exactly. Current11QC has no alternate locations or nonfinite coordinates. Replacement of the predicted domain by the complete native receptor preserves native context without splicing a deformed backbone.
- Existing tests meaningfully cover rigid transformation, omitted-protein overlap, mutation, missing domain residue, patch deformation and preparation coordinate preservation. Add regressions for reference mutation/nonfinite coordinates, unknown binder residues, explicit observed-NAG clash, and output identity preservation for every binder atom and both target chains. The current single binder-Cα assertion alone does not prove whole-binder preservation.
- The positive fixture places a copied protein fragment500Å away. It tests restoration mathematics and clearance rejection, **not** a realistic bound pose, patch acquisition, free-binder folding or native refold survival. Keep that distinction in CPU artifacts.

No scientific reason to relax existing gates emerged. Domain preparation/restoration remains a CPU alternative pending these integrity checks and an explicit native full-context refold path; it is not ready for automatic candidate admission or pH scoring.
