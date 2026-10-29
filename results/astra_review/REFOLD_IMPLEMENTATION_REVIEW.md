# Hard-sequence refold implementation review — 2026-10-22

**Verdict:** no launch-blocking defect found in the reviewed two-sequence diagnostic. This is a bounded test of whether the rejected sequences predict plausible isolated folds, followed conditionally by native-monomer-target complex predictions. It is not a candidate qualification or a validated binding experiment. The current image must reproduce its CPU preflight before allocation; the earlier image preflight does not validate subsequently edited code.

## Reviewed implementation and checks

Reviewed `scripts/hard_sequence_refold.py`, `prepare_hard_refold.py`, `run_hard_refold.py`, `audit_hard_refolds.py`, the `sequence_refold` Modal function, and relevant pinned upstream feature/prediction/ensemble/export source. Upstream revision: `b6fb3ae15ab6a702f1156cffc11bcc3323c2fd06` (with the separately audited frame guard).

The reviewed worker SHA256 is `2db197cb99a918571ca01d2bdf6f071a2ec7994adae5dee68eda456f94e2ea3f`; preparer SHA256 is `19e14fcb49b5f18c6006b6466b753d529cb2c42439ed9f723324319f62d09101`. Later edits require a fresh hash/preflight record.

- The manifest is restricted to exactly the predeclared v9 and seed23 sequence hashes. Sequence contents and lengths are verified before inference. Fresh hard-one-hot proteins have SEQUENCE-only flags, zero coordinates and false atom masks. No discarded designed binder coordinates enter prediction. The native receptor is separately templated, with checked file identity, sequence, length639 and author residue indices122–760.
- `large_complex_initialization=False`, dropout disabled and three recycles match the protocol. The upstream fixed-sequence feature path does not reinterpret these inputs as optimization logits. The empty-template alignment path has an identity fallback; I found no recurrence of the unresolved-frame degeneracy in this route.
- Both declared models run for both monomers before complexes. No stage filters are supplied to `predict_validation_ensemble`, avoiding its early-return branch. The scheduled-job cap is checked before ensemble calls. Eight prediction jobs mean up to32 internal network evaluations with three recycles.
- Binder-first ordering is asserted before extracting binder confidence. Exported monomers use chainA; complexes use targetA/binderB with reordered confidence. Native hotspot indices are checked. Fresh internal-fold RMSD is distinctly named; it does not redefine the historical `Binder_RMSD`. Upstream aligns complex outputs to the target template, supporting the additional template-frame binder RMSD diagnostic.
- Advancement requires both monomer models to pass confidence and geometry plus their mutual-fold agreement. The complex phase has no automatic admission. CPU restoration is invoked only after verified GPU release; its entry point also checks that release record.

I independently ran the five focused `test_hard_refold.py` tests: all passed. Native peptide geometry passes, collapsed/missing/nonfinite or overlapping backbones fail, and rigid-motion removal does not hide internal deformation. I did not run the network, allocate a GPU, or modify implementation files.

## Interpretation and residual limits

The two AF2 checkpoint families are held out from the design models, but remain AF2. Agreement is a useful check, not independent-architecture validation. A monomer failure is a resource-triage decision for this strategy, not proof that a sequence cannot bind or fold upon binding.

The geometry screen catches the previously observed gross collapse and several backbone defects. It is not a complete stereochemical, side-chain, energetic or experimental validation. Confidence and cross-checkpoint agreement can both be confidently wrong. An apparently plausible monomer says nothing by itself about patch acquisition.

Complex inference contains the complete native639-residue target monomer, not the dimer. Later rigid restoration of chainB and checks against observed NAG atoms detect some exclusions; they cannot reproduce dimer-aware prediction, complete glycan shielding or assay biology. Existing native-context and submission criteria remain in force, without automatic pH work or candidate promotion.

One non-admission reporting gap was sent to root: the initially reviewed `audit_hard_refolds.py` marks an empty file glob as completed even if complex exports were expected. Require expected-versus-observed export/audit counts before calling that follow-up complete. Its fixed `qualified_candidates=0` prevents this omission from admitting a candidate. A legitimately empty audit is expected when neither monomer advances.

No new GPU spending from this review. Prior estimates supplied by root remain Modal $9.2964416754, Runpod $35.6652 including historical use, consultation $1.1922999; other costs unreconciled.
