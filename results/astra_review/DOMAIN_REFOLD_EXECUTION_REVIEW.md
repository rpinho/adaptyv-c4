# Domain protocol2 execution review — 2026-10-22

**Verdict: no source-level blocker found for the declared one-sequence diagnostic.** This is a review of execution code, not of a generated sequence or completed protocol2 image. At inspection, HEAD was `cd3f34c` and `results/domain_hard_refold_v17` did not yet exist. Actual source extraction, a frozen manifest and its fresh CPU image preflight remain prerequisites. The active v16 run was not changed by this review.

## Checked

- Protocol1 retains its exact two-sequence whitelist. Protocol2 permits exactly one69-residue `domain_primary` entry from declared campaignv16, with self-consistent sequence SHA, source provenance prefix, two fixed AF2 checkpoints, three recycles and four-job cap. An invalid protocol, expanded sequence list or expanded cap is rejected. The cap is checked before each ensemble call; two monomer jobs precede up to two conditional complex jobs, at most16 internal evaluations.
- The real selected sequence is frozen into the deployed image via `add_local_file(copy=True)`, followed by CPU preflight. Without that image manifest, the remote route forv17 does not exist. No request overrides are accepted. Runtime copies the image manifest and CPU preflight into its output, preserving the actual input evidence. The preflight records worker and manifest hashes.
- `prepared()` reconstructs a fresh binder with hard one-hot sequence, SEQUENCE-only flags, zero atoms and false masks. The reduced-domain pose is not read or passed into inference. The complex target is the separate hash-checked639-residue native monomer, with verified sequence, author indices122–760, TEMPLATE flags and no DESIGN flags. Forced targeting is disabled and large-complex initialization remains false.
- The local launcher requires a prepared manifest, checks all three deployed functions before a fresh allocation, and uses the existing bounded/release-verified driver. Follow-up CPU audit runs only after completed inference and verified release. There is no pH task, redesign successor or automatic admission.
- Extraction remains deterministic: first ranked design, otherwise the sole saved trajectory. It now rejects multiple model identities, alternate conformers, insertions and conflicting residue names; requires native195-residue target author mapping and complete binder1–69. A rejected restored source pose remains explicitly rejected while its sequence can undergo the separate diagnostic.
- Output writers retain monomer chainA and complex targetA/binderB, with binder-only confidence and reordered export confidence. The audit checks the exact expected set of completed complex exports, including the valid zero-complex case when monomer triage fails. Canonical restoration is staging-only; every audit reports zero qualified candidates and no submission recommendation.
- The new checkpoint-scope disclosure correctly avoids calling these models held out when an accepted upstream domain design may already have been selected using them. Both are AF2, and the full-target prediction is monomer-based with later CPU dimer restoration, not dimer-aware inference.

Six focused hard-refold tests passed, including protocol2 list/campaign restrictions and geometry/sequence checks. I did not execute neural inference, deploy an image or allocate a GPU.

## Provenance boundary

Protocol2's generic validator cannot establish authentic campaign origin merely from a source-path prefix and a64-character source hash. That evidence comes from the local extractor and the actual immutable image manifest. Before launch, preserve/compare the extracted manifest SHA against the new image preflight and later mirrored manifest. Do not rerun extraction to replace the selected sequence after observing a refold result or use an older deployed image with a different local manifest. These are checks on the declared workflow, not a new modeling branch or an instruction to change acceptance criteria.

## Earlier domain preflight correction

The corrected code explicitly verifies actual HOTSPOT flags at205,206,208,209,210 and matches the domain input hash to `PREPARED_INPUT.json`. Its CIF roundtrip now follows the actual downstream path, `cif2pdb` using author IDs followed by the PDB reader. This correctly tests preservation of189–383 rather than expecting the upstream direct-CIF reader's label numbering1–195 to preserve author IDs. No global numbering behavior was changed.

## Reviewed file identities

- `scripts/hard_sequence_refold.py`: `5a759c8a9a032407d31bed46ddaafab52194c1ad8f296c7568916df0aa921714`
- `scripts/modal_campaign.py`: `a01d1b4b5e6c2a280f5a76e2954baa4ad4b1d4e6e8e63971bd915f7381d6b3fd`
- `scripts/run_hard_refold.py`: `5fb840b419647320f2ea9c1dd52257635966b7e031619b9ae5b289ac4fa4d0e9`
- `scripts/audit_hard_refolds.py`: `0fe26dca5ecdac68289b0e2b75f6758e7948f59c9e25131e4b31670dfb6c350e`
- `scripts/prepare_domain_refold.py`: `87bde160604fb8a129b357deb167561c9970ff499a4284b904aba0bf2ff86620`

No GPU spending from this review. Prior completed Modal estimate supplied by root is$9.5172782; activev16 accrual is separate. Runpod$35.6652 and consultation$1.1923 remain as reported, with extras unreconciled.
