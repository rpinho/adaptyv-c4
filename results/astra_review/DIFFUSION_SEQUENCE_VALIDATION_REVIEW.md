# RFdiffusion → MPNN → protocol3 review — 2026-10-22

**Verdict: no launch-blocking source defect found in the reviewed conditional validation path.** This does not assert that an RF backbone passed staging, that two sequences exist, or that a protocol3 image is ready. Actual backbone staging, successful CPU MPNN extraction, frozen manifest and fresh manifest-specific image preflight remain prerequisites. No sampling, deployment or GPU allocation was performed by this reviewer.

## Focused findings

- Sequence preparation requires the successful backbone audit and its exact restored-output hash, plus pinned MPNN revision/source/checkpoint hashes. Existing sampling/output folders prevent silent resampling or replacement. It requests two samples from one backbone, with binderC designed and nativeA/B fixed.
- I inspected the locally pinned official ProteinMPNN runner. Its `save_probs` NPZ contains concatenated full sampled sequences, chain ordering and `mask_for_loss`. Thus the current NPZ checks correctly verify the actual design mask and every fixed target position against the restored input, rather than trusting the FASTA header alone. FASTA and sampled binder tensors must agree. Batch size1 and one temperature yield exactly two samples in this code path.
- Duplicate sequence draws are possible and correctly reported through `unique_sequences`; they do not become independent backbones. The manifest retains `independent_backbones=1` and two distinct sample names. No MPNN score is treated as binding evidence.
- Protocol3 requires exactly `rfd_mpnn_1` and `rfd_mpnn_2`, both69 residues, fixed model pair/recycles and an eight-job cap. It requires declared native-target verification and the v18 source provenance. Generic path-prefix/hash-shape checks alone cannot establish authentic provenance; the local extraction and source-hash checks supply that evidence.
- Coordinate-free preparation remains shared with the reviewed refold path: fresh one-hot binder sequence, zero binder atoms, false binder masks, SEQUENCE-only flags. The complex target comes from the separate verified639-residue native monomer, with no forced surrogate or supplied RF binder coordinates. Both monomers run before conditional complexes. The cap is enforced before prediction calls.
- The image route exists only after a manifest is present and its CPU preflight succeeds. Local launch requires matching manifest-specific CPU proof and unchanged source artifacts. The remote function checks image proof against both manifest and worker hashes, copies its actual input evidence, and retains finite supervision. The launcher checks all four relevant functions before a new allocation.
- Actual release precedes CPU audit. The expected complex-export set must match completed jobs, including the legitimate empty set when monomer triage rejects both. Output labels, chain arrangement and confidence slicing retain the reviewed conventions. Native restoration remains staging-only, with zero qualified candidates and no automatic pH/submission step.

Seven focused `test_hard_refold.py` tests passed. This is source/CPU review, not a test of the full-length MPNN run or neural refolding. The two AF2 validators differ from the RFdiffusion/RoseTTAFold generator but share AF2 architecture with each other. Complex validation still uses a639-residue native monomer with later CPU dimer restoration, not dimer-aware inference.

## Reviewed identities

- `scripts/prepare_diffusion_sequences.py`: `17795dd32ce94a899fb4e0b04c1f6ebe0bda33949d795df74ccc2557a03024ea`
- `scripts/hard_sequence_refold.py`: `0302e037be4b00144ff1b9491695162878548c89703d2e67084645408a950d46`
- `scripts/run_hard_refold.py`: `2ec2f4cd586500adaafafbe5e192b6335bb22de1eed7f6cd253f8dfdcdf1232f`
- `scripts/audit_hard_refolds.py`: `b597d463308fb02f922860d3a3563a7775ae688846b671f0a8f11f77f7987fc0`
- `scripts/modal_campaign.py`: `9ab344282f1bddfda3a41034b714a974cc63b4f000a9b79b09542ef9780cd1de`

Spending provenance supplied by root: latest posted Runpod billing$56.589988 across five pods, with known three-current-pod subset$51.9677576; do not add the old$1.7952 again. Modal runtime upper estimate approximately$10.5737 including failed starts; consultations$1.87865, other costs unreconciled. This review allocated no compute.
