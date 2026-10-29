# RFdiffusion pilot implementation review — 2026-10-22

**Initial verdict: fix/verify the following preparation and configuration issues before launch.** No deployment, installation or GPU allocation was performed by this reviewer. The issues were sent to root immediately rather than waiting for this report.

## Launch/preflight issues

1. **Confirmed hotspot indexing error in the CPU assertion.** `sampler.hotspot_0idx` is returned by `get_idx0_hotspots` from `receptor_con_hal_idx0`, which is receptor-local. For this exact contig its values are16,17,19,20,21. The expected full concatenated feature indices are85,86,88,89,90 because69 binder residues precede the195-residue target. The preflight currently compares the first representation with the second and will fail. Assert the receptor-local list separately; retain the global-index assertion on the actual `_preprocess` hotspot feature. Do not modify upstream indexing to make them equal.

   I independently executed the pinned official `ContigMap` on CPU, with the actual contiguous target IDs and declared contig, without installing RFdiffusion. It produced264 residues and exactly the two index lists above. Source SHA256: `5df8d63a3b6efa3bd309585982864cd2af8cf4d1a6f681636be6da35ab77f1c0`.

2. **The50-step envelope is checked before effective configuration is known.** `diffusion_worker.py` checks the loaded base config before calling upstream inference. Pinned `Sampler.assemble_config_from_chk` overwrites model/diffuser/preprocess values from the checkpoint; it reapplies only Hydra task overrides. Programmatic `inference(conf)` does not establish those overrides. The CPU preflight records the resolved sampler config but initially does not assert effective `T`, `t_step_input`, absence of partial diffusion, or final step. Enforce the declared50-step trajectory on the post-checkpoint configuration before any network forward, and save that effective configuration. I have not established that this checkpoint actually changesT; the defect is the missing guarantee and potentially misleading pre-resolution record.

3. **Exact seeded initialization needs the same seed timing.** The CPU preflight seeds before constructing the sampler. Official inference constructs the sampler and then reseeds with the design index immediately before `sample_init`. Model construction can consume RNG state. Reseed the CPU check immediately before `sample_init` to match that operation order. Masks/identity assertions are still useful without this change, but it should not be described as auditing the exact initialized random sample.

Source basis: official pinned [sampler](https://github.com/RosettaCommons/RFdiffusion/blob/86507b6538f51fce57b5a72477165f03999ed7ae/rfdiffusion/inference/model_runners.py), [hotspot helper](https://github.com/RosettaCommons/RFdiffusion/blob/86507b6538f51fce57b5a72477165f03999ed7ae/rfdiffusion/inference/utils.py), [contig mapper](https://github.com/RosettaCommons/RFdiffusion/blob/86507b6538f51fce57b5a72477165f03999ed7ae/rfdiffusion/contigs.py), and [inference entry point](https://github.com/RosettaCommons/RFdiffusion/blob/86507b6538f51fce57b5a72477165f03999ed7ae/scripts/run_inference.py).

## Sound parts of the reviewed path

- Separate digest-pinned official-derived image, pinned source revision, checkpoint checksum verification and CPU construction precede GPU use. The actual image still needs to pass its build; this review does not establish dependency compatibility.
- The target file is hash-checked. The preflight inspects fixed target sequence/structure masks, generated binder masks, actual target sequence, model hotspot feature, and actual writer output. TargetA author189–383 and binderB1–69 are checked directly, without substituting legacy TRB chain conventions.
- The worker refuses CPU fallback and non-A100 GPUs. The remote function requests A100-80GB only, accepts one fixed run identity with no overrides, requires a fresh output directory, and enforces a finite wall-time stop. Its finalizer kills a surviving worker. The local controller uses actual release verification and checks four deployed functions before a fresh allocation. These bounds do not guarantee release timing during a cloud/API outage.
- Exactly one backbone output is required and hashed. The driver mirrors `.trb` metadata as well as the PDB and logs. The post-release CPU audit checks source identity, one model, exact chain/residue layout, poly-Gly binder, backbone geometry, declared patch contact and native full-context restoration. Only successful staging permits later MPNN; no sequence generation or successor is launched by this code.
- All output qualification fields remain zero/false. Poly-Gly geometry cannot be treated as a designed binding sequence. Patch contact is a pilot triage rule, not evidence of affinity. Target masks being fixed during preparation do not replace the final restoration audit of actual predicted target geometry.

## Reviewed file hashes

- `scripts/modal_diffusion_image.py`: `841ea5656e40fad7023f16b2287d78db727c84df1edebb877e081337b93a6c3b`
- `scripts/rfdiffusion_preflight.py`: `4c722310299325b64bd834e185fc434cb5165fd6e9c9bedfbcacf835c0693f77`
- `scripts/diffusion_worker.py`: `fc553c5a2747a0b331c2b5a45ed8934b4157c5d04f9a5727c415e5273b3740a8`
- `scripts/run_diffusion_probe.py`: `fc6f1b38c7e6faec07c26606cce8cf001f307f9e0d4944c7ecc42618a0633992`
- `scripts/audit_diffusion_backbone.py`: `87012a1ce01effe9016670c6fb16bf0936db81bf33e67966a832df5bede46cd0`

These identify the inspected working files; corrections require fresh CPU preflight evidence. Reported completed Modal estimate$9.9065633, other spending as supplied by root; no new GPU spending from this review.

## Focused correction/custom-loop review

The revised code fixes the original receptor-local/global hotspot comparison and reseeds immediately before both CPU and GPU `sample_init`. It now inspects the actual checkpoint-resolved sampler limits before any neural forward and writes the resolved configuration. These address the original three findings at source level, subject to the new CPU image actually passing.

**New blocker sent to root:** the revised Hydra task list includes `diffuser.partial_T=null`. Pinned `assemble_config_from_chk` replays each diffuser override by calling the current value's type on its string. With a null current value this becomes `NoneType('null')`, which raises `TypeError`. Leave partial diffusion at its existing defaultNone, omit that null override, and keep the post-resolution `is None` assertion. The real `diffuser.T=50` override is still required. This is a parser incompatibility, not a reason to alter the scientific protocol.

The custom numerical loop was compared against official `scripts/run_inference.py` and `SelfConditioning.sample_step` at the pinned revision:

- With SelfConditioning and symmetry disabled, each returned sequence is an unchanged clone. Keeping the initial sequence through all50 calls is equivalent for this declared mode. This equivalence must not be generalized to a different sampler/symmetry mode.
- The range50..1 with final_step1 matches upstream. The final returned `xyz` equals the latest state selected by upstream's reversed trajectory stack; writing its first four atoms is equivalent to the official backbone output.
- Self-conditioning state is updated inside `sample_step`; the wrapper does not omit or replace it. Sampling/noising is still handled by the pinned sampler. The extra hotspot-feature preparation is deterministic, with no new stochastic draw identified.
- Metadata differs: the custom loop writes step JSON and an explicit tensor/input/output map instead of an upstreamTRB/full trajectory. This is acceptable for the declared bounded pilot but must be described accurately. The PDB writer's default B factors are zero for all residues; they are not pLDDT. Logged mean pLDDT spans the whole modeled system and is not a binder-only confidence or interface score.

The removed MD5 assumption was incorrect: an opaque directory name in an official download URL is not a published checksum assertion. Downloading from the official HTTPS endpoint, checking expected size, loading the checkpoint onCPU and freezing its observed SHA256 gives a reproducible runtime identity. It should be described as the downloaded artifact's measured hash, not as verification against an independently published official digest. I found no reason to block this bounded pilot solely because no separate published SHA was available.

Reviewed corrected worker SHA256: `fc553c5a2747a0b331c2b5a45ed8934b4157c5d04f9a5727c415e5273b3740a8`.
Reviewed corrected preflight SHA256: `84bda1f30be41bc0d9526b658050cff1ae4bc5d764dd0722b26630cffbb951e2`.

No GPU allocation or implementation edits by the reviewer. The null-override fix and successful fresh CPU build remain outstanding at this addendum's writing.
