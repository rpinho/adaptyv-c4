# Completed frame-guard integration diagnostic

One actual, unoptimized AlphaFold forward and sequence-gradient call completed
on CUDA with the repaired initialization path. This is software integration
evidence, not a usable binder or a completed optimization trajectory.

Raw artifacts: `../modal_frame_diagnostic_v8/mirrored/diagnostic.json`,
`first_forward.cif`, worker logs and `../modal_frame_diagnostic_v8/supervisor.json`.
Call `fc-01M5GG1MCW8TB2BP8RN6B9B6FP`, seed 20261021, trajectory 1, 69 residues,
model `model_5_multimer_v3`. Input is the complete native 639-residue target in
the audited local coordinate frame. Loss, epitope and acceptance gates were not
changed to obtain this diagnostic.

| Measurement | Result |
|---|---:|
| Loss | 15.7494259 |
| All sequence gradients finite | yes |
| Binder sequence-gradient norm, float64 accumulation | 6.2586025 |
| Binder CA positions / unique positions | 69 / 69 |
| CA span, xyz (A) | 24.219 / 19.584 / 58.246 |
| Adjacent CA minimum / median / maximum (A) | 1.562 / 3.628 / 4.173 |
| Minimum requested-hotspot to binder CA distance (A) | 38.541 |

The exact collapse is absent and the real sequence gradient is finite/nonzero.
The geometry does **not** qualify as a physical backbone: even the adjacent-CA
minimum is abnormal. The patch remains far away. No affinity, pH release,
mouse cross-reactivity, or experimental validation follows from these values.

The guard-off control was not run in this diagnostic. Comparing this initial
forward with the old v5 50-step endpoint is not a matched causal experiment.
The source/CPU evidence establishes the singular-frame mechanism; the GPU
result establishes that the repaired path runs without that output collapse.
V3's separate off-patch behavior remains unresolved.

Remote elapsed time was 102.493 seconds, return code 0. The single-use GPU
resource was released, independently verified with zero runners, zero inputs,
and zero backlog at release time. No successor design campaign is queued.
Requested H100 runtime upper rate: $4.649328/hour including CPU/RAM. Estimated
diagnostic runtime: **$0.1323675**; cumulative Modal runtime **$4.0531838**.
Earlier newly supervised Runpod runtime **$1.7952**. These are estimates;
setup/storage, older campaigns, and consultation usage are separate.

The old adaptive controller remains cancelled. Native-context error handling,
unique-trajectory accounting and physical patch-acquisition routing are still
required before its search can be treated as a qualifying-candidate pipeline.
