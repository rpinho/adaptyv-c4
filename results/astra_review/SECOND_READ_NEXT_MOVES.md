# Second read of next moves — 2026-10-22

## Independent recommendation

**Change the original ordering. A small sequence-only refold diagnostic should precede a therapeutic structural control and another generation campaign.** Stop additional same-recipe native seed repeats. Keep the contiguous apical-domain route conditional and time-boxed; do not turn its validation infrastructure into another multi-day prerequisite before testing a useful hypothesis.

The original plan correctly separated eligibility from internal metrics and retired the invalid crop. Its weakness was making a mismatched positive control the first gate. A6W3H restoration success cannot tell us whether the minibinder generator/runtime works, and its227-residue Fc plus rearranged receptor fragment require special handling outside the current30–150-residue binder adapter. That work can consume time without improving candidate generation.

Root reports updating NEXT_MOVES.md to incorporate this revised priority during the consultation. This review explains the change rather than automatically endorsing the initial plan. No new GPU, account, submission or code mutation occurred. Recorded spend remains Modal approximately$9.30, Runpod$35.6652 including history, consultations$1.1922999; extras unreconciled.

## What the actual evidence supports

I read the completed v14 series ledger and six saved native-run CIFs: original repaired v9 plus its five independent-seed replicates. V14 completed five calls, zero generator acceptances, all scientific follow-ups cleanly reporting no accepted scaffold; cumulative recorded Modal runtime estimate is$9.2964416754. Actual resource release/fresh zero stats are reported by root; I did not query cloud state in this consultation.

There are nine recent failed attempts, **not nine independent identical trials**: four morning ablations reused the original initialization, while the five new seeds extend the native v9 recipe. Six distinct native-recipe seeds provide evidence that this setup is repeatedly unsuccessful. They do not establish target impossibility, a precise hit rate, or the cause of failure.

My saved-structure measurements:

| Run | Length | Minimum adjacent binder Cα distance | Minimum requested-patch Cα distance |
|---|---:|---:|---:|
| Native v9 |69|3.392Å|52.111Å|
| Replicate2, seed20261022 |69|0.814Å|60.607Å|
| Replicate3, seed20261023 |67|2.647Å|55.228Å|
| Replicate4, seed20261024 |68|1.751Å|60.171Å|
| Replicate5, seed20261025 |74|1.164Å|56.479Å|
| Replicate6, seed20261026 |75|1.265Å|45.880Å|

These are the saved selected predictions, not necessarily final-round structures. Even v9's less-compressed binder has consecutive C–N distances0.952–1.979Å; replicate3 spans0.806–2.058Å. **There is currently no plausible native interface being withheld solely by a strict confidence cutoff.** Lowering a screen threshold would not repair these off-patch poses. Conversely, a failed soft-sequence optimization prediction does not experimentally disprove the corresponding discrete sequence.

## First bounded experiment: two fixed sequences, monomer first

Predeclare these two inputs before any refolding result; do not replace a poor result with a newly cherry-picked sequence:

1. `results/modal_repaired_native_v9/mirrored/campaign/1_Trajectories/c4_repaired_native_v9_denovo_l69_487e943cff8955d6/c4_repaired_native_v9_denovo_l69_487e943cff8955d6_trajectory.cif`
   - Artifact SHA256: `8311b87eff4d30aa7c4323d4ed4051cf4ef7d50d6f631af9a3d609fd85949a5d`.
   - BinderB sequence,69 residues, SHA256 of uppercase one-letter sequence without newline: `2699888c00c9c6af76b61c8547687cce87ccb0152754ae0d0a32c10e5f42b85d`.
2. `results/modal_native_replicates_v14/replicate_3/mirrored/campaign/1_Trajectories/c4_repaired_native_v9_replicate_3_denovo_l67_a207b66be0c44216/c4_repaired_native_v9_replicate_3_denovo_l67_a207b66be0c44216_trajectory.cif`
   - Artifact SHA256: `a5e467a235cfa2237c36740e0d036c6bc8d6985332a72fe0ea4ddb2de5beaf69`.
   - BinderB sequence,67 residues, SHA256: `8b0e2189b54820dd47b3d25e26969dd1a519f6c9a7267aee9d6f50a67fefa2b3`.

V9 is the less-geometrically-distorted original native reject; replicate3 has the least adjacent-Cα compression among the five saved independent-seed outputs. This selection is for a diagnostic, not an assertion of likely binding. Extract the discrete sequence exactly once and record its provenance. Do not redesign it, inject its saved backbone, initialize it from the saved coordinates, or provide the designed binder as a template.

**Bound:** two sequences × the two existing validation checkpoints for binder-alone folding; only plausible monomers proceed to at most two native-target complex predictions each. At most eight prediction jobs, one explicitly capped allocation, no automatic sequence redesign or successor. A modest wall-time cap may leave this partial; preserve that result rather than expanding automatically. The prediction adapter should reuse the pinned validation interface, with CPU checks for sequence/flags/template absence, rather than build a new general framework. The source already exposes validation ensembles through a predictor interface. [Pinned validation source](https://raw.githubusercontent.com/PacesaLab/BindCraft2/b6fb3ae15ab6a702f1156cffc11bcc3323c2fd06/bindcraft/MPNN_stage.py).

Actual archived metadata lists all five `model_*_multimer_v3` checkpoints for design and `model_1_ptm`/`model_2_ptm` for validation. These are **held-out AF2 checkpoint families within the same model lineage**, not genuinely different-architecture confirmation. Use that terminology. A later independent architecture would be stronger additional evidence; installing one should not block this small diagnostic.

Monomer advancement should require physically credible peptide/backbone geometry, absence of severe self-overlap, and confidence consistent with the existing0.70 unbound-binder criterion across the declared predictions. Compare the two fresh folds for gross disagreement. Predeclare the exact deterministic checks before allocation; do not select whichever checkpoint looks best. Borderline/disagreeing folds are inconclusive and do not automatically advance. This is resource triage for a stable-minibinder design strategy: failure does not rule out induced folding or all biological binding, but provides little reason to spend more on these sequences now.

Complex prediction must use native target chemistry and no designed-binder coordinates. Prefer full native receptor context where feasible. If memory requires the already-declared complete639-residue monomer approximation, retain that explicit limitation and restore/check the second11QC subunit under unchanged full-context admission afterward; do not call that a dimer-aware refold. A good free fold alone is not evidence of TfR1 affinity. A newly plausible complex could justify further review even though the original trajectory failed screening, without rewriting its historical status or weakening native-context checks.

## Decisions after this diagnostic

- **Monomers fold plausibly; complexes fail:** the next bottleneck is acquisition/interface generation. One bounded native contiguous195-residue apical-domain search becomes the next useful generation experiment, followed by native-context refolding/restoration. The refold adapter just exercised is directly reusable. No domain-score promotion.
- **Both monomers fail:** do not run a large rescue/refold sweep. One small domain-generation test is still defensible because context/optimization changes, but more weight or length tweaking has weak support. In parallel planning, evaluate one documented alternative generator/reference runtime as a real fallback, not another speculative installation project.
- **A native complex becomes plausible:** stop expanding methods and spend the next effort on its robustness, novelty and submission evidence, then pH/cross-species variants.

Changing binder length alone is not the best next test:65–75 already sampled and the failure is not demonstrated to be a size deficit. Do not simultaneously change target representation, length, weights and runtime. Arbitrary dependency upgrades would destroy attribution. A different generator could ultimately outperform this configuration, but access, installation, target support and validation costs are not yet established. The experimentally published original BindCraft is a reasonable reference implementation to assess; its published successes do not establish TfR1 performance or validate this modified BindCraft2 setup. [Primary BindCraft paper](https://www.nature.com/articles/s41586-025-09429-6).

6W3H remains an optional, short CPU geometry-calibration task: can mapping/restoration correctly handle an experimental therapeutic complex? Its outcome validates only that question. Do not force it through minibinder architecture limits, use its sequence as a candidate seed, or claim that a passed control validates generation. A generation control would require a separate published target/design benchmark with the relevant runtime and preset.

## Competitive submission and milestones

The official C4 page allows de novo formats, up to20 entries for Tracks2/3, and does not prescribe our pLDDT/PROPKA thresholds. The winner must satisfy human affinity, pH selectivity and mouse recognition; an early human-only prospect is incomplete evidence, not a competitive winner yet. The recommended epitope is optional. The current apical patch has substantial human/mouse differences, so perform the cheap conservation/epitope check now even while deferring expensive mouse/pH modeling. Do not discover this incompatibility only after optimizing a human-specific interface. [Official challenge](https://proteinbase.com/competitions/anthropic-adaptyv-2026/challenges/tfr1).

Selection considers multiple metrics, novelty and diversity; no single confidence metric decides. This supports ranking credible incomplete designs honestly, not submitting physically contradictory models or padding the list to20. [Official selection description](https://proteinbase.com/competitions/anthropic-adaptyv-2026/selection).

Early upload reduces logistics risk; no verified early-timestamp scoring advantage or daily experimental feedback exists. Public FAQ estimates validation after the window, with first data in December; official terms give a broader December/January schedule, subject to changes. Replacement/editing rights remain **unverified**, including after reading the public terms. Do not spend scarce slots assuming later uploads replace earlier ones. Authenticated submission administration remains Ricardo's task. [Official terms, sections5.3–5.6](https://proteinbase.com/competitions/anthropic-adaptyv-2026/terms).

October24 is credible as an **evidence/package decision review**, not a promise of eligible strong candidates. October31 is a sensible internal finalization target ahead of C4's November1 cutoff. Set nearer decision milestones: finish the small refold diagnostic today if implementation is ready; decide the next generation lane from that result; by October24 either present a defensible shortlist or explicitly record zero recommendations and the remaining bounded experiment. Do not let an empty shortlist become an excuse for endless pipeline work, and do not call weak sequences credible to satisfy a date.

`SUBMISSION_STATUS.md` still says the seed queue is active and quotes Modal$6.50; those end paragraphs are stale relative to the completed ledger. Correct that status/spend, keeping historical claims visibly dated. The underlying distinction between software checks, computational candidates, submission recommendations and experiments is sound.
