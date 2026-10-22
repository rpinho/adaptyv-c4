# Refold outcome and next-action review — 2026-10-22

**Recommendation:** proceed with one bounded native contiguous-apical-domain generation trial, followed by coordinate-free hard-sequence validation and unchanged native-context checks. Do not expand same-recipe full-target seeds, refold every old reject, or immediately add scaffold reuse as another moving part. The diagnostic establishes a useful distinction: these two generated sequences can have plausible predicted isolated folds, while their current predicted complexes provide no credible evidence of intended binding.

## Observed results

Source: `results/hard_sequence_refold_v15/mirrored/refold_results.json`, completed with eight prediction jobs; accompanying `native_context_audits.json` and `supervisor.json` were read after completion.

| Sequence | Monomer pLDDT, two checkpoints | Monomer mutual Cα RMSD | Complex iPTM | Complex minimum hotspot Cα distances | Complex template-frame binder disagreement |
|---|---|---|---|---|---|
| v9,69aa | .829/.935 | .981Å | .194/.148 | 57.47/70.50Å | 15.91Å |
| seed23,67aa | .796/.844 | 1.280Å | .197/.150 | 52.38/66.81Å | 18.09Å |

All eight binder backbones passed the declared geometry screen. The v9 complex binder confidence remains .804/.872; seed23 drops to .581/.614. Complex pTM≈.87 includes the large confident receptor and is not evidence of reliable docking. Internal binder-fold agreement and target-relative pose agreement answer different questions.

Native-context staging rejected three complexes for backbone overlap with the complete restored receptor: minima .340,.343,.635Å. Only seed23/model2 passed this staging check; its low interface confidence, off-patch location and pairwise pose inconsistency remain. Staging is not qualification. The audits correctly report zero qualified candidates and four expected complex jobs.

The supervisor records actual GPU release and zero runners/inputs/backlog for this function. I did not independently query cloud state. No GPU was launched by this review.

## What this changes—and what remains unresolved

The strong hard-sequence monomer predictions argue against treating every low-confidence search snapshot as evidence that its encoded sequence cannot fold. They do not establish an erroneous binder rejection: the complex checks still fail. They also do not isolate a single software cause. Checkpoint family, input representation, coordinates, recycling and dropout differ from the earlier coupled optimization. A same-checkpoint hard/soft comparison would study that mechanism, but it is lower priority than obtaining a useful interface now.

Both validators are AF2 checkpoints. Neither experimental stability nor different-architecture corroboration has been demonstrated. High confidence, plausible geometry and a shared fold can still be wrong. These two preselected sequences do not estimate the foldability rate of the whole campaign.

## One next bounded experiment

Use the already tested contiguous native-domain preparation/restoration path for one generation trajectory. Keep the native domain chemistry, existing seed/length draw, loss weights and repaired runtime; do not combine the target change with forced lysination, new length, stronger confidence weight or own-scaffold conditioning. This changes target representation substantially, but avoids gratuitous confounding.

Predeclare one saved hard-sequence snapshot using the existing trajectory saving/selection convention before examining the outcome. Preserve and diagnose that snapshot even if the coupled soft screen rejects it; label it a diagnostic sequence, never an accepted candidate. If the trajectory produces no usable artifact, retain that failure rather than switching selection rules. This prevents a repeat of discarding foldable hard sequences solely because of their search-state confidence.

Bound follow-up to the existing two-checkpoint monomer triage and, only if it advances, two coordinate-free full-native-target complex jobs. Keep generation plus validation under an explicit finite allocation/runtime envelope with no automatic successor. CPU restore/audit any reduced-domain pose before interpreting its local interface. Full-target complex inference must retain the accurately disclosed monomer approximation unless dimer-aware inference is actually implemented and verified; rigid chainB restoration is not equivalent.

Success for this experiment requires a physically credible interface that survives native context and fresh sequence-only prediction. A high isolated-fold score, a reduced-domain-only score or mere proximity/overlap is insufficient. If this one trial again yields foldable but off-target sequences, stop local seed/weight rescue and assess a different interface-generating method against the established validation path.

## Why own-scaffold reuse is second priority

The v9 fold is the better current scaffold prospect because it retains high confidence in both contexts, but it has no demonstrated TfR-binding surface or reliable pose. Conditioning on it could reduce folding search difficulty; it does not solve interface placement by itself.

Pinned source supports scaffold inputs, but the semantics matter: `prepare_binder_chains` calls `Protein.from_structure(...).with_scaffold(...)`. Unedited residues are template/sequence-held; DESIGN edits replace selected residues with random logits and unresolved coordinates. Simply supplying a folded CIF does not perform a validated continuation of the original fully designable sequence. Choosing protected regions, redesigned regions and placement would require a separate declared experiment. [Pinned scaffold implementation](https://github.com/PacesaLab/BindCraft2/blob/b6fb3ae15ab6a702f1156cffc11bcc3323c2fd06/bindcraft/protein.py).

The official FAQ requires designs produced from scratch and prohibits modifying an existing binder; it also requires diversity from known proteins. These current project intermediates have documented from-scratch generation and no demonstrated binding. My interpretation is that continuing their computational design can preserve de novo provenance, rather than borrowing a known binder. The FAQ does not explicitly adjudicate this particular reuse scenario, so eligibility is not guaranteed. Preserve parent sequence hashes, source runs and all transformations, disclose reuse, and do not treat related variants as independent novel families. No known therapeutic or pathogen-derived binder is proposed as a seed. [Official challenge FAQ8](https://proteinbase.com/competitions/anthropic-adaptyv-2026/challenges/tfr1).

Original provenance alone does not establish novelty: the organizer's classification considers sequence and structural similarity to known proteins. A freshly generated common fold can still score poorly. No novelty search or organizer eligibility decision for these sequences has been obtained here. [Organizer novelty methodology](https://www.adaptyvbio.com/blog/novelty).

## Priority versus causal purity

Getting a credible submission need not wait for a factorial causal study. A coherent multi-component method change can be justified when required to make the method work, but must be declared as such rather than attributed to one factor. Here the domain adapter and refold path already exist, so one domain trial is the smallest useful next production hypothesis. Adding scaffold conditioning now introduces more choices without binding evidence and risks another pipeline detour. Controls and full mechanistic diagnosis should not block this bounded test.

Recorded runtime estimate for the completed diagnostic is $0.19353935073, bringing the supplied prior Modal total to approximately $9.48998102613. Runpod remains $35.6652 including historical use; consultation $1.1922999; other costs unreconciled. No additional GPU spend from this consultation.
