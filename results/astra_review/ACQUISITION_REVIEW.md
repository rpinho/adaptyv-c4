# Target acquisition review — 2026-10-22

## Decision

**No additional concrete hotspot-mapping, loss-sign, padding or recycling bug was established. Run one bounded optimization trajectory with the repaired native v5 configuration before changing targeting settings.** The repaired v8 path has produced one unoptimized prediction, not an optimization trajectory. Its off-patch result cannot establish that repaired optimization fails.

The source and CPU checks support the intended hotspot guidance. They also expose why guidance and initialization are insufficient guarantees: hotspots affect the optimization loss, not a direct neural-network docking instruction; the current contact objective rewards predicted distances below20Å; and initial placement is not a coordinate restraint. These are limitations to measure, not proof that a particular setting change will solve acquisition.

This consultation used existing artifacts and CPU/read-only source checks. **No new GPU spend.** Prior runtime estimates supplied for this checkpoint: Runpod **$1.7952**, Modal **$4.0531838**, plus separate/unreconciled extras; consultations **$1.1922999**. No models ran during the reported overnight gap after v8 completed. No worker, campaign or account was created by this review.

## What was checked numerically

I executed selected, unmodified function definitions from the pinned source with JAX CPU; I did not import the full inference stack or load model weights. The local `protein.py`, `protein_preparation.py`, `loss.py` and `trajectory.py` bytes were independently compared with GitHub revision `b6fb3ae15ab6a702f1156cffc11bcc3323c2fd06`; all matched. AF2's intentional frame-guard patch is tracked separately by its existing manifest.

### Hotspot mapping is correct in the reconstructed input path

`target_binding_site_flags()` generates chain-qualified HOTSPOT annotations; `annotated_residue_flags()` matches them against original residue numbers, not zero-based positions or CIF label numbers. On the actual local-frame target:

| Original A residue | Identity | Zero-based target row | Flags | Global row with96-residue padded binder |
|---|---|---:|---:|---:|
| 205 | LYS | 83 | 22 | 179 |
| 206 | ASN | 84 | 22 | 180 |
| 208 | ARG | 86 | 22 | 182 |
| 209 | LEU | 87 | 22 | 183 |
| 210 | VAL | 88 | 22 | 184 |

Flags22 are SEQUENCE4 + TEMPLATE2 + HOTSPOT16. All639 target residues retain TEMPLATE. Padding adds27 positions to the69-residue binder; the target offset is therefore96, not69. Canonical state/chain renaming preserves the Protein objects and flags, and `chain_residue_slices()` uses the same sorted chain ordering as model array construction. No lost HOTSPOT bit or off-by121 numbering error was found.

Limitation: this reconstructs the exact parser behavior on the saved input/configuration. The actual GPU-side hotspot mask and padded global indices were not archived, so direct runtime-mask verification remains a cheap useful preflight assertion.

### The contact loss selects the intended pairs

I executed the actual `interface_contacts_loss()` and its mask/contact helper functions on a synthetic distogram with the real hotspot identities, real binder residues, padding and non-hotspot target decoys. This isolates routing/sign from model behavior:

| Synthetic case | Loss |
|---|---:|
| All pairs predict far distance | 11.99955 |
| Two real binder partners per hotspot predict a2.94Å bin | 0.00000625 |
| Only non-hotspot decoys get favorable contacts | 11.99955 |
| Only binder padding gets favorable contacts | 11.99955 |
| Intended partners predict a19.19Å bin | 0.00000625 |

The distogram-gradient norm was0.31880, entirely on hotspot-to-real-binder pairs; decoy and padded blocks were zero. Thus the kernel has the expected direction and does not reward padding or unrelated receptor surfaces in this test.

The equal2.94Å/19.19Å scores are an important **objective limitation**, not a new implementation bug. `cutoff=20.0`, `contacts_per_residue=2` means broad predicted proximity can satisfy this term. It does not measure Cartesian atom contact, steric feasibility or final hotspot coverage. `distogram_pair_loss(binary=False)` is an entropy-containing conditional cross-entropy, so its raw magnitude is not simply a probability of contact. Comparing scalar magnitudes across loss terms does not establish which sequence gradient dominates. The available v8 aggregate gradient proves only that the total objective has a finite nonzero sequence gradient; it does not isolate the hotspot contribution through the network.

Source: [pinned contact-loss implementation](https://raw.githubusercontent.com/PacesaLab/BindCraft2/b6fb3ae15ab6a702f1156cffc11bcc3323c2fd06/bindcraft/loss.py), especially `best_contact_mean`, `distogram_pair_loss`, `interface_contacts_loss`.

## Geometry: distinguish three different experiments

All measurements below were recalculated from saved CIF atom tables using `auth_asym_id` and original `auth_seq_id`, with chain A as target and B as binder. Cα distances are geometric diagnostics, not substitute acceptance criteria.

| Artifact | Optimization exposure | Binder Cα uniqueness | Minimum hotspot Cα distance | Minimum hotspot heavy-atom distance |
|---|---|---:|---:|---:|
| Native v3, six saved rejected trajectories | Five50-step screens; one75-step screen/refine | Every binder Cα distinct | 46.68–57.74Å | 42.25–52.13Å |
| v5 first saved rejected trajectory | 50-step screen with singular frames | 1/69 | 11.48Å | 5.45Å |
| Repaired v8 `first_forward.cif` | One unoptimized sequence-gradient invocation, including one recycle | 69/69 | 38.54Å | 33.30Å |

V5's apparent proximity is meaningless: it is a collapsed binder. Its contact loss nevertheless fell from6.8 to4.6 across recorded steps, directly illustrating that this distogram objective can improve while Cartesian geometry remains impossible.

The six v3 snapshots are off-patch and often overlap another receptor surface: minimum any-target heavy-atom distances range0.26–1.32Å. Their nearest receptor residues cluster around the opposite surface, including318–322,688–691 and752–754; the first69-residue snapshot instead overlaps the627–644/759–760 region. None acquired the requested205–210 patch. Recorded contact-loss minima were5.2–6.2 and final logged values5.9–7.0, showing no strong contact-objective solution in those logs. These are model-dependent, rounded scalar outputs; they do not diagnose the cause.

The repaired v8 binder is noncollapsed, with minimum any-target Cα/heavy-atom distances9.89/6.05Å and nearby receptor residues689,753,740,472. It lies on a similar broad face to several v3 failures. Its adjacent-Cα minimum1.56Å remains abnormal. This is an observation from one unoptimized sample, not evidence that it is trapped in that pose after optimization.

Do not pair a saved trajectory structure blindly with the last CSV row: `run_gradient_design_stage()` can retain an earlier best-loss prediction. The files are selected stage outputs, and per-round structures were not retained here. Their temporal path toward or away from the patch cannot be recovered from scalar losses alone.

## Initialization, features and recycling

The raw v3 target coordinates were far from the origin, but v3 had `large_complex_initialization=false`; that distance is not an explanation of its failure. V5 enabled coordinate initialization and suffered the now-demonstrated zero-rotation defect. The local-frame translation preserves the receptor and initially places the unresolved binder origin15Å from the patch centroid. V8 repairs those frames and exercises the guarded initial and recycled paths.

`alphafold_input_features()` does not encode HOTSPOT as a special input channel. TEMPLATE, DESIGN, PADDING and CYCLIC affect features; hotspot flags are consumed by losses. Holding sequence/coordinates fixed, changing only hotspot annotation is therefore not expected to direct the **first** predicted structure toward the newly labeled patch. It changes optimization gradients. The corrected placement supplies an initial relative position, but the multimer structure module calls every fold iteration with `update_rigid=True`; it does not freeze target frames or tether binder frames to the patch. Recycling then uses predicted positions, with stop-gradient through preceding passes. No origin-reset or wrong-template-mask bug was found in the repaired recycling path.

Re-alignment to the target template applies one rigid transform to the whole predicted complex; it does not change binder-to-hotspot distances. It cannot itself explain the38–58Å separation. Neither the first pre-recycle prediction nor raw model frames were saved, so the existing artifacts cannot say whether the first structure-module pass or recycling displaced the binder most. That is a falsifiable future diagnostic, not an established root cause.

## Immediate bounded test

Proceed with the parent's proposed **one repaired-native v9 trajectory**, seed20261021, expected69-residue draw from the same65–75 range, full local-frame639-residue target, `forced_targeting=false`, unchanged objectives and acceptance gates. Only the frame repair changes scientifically from v5. Use a finite runtime ceiling, no automatic successor and verified GPU release before CPU qualification.

Before allocation, assert the runtime hotspot mapping/flags/global indices above and the target/guard hashes. During the trajectory, record minimum patch distances, per-hotspot physical contact counts, binder geometry validity, confidence, model identity and objective terms alongside the round number. Save the selected structure's actual round. Halt on collapse or nonfinite geometry; ordinary rejection remains a rejection.

Interpret outcomes separately:

- Noncollapsed geometry with measurable movement toward the patch validates useful optimization progress, not a binder.
- An accepted generator output must still pass native-context admission, unique-parent accounting and unchanged final gates.
- A normal rejected trajectory that remains off-patch is evidence for an acquisition-focused next test. It does not establish target impossibility.
- Scalar contact improvement without physical approach motivates checking per-term sequence gradients and distogram/Cartesian agreement before simply increasing a weight.

A further CPU source inspection cannot answer whether the repaired model optimizes successfully. The missing experiment is one actual optimization trajectory; there is no reason to await another unrelated consultation before that bounded test.

## Forced targeting versus changing cutoff/weight

**Forced targeting is a defensible later paired hypothesis, but not the cleanest immediate test.** The pinned implementation identifies a10Å shell around the hotspot atoms, protects it, and converts exposed residues outside it to lysine while retaining the backbone. It is a broad change to sequence/template chemistry; an intact input file does not make the optimization target native. The previous full-monomer forced pilot reported371/639 substitutions. Upstream restores the native target at `harden`; screen/refine/anneal failures can still contain the surrogate sequence. Such outputs cannot enter native binding or pH qualification.

Its falsifiable hypothesis is that removing attractive off-patch surface chemistry helps acquisition while the complete receptor backbone supplies steric context. Compare repaired forced=false and forced=true arms with matched initial sequence, seed, model schedule, lengths and losses; then require native restoration and all existing gates. Patch acquisition on the surrogate alone would not establish a native binder. Do not present toggling this flag as merely strengthening a native hotspot constraint.

Sources: [pinned forced-targeting implementation](https://raw.githubusercontent.com/PacesaLab/BindCraft2/b6fb3ae15ab6a702f1156cffc11bcc3323c2fd06/bindcraft/epitope_targeting.py), [restoration in trajectory stages](https://raw.githubusercontent.com/PacesaLab/BindCraft2/b6fb3ae15ab6a702f1156cffc11bcc3323c2fd06/bindcraft/trajectory.py).

Changing the contact cutoff is narrower with respect to target representation, but changes what the contact objective rewards. A tighter cutoff may improve locality or weaken useful long-range optimization gradients; the present evidence does not decide. Raising its weight preserves its broad distance definition and can emphasize the same surrogate without fixing physical docking. Neither should be combined with a length/epitope change in the first diagnostic pair. Choose between these only after the repaired-native trajectory establishes which problem persists, using physical acquisition and native-context outcomes rather than stage labels.

Only this new review file was written. No pipeline settings, running processes or acceptance gates were changed by the reviewer.

## Addendum — pre-deployment v9 monitor check

Reviewed `scripts/monitored_design.py`, the repaired-v9 policy, Modal invocation and one-run driver. No additional critical execution or sampling blocker was found. The wrapper calls the original `sequence_gradients` once with unchanged arguments and returns its original values; observation does not consume RNG or change objectives. Default binder-alone validation uses `predict`, so it does not hit the complex-only gradient observer. The repaired profile retains the native local-frame target, seed, length distribution, losses and acceptance gates, with one trajectory and no automatic successor.

The necessary `compile_only` passthrough and Python3.12-safe `next(iter(entries))` entry-point access were verified after correction. The first true step now asserts69 residues and records initial sequence-array provenance; imported AF2 and monitor hashes are recorded. I also flagged structured `.jsonl` mirroring and unguarded selection of a complex from future binder-alone gradient calls; the parent reports the former added and the latter receiving explicit non-complex handling. These do not justify delaying the current native one-target experiment. This was a CPU/source review, not a full installed CLI integration execution or GPU launch by the reviewer.

## Conditional follow-up — contact cutoff20→8Å

At the reviewed partial v9 checkpoint,32 true gradient calls had been mirrored; the trajectory was incomplete. The parent reported total loss falling15.68→best10.24 while no step approached the patch more closely than the initial38.512Å. This is not yet a completed rejection. Parent's preceding runtime checkpoint was cumulative Modal upper estimate$4.389; this CPU review adds no GPU usage.

**If v9 finishes off-patch, one matched cutoff8 trajectory is a defensible objective ablation. It is not an established remedy or demonstrably a better hit-rate bet than additional seeds.** Preserve native target, frames, initial sequence/model sampling, seed, length, contacts-per-residue2, every other weight, final acceptance gates and finite stop limits. Its advantage over another unchanged seed is diagnostic information about the broad contact objective. More seeds instead probe stochastic variability; neither has a proven success advantage from this sample size.

I executed the exact pinned `distogram_pair_loss(binary=False)` and its JAX gradient on synthetic64-bin logits. A20Å cutoff selects58 bins;8Å selects20. These are single-pair synthetic values, not real v9 network-gradient measurements:

| Predicted logits | Loss at20 /8Å | Logit-gradient norm at20 /8Å |
|---|---:|---:|
| Uniform | 4.15888 /4.15888 | 0.04020 /0.18541 |
| Strong peak in far/overflow bin | 12.00039 /12.00039 | 1.00819 /1.02431 |
| Strong peak at19.19Å | 0.00459 /12.00039 | 0.00427 /1.02431 |
| Strong peak near5Å | 0.00459 /0.00179 | 0.00427 /0.00170 |

The narrower cutoff closes the19Å proximity loophole and redistributes encouragement into shorter bins. It **does not fix a demonstrated vanishing far-distance gradient**: both cutoffs already have nonzero, similar logit gradients when the distribution strongly predicts far separation. For this loss, the sum of selected-bin logit derivatives is `P(selected bins)−1`, so small selected probability does not by itself imply saturation. Gradients through the actual network to sequence can still be weak, conflicting or unhelpful. The distogram's far bin also does not distinguish a38Å Cartesian miss from a63Å miss as a smooth physical distance restraint would. Changing the cutoff may make optimization harder; predicted short-range mass need not produce a physical interface.

Judge the conditional experiment by actual approach/contact and retained geometry/confidence, followed by unchanged native qualification. Do not compare raw total losses between the20Å and8Å runs as if they were the same objective. A cutoff8 result that lowers its objective while remaining off-patch would argue against further blind weight/cutoff escalation.

**Logging weighted scalar loss components is useful, but does not establish gradient dominance.** In `af2.py`, `weighted_losses` are summed into the returned loss and individually inserted into `prediction.metrics`; these values are already weighted. Record the exact active loss keys, without multiplying weights again or mixing confidence metrics such as iPTM into the sum. For this one-state run, their sum should reproduce the total within numerical tolerance. This identifies which numeric terms account for improvement while the patch remains distant. It cannot reveal individual sequence-gradient norms, gradient cancellation or causal control of geometry. Those require a separate derivative decomposition; an aggregate nonzero gradient is insufficient. In future multi-state use, avoid summing a duplicated shared loss from multiple prediction dictionaries.

No root code or running computation was changed for this follow-up. Recommendation remains conditional on v9 completion, with at most one declared successor experiment rather than an automatic sweep.

## Conditional forced-targeting preparation — native restoration and admission

V9 subsequently completed50 rounds and failed screen (pLDDT0.54/iPTM0.17), with no patch acquisition. Its GPU release was verified by the parent. V10 is a separate cutoff8 test underway; **this note neither predicts its result nor authorizes a forced-targeting successor**. Prior completed Modal runtime upper estimate is$4.6368557, with v10 accruing separately; this consultation adds no GPU usage. Runpod$1.7952 and other/consultation costs remain separate.

For the declared single-target, ordinary binder configuration, the pinned source's normal success path restores native target chemistry before later qualification:

1. `run_trajectory()` saves `wild_type_states` before constructing the lysinated representation. Screen, refine and anneal therefore optimize and filter the artificial target; failures saved there are **not native candidates**.
2. Restoration is the **prepare operation at the start of harden**, before its gradient calls and stage filters (`trajectory.py:170–173,220–221,308`). It replaces the complete target Protein objects and the schedule's target-state reference with the originals, including sequence, atom positions, masks and flags. It retains the designed binder.
3. `build_design_schedule()` calls `restart_for_stage()` with those restored targets. The schedule replaces its target entries and clears cached predictions (`target_schedule.py:20–24,155–160`). I found no stale lysinated-target cache feeding harden in this path.
4. Harden predictions are newly generated; earlier best-round predictions are not reused across the new gradient-stage call. Optional mutation polish and the trajectory final gate follow native restoration. A failure in an earlier stage blocks the MPNN branch (`campaign.py:203`).
5. MPNN preparation explicitly clears target DESIGN flags. Candidate validation inserts only the redesigned binder into the supplied restored target states. The complete PDB target does not take the FASTA-cropping path. Complex refolds and their final filters therefore use the restored target; binder-alone validation has no receptor. Only candidates without failed filters are returned, and accepted structures are written from those validated predictions.

Sources: pinned [trajectory flow](https://raw.githubusercontent.com/PacesaLab/BindCraft2/b6fb3ae15ab6a702f1156cffc11bcc3323c2fd06/bindcraft/trajectory.py), [schedule restart](https://raw.githubusercontent.com/PacesaLab/BindCraft2/b6fb3ae15ab6a702f1156cffc11bcc3323c2fd06/bindcraft/target_schedule.py), [MPNN preparation and refold validation](https://raw.githubusercontent.com/PacesaLab/BindCraft2/b6fb3ae15ab6a702f1156cffc11bcc3323c2fd06/bindcraft/MPNN_stage.py), and [accepted-output writing](https://raw.githubusercontent.com/PacesaLab/BindCraft2/b6fb3ae15ab6a702f1156cffc11bcc3323c2fd06/bindcraft/campaign.py).

**This is a source-backed success-path expectation, not an independent native-identity guarantee at every output boundary.** Upstream treats its supplied original target as the reference and does not explicitly compare the final receptor sequence with the expected UniProt122–760 sequence immediately before acceptance. Early saved files can legitimately remain lysinated. Keep the explicit local `canonicalize()` sequence/numbering/connectivity/full-context admission guard authoritative, and require native identity in the actual refold/output before counting a qualifying parent or running pH scoring. A native-sequence/hash assertion at restoration and on final prediction/output would make a future representation test easier to audit; physical backbone/dimer/glycan checks remain necessary regardless of identity.

The conditional hypothesis is valid only with accurate labeling: **artificial off-patch chemistry during search, complete receptor backbone retained, native restoration and qualification afterward**. It satisfies a native-context admission requirement only when the actual restored output passes that admission. It does not preserve native chemistry throughout optimization, and a surrogate-only confidence or acquisition improvement cannot be counted as native binding. No code, setting or GPU state was changed by this source review.

## Existing COM loss is not local-patch guidance

The pinned `com_distance_loss` (`loss.py:457–466`) uses **all resolved predicted target Cα atoms**, without a HOTSPOT mask or a patch-selector parameter. Its expression is `ReLU(distance(binder_CA_centroid, whole_target_CA_centroid) − radius_ratio × binder_radius_of_gyration)`, with default `radius_ratio=0.9`. Despite the name, these are geometric Cα centroids, not atom-mass-weighted centers.

Do not enable this term as though it attracts the binder to the apical hotspot. It can reward movement toward the whole receptor's interior or another surface; a larger binder radius can also reduce the penalty. A hotspot-specific Cartesian guidance term would be a new objective requiring separate numerical checks for mapping/padding, rigid-transform invariance, finite gradients, and size/clash loopholes. It should not be added simultaneously with the proposed single-variable forced-targeting experiment.

The parent's partial v10 component comparison (first→step26: compactness3.961→−0.495, binder_contacts4.783→3.706, interface_contacts8.085→8.167) supports attributing much of the **numerical** objective reduction to compactness, without establishing sequence-gradient dominance. V10 was incomplete at this checkpoint. This COM assessment used source/CPU only and incurred no new GPU usage or code changes.

## V11 follow-up: confidence-weight sensitivity test

V11 completed50 gradients and failed screen (parent-reported pLDDT0.36/iPTM0.13); none of its50 predictions reached native restoration. All retained375 surrogate target substitutions. The closest hotspot approach was an overlap (Cα1.255Å), and the best-loss round also had malformed adjacent binder Cα spacing0.769Å. This is neither productive native acquisition nor evidence of a valid fold. GPU release was verified by the parent; this review used no GPU.

**A single matched0.1→0.5 pLDDT-weight test is defensible, with a configuration correction: use top-level `weights_plddt_loss` for the pinned revision.** It is not obsolete. `settings.py` accepts only `params`/`prediction_state` within a loss entry; `build_losses()` reads `settings['weights_plddt_loss']`. V11's actual campaign metadata records that key at0.1 and binds the loss to `prediction_state='complex', chain='binder'`. CPU execution of the exact pinned schema-validation functions produced:

- `losses.plddt_loss.weight=0.5`: `ValueError: unrecognized campaign settings: losses.plddt_loss.weight`.
- Top-level `weights_plddt_loss=0.5` with existing loss parameters: accepted by the schema check.

This CPU extraction is not a full installed-stack settings load. Before allocation, run the following in the actual pinned image, passing the archived V11 campaign path; it checks the effective objective rather than merely the JSON spelling:

```python
import copy, json
from pathlib import Path
from bindcraft.settings import load_settings
from bindcraft.loss import build_losses

raw = json.loads(Path(V11_CAMPAIGN_PATH).read_text())
changed = copy.deepcopy(raw)
changed['weights_plddt_loss'] = 0.5
a = build_losses(load_settings(raw), seed=20261021)
b = build_losses(load_settings(changed), seed=20261021)
assert a.keys() == b.keys()
assert a['plddt_loss'].weight == 0.1
assert b['plddt_loss'].weight == 0.5
for name in a:
    assert a[name].required_states == b[name].required_states
    assert a[name].function.keywords == b[name].function.keywords
    if name != 'plddt_loss':
        assert a[name].weight == b[name].weight
assert b['plddt_loss'].function.keywords['prediction_state'] == 'complex'
assert b['plddt_loss'].function.keywords['chain'] == 'binder'
print({name: {'before': a[name].weight, 'after': b[name].weight,
              'params': b[name].function.keywords} for name in a})
```

The pinned loss minimizes the DESIGN-residue mean of `1−predicted_pLDDT`; in this configuration it evaluates the binder in the complex, not an independent free-binder refold. Increasing the coefficient multiplies **this term's contribution to the sequence derivative by five at an identical model/sequence state**. It does not establish that the term dominates the aggregate gradient or that later updates will be five times larger. Scalar loss magnitudes, including compactness accounting for most numerical improvement, cannot establish gradient dominance, cancellation or relative control of geometry. [Pinned loss implementation](https://raw.githubusercontent.com/PacesaLab/BindCraft2/b6fb3ae15ab6a702f1156cffc11bcc3323c2fd06/bindcraft/loss.py); [pinned settings parser](https://raw.githubusercontent.com/PacesaLab/BindCraft2/b6fb3ae15ab6a702f1156cffc11bcc3323c2fd06/bindcraft/settings.py).

This is an **indirect fold-confidence hypothesis**, not a direct geometry-repair objective: there is no peptide-bond, adjacent-Cα, nonbonded-clash or hotspot-distance term inside `plddt_loss`. It could favor a better folded binder, merely raise predicted confidence in an incorrect pose, or trade away approach to the patch. The value0.5 is a declared sensitivity level, not a calibrated optimum.

Keep V11's complete scientific profile, initial logit draw/length69, seed, forced representation, target, other weights/cutoff, sampling, gates and finite limits unchanged. The archived V11 interface-contact cutoff is20Å; do not silently substitute V10's8Å while calling this a one-variable test. Compare unweighted pLDDT, actual geometry and valid patch contacts, not raw totals across differently weighted objectives. Confidence improvement with continued geometry defects or surrogate-only overlap does not support advancement. Any later output must actually restore native chemistry and pass existing full-context qualification before admission. One trajectory can reveal sensitivity for this initialization; it cannot estimate hit rate or justify an automatic weight ladder.
