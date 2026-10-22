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
