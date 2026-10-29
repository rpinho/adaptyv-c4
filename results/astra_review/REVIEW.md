# Independent review — 2026-10-22 UTC

## Verdict and urgent finding

There is no qualifying native-context binder. The original crop produced physically inadmissible poses; the native v3 pilot produced low-confidence folds; v5 introduced a separate, deterministic coordinate-initialization defect. Treating these as one problem and responding with more length/weight changes is not justified.

**P0: v5's unresolved binder is trapped by singular residue frames.** I independently read its saved trajectory with `Bio.PDB.MMCIF2Dict`: target A has 639 distinct Cα coordinates; binder B has 69 Cα atoms but **one unique coordinate, zero span on every axis, and 68 consecutive Cα distances equal to zero**. All six saved v3 binders have distinct Cα coordinates and nonzero spans. Thus this defect explains the collapsed v5 trajectory, not the earlier v3 failures.

The pinned implementation gives a complete causal mechanism:

1. `Protein.empty()` uses zero atom coordinates and false atom masks (`bindcraft/protein.py:189–195`).
2. `alphafold_input_features()` passes these coordinates directly as `initial_atom_pos` when `large_complex_initialization` is enabled (`bindcraft/af2.py:186–192`).
3. The multimer structure module builds a residue frame from N, CA and C without checking whether those atoms exist (`folding_multimer.py:508–516`).
4. `make_transform_from_reference()` delegates to `Rot3Array.from_two_vectors()`. Coincident N/CA/C give zero axes and hence the **zero matrix**, not an identity rotation. Epsilon-protected normalization keeps the result finite, so a NaN check alone misses it.
5. `FoldIteration` composes `rigid @ rigid_update`. Composition uses `R_new = R R_update`, `t_new = R t_update + t`. With R=0 and t=0, every subsequent update preserves the collapse.

The biotite divide warning is compatible with this collapse but is not the proof. The coordinate audit plus frame construction/composition establish the defect. Translating the receptor cannot repair a singular binder rotation. Increasing binder length or pLDDT weight cannot repair it either.

Primary pinned source: [frame construction](https://raw.githubusercontent.com/PacesaLab/BindCraft2/b6fb3ae15ab6a702f1156cffc11bcc3323c2fd06/bindcraft/af/alphafold/model/all_atom_multimer.py), [rotation construction](https://raw.githubusercontent.com/PacesaLab/BindCraft2/b6fb3ae15ab6a702f1156cffc11bcc3323c2fd06/bindcraft/af/alphafold/model/geometry/rotation_matrix.py), [rigid composition](https://raw.githubusercontent.com/PacesaLab/BindCraft2/b6fb3ae15ab6a702f1156cffc11bcc3323c2fd06/bindcraft/af/alphafold/model/geometry/rigid_matrix_vector.py). Local source inspected under `/tmp/c4-upstream` at revision `b6fb3ae15ab6a702f1156cffc11bcc3323c2fd06`.

**Operational snapshot:** the parent controller cancelled v5 during review. `results/modal_campaign_v5/search.json` and `batch_1/supervisor.json` now record cancellation and verified zero runners/inputs/backlog. This was not a completed two-trajectory biological comparison. Cumulative Modal runtime upper estimate **$3.920816**; newly supervised Runpod runtime **$1.7952**; older campaigns, image builds and storage remain separate/unreconciled. This review launched no GPU, changed no running code and created no account or submission.

## Minimal repair and falsifiable calibration

Repair unresolved **initialization frames**, preserving the actual target and atom masks. A centralized feature-construction repair is reasonable: create a separate initialization-coordinate array, substitute a canonical noncollinear N/CA/C triplet for residues without a valid backbone frame, and leave `all_atom_positions`, template features and all atom masks untouched. Do not change the actual protein coordinates to pretend unresolved atoms are observed. Because monomer and multimer frame conventions differ, validate both actual upstream conversion functions; do not assume one triplet gives identity in both.

Validity must include N/CA/C presence, finite coordinates, nonzero axis lengths and noncollinearity. A mask-only check misses malformed observed frames. Preserve valid resolved frames exactly. For a residue with a valid CA but missing other backbone atoms, preserve a deliberate translation rather than silently teleporting it; distinguish this case from a fully unresolved de novo binder. Check padding as well. Verify recycling, because `af2.py:199–200` replaces initial positions with the previous prediction: a feature-only repair does not intercept that second path.

Before another GPU allocation, run a small CPU numerical regression using the **actual pinned frame functions**, covering zero, partial, collinear, finite valid and nonfinite inputs. Require finite matrices, `RᵀR≈I`, determinant≈+1, preserved valid target frames, unchanged true atom/mask features, and a nonzero test translation surviving composition. The old behavior must fail this test. Track patch/source hashes in runtime metadata and validate the deployed function, not just a local NumPy reimplementation.

Then use one bounded GPU diagnostic with the same local-frame target, length, random keys and loss schedule as v5, with only the frame repair changed. Save the first forward prediction and stop immediately if the binder remains collapsed or the frames become invalid. The falsifiable success criterion is ordinary nondegenerate binder geometry and finite losses/gradients; **passing the pLDDT gate is a separate outcome**. Use the already archived failed v5 trajectory as the negative control rather than paying to reproduce the known defect. A repaired noncollapsed low-pLDDT result validates only the software repair.

## What the earlier evidence does and does not establish

**Old crop:** the restored-context audit and buried full-dimer SASA values justify retiring these poses. Apparent acceptance, favorable pKa scores and sequence variants of the same backbone cannot overcome missing receptor sterics. This is a model-context failure, not evidence that the whole target is undesignable.

**Native v3:** five of six trajectories failed screen and one failed refine; none reached sequence redesign. Six attempts are a small pilot. Even at a hypothetical 14% per-trajectory success rate, zero successes in six has probability about 40%; the crop's historical 14% is not an applicable estimate for the native target. Neither impossibility nor adequate sampling is established.

My additional coordinate audit shows a neglected failure mode: **none of the saved v3 binders acquired the requested patch**. I selected chain A hotspot Cα atoms by original `auth_seq_id` (205/206/208/209/210; verified the target starts at auth residue122), and chain B binder Cα atoms. Minimum binder-to-hotspot Cα distances were:

| v3 batch / length | Minimum hotspot Cα distance (Å) | Minimum any-target Cα distance (Å) |
|---|---:|---:|
| 1 / 69 | 52.02 | 1.38 |
| 1 / 73 | 49.92 | 4.11 |
| 2 / 69 | 48.44 | 6.69 |
| 2 / 70 | 57.74 | 5.27 |
| 3 / 67 | 48.63 | 3.97 |
| 3 / 65 | 46.68 | 7.40 |

Binder adjacent-Cα median distances are 3.68–3.75Å; these structures are noncollapsed. They can contact other receptor regions while missing the chosen patch by tens of ångströms. The first model's 1.38Å target/binder Cα separation is an additional severe clash. These are rejected trajectory snapshots, not final validated structures; nevertheless they show that “fold-confidence failure” is an incomplete diagnosis. No causal explanation of failed patch acquisition is proven here.

The actual upstream `interface_contacts_loss()` (`loss.py:545–555`) uses hotspot flags when present, with a distogram-derived contact objective. The current config uses cutoff20Å, not a direct atomic docking constraint. Thus disabling `forced_targeting` does **not** remove all hotspot guidance. Conversely, a final hotspot filter and a low distogram loss cannot be assumed to guarantee physical patch acquisition. Inspect resolved flags and Cartesian patch distances during optimization. The archived weighted contact losses remain about5.9–7.0 at the end; increasing pLDDT alone could reward a confident fold in the wrong location.

**Epitope:** reopening the accessible apical surface after the buried-groove audit was justified. Accessible Leu209/Val210 and a remote dimer interface make it a reasonable hypothesis. That does not establish a productive binder interface, mouse conservation sufficient for cross-reactivity, or a pH switch. Acid counts and one target histidine census cannot decide those questions. Keep full-context sterics and observed-glycan checks; the present monomer generation remains an approximation and the unresolved glycan envelope remains untested. Do not move the epitope merely because two trajectories fail a fold screen.

**pH method:** the linkage equation is a useful model; the favorable ranking is partly a consequence of that equation, not independent validation of predicted pKas. Repeated sequence variants at one residue are not independent biological replicates. The independent free-state preparation already eliminated an attractive high-pKa claim. Preserve this negative result. Establish plausible native-context binding first, then require matched free/bound preparation, complete group accounting and perturbation controls. A Rosetta ddG cutoff is neither affinity nor the protonation free-energy ensemble. The internal −2 threshold does not equal the challenge assay endpoint. Several historical METHOD.md paragraphs still use “confirmed,” “measured” and “independent APBS” too strongly despite later corrections; consolidate an authoritative current section so those claims cannot re-enter decisions.

## Two controller bugs to fix before calling a set ready

1. **Raw ranked rows are counted as qualifying scaffolds.** `drive_modal_adaptive.outcome()` counts `3_Ranked/!_Ranked.csv`, while the driver checks only the CPU subprocess return code. `postprocess_gpu_campaign.py` catches per-design canonicalization/context errors into JSON, then can finish normally with exit0 and `status='completed'`. The controller can therefore increment acceptance and choose more of a configuration even when every candidate failed native-context admission. Read and enforce the scientific follow-up outcome explicitly. Distinguish generator acceptance, native-context qualification and switch evidence; none implies the next.
2. **Rows are not independent scaffolds.** Upstream `campaign.py:209–236` writes multiple `_seqN` candidates from one trajectory under the same hash. With `kept_sequences=2`, three ranked rows can represent fewer than three scaffolds. Count unique backbone/trajectory hashes and apply diversity checks after native-context qualification. The present name `scaffold_set_ready_for_review` overstates the raw-row count.

The controller also treats every `screen`/`refine` termination as a fold-confidence diagnosis; v5 demonstrates why that is unsafe. Numeric/geometry failure must stop allocation rather than trigger longer sequences. Later-stage/empty acceptance warrants inspection of explicit rejection reasons, not just the stage label. The existing single-GPU limit, bounded spend reserve, exact-call cancellation and verified release are valuable and should remain.

## Prioritized next work

| Priority | Action and stop criterion | Time/compute | Scientific value |
|---|---|---|---|
| P0 | Keep v5 and its queued length/weight profiles stopped. Repair and test initial frames as above; fix scientific acceptance accounting. | CPU only, tens of minutes; no GPU needed | Eliminates a proven defect and false success reporting. |
| P1 | Run the one repaired-frame diagnostic above. Inspect first-forward geometry, then the normal screen. Stop on collapse/nonfinite frames. | One finite H100 diagnostic; first-forward check can precede a full50-step screen | Tests the repair directly with minimal confounding. |
| P1 | Record physical patch distances, actual hotspot flags, contact progress, binder geometry and confidence alongside losses. Compare with saved v3/v5 using matched random keys where feasible. | Mostly CPU; reuse generated artifacts | Separates failed folding from failed docking and broken representation. |
| P2 | If repaired trajectories are noncollapsed but still remote from the patch, test the positioning/contact hypothesis first; if they reach the patch but fail fold confidence, test one fold objective or one length change, not both at once. Use matched seeds and retain identical final acceptance gates. | Small paired pilot, bounded by remaining budget; scale only after a meaningful intermediate improvement | Every outcome changes the next decision rather than cycling through presets. |
| P2 | Run one small official upstream example or documented positive structural control in the same pinned runtime if implementation behavior remains suspect. Keep any known sequence strictly a calibration control, outside the candidate/submission pool. | One short diagnostic; do not turn it into another open-ended campaign | Distinguishes general runtime/pipeline failure from target-specific difficulty. |
| P3 | After a native-context scaffold qualifies, add independent complex prediction and mouse/context checks; only then run matched pH-switch comparisons. | Deferred until a real substrate exists | Avoids expensive polishing of impossible poses or unqualified folds. |

Do not promise a binder from these repairs. The immediate deliverable is a functioning, auditable generation experiment with an interpretable failure mode. Sampling can then be justified by observed progress rather than by the fact that a GPU is available.

## Review scope and verification

Read HANDOFF_CODEX.md, METHOD.md, DECISIONS.md, the current v5 decision/status/search/supervisor artifacts, v3 outcomes, full-context and free-network reports, local policy/driver/postprocess/frame-preparation code, actual saved CIF/CSV outputs, and the relevant pinned upstream initialization/loss/campaign implementations. Source-level frame reasoning and local coordinate calculations were independent checks. No upstream Python frame function was executed in this review; that execution remains an explicit pre-allocation requirement. Only this review file was written.

## Addendum — review of the implemented repair and one-call diagnostic

Reviewed `scripts/frame_guard.py`, `patch_frame_guard.py`, `validate_frame_guard.py`, `frame_gpu_diagnostic.py`, `run_frame_diagnostic.py`, the new Modal `frame_diagnostic` function and driver changes, and `results/frame_guard/CPU_VALIDATION.json`.

**The core repair is appropriate; I found no blocking correctness error in the initialization guard.** It changes only the separate initialization positions; valid N/CA/C frames pass through, missing/degenerate/nonfinite frames receive a canonical triplet, and an observed finite CA preserves translation. The source patch is SHA-locked and touches both the initial feature and recycled prediction paths. Using sequence-derived `atom37_atom_exists` on recycled predictions is correct: the original false template mask would incorrectly reset a now-resolved binder. Padding receives finite frames while its existing sequence mask remains unchanged.

The saved numerical regression reports actual pinned monomer and multimer conversion functions producing orthonormal determinant+1 matrices, preserved valid target frames and a nonzero translation update; the old zero-frame control fails. This directly addresses the established mechanism. Two precision limits remain: the initial-input/recycling **integration** is checked by source anchors and direct helper calls, not by executing the complete feature/recycle pipeline on CPU; and unchanged NumPy input assertions demonstrate no mutation, not a complete template-feature equivalence test. Source inspection supports the intended narrow modification, and the real GPU diagnostic is the appropriate next integration test.

**Two inexpensive diagnostic correctness fixes were communicated before allocation:**

1. Assert and record a real GPU JAX backend/device. The new diagnostic did not inherit the previous campaign's CPU-fallback check. Allocation of a GPU by Modal does not by itself prove JAX used it.
2. Compute the binder gradient norm after casting to float64 and require `isfinite(norm) and norm > 0`. Upstream gradients can retain float16 dtype; individually finite elements do not guarantee a finite norm because squaring/summing may overflow. The reviewed condition `gradient_norm > 0` would accept infinity.

The one-call diagnostic has an appropriate limited aim: a real sequence-gradient invocation with initial/recycled frames and a noncollapsed binder, not a candidate. Its matching trajectory-key split follows `run_trajectory()` and the explicit69-residue assertion checks the intended first-v5 control length. Record the selected model, input SHA and sequence-logit hash if claiming exact first-forward reproducibility; matching length/seed alone is weaker evidence. `design_recycles=1` means the single gradient invocation includes recycling, so it exercises both repaired paths.

The30-minute Modal function limit,25-minute subprocess limit, zero retries, unique output directory, single-use container, release-checking driver and disabled CPU candidate follow-up keep this bounded. The old adaptive queue must remain disabled; this diagnostic intentionally does not solve the separate qualification-count and policy-routing bugs. A passing noncollapse test still does not show correct peptide geometry, useful confidence, patch acquisition, affinity or pH release. Save all those as separate future questions. I made no code changes or allocations during this addendum review.
