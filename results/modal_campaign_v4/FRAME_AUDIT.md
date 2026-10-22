# Coordinate-frame defect in the new initialization test

The new v4 initialization test was cancelled early and its actual GPU release
verified, rather than running its queued replicates. Estimated v4 runtime upper
cost$0.2559. This is not a completed0/2 design result.

CPU inspection found that the supplied target hotspot centroid is
(92.922,135.386,93.208)A, **188.816A from the origin**, while an unresolved
de novo binder starts at the origin. The nearest full-receptor CA is180.090A
away. Upstream Protein.from_structure and initialize_design_trajectory retain
these input coordinates; folding_multimer directly initializes its residue
frames from initial_atom_pos when the new flag is enabled. There is no target
centering in those paths.

This is an avoidable frame/placement defect in the newly added initialization
experiment. It **does not establish the cause of v3's poor binder folds**:
v3 did not inject these input positions into the structure module. The v4
decision note's claim of a well-controlled initialization comparison was
premature without this geometric preflight.

Corrected v5 input retains all639 native target residues, atom identities,
numbering and shape. A rigid translation places the binder's starting origin
15A outside the hotspot center along the target-center-to-patch direction.
CPU checks against the whole native dimer give7.955A nearest-atom clearance at
that point. This placement is another hypothesis, not proof of better binding.
The exact transform, full atom identity comparison and all CA pair-distance
checks are saved in ../modal_campaign_v5/human_11qc_local_frame.audit.json.

Sources at the pinned upstream revision:
`bindcraft/protein.py` Protein.empty/from_structure;
`bindcraft/protein_preparation.py` initialize_design_trajectory;
`bindcraft/af2.py` alphafold_input_features;
`bindcraft/af/alphafold/model/folding_multimer.py` generate_affines.

Operational correction: three pre-set repeats are sampling, not continuous
autonomous model improvement. Future background status must name the actual
finite policy. The next search changes settings only through a declared
failure-dependent policy, with GPU release and CPU validation between calls.
