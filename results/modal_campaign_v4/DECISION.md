# Large-complex coordinate initialization: bounded next hypothesis

The native-target series finished0/6 accepted, with five failures at the initial
confidence gate and one at refine. No gate was relaxed and no adaptive setting
change was recorded. Saved structures retain a confident receptor (mean CA
pLDDT91.9–92.7/100) alongside weak binder folds (40.3–55.0/100). This diagnoses
poor binder generation, not an experimentally confirmed epitope failure.

The pinned upstream source explicitly supplies a `large_complex` property that
sets `large_complex_initialization=true`. `alphafold_input_features` then passes
the supplied native target coordinates as `initial_atom_pos` rather than
starting the structure module from its default. That is a plausible, testable
change for a full639-residue target. It does not establish that initialization
caused the previous failures, and it may fail in exactly the same way.

Sources, same pinned revision as every compared run:

- [Upstream property](https://github.com/PacesaLab/BindCraft2/blob/b6fb3ae15ab6a702f1156cffc11bcc3323c2fd06/settings/property/large_complex.json)
- [Feature construction and recycling](https://github.com/PacesaLab/BindCraft2/blob/b6fb3ae15ab6a702f1156cffc11bcc3323c2fd06/bindcraft/af2.py)

Repeat the same three seed/length draws and two-trajectory batches with that one
flag changed. Keep the entire native639-residue monomer, same epitope, losses,
terminus constraints, stage confidence thresholds and all final acceptance
gates. Native sequence, restored opposite subunit and clash checks still precede
CPU pH analysis. No old crop or rejected trajectory enters the accepted pool.

The existing Modal image is cached; CPU deployment precedes GPU allocation.
H100 GPU $3.9492/hour; requested CPU/RAM runtime upper rate $4.649328/hour.
Fallback A100-80GB requested runtime upper rate $3.198528/hour. At most three
sequential two-hour functions, one GPU at once, zero retries, single-use
containers and verified zero runners between batches. The six-hour runtime
envelope is $27.90 H100/$19.20 A100 plus setup/teardown/storage; the previous
series actually used an estimated$2.8760 upper runtime. Ricardo's higher-spend
and continued-modeling authorization persists; no new approval wait.

New live pointer: `series.json`; first-batch state `supervisor.json`, later
`replicate_2/` and `replicate_3/`. A `CANCEL` file here cancels the exact active
call including its container. Completion/failure releases the actual GPU.
An execution, release or CPU follow-up failure blocks further allocation for
diagnosis. No account creation/submission, KB main or repo main changes.
