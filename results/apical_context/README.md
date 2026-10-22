# CPU groundwork for a conditional contiguous-domain route

This directory is not a candidate set or a ranked-design input. No domain GPU
search has been launched. The current V12 GPU test still uses the complete639-
residue receptor backbone, with its explicitly tracked forced-search chemistry.

`scripts/restore_apical_context.py` prepares exactly native author residues
A189-383, preserving all195 consecutive residue identities and input coordinates.
The prepared file here comes from the existing full native local-frame input.
Its audited internal peptide C-N range is1.322-1.349A. This reduced structure
still omits interdomain contacts, partner subunit and glycan context; it does not
establish that the native isolated domain is stable.

The separate restoration operation fits all195 domain CA atoms to native11QC,
reports global and local five-hotspot RMSD, refuses excessive deformation, and
transforms the binder into both complete native639-residue receptor chains.
Existing2A binder/full-target backbone and observed-NAG clearance limits are
retained. Mutated/lysinated, disconnected, wrongly numbered, incomplete or
nonfinite domain inputs fail closed. The adapter is deliberately not connected
to `postprocess_gpu_campaign.py` or its scientific admission path.

Tests include a rigidly moved native domain with a deliberately distant
nonbinding30-residue fixture, an omitted-receptor-protein overlap, a local patch
deformation, wrong chemistry, deleted internal residue and atom-preserving
preparation. The positive fixture tests restoration only; it is not an engineered
binder. The expanded project suite currently passes72 tests in `.venv-rosetta`;
the log is in
`../pipeline_repair/regression_tests.log`.

**A restored staging PDB remains unqualified.** It still requires a separate
native full-target refold, unchanged quality/patch/context gates and explicit
provenance before any pH engineering. Existing strict canonicalization still
rejects direct crops. No new search or admission relaxation is authorized by
this software artifact. Published therapeutic binders remain calibration
controls only, never candidate starting points.

Independent review: `../astra_review/DOMAIN_ADAPTER_REVIEW.md`. Its two integrity
gaps were fixed before integration: exact11QC SHA256 pinning, and explicit
unknown/nonprotein residue rejection. Broken binder peptide bonds are also
rejected. Added real regressions for each. Observed NAG atoms are checked for
clearance, then excluded from the protein-only staging PDB; their absence is
explicit in its audit. This work allocates no
GPU. Latest completed Modal runtime estimate before the active V12$5.8023551;
V12 accrues at requested upper$4.649328/hr. Runpod$1.7952 plus historical$33.87;
consultations$1.1922999; setup/storage/older campaigns remain unreconciled.
