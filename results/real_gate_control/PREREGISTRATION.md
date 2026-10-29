# 6W3H real-gate preregistration

Registered at 2026-10-22T15:47:45.545731+00:00 before executing any new real-gate control output.
Prior custom calibration exists in results/therapeutic_control and is known;
this is not a blinded experiment. It did not exercise audit_full_context.py.

## Locked bands supplied by Ricardo

- GATE VINDICATED: target CA RMSD <=2 Å, zero or low-single-digit backbone pairs
  under2 Å, and minimum backbone separation about1.5 Å or greater.
- GATE FALSIFIED: tens to hundreds of backbone pairs under2 Å and minimum
  distance below1 Å, resembling the retired66 candidates.
- INCONCLUSIVE: anything between, or partner_fit >2 Å. Revisit annotation;
  never change thresholds to rescue a result.
- Zero clashes is not required by this interpretation. The production script's
  stricter rigid_context_pass flag is recorded separately, without alteration.

## Locked execution

Independently read every P02786 struct_ref_seq segment and reconcile all166
construct positions through entity_poly_seq / pdbx_poly_seq_scheme, including
unmapped and unobserved positions. No alignment inferred from chain letter.
Use every observed, annotated native-mapped CA; no outlier trimming or local fit.
Test both deposited Fc/receptor copies, identified by deposited assembly and
source interface geometry, against both native11QC subunits A/B. Report allfour
results; neither a bad nor good copy may be hidden. Each copy gets its own band.
A compatible real complex vindicates the gate's ability to distinguish gross
omitted-context occupancy; an incompatible copy remains inconclusive and does
not become a negative control. If any compatible control falsifies, the overall
answer is falsified. If none are compatible, overall inconclusive.

Call the actual fit(mobile,fixed), with source receptor as mobile and native
11QC as fixed; transform the source Fc with the same x@R+t. Save mapped receptor
crop A/B with Fc renamed C for the unchanged audit, plus an independently
restored full-native-dimer coordinate deliverable. Run audit_full_context.main()
unchanged with only its ROOT redirected to an isolated fixture directory;
copy/symlink identical11QC data there. Preserve its raw outputs, source/helper
hashes, fitting metadata, complete mapping, and original66 statistics. The script
restores omitted context itself; feeding a fully restored dimer would leave an
empty omitted-context set, so that is a separate artifact, not audit input.
No changes to audit_full_context.py, run_context_ensemble.py, canonicalize(), or
thresholds. Fc is a crystallographic control, not a submission; its >150-residue
architecture is outside canonicalize's minibinder admission restriction. Do not
truncate Fc to evade that rule. No GPU, no generator, no account or submission.
