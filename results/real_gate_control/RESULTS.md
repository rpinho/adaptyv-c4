# 6W3H control: GATE VINDICATED

The compatible deposited FcB/TfRC complex passes the unchanged omitted-context audit in both native11QC orientations. The alternate FcA/TfRD complex is INCONCLUSIVE because its all-mapped-residue target fit exceeds2 Å. No outliers were removed, no thresholds changed, and no GPU was used.

## Complete independent mapping

Receptor entity2 was selected by the P02786 accession in struct_ref_seq, then joined through struct_ref.ref_id/entity_id. Chain letters were assigned only afterward. All four P02786 alignment rows are:

| align_id | author chain | label chain | label residues | author residues | UniProt residues | count |
|---|---|---|---|---|---|---|
| 2 | C | B | 2–55 | 22–75 | 326–379 | 54 |
| 3 | C | B | 56–158 | 76–178 | 194–296 | 103 |
| 4 | D | C | 2–55 | 22–75 | 326–379 | 54 |
| 5 | D | C | 56–158 | 76–178 | 194–296 | 103 |

Thus each166-residue construct has157 UniProt-mapped residues, not166 native receptor residues. The remaining9 are explicitly annotated in struct_ref_seq_dif: label1/auth21 THR (expression tag), label159–160/auth179–180 ALA-SER (linker), and label161–166/auth181–186 six HIS (expression tag). All9 are unresolved in both copies. There is no missing third P02786 segment. MAPPING.csv enumerates every one of166 positions for each copy, including native mapping, observed status, author/label numbering and unmapped reason. ReceptorC has144 observed mapped CAs; receptorD147. All mapped observed residue identities match11QC; no sequence alignment or extrapolated tag mapping was used.

Deposited assembly1 uses labelA/C/E: Fc authA with receptor authD (labelC). Assembly2 uses labelB/D: receptor authC (labelB) with Fc authB. No symmetry operator beyond identity is required. Raw annotation tables are preserved in ANNOTATION_TABLES.json.

## Actual execution and fit direction

Preregistered bands and copy handling are in PREREGISTRATION.md, with its frozen hash. Earlier custom calibration was already known and explicitly disclosed; this was not a blinded experiment. No new real-gate output was inspected before registration.

run_context_ensemble.fit(mobile,fixed) returns R,t with mobile@R+t in the fixed frame. canonicalize() calls fit(templateA,sourceA) to move the template partner into the source frame. Here the deliberate inverse call fit(source receptor,native11QC) moves the source Fc into the native frame. Both directions were run; equal RMSDs and inverse rotations were asserted. Every observed annotated CA was used. FITS.json records transforms, copy selection and partner_fit.

run_control.py invokes the original audit_full_context.main() unchanged; only its ROOT variable points to an isolated fixture directory containing byte-identical11QC and the four mapped source controls. The original66 inputs/outputs were untouched. Audit inputs contain the observed source receptor crop and complete observed Fc; the real script performs its own omitted-context restoration. Separate *_restored.pdb artifacts contain the full native dimer plus transformed Fc. Feeding those full dimers to the omitted-context audit would leave no omitted atoms and is therefore inappropriate.

The audit itself refits the rounded PDB inputs. Raw outputs are in audit_fixture/results/full_context_audit/. Its automatically generated RESULTS.md contains generic legacy-candidate prose; the per-control interpretation is this report and RESULTS.json. Hashes before/after confirm audit_full_context.py, run_context_ensemble.py and postprocess_gpu_campaign.py are unchanged.

| control | audit target CA RMSD Å | partner_fit Å | restored BB pairs <2 Å | restored BB min Å | full-native BB min Å | band |
|---|---:|---:|---:|---:|---:|---|
| assembly1_FcA_TfRD_nativeA | 3.709 | 3.709 | 0 | 12.985 | 2.672 | INCONCLUSIVE |
| assembly1_FcA_TfRD_nativeB | 3.708 | 3.708 | 0 | 12.991 | 2.588 | INCONCLUSIVE |
| assembly2_FcB_TfRC_nativeA | 0.963 | 0.963 | 0 | 10.717 | 2.884 | GATE VINDICATED |
| assembly2_FcB_TfRC_nativeB | 0.950 | 0.950 | 0 | 10.737 | 2.899 | GATE VINDICATED |

The audit minimum measures only receptor residues omitted from the construct; the smaller full-native minimum includes contacts against mapped receptor residues as well. Supplementary full-native minima use the existing min_dist helper on restored coordinates, not a replacement audit. The compatible controls also exceed2 Å against the entire native backbone. Neither copy was truncated, relaxed or locally refitted.

The Fc control is larger than the30–150 residue architecture allowed by canonicalize(); this test validates its requested rigid omitted-context clash diagnostic, not admission of an Fc as a minibinder. We did not bypass or relax canonicalize(). This is one experimental-control system, not calibration of generator success, pH selectivity, novelty, or every possible false rejection.

## Consequence for the66

Original66 have 51–330 backbone pairs below2 Å, with minimum separations 0.140–0.638 Å. Of these,65 have target fits<=2 Å, with54–330 pairs below2 Å. In comparison, compatible6W3H has zero and a10.7 Å omitted-context minimum (2.88–2.90 Å against all native backbone). The control therefore does not support reinstating the66. Their retirement stands; the one original>2 Å fit remains an inconclusive context transfer rather than a precise clash determination. This is structural evidence, not an experimental binding measurement.

## Release, spend and scope

All five functions on our Modal app had zero runners, inputs and backlog before the first reply. All3 current Runpod pods were EXITED; the prior local driver PID36142 was no longer running. Raw release evidence: RELEASE_VERIFIED.json and runpod_current_pods.json. Our generator lane is superseded and has no authorized successor. The separate approved RFD4 lane is untouched.

Our Modal conservative submission-to-release envelope is$10.718947; Runpod posted pod charges are$56.589988; combined$67.308935 (~$67.31). New control GPU spend$0.00. Other-agent RFD4 spend, consultations, CPU builds and persistent storage are excluded. Modal is an estimate, Runpod posted charges cover Oct19–current bucket across historical pods. SPEND.json and the corrected V18 spend checkpoint preserve the accounting.

Portal facts supplied by Ricardo reaffirmed in records: molecule_class=single_chain, once every24hours, automatic novelty minimum3/4 after upload, Challenge4 closes2026-11-02 11:59UTC. Ricardo creates accounts and submits entries. No upload, account, submission, GPU start, main-branch write or KB write occurred.

Reproduce CPU-only: `.venv-rosetta/bin/python results/real_gate_control/run_control.py`.
