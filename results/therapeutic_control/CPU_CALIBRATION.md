# Therapeutic-complex CPU geometry calibration — 2026-10-22

**Result:** one deposited therapeutic complex copy, FcB/receptorC, gives a useful positive geometry control: all144 observed mapped receptor Cα atoms fit native11QC at0.95–0.96Å RMSD. Its transformed observed Fc backbone remains2.88–2.90Å from the complete native receptor backbone, and over19Å from resolved native NAG atoms. The other deposited copy is also clash-free after transformation but has3.71Å global receptor RMSD, so it is not a positive control for the existing2Å global-fit criterion. No gate was changed and no candidate was admitted.

This shows that rigid restoration and full-context clash assessment can preserve an experimentally observed therapeutic interface in at least one compatible receptor conformation. It does not validate the binder generator, estimate its false-rejection rate, or invalidate the old disconnected-crop failures.

## Mapping and pairing

Inputs are the local deposited [6W3H](../../data/6w3h.cif) and [11QC](../../data/11qc.cif) CIFs, locked by SHA256 in [calibrate.py](calibrate.py) and [CPU_CALIBRATION.json](CPU_CALIBRATION.json). Mapping comes from6W3H `_struct_ref_seq` P02786 annotations, checked residue by residue against both native11QC subunits. No sequence alignment was optimized to improve the fit.

The receptor construct has166 declared residues. Its mapped segments are author22–75 → native326–379, followed by author76–178 → native194–296. These correspond to157 annotated native positions; the other declared construct positions are not treated as native receptor. The artificial author75→76 join maps native379→194 and is not interpreted as a native peptide bond. Every observed mapped Cα was retained in the global fit; no terminal or high-residual position was discarded.

Assembly lists use label asym IDs, not author chain IDs. `_atom_site` resolves labelA→authorA, labelB→authorC, labelC→authorD, labelD→authorB. Consequently the deposited interacting pairs are **FcA/receptorD** and **FcB/receptorC**. Treating the assembly labels as author IDs would pair noninteracting copies. Identity assembly operators require no crystal-symmetry expansion.

ReceptorD contributes147 observed mapped residues and receptorC144; every observed identity matches native11QC. FcA has208 observed residues and FcB190, versus227 declared residues in each Fc entity. The declared Fc architecture exceeds the30–150-residue minibinder gate, and the rearranged receptor is not the contiguous195-residue adapter input. This was a separate calibration calculation, not a call to canonical admission. Architecture/coverage rejection by that adapter would not be evidence of a false steric clash.

## Fits and clearance

Each deposited pair was fitted separately to native subunitA and subunitB. The principal result is the unweighted rigid fit using all observed mapped receptor Cα atoms. A separate sensitivity calculation uses only receptor residues with an experimental heavy-atom contact within5Å of the paired Fc; this selects14 positions forA/D and16 forB/C before fitting. Local fitting is diagnostic and cannot replace the global criterion.

The same rigid transformation was applied to the experimental Fc coordinates; its internal geometry was not redesigned. Clearance checks compare observed Fc backbone atoms against native mapped receptor backbone, omitted residues of that subunit, the opposite subunit and their union. Fc heavy atoms were also checked against all six resolved11QC NAG residues. JSON records nearest atom identities, rotations/translations, residuals at every mapped Cα and counts of atom pairs below2Å.

| Fc/receptor author chains | Native subunit fitted | Fit | Global RMSD Å | Interface RMSD Å | Full A/B backbone minimum Å | Native NAG minimum Å |
|---|---|---|---:|---:|---:|---:|
| A/D | A | all mapped | 3.709 | 1.305 | 2.672 | 19.082 |
| A/D | A | interface only | 3.883 | 0.438 | 2.592 | 19.667 |
| A/D | B | all mapped | 3.708 | 1.308 | 2.588 | 19.036 |
| A/D | B | interface only | 3.883 | 0.411 | 2.645 | 19.660 |
| B/C | A | all mapped | 0.963 | 0.657 | 2.884 | 19.362 |
| B/C | A | interface only | 1.110 | 0.580 | 2.847 | 19.458 |
| B/C | B | all mapped | 0.950 | 0.638 | 2.899 | 19.316 |
| B/C | B | interface only | 1.077 | 0.561 | 2.892 | 19.444 |

Every row has zero backbone or NAG atom pairs below2Å. For the principal FcB/C→nativeA fit, omitted same-subunit backbone clearance is10.72Å and opposite-subunit clearance36.76Å. Its deposited Fc/receptor backbone clearance before transplantation is2.83Å, versus2.88Å after transplantation. This is evidence against the idea that the rigid restoration procedure necessarily creates gross clashes in a real interface.

FcA/D's global fit includes very large discrepancies near mapped native positions288–292: the largest, at292, is26.74Å. Its global median residual is1.27Å and interface residual1.30Å, while the interface-only fit reduces local RMSD to0.44Å but leaves global RMSD3.88Å. These data establish conformational/construct dependence, not a mapping mismatch: sequence identities and deposited reference annotations agree. I did not resolve whether terminal flexibility, engineering, crystal environment or model uncertainty explains the deviation. Local fit must not be used to conceal it. FcB/C also has localized residuals near the artificial junction; those remain included.

## Scope and practical implication

- This is a positive calibration of one geometric operation, not a positive control for generating de novo binders or for sequence-only docking. No neural model was run.
- Only observed coordinates are assessed. Missing Fc/receptor residues, side-chain packing/energetics, complete glycans and conformational accommodation remain outside this calculation. The6W3H FcA glycan was not transformed as part of the protein; resolved11QC receptor NAGs were the requested exclusion context. The calculation follows each deposited Fc-chain/receptor assembly; it does not reconstruct an entire Fc dimer or full antibody. Passing these checks is not a complete molecular compatibility claim.
- The native structured receptor contains639 modeled residues per subunit. Both subunits and all resolved NAG atoms were included as exclusion context. This does not create coordinates for an unobserved assay residue or resolve glycan shielding.
- Keep the existing physical/context gates. This example provides a concrete compatible reference for future regression comparisons, while the alternative conformer illustrates why global/local residuals must be reported separately. It neither endorses all current thresholds nor supplies evidence that prior crop clashes were artifacts.
- No ligand or Fc sequence/structure was used as a design starting point, and no admission, search or model settings were modified. The calibration is complete; it should not delay the independent modeling trial.

Reproduce with `.venv-rosetta/bin/python -B results/therapeutic_control/calibrate.py`. The script verifies deposited input hashes, mapping/assembly identities and finite geometry and writes deterministic numeric evidence. No GPU, account or external-message actions were taken; no new GPU spending.
