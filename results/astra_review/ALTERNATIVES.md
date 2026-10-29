# Alternatives after full-target acquisition failures

Reviewed 2026-10-22. This is groundwork, not a pipeline switch. Only this review file was written; no GPU, worker, account, submission or external message was created. V11 was already running under the parent agent. Parent's last completed-runtime estimates: Runpod $1.7952; Modal $5.2184357 before V11; extras separate. This review adds no GPU spend and does not reconcile consultation billing.

## Verdict

A **native, contiguous apical-domain search representation** is scientifically more defensible than the old disconnected spatial crop, but remains an approximation requiring a separately tested restoration and validation path. It is a reasonable conditional alternative if V11 also fails to acquire the patch. It is not yet ready to run through the existing admission code.

Published therapeutic antibody complexes offer useful positive controls and epitope evidence. **They cannot be modified into competition candidates under the official rules.** There is no evidence here that changing to an antibody format would solve acquisition, cross-species recognition or pH switching.

## Rules: controls are allowed; existing-binder starting points are not

The official FAQ8 requires de novo/zero-shot designs, prohibits starting from or modifying an existing binder, and requires sequence/structural diversity (CDR diversity for antibodies). It expressly permits known binders for evaluating/calibrating filtering metrics or training models. FAQ1 permits common antibody frameworks as starting points. This is permission for a generic framework with newly designed binding regions, not for retaining a known anti-TfR paratope. [Official challenge and FAQ](https://proteinbase.com/competitions/anthropic-adaptyv-2026/challenges/tfr1).

Repository `DECISIONS.md:291–293` additionally says to use native full receptor context and not transplant known binder sequences or scaffolds. Thus the current project should use therapeutic complexes as **controls only**. A generic antibody-framework route would be a separate method decision, not a silent extension of the current minibinder run. No numerical novelty threshold was verified here; do not invent one.

## Native contiguous domain: actual CPU observations

The apical domain is conventionally residues189–383; it is a continuous insertion within the protease-like domain. Interdomain contacts and a calcium site connect it to the rest of the receptor, and its long arm can change conformation with transferrin binding. Domain boundaries do not establish independent native stability. [Primary structural/mutagenesis study](https://pmc.ncbi.nlm.nih.gov/articles/PMC3302921/).

I inspected `data/11qc.cif` (SHA256 `4eaa5db67ed83a896713c4b42d1dc214c778bdac240fdf3a9a82e86a66f79495`) and the existing `results/gpu_campaign/human_11qc_full_dimer.pdb` with Bio.PDB/numpy, without generating a new target:

| Check | Observation |
|---|---|
| Native modeled chains A/B | Each has639 residues,122–760 |
| A189–383 |195 consecutive residues; no missing N/CA/C/O |
| Internal peptide C–N distances |1.322–1.349Å |
| Prepared-domain Cα versus original CIF |Identical for both A and B |
| Terminal189/383 Cα to nearest current205/206/208/209/210 hotspot Cα |34.85/35.51Å |
| Omitted same-chain residues with a heavy atom within4Å of the domain |27, including flanking residues and remote protease/helical regions |

The omitted-contact set is184–188,384,412,457–461,464–465,468–469,537–538,543,562,564,567–570,693,697. This is a geometric observation, not a list of residues to reconnect or mutate. It demonstrates that cropping exposes a new surface even when the internal backbone is continuous. The second subunit and unresolved glycan volume add further context not captured by this count.

The original73-residue crop failed because favorable poses occupied omitted receptor protein;65 reliable restorations all clashed. A contiguous domain removes fictitious internal peptide joins, but does not by itself prevent this omission artifact.

Published AP01/6Y76 is evidence that computational redesign can make an isolated apical domain soluble. It is **not** evidence that an unchanged native domain is autonomously stable, and its redesigned sequence/backbone must not replace the assay receptor. [Primary paper](https://doi.org/10.1002/prot.25974); [deposited structure](https://www.rcsb.org/structure/6Y76).

### What a defensible future domain route would require

1. Preserve the native contiguous sequence, original residue mapping and coordinates; audit connectivity, atom completeness, frame validity and hotspot mapping. Track it explicitly as a reduced search representation. Do not concatenate spatially separated pieces.
2. Before interpreting any candidate score, align the entire retained domain to native11QC, report both global-domain and local-patch fit, transform the binder, and restore **both complete native A/B chains**. Use all retained domain residues for fitting rather than choosing only favorable hotspot atoms. Reject severe domain deformation, newly exposed-face binding, and full-context backbone/NAG overlap.
3. Carry the restored binder into a **native full-target validation/refold** and apply unchanged quality and context gates. Rigid restoration plus absence of clashes is necessary but insufficient: it cannot prove that the sequence retains its pose when the omitted receptor is present. Domain confidence/loss values remain search diagnostics; do not promote them to native-complex confidence or affinity.
4. Perform pH-related scoring only after native full-context admission. Mouse binding remains a separate untested requirement. Modelled122–760 also lacks assay residue121; retain that existing limitation explicitly.

**Immediate software constraint:** `scripts/postprocess_gpu_campaign.py:canonicalize()` deliberately requires complete native639-residue A before restoring B; it will reject a domain input. A separate, tested partial-domain restoration adapter is needed upstream. Keep that strict admission check. Its existing≤2Å full-monomer fit cannot simply be renamed a full-context fit when only195 residues were observed. Include exact native self-restoration and a deliberately omitted-protein-clashing pose as positive/negative CPU fixtures before any domain GPU test.

## Noninfectious therapeutic structural controls

**Best verified deposited example: Denali TV35–TfR, PDB6W3H.** This is an engineered antibody Fc transport vehicle developed for therapeutic brain delivery, with a3.38Å experimental complex and a binding site away from transferrin/FcRn sites. It is a useful native-context restoration/epitope control, not a minibinder seed or an eligible candidate. [Primary publication](https://doi.org/10.1126/scitranslmed.aay1359); [PDB6W3H](https://www.rcsb.org/structure/6W3H).

Caution verified directly from the deposited mmCIF: the receptor entity has166 residues, and its reference mappings run native326–379 followed by194–296. This is a rearranged fragment, not native189–383 or a complete monomer. Map by sequence/reference annotation before superposition; never infer native numbering from its author residue IDs or the PDB page's “Mutation(s): No” label. Restoring native11QC around the experimental bound Fc would provide a concrete test of the context checker. This review did not execute that superposition or establish that it passes our strict clash thresholds.

**Useful complementary therapeutic evidence:** the TED-TfR study reports two experimentally studied antibody classes: one contacting apical plus protease-like domains, another contacting the apical domain alone, with different pH responses. It supports the importance of epitope/context and cautions against assuming all apical binders share a pH mechanism. Its acidic-condition observations do not establish the challenge's exact pH6.0 endpoint or mouse cross-reactivity. I did not verify a downloadable coordinate accession for these complexes, so they are literature controls rather than a ready structural benchmark. [Primary study](https://pmc.ncbi.nlm.nih.gov/articles/PMC12952256/).

No viral/pathogen receptor ligand is proposed or used as a design starting point in this review.

## Priority and falsifiable decisions

- **First, CPU control/calibration:** verify the restoration adapter on native self-restoration and an intentionally invalid pose, then assess whether the therapeutic experimental complex can be correctly mapped/restored. High information value, no new modeling GPU. A failed known-complex restoration must be diagnosed before candidate rejection thresholds are interpreted.
- **Conditional next modeling hypothesis:** a single bounded contiguous-domain search, only after the adapter/control passes and V11's native-restoration outcome is known. Its falsifiable endpoint is actual patch acquisition followed by successful native full-context validation at unchanged gates. Acquisition that vanishes on restoration/refolding is a failed approximation, not progress toward an admitted binder.
- **Defer antibody-format migration:** permissible de novo framework design is a larger change in model/architecture and calibration. Known-binder reuse is disallowed, and the current repository pipeline imposes30–150-residue single-binder architecture checks; a Fab/scFv would require deliberate support rather than bypassing those checks.

These checks cannot establish hit rates from the small number of trajectories, nor prove that target size caused the failures. A smaller domain changes receptor context, tensor size and network behavior simultaneously. Treat improved acquisition as evidence about this approximation, not a unique causal explanation for the full-target failures.
