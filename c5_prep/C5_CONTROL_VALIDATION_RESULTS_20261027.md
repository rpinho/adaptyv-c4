# C5 positive-control validation results

Run: 2026-10-27T22:16:04Z
Pod: glsqvthcvcj65k (A40 48GB, EU-SE-1)
Pipeline: ColabFold AF2-Multimer v3, 5 models, 3 recycles, mmseqs2_uniref_env

## Verdict: FAIL

ipTM(G12D) = ipTM(WT) = 0.84. The pipeline cannot discriminate the known
binder from the known non-binder at the best-model level.

## Scores (all 5 models per state)

### G12D pMHC — should bind (KD 63 uM experimentally)

| rank | model   | ipTM | pTM  | pLDDT |
|------|---------|------|------|-------|
| 1    | model_2 | 0.84 | 0.86 | 91.7  |
| 2    | model_1 | 0.84 | 0.85 | 91.1  |
| 3    | model_4 | 0.83 | 0.84 | 91.2  |
| 4    | model_5 | 0.80 | 0.82 | 90.2  |
| 5    | model_3 | 0.78 | 0.81 | 90.5  |

**Mean ipTM: 0.818** | Wall: 449s

### WT pMHC — should NOT bind (NB experimentally)

| rank | model   | ipTM | pTM  | pLDDT |
|------|---------|------|------|-------|
| 1    | model_1 | 0.84 | 0.85 | 90.9  |
| 2    | model_2 | 0.83 | 0.85 | 91.6  |
| 3    | model_5 | 0.79 | 0.81 | 89.4  |
| 4    | model_3 | 0.77 | 0.80 | 90.3  |
| 5    | model_4 | 0.63 | 0.68 | 88.5  |

**Mean ipTM: 0.772** | Wall: 694s

### Peptide-free HLA+b2m — should NOT bind

| rank | model   | ipTM | pTM  | pLDDT |
|------|---------|------|------|-------|
| 1    | model_1 | 0.77 | 0.80 | 90.9  |
| 2    | model_2 | 0.74 | 0.77 | 90.6  |
| 3    | model_5 | 0.70 | 0.74 | 89.1  |
| 4    | model_4 | 0.68 | 0.73 | 90.5  |
| 5    | model_3 | 0.66 | 0.71 | 89.7  |

**Mean ipTM: 0.710** | Wall: 457s

## Analysis

### What the pipeline CAN distinguish

Peptide presence: G12D best ipTM (0.84) > pepfree best ipTM (0.77), delta 0.07.
The pipeline sees the TCR reading the peptide groove. This is expected -- it's
a large structural difference (10 aa present vs absent).

### What the pipeline CANNOT distinguish

G12D vs WT: both score ipTM=0.84 at best-model. The one-residue difference
(Asp at position 12 vs Gly) does not produce a detectable ipTM signal.

### Why (mechanism)

Poole et al. 2022 showed the JDI TCR's G12D selectivity is driven by
**indirect electrostatic interactions**, not by direct structural contacts.
Crystal structures of TCR bound to G12D vs WT pMHC show "few discernible
differences in TCR interactions." The selectivity comes from:
- The charge of Asp (G12D) vs the neutral Gly (WT) at position 12
- This charge difference propagates through an electrostatic network
- Visible in thermodynamics and MD simulations, NOT in static structures

AF2-Multimer predicts static structures. It cannot see indirect electrostatic
effects. This is a fundamental limitation for this target, not a pipeline bug.

### Mean-based signal (suggestive but unreliable)

Mean ipTM across all 5 models: G12D=0.818, WT=0.772, pepfree=0.710.
The ordering is correct (G12D > WT > pepfree). However:
- The delta (0.046) is driven largely by WT model_4 being an outlier at 0.63
- The top-2 models for both G12D and WT are essentially identical
- Using mean ipTM as a ranking criterion would be fragile and unjustified

### Implications for design

1. **AF2-Multimer ipTM alone cannot rank designs for G12D selectivity.**
   A design scoring 0.84 against G12D might score 0.84 against WT too.

2. **BUT: designed binders may use different discrimination mechanisms.**
   The TCR's selectivity works through an unusual indirect path. A designed
   binder could use direct hydrogen bonds to the Asp side chain, which AF2
   might detect. The FAIL on this control does not prove the pipeline is
   useless for ALL binders, only that it cannot see THIS mechanism.

3. **The pipeline IS useful for:** binding vs non-binding (ipTM difference
   for peptide-present vs peptide-free is clear), structural quality (pLDDT),
   and fold confidence (pTM).

## Options for scoring selectivity

### A. Physics-based rescoring (complement AF2)

Use AF2-Multimer for fold quality + binding, then rescore with a physics-based
method that captures electrostatics:

- **FoldX ΔΔG**: fast, captures electrostatic contributions. Available via
  the FoldX suite. Can score the G12D→WT mutation directly.
- **Rosetta interface ΔΔG**: more thorough, repacking + minimization around
  the mutation site. Requires PyRosetta.
- **APBS electrostatic calculations**: directly compute the electrostatic
  potential change from Asp→Gly. Targets the known mechanism.

### B. Composite scoring

Score = α·ipTM(G12D) - β·ipTM(WT) + γ·ΔΔG(FoldX, G12D→WT)

Use AF2 for fold quality and FoldX/Rosetta for selectivity.

### C. Alternative structure prediction

- **AF2-Multimer with reduced MSA depth**: force the model to rely on
  structural features rather than co-evolutionary signal. May increase
  sensitivity to point mutations.
- **ESMFold**: single-sequence structure prediction. No MSA dependency.
  May or may not be better for point-mutation discrimination.

### D. Accept and mitigate

Design binders that use direct contacts to the Asp side chain (position 12).
If the design mechanism is direct, AF2-Multimer may detect it even though
it cannot detect the TCR's indirect mechanism. Screen designs for direct
Asp contacts in the predicted structures.

## Decision needed

Which scoring strategy to pursue. This is a strategic decision per the C5
constraints table.

## Files

- `control_results/control_validation_summary.json` — machine-readable verdict
- `control_results/scores_g12d.json` — full G12D per-model scores
- `control_results/scores_wt.json` — full WT per-model scores
- `control_results/scores_pepfree.json` — full pepfree per-model scores
- `control_results/run.log` — runner stdout
