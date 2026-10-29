# C5 positive control — TCR from 7OW6

Pre-GPU gate step 3: validate that a differential scoring readout can recover a known
binder's preference for G12D over wild-type pMHC.

## The control: JDI TCR (PDB 7OW6, chains D+E)

The TCR co-crystallised in 7OW6 is from Poole et al., "Therapeutic high affinity T cell
receptor targeting a KRASG12D cancer neoantigen," *Nature Communications* 13, 5333 (2022).
DOI: 10.1038/s41467-022-32811-1. PMC9464187.

Source: Immunocore. The TCR was isolated from an HLA-A*11:01+ healthy donor and is the
parent of the affinity-matured IMC-KRASG12D bispecific.

### Published SPR selectivity (parent JDI TCR)

| pMHC | KD | source |
|---|---|---|
| HLA-A*11:01 + KRAS G12D (VVVGADGVGK) | **63 µM** | Poole et al. Table 1 / Fig 1c |
| HLA-A*11:01 + KRAS WT (VVVGAGGVGK) | **NB** (no binding) | Poole et al. SPR panel |
| HLA-A*11:01 + KRAS G12C | NB | Poole et al. SPR panel |
| HLA-A*11:01 + ERAS | NB | Poole et al. SPR panel |
| HLA-A*11:01 + RRAS | NB | Poole et al. SPR panel |
| HLA-A*11:01 + MRAS | NB | Poole et al. SPR panel |
| HLA-A*11:01 + DIRAS2 | NB | Poole et al. SPR panel |
| HLA-A*11:01 + RAP1A | NB | Poole et al. SPR panel |
| HLA-A*11:01 + self-peptide pool | NB | Poole et al. SPR panel |

**This is an ideal positive control.** The TCR:
- Binds G12D with measurable affinity (63 µM)
- Shows NO binding to wild-type (the one-residue-different comparator)
- Shows NO binding to peptide-free HLA (not explicitly listed; HLA+b2m alone is target 3)
- Was co-crystallised with the exact complex we're designing against

### What the control validation must show

A scoring pipeline that cannot recover the following ranking is broken:

1. TCR + state1 (G12D complex): **should score as binding** (KD 63 µM)
2. TCR + state2 (WT complex): **should score as non-binding** (NB experimentally)
3. TCR + state3 (peptide-free): **should score as non-binding** (expected — TCR reads
   the peptide in the groove)

The minimum passing criterion: the predicted affinity or interaction score for
TCR+state1 must be **unambiguously better** than TCR+state2 and TCR+state3, with
state2 and state3 comparable. A pipeline that ranks state2 ≥ state1 cannot discriminate
and will not rank designs either.

### Mechanism note (for the design session)

Poole et al. report that while the crystal structures of TCR bound to G12D vs WT pMHC
show "few discernible differences in TCR interactions," the specificity is driven by
**indirect electrostatic interactions** revealed by thermodynamic analysis and MD
simulations. This means:
- Static scoring (docking scores, interface energy) may underestimate the discrimination
- The Asp (G12D) vs Gly (WT) charge difference is the key, but its effect is indirect
- Dynamics-aware methods or electrostatic calculations may be needed to recover the signal

This is a harder discrimination problem than it looks from the structure alone.

## Files

The TCR chains are in `targets/c5_pmhc/control_tcr_only.pdb`:
- Chain D: TCR alpha, 181 residues (TRAV19/TRAJ28)
- Chain E: TCR beta, 245 residues (TRBV6/TRBJ1)

To reconstruct the TCR-pMHC complex for scoring, combine the TCR chains (D, E) with
each of the three target states (chains A, B, C from state1/state2; chains A, B from
state3). The chain IDs are already distinct and non-overlapping.
