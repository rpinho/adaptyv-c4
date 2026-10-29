# STATUS — read before using anything in this directory

Written 2026-10-27 ~16:11 EDT. **Validated 2026-10-27 ~20:30 EDT.**

These are chain selections and one side-chain truncation derived from **PDB 7OW6** (public;
`curl https://files.rcsb.org/download/7OW6.cif`, source file sha256
`8ce7adf9495d3d8c7d79da690bba9be08f1aee83ba63b409e8b2a0cb35957b0e`). Anyone can regenerate
them from that entry.

## Validation results (all passing)

**Critical assertion 1: State1 vs State2 differ ONLY at peptide position 6.**
- Chains A+B: 3,079 atoms, ALL coordinates identical between states. PASS.
- Chain C: differs at exactly one position — residue 6, ASP→GLY (D→G). PASS.
  - State1 (G12D): VVVGADGVGK — 63 chain-C atoms (ASP has 8 atoms at pos 6: N,CA,C,O,CB,CG,OD1,OD2)
  - State2 (WT):   VVVGAGGVGK — 59 chain-C atoms (GLY has 4 atoms at pos 6: N,CA,C,O)
- No other differences at any level. PASS.

**Critical assertion 2: State1 vs State3 differ ONLY by absence of chain C (peptide).**
- State3 has chains A+B only — chain C absent. PASS.
- Chains A+B: 3,079 atoms, ALL coordinates identical to state1. PASS.

**Control (TCR-only):** chains D (181 res) and E (245 res) — the TCR alpha and beta chains
from 7OW6, separated from the MHC complex.

## Sequence comparison with competition page

| chain | PDB | competition | match? |
|---|---|---|---|
| A (heavy) | 276 aa | 281 aa | First 275 identical; pos 276 PDB=P vs comp=L; comp has 5 more C-terminal residues (SSQPT) |
| B (β2m) | 99 aa | 99 aa | **EXACT MATCH** |
| C (G12D peptide) | VVVGADGVGK | VVVGADGVGK | **EXACT MATCH** |
| C (WT peptide) | VVVGAGGVGK | VVVGAGGVGK | **EXACT MATCH** |

### Heavy chain detail
PDB auth_seq_id 1–276. Residues 1–275 match the competition sequence exactly. Residue 276:
PDB has PRO, competition has LEU. The competition sequence continues for 5 more residues
(positions 276–281: LSSQPT). This is the α3 domain C-terminal tail approaching the
transmembrane anchor — far from the peptide-binding groove. The PRO at position 276 in the
crystal is likely an expression-construct artifact (not the natural sequence).

### β2m — earlier STATUS.md was WRONG
The earlier version of this file claimed a "100 aa with initiator Met" offset. **That was
incorrect.** The PDB chain B is 99 residues and matches the competition's 99-aa mature β2m
sequence exactly. No offset.

## What has NOT been done

- No geometry check, no clash check, no protonation, no minimisation.
- No positive-control validation (TCR selectivity readout — step 3 of the pre-GPU gate).
- These are crystal coordinates, not modelled. The missing 5 C-terminal heavy-chain residues
  are not present and the P→L at position 276 is not corrected. The executing session must
  decide whether to use as-is (the groove region is complete) or model the missing tail.

## File inventory

| file | atoms | chains | description |
|---|---|---|---|
| state1_g12d_complex.pdb | 3,142 | A(276), B(99), C(10) | Target 1: HLA + β2m + G12D peptide |
| state2_wt_complex.pdb | 3,138 | A(276), B(99), C(10) | Target 2: HLA + β2m + WT peptide |
| state3_peptide_free.pdb | 3,079 | A(276), B(99) | Target 3: HLA + β2m, no peptide |
| control_tcr_only.pdb | 3,366 | D(181), E(245) | TCR alpha+beta from 7OW6 (positive control) |
