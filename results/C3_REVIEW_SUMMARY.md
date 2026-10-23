# C3 (GLP1R Agonist) — Submission Review

**Challenge**: Design a GLP1R agonist (binding affinity + receptor activation)
**Deadline**: Oct 25, 23:59 AoE (~2 days)
**Scoring**: (1) lowest KD among all binders, (2) agonist confirmation in cAMP assay. A plain binder is scoreable — agonism is upside.

## Pipeline

1. **BindCraft** forced-targeting campaign on GLP1R ECD groove (PDB 6X18, chain A)
   - Hotspots: A67–A69, A88–A91, A121, A123, A128, A137–A145
   - Key groove residues: W33, W39, R40, E68, Y69, W91
   - Binder lengths 45–100, max 3000 trajectories, 20 final designs requested
   - 10 designs ranked
2. **ESMFold** independent refold (Modal L40S GPU, `facebook/esmfold_v1`)
3. **Quality gates**: refold FAIL → drop; >10 all-atom clashes → drop

## Results: 9 of 10 pass

| Rank | Design | Len | iPTM | pLDDT | Refold | Hotspot | Clashes | Charge | Note |
|------|--------|-----|------|-------|--------|---------|---------|--------|------|
| 3 | l57_fbcf76 | 57 | 0.85 | 0.93 | 87.2 | 0.53 | 1 | -1.9 | **TOP** — best iPTM+pLDDT, near-neutral |
| 5 | l100_c674f8 | 100 | 0.81 | 0.91 | 89.6 | 0.67 | 6 | -6.9 | Highest refold pLDDT |
| 2 | l65_8f7b63 | 65 | 0.78 | 0.87 | 84.8 | 0.60 | 2 | -3.9 | Clean, minimal clashes |
| 10 | l78_3d63fd | 78 | 0.78 | 0.92 | 87.3 | 0.73 | 5 | -7.9 | Strong refold, good hotspot |
| 1 | l70_7781fe | 70 | 0.78 | 0.87 | 84.0 | 0.67 | 10 | -2.9 | At clash limit, OK charge |
| 8 | l52_6950ff | 52 | 0.78 | 0.88 | 85.4 | 0.60 | 2 | -6.0 | Compact, clean |
| 6 | l99_e26aff | 99 | 0.79 | 0.88 | 82.1 | 0.67 | 4 | -4.9 | Solid all-round |
| 9 | l100_4bd789 | 100 | 0.76 | 0.83 | 85.6 | 0.80 | 4 | -9.9 | Lowest iPTM, extreme charge |
| 7 | l52_fff0d0 | 52 | 0.81 | 0.89 | **68.4** | 0.73 | 4 | -5.8 | MARGINAL refold (38% >70) |

**Dropped**: Rank #4 (l53_aaac96) — 13 all-atom clashes

## Flags

- **Rank #7** (l52_fff0d0): MARGINAL ESMFold refold (68.4, only 38.5% residues >pLDDT 70). BindCraft scores are decent (iPTM 0.81, pLDDT 0.89) but the independent fold doesn't confirm the structure well. Consider dropping.
- **Rank #9** (l100_4bd789): Lowest iPTM (0.76), lowest BindCraft pLDDT (0.83), extreme negative charge (-9.9). Still passes ESMFold (85.6). Codex analysis recommended dropping if constrained.
- All 10 designs contact all 6 key groove residues (W33, W39, R40, E68, Y69, W91).

## V2: Point-mutation variants (Codex analysis)

Codex did structural analysis (BioPython SASA, distance-to-groove) on the top 3 parents and proposed 9 interface-edge mutations — all E→Q or Q→K/F→Y to improve charge neutrality and contacts. All 9 variants pass ESMFold (pLDDT 84.2–89.6, matching parents).

| Variant | Parent | Mutations | Len | Refold pLDDT | Charge |
|---------|--------|-----------|-----|-------------|--------|
| c3_v_l57_Q18K | #3 | Q18K | 57 | 86.6 | -1 |
| c3_v_l57_F38Y | #3 | F38Y | 57 | 87.1 | -2 |
| c3_v_l57_Q18K_F38Y | #3 | Q18K+F38Y | 57 | 86.6 | -1 |
| c3_v_l100_E93Q | #5 | E93Q | 100 | 89.6 | -6 |
| c3_v_l100_E79Q | #5 | E79Q | 100 | 89.6 | -6 |
| c3_v_l100_E79Q_E93Q | #5 | E79Q+E93Q | 100 | 89.5 | -5 |
| c3_v_l65_E53Q | #2 | E53Q | 65 | 84.7 | -3 |
| c3_v_l65_E5Q | #2 | E5Q | 65 | 84.2 | -3 |
| c3_v_l65_E5Q_E53Q | #2 | E5Q+E53Q | 65 | 84.8 | -2 |

## RFD4 designs (2 selected)

200 de novo binders generated via RFD4-Proteina on Modal (PDB 6X18, chain F/label numbering, same groove hotspots). 104 extracted, filtered (no Cys, 45-100 residues), ranked by charge neutrality + length. Top picks validated via ESMFold:

| Design | Len | Net Charge | Refold pLDDT | Status |
|--------|-----|-----------|-------------|--------|
| c3_glp1r_ecd_96 | 72 | 0 | 83.8 | PASS |
| c3_glp1r_ecd_22 | 65 | -1 | 89.6 | PASS |

Rejected: c3_glp1r_ecd_61 (charge 0, pLDDT 64.7 MARGINAL), c3_glp1r_ecd_94 (pLDDT 61.5 MARGINAL).

## Submission files

- `results/c3_submission_v3_20261023.csv` — **20 designs** (9 BindCraft + 9 variants + 2 RFD4)
- `results/c3_methodology.txt` — for Proteinbase submission form
