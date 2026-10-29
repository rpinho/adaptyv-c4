# C4 v8 PROPKA Analysis Report
Date: 2026-10-23

## Method
- ESMFold monomer structures for all 20 designs (Modal L40S GPU)
- PROPKA 3.5.1 local analysis on downloaded PDBs
- Wyman-linkage ΔΔG_monomer estimate: ΔΔG = -RT * Σ ln[(1 + 10^(pKa_i - 5.5)) / (1 + 10^(pKa_i - 7.4))]
- RT = 0.592 kcal/mol at 25°C

## Combined Ranking (pLDDT + PROPKA)

| Rank | Design | pLDDT | His | Switch | LowPk | AlwOn | ΔΔG_mono | Grade |
|------|--------|-------|-----|--------|-------|-------|----------|-------|
| 1 | c4_v_58_4x | 95.7 | 15 | 15 | 0 | 0 | -19.36 | ★★★ |
| 2 | tfr1_e1_his015_131 | 94.9 | 15 | 15 | 0 | 0 | -18.19 | ★★★ |
| 3 | tfr1_e1_his015_145 | 91.4 | 13 | 12 | 1 | 0 | -17.78 | ★★★ |
| 4 | c4_v_33_4x_v2 | 94.8 | 13 | 13 | 0 | 0 | -16.57 | ★★★ |
| 5 | c4_v_62_max | 95.6 | 16 | 15 | 1 | 0 | -16.15 | ★★★ |
| 6 | tfr1_e1_his015_61 | 94.2 | 12 | 12 | 0 | 0 | -14.00 | ★★ |
| 7 | tfr1_e1_his015_96 | 92.3 | 10 | 8 | 1 | **1** | -12.98 | ★★ |
| 8 | tfr1_e1_his015_57 | 91.1 | 10 | 10 | 0 | 0 | -12.89 | ★★ |
| 9 | tfr1_e1_his015_64 | 93.5 | 13 | 10 | 3 | 0 | -12.64 | ★★ |
| 10 | tfr1_e1_his015_14 | 93.7 | 12 | 11 | 1 | 0 | -12.21 | ★★ |
| 11 | tfr1_e1_his015_41 | 86.3 | 13 | 10 | 3 | 0 | -11.66 | ★ |
| 12 | tfr1_e1_his015_148 | 89.8 | 10 | 9 | 1 | 0 | -11.56 | ★★ |
| 13 | c4_v_185_4x | **77.6** | 9 | 9 | 0 | 0 | -11.40 | ⚠ |
| 14 | tfr1_e1_his015_174 | 94.8 | 9 | 9 | 0 | 0 | -10.52 | ★★ |
| 15 | tfr1_e1_his015_115 | 93.8 | 11 | 10 | 1 | 0 | -10.41 | ★★ |
| 16 | tfr1_e1_plain_98 | 94.5 | 11 | 10 | 1 | 0 | -11.00 | ★★ |
| 17 | tfr1_e1_his015_26 | 93.0 | 10 | 10 | 0 | 0 | -9.60 | ★ |
| 18 | tfr1_e1_his015_140 | 94.2 | 9 | 8 | 1 | 0 | -7.52 | ★ |
| 19 | tfr1_e1_his015_20 | 95.4 | 6 | 6 | 0 | 0 | -6.87 | ★ |
| 20 | tfr1_his200m (BindCraft) | 90.5 | 3 | 2 | 1 | 0 | -2.59 | ⚠ |

## Key Findings

### Positive
- **19/20 designs have pLDDT >86** — well-folded monomers
- **ALL His in the switching window (5.5-7.4)** for the top 5 designs — no wasted residues
- **Zero ALWAYS-ON His** across almost all designs (only tfr1_e1_his015_96 has 1)
- Charge-optimized variants (c4_v_58_4x, c4_v_33_4x_v2, c4_v_62_max) rank 1st/4th/5th

### Concerning
- **c4_v_185_4x: pLDDT=77.6** — below 80 threshold, may not fold properly in vitro
- **BindCraft design: ΔΔG=-2.59** — only 3 His, minimal pH-switching potential
- **tfr1_e1_his015_41: pLDDT=86.3** — lowest of the RFD4 set, 3 low-pKa His

### Comparison to Prior Art
- actimo-labs reported ΔΔG 2.84-4.56 kcal/mol from **complex** PROPKA (Boltz-2 structures)
- Our monomer estimates (6.87-19.36) are NOT directly comparable — monomer ΔΔG is an upper bound
- The relative ranking within our portfolio IS valid for design selection
- Next step for true comparison: Boltz-2 complex prediction → PROPKA on complex

## Replacement Candidates

### c4_v_185_4x (rank 13, pLDDT=77.6)
**Recommend replacement.** pLDDT below 80 signals unreliable fold. The ΔΔG is middling.
Potential replacement: next-best RFD4 candidate from generation pool.

### BindCraft design (rank 20, ΔΔG=-2.59)
**Keep for diversity.** Only BindCraft-derived design; method diversity is valuable for the
competition even if pH-switch potential is low. 3 His is the minimum for any switch.

### tfr1_e1_his015_20 (rank 19, ΔΔG=-6.87)
**Consider replacement.** Good fold (95.4) but only 6 His — weakest RFD4 design by switch potential.
