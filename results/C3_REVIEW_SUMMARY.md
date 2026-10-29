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
4. **Fable review** (v4): 2 designs dropped, 2 replaced from RFD4 pool

## Results: 9 of 10 pass (BindCraft)

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
| 7 | l52_fff0d0 | 52 | 0.81 | 0.89 | **68.4** | 0.73 | 4 | -5.8 | **DROPPED v4** — MARGINAL refold |

**Dropped**: Rank #4 (l53_aaac96) — 13 all-atom clashes

## V2: Point-mutation variants (Codex analysis)

9 interface-edge mutations (E→Q, Q→K, F→Y) on top 3 parents. All pass ESMFold.

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

## RFD4 designs (3 selected after Fable review)

200 de novo binders generated via RFD4-Proteina on Modal. Pool screened for TM-like artifacts, charged fraction, diversity.

| Design | Len | Charge | GRAVY | Charged% | Refold pLDDT | Status |
|--------|-----|--------|-------|----------|-------------|--------|
| c3_glp1r_ecd_22 | 65 | -1 | -0.65 | 26% | 89.6 | PASS |
| c3_glp1r_ecd_172 | 67 | 0 | -0.47 | 28% | 89.3 | PASS (v4 NEW) |
| c3_glp1r_ecd_171 | 68 | -4 | -0.49 | 41% | 88.7 | PASS (v4 NEW) |

Dropped: ecd_96 (TM-like: GRAVY +0.72, 55.6% hydrophobic, max 19-mer window +2.32, 11% charged — will partition to membranes), ecd_61 (MARGINAL pLDDT 64.7), ecd_94 (MARGINAL pLDDT 61.5).

## V4 Fable assessment

**Strengths**: No poly-AA ≥4, no exact ≥6-res repeats, zero Cys, good E. coli expression predicted.
**Weaknesses**: (1) Effective diversity = 8 independent backbones (point mutants correlated); (2) Missing Trp/Tyr-enriched groove complementarity; (3) No agonism potential — all ECD binders without TMD-reaching extension.

## V4 → V5: aromatic complementarity upgrade

Fable flagged missing Trp/Tyr for W33/W39/Y69 groove stacking, and recommended trading redundant l100 variants for fresh backbones. Dropped 2 of 3 single l100 variants (kept double E79Q_E93Q), added 2 WY-rich RFD4 designs:

| Slot | V4 design | WY% | V5 replacement | WY% | pLDDT | Charge |
|------|-----------|-----|----------------|-----|-------|--------|
| 9 | c3_v_l100_E93Q (single variant) | 5.0% | **c3_glp1r_ecd_29** | **9.4%** | **93.5** | −3 |
| 10 | c3_v_l100_E79Q (single variant) | 5.0% | **c3_glp1r_ecd_12** | **12.0%** | **84.2** | −4 |

ecd_12 has 2W+8Y (12% aromatic) — the highest groove-complementarity in the set.
Mean portfolio WY improved from 5.5% → 6.3%. Independent backbones: 10 (up from 8).

## V5 → V6: replace low-WY ecd_22

ecd_22 had only 1.5% WY (1 Tyr in 65 res) — poor aromatic complement for the W33/W39/Y69-lined groove. Replaced with ecd_151 (11.5% WY, 2W+4Y, 52 res, pLDDT 80.3, charge -5).

| Slot | V5 design | WY% | V6 replacement | WY% | pLDDT | Charge |
|------|-----------|-----|----------------|-----|-------|--------|
| 16 | c3_glp1r_ecd_22 | **1.5%** | **c3_glp1r_ecd_151** | **11.5%** | **80.3** | −5 |

Mean portfolio WY improved from 6.3% → 6.7%. Now 3 designs above 9% WY.
6 RFD4 designs (up from 5), still 10 independent backbones.

## Submission files

- **`results/c3_submission_v6_20261023.csv`** — 20 designs (8 BindCraft + 7 variants + 6 RFD4)
- `results/c3_methodology.txt` — for Proteinbase submission form
