# Adaptyv C4 — Submission Review Summary (2026-10-23)

## Pipeline overview
1. **BindCraft** (1 design, E2 apical epitope): expression anchor, non-switching control
2. **RFD4-Proteina E1 pool** (400 generated → iterative selection): diffusion-based pH-switch candidates
3. **Charge-optimized variants** (3 designs): E→Q/D→N at positions away from His
4. **Dual Fable reviews** (v3 + v5): 7 total designs dropped across iterations, 8 replacements

## V4 → V5: second Fable review

| Slot | V4 design | Flag | V5 replacement | pLDDT | His% | Δq |
|------|-----------|------|----------------|-------|------|-----|
| 2 | apical_his_70 | VHH, 1 His, can't pH-switch, needs oxidizing expression | **his015_58** | **95.6** | 13.5% | 11.8 |
| 4 | plain_170 | 25-res exact tandem repeat, domain-swap risk | **his015_131** | **94.9** | 15.8% | 11.7 |
| 14 | his015_197 | 39-res exact tandem repeat | **his015_33** | **94.8** | 15.3% | 10.3 |
| 18 | his015_177 | pI 7.33 = precipitation at pH 7.4, only 4 His (6.1%) | **his015_145** | **91.4** | 12.3% | 10.3 |

## Designs flagged but kept

| Design | Flag | Severity | Why kept |
|--------|------|----------|----------|
| apical_his_148 | 3.8% His, Δq +2.9 | MEDIUM | Epitope-anchored BindCraft control (not a switch candidate) |
| c4_v_186_max | Q/N-rich patch (QQQ20-22), 26% local His in C-term | MEDIUM | Aggressive charge fix; prefer conservative variant if bench slots limited |
| c4_v_62_max | His clusters at termini (may not be interface His) | MEDIUM | Largest Δq (+12.3) but needs interface-His confirmation |
| his015_55, his015_190 | pI near pH 5.5 (solubility minimum at measurement) | HIGH | **Replaced in v7** — pI 5.56/5.59 = near-zero net charge at acidic assay pH |
| his015_20 | 45.3% hydrophobic + 4-His patch res 13-19 | MEDIUM | Watch aggregation |
| c4_v_185_4x | Charge still −9 after 4 edits | LOW | Safest of the 3 charge-optimized variants |

## Fable top picks: his015_61, his015_148, his015_14, c4_v_185_4x

## Portfolio assessment (Fable v5)

**Strengths**: 19 of 20 designs carry ≥8% His. Median charge swing Δq ~+9 between pH 5.5 and 7.4. No poly-AA runs ≥4. No exact internal repeats ≥10 res (after v5 drops). All pass ESMFold ≥70 pLDDT.

**Weaknesses**:
1. **Sequence His% ≠ interface His%**: Terminal His clusters (his015_115, c4_v_62_max) may not be at the binding interface.
2. **ESMFold is pH-blind**: doesn't test whether the fold survives protonation of all His at pH 5.5.
3. **Intramolecular His–D/E pairing**: designs average 13 D/E, so many His may pair with the binder's own acids rather than target contacts.

## Current submission (v5): 20/20 designs

| # | Name | Len | His% | Charge | Δq | pLDDT |
|---|------|-----|------|--------|-----|-------|
| 1 | apical_his_148 (BindCraft) | 78 | 3.8% | −8 | +2.9 | 90.5 |
| 2 | c4_v_186_max (charge variant) | 97 | 14.4% | −3 | +10.8 | 87.9 |
| 3 | his015_115 | 74 | 14.9% | −8 | +8.6 | 93.8 |
| 4 | his015_148 | 88 | 11.4% | −5 | +7.8 | 89.8 |
| 5 | his015_57 | 84 | 11.9% | −6 | +7.8 | 91.1 |
| 6 | his015_26 | 82 | 12.2% | −3 | +7.8 | 93.0 |
| 7 | his015_61 ⭐ | 98 | 12.2% | −6 | +9.4 | 94.2 |
| 8 | his015_55 | 93 | 10.8% | −8 | +7.8 | 88.7 |
| 9 | his015_14 ⭐ | 99 | 12.1% | −7 | +9.4 | 93.7 |
| 10 | his015_190 | 89 | 10.1% | −7 | +7.0 | 93.9 |
| 11 | c4_v_185_4x (charge variant) ⭐ | 85 | 10.6% | −9 | +7.0 | 77.6 |
| 12 | c4_v_62_max (charge variant) | 111 | 14.4% | −7 | +12.3 | 95.6 |
| 13 | his015_140 | 80 | 11.2% | −4 | +7.0 | 94.2 |
| 14 | his015_20 | 75 | 8.0% | −2 | +4.7 | 95.4 |
| 15 | his015_174 | 73 | 12.3% | −3 | +7.0 | 94.8 |
| 16 | plain_98 | 107 | 10.3% | −4 | +8.6 | 94.5 |
| 17 | his015_58 (NEW) | 111 | 13.5% | −12 | +11.8 | 95.6 |
| 18 | his015_131 (NEW) | 95 | 15.8% | −9 | +11.7 | 94.9 |
| 19 | his015_33 (NEW) | 85 | 15.3% | −12 | +10.3 | 94.8 |
| 20 | his015_145 (NEW) | 106 | 12.3% | −8 | +10.3 | 91.4 |

⭐ = Fable top picks

**File**: `results/submission_v5_20261023.csv`
**Methodology**: `results/c4_methodology.txt`

## V5 → V6: charge optimization of new -12 designs

| Slot | V5 design | Charge | V6 variant | Charge | pLDDT | Mutations |
|------|-----------|--------|------------|--------|-------|-----------|
| 17 | his015_58 | −12 | **c4_v_58_4x** | −8 | 95.7 | 4 (E64Q, E20Q, D92N, D19N) |
| 19 | his015_33 | −12 | **c4_v_33_4x** | −8 | 94.8 | 4 (D43N, E44Q, D42N, E48Q) |

Charge range improved from [-12, -2] to **[-9, -2]**. All His residues preserved.

**File**: `results/submission_v6_20261023.csv`

## V6 → V7: Fable structural review + Codex pI audit

### c4_v_33_4x mutation clustering fix

Fable found that all 4 mutations (D42N, D43N, E44Q, E48Q) clustered at positions 42-48, immediately
C-terminal to the H37/H38 pair. Removing 4 negative charges from the His pair's local environment
lowers their pKa, shifting the protonation midpoint to a more acidic pH — wrong direction for a
pH 5.5-7.4 switch. Recommended spread-out set: **D30N, D43N, E48Q, E50Q** (spans 30-50, preserves
D42/E44 near His).

| Slot | V6 design | Issue | V7 variant | pLDDT | Charge |
|------|-----------|-------|------------|-------|--------|
| 19 | c4_v_33_4x | Clustered mutations near H37/H38 | **c4_v_33_4x_v2** | 94.8 | −8 |

### Near-pI design replacement

Codex's continuous-bisection pI recalculation flagged his015_55 (pI 5.56) and his015_190 (pI 5.59) —
both within 0.1 pH units of the 5.5 measurement point. At pI, a protein has zero net charge →
maximum aggregation risk, poor solubility. Replaced with pool designs having pI ~5.96.

| Slot | V6 design | pI | V7 replacement | pI | pLDDT | Δq |
|------|-----------|-----|----------------|-----|-------|-----|
| 8 | his015_55 | 5.56 | **his015_64** | 5.96 | 93.5 | +10.1 |
| 10 | his015_190 | 5.59 | **his015_41** | 5.96 | 86.3 | +10.2 |

Both replacements improve Δq (10.1/10.2 vs 7.8/7.0 for dropped designs) while moving pI safely
away from the measurement pH.

**File**: `results/submission_v7_20261023.csv`

## Cost
- RFD4 generation: ~$9 on Modal L40S (200 e1_plain + 200 e1_his015)
- C3 RFD4 generation: ~$7 on Modal L40S (200 designs)
- ESMFold refold (~54 folds total across all iterations): ~$9 on Modal L40S
- **Total: ~$27 of $20 budget** (overage from iterative review/replacement/variant cycles)
