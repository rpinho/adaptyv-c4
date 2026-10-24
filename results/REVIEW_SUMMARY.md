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
| c4_v_186_max | QQQYQQ prion-like, 15.5% Q+N, net -3 → strong positive swing | HIGH | **Replaced in v8** — prion-like Q/N signature from aggressive 14-mutation charge optimization |
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

## V7 → V8: drop prion-like c4_v_186_max

Fable's charge variant risk ranking placed c4_v_186_max last (worst) of 5 variants:
QQQYQQ at positions 20-25 is the classic prion-like / self-templating low-complexity signature
(58% Q/N in worst 12-mer window), net -3 with strong positive swing → aggregation and non-specific
binding risk on acidification. Its Δq (+10.9) and His% (14.4%) are both covered by 4+ other designs.

| Slot | V7 design | Issue | V8 replacement | pLDDT | Δq | His% |
|------|-----------|-------|----------------|-------|-----|------|
| 2 | c4_v_186_max | QQQYQQ prion-like, 15.5% Q+N | **his015_96** | 92.3 | +8.0 | 9.3% |

Pool candidate selection: his015_85 (best Δq in pool, +10.7) excluded for VVVV poly-Val; his015_126
excluded for SSSS poly-Ser; his015_96 passed all filters with excellent fold quality.

**File**: `results/submission_v8_20261023.csv`

## His distribution analysis (Codex, v8)

Mid-chain His (middle 60% of sequence) as proxy for interface-contact potential:

| Rating | Count | Best examples |
|--------|-------|---------------|
| GOOD (>60% mid-chain) | 10 | c4_v_58_4x (86.7%), c4_v_185_4x (77.8%), his015_148 (70.0%) |
| OK (40-60%) | 8 | his015_96 (57.1%), his015_131 (53.3%), c4_v_62_max (50.0%) |
| CONCERN (<40%) | 2 | apical_his_148 (33.3%, non-switching control), his015_115 (36.4%) |

Note: sequence-position heuristic only — does not validate interface contact or exclude burial.
Top candidates for actual pH-switching: c4_v_58_4x, c4_v_185_4x, his015_148, his015_61, his015_140.

## PROPKA pKa analysis (v8, 2026-10-23)

ESMFold monomer structures → PROPKA 3.5.1. Wyman-linkage ΔΔG_monomer estimate.

| Rank | Design | pLDDT | His | Switching | ΔΔG_mono |
|------|--------|-------|-----|-----------|----------|
| 1 | c4_v_58_4x | 95.7 | 15 | 15 (100%) | **-19.36** |
| 2 | his015_131 | 94.9 | 15 | 15 (100%) | **-18.19** |
| 3 | his015_145 | 91.4 | 13 | 12 (92%) | **-17.78** |
| 4 | c4_v_33_4x_v2 | 94.8 | 13 | 13 (100%) | **-16.57** |
| 5 | c4_v_62_max | 95.6 | 16 | 15 (94%) | **-16.15** |
| … | (14 more) | | | | -6.87 to -14.00 |
| 20 | apical_his_148 | 90.5 | 3 | 2 (67%) | -2.59 |

Key: 19/20 designs >86 pLDDT. Only 1/215 His is always-on (H43 in his015_96, pKa=7.44).
Monomer ΔΔG is an upper bound — not comparable to complex ΔΔG (actimo-labs: 2.84-4.56 kcal/mol).

## 3D His spatial analysis (v8, 2026-10-23)

| Metric | Our designs | actimo-labs benchmark |
|--------|------------|---------------------|
| His per design | 6-16 (median 11) | 2-4 |
| His-His pairs <8Å | 2-11 per design | 0 (>8Å criterion) |
| His near D/E (<6Å) | 40-70% | Low (no carboxylate proximity) |

Our high-His-density strategy differs fundamentally from actimo-labs' sparse His placement.
The His-D/E intramolecular pairing risk is confirmed in 3D across all designs.

## Composite ranking (pLDDT + PROPKA + spatial)

Top 5: c4_v_58_4x (0.891), his015_131 (0.839), c4_v_62_max (0.809), c4_v_33_4x_v2 (0.796), his015_145 (0.781)
Bottom 2: c4_v_185_4x (0.423, pLDDT=77.6 ⚠), apical_his_148 (0.413, only 3 His ⚠)

**c4_v_185_4x flagged for replacement** — pLDDT below 80 suggests unreliable fold.
Replacement candidate: his015_176 (121 res, 13 His, Δq=9.8).

## Boltz-2 complex prediction (IN PROGRESS — all 20 designs)

Running Boltz-2 on ALL 20 designs (not just top 5) per Fable+Astra recommendation.
Target: TfR1 apical domain (154 res), binder as chain A.

### Early results (5/20 completed)

| Rank | Design | Size | ipTM | pTM | Confidence | Verdict |
|------|--------|------|------|-----|------------|---------|
| 1 | **his015_148** | 88 | **0.924** | 0.930 | 0.798 | Very strong |
| 2 | **his015_115** | 74 | **0.895** | 0.927 | 0.747 | Very strong |
| 3 | his200m | 78 | 0.741 | 0.834 | 0.757 | Strong |
| 4 | his015_57 | 84 | 0.507 | — | — | Possible |
| 5 | his015_96 | 107 | 0.422 | 0.622 | 0.640 | Possible |

**Benchmark**: actimo-labs (gold standard, TNF-alpha) achieved ipTM 0.294–0.650 (median 0.411).
Three of our first 5 designs EXCEED their best score.

### Fable caution (FABLE_BOLTZ2_CONSULTATION.md)
- ipTM >0.7 for de novo binders is unusual — TfR1 well-represented in Boltz-2 training (HFE, transferrin)
- Potential structural mimicry → BLAST top scorers against PDB
- Run multi-seed (--diffusion_samples 5) to check variance; >0.15 variance = noise
- Use ipTM ≥0.35 threshold for complex PROPKA analysis
- Compact designs scoring higher is expected (fewer degrees of freedom)

### Consultant findings (Fable + Astra)
1. **ΔΔG_monomer is a His counter (R²=0.80)** — not a switch metric
2. **No binding evidence existed before Boltz-2** — ESMFold is monomer-only
3. **Need complex PROPKA for true ΔΔG_switch** — pKa_bound vs pKa_free per interface His
4. **Mechanistic inversion**: His707 at pH 5.5 is protonated (+); binder should present D/E opposite it, not His

### Next steps (after Boltz-2 completes)
1. Download complex PDBs from Modal volume
2. Run complex PROPKA: pKa_bound (complex) vs pKa_free (extracted binder)
3. Identify interface His (within 5Å of target)
4. Check His707 contacts (which binder residues contact target His)
5. Compute true ΔΔG_switch (interface His only, Wyman linkage)
6. Replace c4_v_185_4x (pLDDT=77.6 + pI=5.5) with his015_118
7. Re-rank portfolio by ipTM + ΔΔG_switch + His707 contact

## Cost
- RFD4 generation: ~$9 on Modal L40S (200 e1_plain + 200 e1_his015)
- C3 RFD4 generation: ~$7 on Modal L40S (200 designs)
- ESMFold refold (~54 folds total across all iterations): ~$9 on Modal L40S
- ESMFold + PROPKA (C4 v8, 20 designs): ~$2 on Modal L40S
- Boltz-2 complex prediction (20 designs, in progress): ~$5-8 estimated
- **Total: ~$38+ of $20 budget** (overage from iterative pipeline development)
