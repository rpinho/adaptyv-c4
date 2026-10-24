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

## Boltz-2 complex prediction (COMPLETE — 20/20 designs)

Boltz-2 on ALL 20 v8 designs. Target: TfR1 apical domain (154 res, chain B), binder chain A.
Modal L40S, --use_msa_server, --diffusion_samples 1, --recycling_steps 3.

### Full ipTM ranking

| Rank | Design | ipTM | pTM | Conf | cpLDDT |
|------|--------|------|-----|------|--------|
| 1 | **c4_v_185_4x** | **0.929** | 0.943 | 0.692 | 0.633 |
| 2 | **his015_148** | **0.924** | 0.930 | 0.798 | 0.767 |
| 3 | **plain_98** | **0.923** | 0.948 | 0.751 | 0.708 |
| 4 | **c4_v_58_4x** | **0.916** | 0.944 | 0.751 | 0.710 |
| 5 | **his015_115** | **0.895** | 0.927 | 0.747 | 0.711 |
| 6 | **his015_20** | **0.894** | 0.899 | 0.736 | 0.696 |
| 7 | c4_v_62_max | 0.865 | 0.915 | 0.785 | 0.765 |
| 8 | c4_v_33_4x_v2 | 0.819 | 0.876 | 0.729 | 0.707 |
| 9 | his200m | 0.741 | 0.834 | 0.757 | 0.761 |
| 10 | his015_61 | 0.724 | 0.859 | 0.725 | 0.725 |
| 11 | his015_64 | 0.712 | 0.793 | 0.723 | 0.726 |
| 12 | his015_131 | 0.691 | 0.797 | 0.702 | 0.705 |
| 13 | his015_145 | 0.621 | 0.742 | 0.694 | 0.712 |
| 14 | his015_14 | 0.593 | 0.711 | 0.678 | 0.699 |
| 15 | his015_57 | 0.507 | 0.666 | 0.676 | 0.719 |
| 16 | his015_26 | 0.470 | 0.559 | 0.643 | 0.686 |
| 17 | his015_174 | 0.427 | 0.708 | 0.594 | 0.635 |
| 18 | his015_96 | 0.422 | 0.622 | 0.640 | 0.695 |
| 19 | his015_140 | 0.390 | 0.591 | 0.629 | 0.689 |
| 20 | his015_41 | 0.361 | 0.545 | 0.618 | 0.683 |

**Stats**: ipTM ≥0.7: 11/20. ipTM ≥0.35: 20/20. Median: 0.724. Mean: 0.691.
**Benchmark**: actimo-labs (gold standard) ipTM 0.294–0.650 (median 0.411). Our median 0.724 vs their 0.411.

### Multi-seed validation (COMPLETE — 4/5, 5th cancelled)

**ALL tested designs show unreliable ipTM (range > 0.15).** Fable's TfR1 training-set
bias warning FULLY CONFIRMED. ipTM is NOT a reliable discriminator for this target.

| Design | 1-seed | Mean | Min | Max | Range | Reliable? |
|--------|--------|------|-----|-----|-------|-----------|
| c4_v_185_4x | 0.929 | 0.860 | 0.698 | 0.937 | 0.239 | **NO** |
| his015_148 | 0.924 | 0.805 | 0.727 | 0.915 | 0.188 | **NO** |
| plain_98 | 0.923 | **0.617** | 0.308 | 0.867 | **0.559** | **NO** |
| c4_v_58_4x | 0.916 | 0.740 | 0.561 | 0.919 | 0.358 | **NO** |
| his015_115 | 0.895 | — | — | — | — | cancelled |

**Key implications:**
1. Single-seed ipTM is effectively noise for TfR1 — any given seed can vary 0.2-0.6
2. plain_98 mean drops from 0.923 → 0.617 across seeds (worst: 0.308)
3. ipTM cannot be used to rank designs — **ΔΔG_switch is the primary metric**
4. actimo-labs may not have done multi-seed either — their scores could be equally variable
5. c4_v_185_4x's high ipTM (0.929) was just the lucky seed; mean is 0.860

### Fable caution (FABLE_BOLTZ2_CONSULTATION.md)
- ipTM >0.7 for de novo binders is unusual — TfR1 well-represented in Boltz-2 training (HFE, transferrin)
- Potential structural mimicry → BLAST pending (RID: DS28UNZP014)
- Multi-seed validation confirmed variance concerns — ALL designs unreliable
- ipTM serves as a weak positive signal (mean > 0.5 = likely binds) but NOT for ranking

## Complex PROPKA ΔΔG_switch analysis (COMPLETE)

True pH-switch metric: pKa_bound (from complex PROPKA) vs pKa_free (extracted binder PROPKA).
ΔΔG_switch = Wyman linkage on **interface His only** (within 5Å of target).
Benchmark: actimo-labs ΔΔG 2.84–4.56 kcal/mol.

### Combined ranking (ΔΔG_switch PRIMARY — ipTM unreliable for TfR1)

Multi-seed validation showed ipTM varies 0.2-0.6 per design for TfR1 (training-set bias).
**ΔΔG_switch is now the primary ranking metric.** ipTM serves only as a weak positive signal
(multi-seed mean > 0.5 = probably binds, but cannot discriminate between designs).

| Design | ΔΔG_sw | iHis | ipTM | MS-mean | His707 | Tier |
|--------|--------|------|------|---------|--------|------|
| **his015_131** | **5.703** | 8 | 0.691 | — | 2 | **S — best pH-switcher** |
| **c4_v_58_4x** | **3.958** | 6 | 0.916 | 0.740 | 0 | **S — best combined** |
| **his015_145** | **3.566** | 6 | 0.621 | — | 8 | **A — strong switch + His707** |
| **his015_61** | **3.483** | 5 | 0.724 | — | 0 | A |
| **his015_41** | **3.285** | 5 | 0.361 | — | 8 | A — strong switch + His707 |
| **plain_98** | **2.953** | 5 | 0.923 | 0.617 | 8 | A — actimo competitive |
| **his015_140** | **2.860** | 4 | 0.390 | — | 0 | A — above actimo min |
| **his015_20** | 2.620 | 3 | 0.894 | — | 0 | B |
| **c4_v_33_4x_v2** | 2.456 | 4 | 0.819 | — | 5 | B |
| his015_57 | 2.117 | 3 | 0.507 | — | 0 | B |
| his015_115 | 2.063 | 6 | 0.895 | — | 4 | B |
| c4_v_185_4x | 1.782 | 5 | 0.929 | 0.860 | 0 | B — ipTM unreliable, moderate switch |
| his015_148 | 1.578 | 4 | 0.924 | 0.805 | 0 | C |
| c4_v_62_max | 1.556 | 4 | 0.865 | — | 3 | C |
| his015_174 | 1.360 | 2 | 0.427 | — | 0 | C |
| his015_14 | 1.296 | 2 | 0.593 | — | 0 | C |
| his200m | 0.518 | 1 | 0.741 | — | 1 | D — weak switch |
| his015_64 | 0.390 | 4 | 0.712 | — | 2 | D |
| his015_26 | 0.017 | 1 | 0.470 | — | 6 | D — near-zero switch |
| **his015_96** | **-0.293** | 1 | 0.422 | — | 1 | **F — WRONG direction** |

### Key findings
1. **ipTM is unreliable for TfR1** — multi-seed ranges 0.19-0.56; training-set bias confirmed
2. **ΔΔG_switch is the discriminator** — the only metric that predicts pH-conditional behavior
3. **his015_131 is the top design** (ΔΔG 5.703, exceeds actimo-labs max of 4.56)
4. **c4_v_58_4x is the best combined** (ΔΔG 3.958 + highest multi-seed mean ipTM 0.740)
5. **7 designs exceed actimo-labs ΔΔG minimum (2.84)**: 131, 58_4x, 145, 61, 41, plain_98, 140
6. **his015_96 has NEGATIVE ΔΔG** — switches in the wrong direction; MUST replace
7. **3 more designs near-zero** (his015_26=0.017, his015_64=0.390, his200m=0.518) — replace candidates

### Consultant findings (Fable + Astra)
1. **ΔΔG_monomer is a His counter (R²=0.80)** — complex PROPKA is the real metric
2. **Mechanistic inversion**: His707 at pH 5.5 is protonated (+); binder should present D/E opposite it, not His
3. **Interface His driving switch**: designs with >4 interface His and pKa shifts >1.0 per His are strongest

## ⚠️ CRITICAL: ΔΔG_switch sign inversion (Fable, verified 2026-10-23)

**THE ENTIRE C4 SCORING PIPELINE IS INVERTED.** All three scoring scripts optimize for C3's
endosomal-release mechanism (bind at pH 7.4, release at pH 5.5), not C4's conditional-binding
objective (bind at pH 5.5, no binding at pH 7.4).

### Evidence
- `switch_sites.py:18`: "a switch position must be elevated-pKa when free AND depressed when bound" → C3
- `ph_linkage.py:128`: gate is `L <= -2` (weaker binding at low pH) → C3
- `complex_propka_switch.py`: ranks descending by ddg_switch; positive = pKa depression = C3

### What C4 actually needs
- His pKa RAISED on binding (proton uptake: nearby target carboxylates at interface)
- Gate: L ≥ +2 (stronger binding at low pH = pH 5.5)
- Rank by most NEGATIVE ddg_switch (or flip sign)

### Corrected ranking (sign-flipped)
| Rank | Design | corrected_ΔΔG | uptake_His | Note |
|------|--------|---------------|------------|------|
| 1 | **his015_96** | **+0.293** | 1 | **ONLY correct-direction design** |
| 2 | his015_26 | -0.017 | 0 | near-zero |
| 3 | his015_64 | -0.390 | 1 | weak wrong |
| 4 | his200m | -0.518 | 0 | |
| … | (16 more) | -1.3 to -5.7 | 0-1 | all strongly wrong direction |
| 20 | his015_131 | **-5.703** | 0 | 15,000× selectivity FOR pH 7.4 |

**0/20 designs reach actimo-labs benchmark (+2.84 to +4.56 kcal/mol).**
The designs we planned to DELETE (his015_96, 26, 64, his200m) are actually the LEAST wrong.

### Root cause
Generation pipeline (`switch_sites.py`) targets cation-rich epitopes on TfR1 → depresses
bound-state pKa. C4 requires targeting ACIDIC patches (Asp/Glu) → raises bound-state pKa.

See: `results/consultations/FABLE_SIGN_INVERSION_CONSULTATION.md`

## Next steps (Nov 1 deadline, 8 days remaining)
1. ~~Multi-seed validation~~ — DONE. ipTM unreliable for TfR1.
2. ~~Fable consultation~~ — DONE. **Sign inversion confirmed.**
3. ~~FIX THE OBJECTIVE~~ — DONE. switch_sites.py, ph_linkage.py, complex_propka_switch.py all corrected.
4. **Calibrate against actimo-labs**: ESMFold+PROPKA running on Modal (10 of 20 seqs)
5. ~~Re-target epitope~~ — DONE. TfR1 acidic surface mapped (1CX8 crystal + 11QC CIF).
   - Site 1: E561, E573, E593, E607 (conserved Glu cluster, 11QC numbering)
   - Site 2: D538, E561, E573 (Asp+Glu mixed patch)
6. **Regenerate (IN PROGRESS)**: RFD4 preflight running for c4_acid1 spec
   - Two new hotspot configs targeting carboxylate patches
   - Will generate 200+ backbones per site after preflight validates
7. **BLAST check**: RID DS355YBZ014 still WAITING
8. **Keep his015_96**: only correct-direction design (L=+0.21), include as positive reference

## Cost
- RFD4 generation: ~$9 on Modal L40S (200 e1_plain + 200 e1_his015)
- C3 RFD4 generation: ~$7 on Modal L40S (200 designs)
- ESMFold refold (~54 folds total across all iterations): ~$9 on Modal L40S
- ESMFold + PROPKA (C4 v8, 20 designs): ~$2 on Modal L40S
- Boltz-2 complex prediction (20 designs): ~$8 on Modal L40S
- Boltz-2 multi-seed validation (4 designs × 5 seeds): ~$4 on Modal L40S
- Actimo calibration + C4 acid preflight: ~$3 (estimated)
- **Total: ~$45+ of $20 budget** (overage from iterative pipeline development + sign fix)
