# Phase Execution Status — 2026-10-23

## Phase 1: TRIM — DONE
Trimmed 20 → 9 designs per independent reviews (Fable + Astra).
Output: `results/trimmed_shortlist.csv`

## Phase 2: REFOLD — BLOCKED
ESMFold/AF2 single-sequence refold not available locally (no GPU, no ESM installed).
Requires Modal or cloud GPU environment.
Fable predicted this would eliminate the His-rich designs (low pLDDT / disordered).

## Phase 3: E1 GENERATION — RUNNING ON MODAL

### Configs created, validated, and deployed:
- `configs/rfd4/tfr1_e1_plain.yaml` — plain E1 binder, no His enrichment
- `configs/rfd4/tfr1_e1_his015.yaml` — E1 binder with His prior 0.15 on the designed chain

### Critical fixes applied during validation:
1. **Hotspot atoms corrected**: Original hotspots copied apical atom names onto E1 residues with different types (SER has no CD/NZ, GLY has no sidechain). Fixed to match actual E1' surface atoms.
2. **Chain B dimer context**: TARGET_PREP §6 warns E1 runs must include chain B. Only B565-575 (11 residues, the E1' dimer mirror directly facing the epitope) is included. B196-205 was dropped because RFD4 forbids two separate output chains with the same source chain letter, and the 360-residue structural gap between B205 and B565 prevents merging them into one chain.
3. **His707 (label 586) as primary hotspot**: The native histidine supplying half the proton linkage for pH switching, per TARGET_PREP mechanism family 6.
4. **Binder is chain C**: RFD4 assigns A (target 50 res), B (dimer 11 res), C (binder 60-130 res). His015 uses `sel('C')` for C_SEQ_PROB_HIS/CYS.
5. **`segments` field removed from his015**: Not recognized by RFD4 condition spec (only: conditions, contig, input, input_as_annotated_cif).

### Validated hotspots (label numbering → author numbering):
| Label | Author | Residue | Atoms | Role |
|-------|--------|---------|-------|------|
| 586   | 707    | HIS     | CE1,NE2 | **Primary** — native His for pH switch |
| 583   | 704    | SER     | CB,OG   | E1' surface |
| 585   | 706    | SER     | CB,OG   | E1' surface |
| 593   | 714    | GLU     | CD,OE1  | E1' surface |
| 609   | 730    | LEU     | CB,CG   | E1' surface |

### Generation runs (detached on Modal L40S):
- **e1_plain**: 200 designs, seed 20261023 — app ap-nc0cMb4yn7tkBz5niNx9dA
- **e1_his015**: 200 designs, seed 20261023 — launching
- Estimated cost: ~$4.50/spec × 2 = ~$9 total
- Results go to `rfd4-runs` volume: `tfr1_e1_plain/`, `tfr1_e1_his015/`

## Phase 4: CLASH RECHECK — DONE

### Method: vdW overlap > 0.4 Å (standard threshold)
The existing `clash_fraction_1p5A` column used a 1.5 Å threshold, which is shorter than a C–C bond (1.52 Å) and thus physically impossible to fail. This recheck uses the standard 0.4 Å overlap criterion.

### Also checked: His anti-switch contacts
Any designed His whose imidazole nitrogen is closer to a target carboxylate (Asp/Glu) than to a basic nitrogen (Arg/Lys) would STRENGTHEN binding at low pH when protonated — the opposite of pH-conditional release.

### Results:
| Rank | Design | Clashes | Max overlap | Anti-switch | Verdict |
|------|--------|---------|-------------|-------------|---------|
| 1 | his020_his_86 | 20 | 2.41 Å | 1 | **FAIL** severe |
| 2 | c200m_binder_126 | 6 | 2.18 Å | 0 | **FAIL** severe |
| 3 | his030_his_133 | 0 | 0.0 | 1 | MARGINAL (borderline anti-switch) |
| 4 | his030_his_115 | 4 | 0.96 Å | 1 | **FAIL** moderate+anti-switch |
| 5 | his020_his_152 | 18 | 2.38 Å | 2 | **FAIL** severe |
| 6 | his200m_his_148 | 1 | 0.49 Å | 0 | MARGINAL (single minor clash) |
| 7 | his200m_his_68 | 2 | 0.87 Å | 0 | **FAIL** moderate |
| 8 | his020_his_13 | 5 | 0.59 Å | 1 | **FAIL** moderate+anti-switch |
| 9 | his030_his_70 | 0 | 0.0 | 0 | **PASS** |

**7 of 9 designs have real clashes. 5 of 9 have His anti-switch contacts.**
**Only 3 survive: ranks 3, 6, 9 (one clean pass, two marginal).**

Common clash partner: A/208 (Arg329 author) — its sidechain collides in 6/7 clashing designs.

Output: `results/clash_recheck_04A.csv`, `results/clash_filtered_shortlist.csv`

## Phase 5: REBUILD — DONE (dry-run submission ready)

### E1 generation complete (400 designs):
- e1_plain: 200 designs, His% 0-11.4% (mean 2.7%)
- e1_his015: 200 designs, His% 0-35.1% (mean 9.5%)
- Cost: ~$9 on Modal L40S (within $20 budget)

### E1 filtering pipeline:
- 400 total → 375 (His cap 15%) → 374 (no poly-His ≥4) → 310 (no Cys)
- Diversity filter (≥30% Hamming): 17 diverse E1 candidates selected
- 16 from his015, 1 from plain (His prior worked as intended)

### E2 survivors (post-clash, post-His-cap):
- his030_his_133: 46.6% His → **EXCLUDED** (poly-His)
- his200m_his_148: 3.8% His, 1 minor clash → INCLUDED
- his030_his_70: 0.9% His, clean pass → INCLUDED

### Combined submission: 19 designs (2 E2 + 17 E1)
- Dry-run validated: `results/submission_dry_run_20261023.csv`
- Selection JSON: `results/combined_selection_20261023.json`

### Independent reviews (Fable + Astra):
Both reviewers converged on dropping 3 designs:
- `his015_193`: net charge -39, 38.6% D/E, likely disordered — **DROPPED**
- `his015_117`: internal tandem repeat, 45% hydrophobic — **DROPPED**
- `his015_168`: 47% hydrophobic, aggregation risk — **DROPPED**

Additional flags (kept with caveats):
- `plain_170`: 25-residue repeat, but only plain-arm design (control)
- `his015_62`: charge -16, borderline
- `apical_his_70`: VHH with 2 Cys (needs oxidizing expression), only 1 His

**Post-review: 19 → 16 designs (2 E2 + 14 E1)**

## Phase 6: REFOLD QUALITY GATE — DONE ✅

ESMFold single-sequence refold on Modal GPU (L40S):
- Model: `facebook/esmfold_v1` via HuggingFace transformers
- **ALL 16 designs PASS** (pLDDT ≥ 70)

| Design | pLDDT | Min | >70% |
|--------|-------|-----|------|
| his015_62 | **95.7** | 85 | 100% |
| plain_170 | **95.3** | 84 | 100% |
| his015_61 | **94.2** | 76 | 100% |
| his015_140 | **94.2** | 81 | 100% |
| his015_190 | **93.9** | 83 | 100% |
| his015_115 | **93.8** | 87 | 100% |
| his015_14 | **93.7** | 78 | 100% |
| his015_197 | **93.7** | 75 | 100% |
| his015_26 | **93.0** | 81 | 100% |
| his015_186 | **92.0** | 79 | 100% |
| apical_his_70 | **91.2** | 55 | 97% |
| his015_57 | **91.1** | 74 | 100% |
| apical_his_148 | **90.5** | 44 | 96% |
| his015_148 | **89.8** | 65 | 99% |
| his015_55 | **88.7** | 36 | 93% |
| his015_185 | **73.2** | 41 | 67% |

Total cost: ~$11 ($9 generation + ~$2 refold) of $20 budget.

### Remaining before final submission:
1. ~~Deploy E1 Modal image~~ DONE
2. ~~Independent reviews (Fable + Astra)~~ DONE — 3 designs dropped
3. ~~ESMFold refold quality gate~~ DONE — 16/16 PASS
4. **Ricardo review + approval** ← NEXT
5. Final `build_submission.py` run (non-dry-run)
