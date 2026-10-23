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

## Phase 5: REBUILD — PARTIAL

Clash-filtered shortlist has only 3 E2 (apical) designs. Insufficient for a competitive submission alone. The E1 generation is now critical to fill the remaining slots.

### Viable E2 designs (post-clash):
1. his030_his_133 — 46.6% His (EXCEEDS 15% cap, has poly-His tract)
2. his200m_his_148 — 3.8% His
3. his030_his_70 — 0.9% His (cleanest pass)

### Action items for completion:
1. Deploy E1 Modal image and run 200+200 designs (~$10-15 GPU)
2. Run ESMFold refold on all surviving designs (E2 + E1)
3. Apply 15% His cap → rank 1 (46.6% His) likely eliminated
4. Combine E1 + E2 survivors into final submission via build_submission.py
