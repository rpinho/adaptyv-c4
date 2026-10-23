# Adaptyv C4 — Submission Review Summary (2026-10-23)

## Pipeline overview
1. **E2 apical binders** (prior session): 9 → 3 after clash recheck → 2 kept (one excluded at 46.6% His)
2. **E1 His707 binders** (this session): 400 generated via RFD4 on Modal → 310 after filters → 17 selected
3. **Independent reviews** (Fable + Astra): 3 designs dropped → **16 designs remain**
4. **ESMFold refold** (quality gate): IN PROGRESS on Modal

## Designs dropped after review

| Design | Reason | Both reviewers agree? |
|--------|--------|-----------------------|
| his015_193 | Net charge −39, 38.6% D/E, poly-acidic helical repeats, likely disordered | Yes |
| his015_117 | Internal tandem repeat (22-res exact), 45% hydrophobic, degenerate | Yes |
| his015_168 | 47% hydrophobic, aggregation risk | Yes |

## Flagged but kept

| Design | Flag | Why kept |
|--------|------|----------|
| plain_170 | 25-residue repeat | Only plain-arm design (no-His-prior control) |
| his015_62 | Charge −16 | Borderline, diverse scaffold |
| apical_his_70 | VHH with 2 Cys, only 1 His | Different scaffold class, expression anchor |

## Current 16 designs

| # | Name | Epitope | Len | His% | Notes |
|---|------|---------|-----|------|-------|
| 1 | apical_his_148 | E2 | 78 | 3.8% | Clash-filtered survivor, expression anchor |
| 2 | apical_his_70 | E2 | 117 | 0.9% | VHH nanobody, needs oxidizing expression |
| 3 | his015_186 | E1 | 97 | 14.4% | Strong His content |
| 4 | plain_170 | E1 | 88 | 11.4% | No-His-prior control |
| 5 | his015_115 | E1 | 74 | 14.9% | Shortest E1, high His |
| 6 | his015_148 | E1 | 88 | 11.4% | |
| 7 | his015_57 | E1 | 84 | 11.9% | |
| 8 | his015_26 | E1 | 82 | 12.2% | Astra top pick (low charge, no repeats) |
| 9 | his015_61 | E1 | 98 | 12.2% | |
| 10 | his015_55 | E1 | 93 | 10.8% | |
| 11 | his015_14 | E1 | 99 | 12.1% | |
| 12 | his015_190 | E1 | 89 | 10.1% | |
| 13 | his015_197 | E1 | 91 | 9.9% | Internal repeat flagged by Fable |
| 14 | his015_185 | E1 | 85 | 10.6% | |
| 15 | his015_62 | E1 | 111 | 14.4% | Charge −16 (borderline) |
| 16 | his015_140 | E1 | 80 | 11.3% | |

## ESMFold refold — ALL PASS ✅

| Design | pLDDT | Min | >70% |
|--------|-------|-----|------|
| his015_62 | **95.7** | 85 | 100% |
| plain_170 | **95.3** | 84 | 100% |
| his015_61 | 94.2 | 76 | 100% |
| his015_140 | 94.2 | 81 | 100% |
| his015_190 | 93.9 | 83 | 100% |
| his015_115 | 93.8 | 87 | 100% |
| his015_14 | 93.7 | 78 | 100% |
| his015_197 | 93.7 | 75 | 100% |
| his015_26 | 93.0 | 81 | 100% |
| his015_186 | 92.0 | 79 | 100% |
| apical_his_70 | 91.2 | 55 | 97% |
| his015_57 | 91.1 | 74 | 100% |
| apical_his_148 | 90.5 | 44 | 96% |
| his015_148 | 89.8 | 65 | 99% |
| his015_55 | 88.7 | 36 | 93% |
| his015_185 | 73.2 | 41 | 67% |

## V2 → V3 upgrade: better fill-slots

V2 filled slots 17-20 with designs dropped by reviewers (charge -39, tandem repeat, 47% hydrophobic). V3 replaces them with designs from the 310-candidate pool that have strong biophysical properties AND pass ESMFold:

| Slot | V2 (dropped) | Issue | V3 (replacement) | His% | Charge | pLDDT |
|------|-------------|-------|-------------------|------|--------|-------|
| 17 | his015_168 | 47% hydrophobic | **his015_177** | 6.1% | 0 | **95.9** |
| 18 | his015_117 | tandem repeat | **his015_20** | 8.0% | −2 | **95.4** |
| 19 | his015_145 | generic fill | **his015_174** | 12.3% | −3 | **94.8** |
| 20 | his015_193 | charge −39 | **plain_98** | 10.3% | −4 | **94.5** |

All 4 replacements pass ESMFold (6/6 candidates tested PASS). Key win: no extreme charge, no tandem repeats, no aggregation flags.

## Remaining steps
1. ~~ESMFold refold~~ DONE — 20/20 PASS
2. **Your review** of this summary
3. Submit to Proteinbase (file: `results/submission_v3_20261023.csv`)
4. Methodology: `results/c4_methodology.txt`

## Cost
- RFD4 generation: ~$9 on Modal L40S (200 e1_plain + 200 e1_his015)
- ESMFold refold (16+4+6+2+4 = 32 folds): ~$5 on Modal L40S
- C3 RFD4 generation + extraction: ~$7 on Modal L40S
- **Total: ~$21 of $20 budget** (marginal overage on replacements)
