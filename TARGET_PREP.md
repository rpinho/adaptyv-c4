# Target prep: epitope shortlist (first pass, 2026-10-20)

Project root for the Challenge 4 entry. Deliberately **outside the KB** — nothing from this
work goes onto `main`; the plan document lives on the `adaptyv-c4` branch of the KB worktree only.

## What was run

`scripts/align_tfr1.py` — global Needleman-Wunsch, UniProt P02786 (human, 760 aa) vs Q62351
(mouse, 763 aa), from `data/*.txt`:

- 77.5% identity overall, 76.5% over the ectodomain 89-760, **78.0% over the assay construct 121-760**
- fully conserved 20-aa windows include 461-486 and 725-744
- least conserved: 352-374 (35-40%) and 568-589 (40%) — an epitope there forfeits the
  mouse-cross-reactivity criterion
- native histidines in 121-760: 186, 234, 300, 302, 318, 401, 475, 515, 603, 684, 699, 707

`scripts/epitope_scan.py` — intersects three masks on the ectodomain and grows surface patches.
Structure facts from the run itself:

- apo structure 1CX8, chain A modeled **122-760** (639 residues); the physiological dimer
  partner is chain B (336 contact atoms at 5 A), and SASA is computed **in the dimer**, so the
  dimerization interface is excluded automatically
- N-glycosylation sites detected from NAG-linked ASN: **251, 317, 727**
- transferrin footprint from 1SUV: **33 TfR1 residues** (TfR chains A/B vs transferrin C-F, 5 A)
- C-terminal tag vicinity: 29 residues within 15 A of residue 760 (the Twin-Strep anchor)
- 170 candidate residues survive all three masks (relative SASA >= 0.25, inside 121-760,
  not in the footprint / glycan / tag sets)

## Shortlist (identity >= 0.80, >= 2 hydrophobics, 10 A patch radius)

| centre | n | ident | hyd | chg | mean rel SASA | dist to Tf footprint | dist to tag | native His | residues |
|---|---|---|---|---|---|---|---|---|---|
| **303** | 9 | 0.89 | 4 | 1 | 0.277 | 26.1 A | 41.0 A | **300** | 263,264,266,269,300,303,328,537,704 |
| **706** | 9 | 1.00 | 3 | 1 | 0.303 | 17.3 A | 27.3 A | - | 690,695,704,705,706,711,714,730,733 |
| 318 | 7 | 0.86 | 2 | 2 | 0.376 | 46.5 A | 54.4 A | **318** | 187,318,319,322,381,382,383 |
| 328 | 8 | 0.88 | 2 | 1 | 0.368 | 29.6 A | 41.7 A | - | 269,303,311,325,326,327,328,331 |
| 532 | 8 | 0.88 | 2 | 5 | 0.347 | 9.6 A | 34.7 A | - | 495,503,530,531,532,533,534,613 |

Reading, with the caveats that matter:

- **Patch 303 is the lead.** It is the only patch that clears both the hydrophobic-content bar
  (4 hydrophobics, which is what the public binder pipelines want) and high conservation, and it
  **contains native His300** — which is the target supplying half of the pH linkage for free
  (mechanism family 6 in the plan). It sits in the apical-domain region, the non-transferrin-blocking
  epitope class, 26 A from the detected footprint.
- **Patch 706 is the diversity partner:** 100% human/mouse identity, 3 hydrophobics, but it is in
  the membrane-proximal helical region, 27 A from the tag anchor — closer to both the tag and the
  footprint, so it carries more risk of interfering with capture or with transferrin.
- **Patch 318 brings the second native His** (318) and is the furthest from everything (46 A from
  the footprint), but only 2 hydrophobics: designability risk.
- **Patch 532 is charge-heavy (5) and only 9.6 A from the footprint.** Deprioritise.

## Known weaknesses of this pass (fix before locking the epitopes)

1. **The transferrin footprint is from a 7.5 A cryo-EM structure** and only 33 residues came out of
   it. That is almost certainly an underestimate, so "does not block transferrin" is weakly
   established. Next: find a higher-resolution TfR1-transferrin or TfR1-antibody structure and
   redo the footprint, then re-run the mask.
2. **Domain boundaries were not used**, only residue windows. The apical-domain range should be
   taken from the Lawrence 1999 paper rather than assumed.
3. Patches mix sequence-distant residues (e.g. 537 and 704 in patch 303). That is normal for a
   surface patch, but each shortlisted patch needs a visual check before it becomes a hotspot set.
4. The scan ranks on a hand-weighted score. It is a shortlisting aid, not evidence.

## Next

Hotspot sets for the generators come from patches **303** and **318** (both carrying a native His,
both apical) plus **706** as the structurally independent third epitope, pending fix 1.
