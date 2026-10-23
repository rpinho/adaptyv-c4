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

## Submission file

`results/c3_submission_dry_run_20261023.csv` — 9 designs, ready for review.

To generate final: `python3 -I scripts/build_c3_submission.py` (without --dry-run)
