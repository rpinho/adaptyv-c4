# C4 structural screen used the WRONG receptor target — measured, reproducible, $0

**Read-only. No candidate file changed. Nothing submitted. No GPU, no pod, no billing.**

## The finding

The shipped C4 ranked file (`3370909a…`) ranks 20 designs by **apical-context** iPTM and L6
(195-residue receptor, `P02786` 189–383), but its **structural screen — d_min and the 2.0 Å clash gate —
was computed on a different receptor**: the 154-residue stitched surrogate in
`results/c4_rebuild_20261024/remote_pdbs/`.

Every one of those 343 local PDBs has **chain B = 154 CA**. The apical predictions on the Modal volume
`rfd4-runs/boltz2_c4_apical/` have **chain B = 195 CA**. Same binder (chain A sequence identical), a
different receptor construct.

`scripts/c4_rank_acid_direction_20261025.py:69` reads the canonical path
`results/c4_rebuild_20261024/remote_pdbs/{n}_complex.pdb` — the 154-residue files — for the clash
screen, while the iPTM/L6 values it ranks come from the 195-residue apical run.

## Measured consequence on the same 20 designs

| screen used | passes (d_min ≥ 2.0 Å, 0 pairs < 2.0 Å) |
|---|---|
| 154-residue surrogate (shipped) | **20/20** |
| 195-residue apical (correct) | **14/20** |

Six ranked designs gain at least one sub-2.0 Å contact under the correct construct:

| rank | design | d_min 154 → 195 | pairs <2.0 Å (195) |
|---:|---|---|---|
| 8 | `tfr1_c4_acid_site2_171` | 2.008 → 1.925 | 1 |
| 9 | `tfr1_c4_acid_site2_14` | 2.110 → 1.711 | 1 |
| 10 | `tfr1_c4_acid_site1_6` | 2.115 → 1.935 | 1 |
| 13 | `tfr1_c4_acid_site1_35` | 2.274 → 1.797 | 2 |
| 17 | `tfr1_c4_acid_site1_20` | 2.093 → 1.799 | 3 |
| 18 | `tfr1_c4_acid_site1_158` | 2.252 → 1.654 | 1 |

d_min shifts by up to **0.60 Å** in either direction (rank 18 −0.598, rank 19 +0.432). The two constructs
are not interchangeable for interface geometry, which is expected: the surrogate draws fragments from
P02786 positions spread across 259–677 and shares only a 16-residue longest common substring with the
apical domain.

## Controls run before recording this

1. **Reproduction control.** For six ranked designs, the shipped ranked `d_min_A` reproduces **exactly**
   (6/6, <0.02 Å) from the 154-residue PDBs — confirming the shipped screen really used those files,
   rather than my re-measurement being off.
2. **Identity control.** For `site2_1` and `site1_35`, binder chain A sequence is **identical** between
   the apical and rebuilt PDBs, while receptor chain B differs (195 vs 154). So this is a
   receptor-construct substitution, not a binder mismatch or a different design.
3. **Independent parser.** Measured with the same standalone heavy-atom parser used by the C3 clash
   audit; chain counts read directly from PDB `ATOM` records.

Artifacts: `structural_recheck.json` (all 20 rows), the 20 apical PDBs in this directory,
script `scripts/c4_apical_pdb_audit_20261025.py`.

## CORRECTION — 11:40 ET Oct 25: L6 was also surrogate-based

This finding originally stated that "iPTM and L6 are apical-context quantities and are ranked on
apical-context values consistently" and therefore the ranking was unaffected. **That is false for L6.**
`scripts/c4_selection_audit_20261025.py:27` computes L6 from the same 154-residue surrogate PDBs used for
the clash screen. Only iPTM and the pKa summary statistics came from the apical 195-residue run.

Re-scored on the correct apical complexes (`results/c4_linkage_apical_20261025/FINDING.md`, sha
`fcb8866c…`), the shipped ranked 20 are **8/20 acid-ON**, not 20/20; the six matched pilot pairs reverse
(1/6 candidates vs 3/6 controls). The structural consequence below still stands, but the ranking criterion
was mixed-construct too. **Treat the "does not invalidate the ranking" clause as superseded.**

## What this does and does not mean

- **It does not invalidate the ranking.** iPTM and L6 are apical-context quantities and are ranked on
  apical-context values consistently.
- **It does invalidate the structural-cleanliness claim for the ranked file.** "20/20 pass the 2.0 Å
  screen" is true only of the surrogate, not of the construct the designs were ranked on. The honest
  statement for the correct construct is **14/20 clean, 6 flagged**.
- **It changes one pilot pair.** `site1_35` (rank 13) is one of the six pilot candidates and one of the
  six flagged designs under the correct screen. That does not remove it from a frozen roster, but its
  structural QC expectation must be recorded in advance.
- **It is a distance screen, not physical plausibility.** 14/20 is not "good", and 20/20 was never
  evidence of anything (the atom-type-blind caveat published with the ranked file still applies).
- The Probe calibration finding is unaffected: that metric does not adjudicate this screen either way.

## Recommended action (Ricardo only; $0 follow-up already done here)

1. Any statement citing "20/20 structurally clean" for `3370909a…` must be corrected to
   "20/20 on the 154-residue surrogate; **14/20** on the 195-residue apical construct actually ranked".
2. The selection script's PDB path should be fixed so the screen and the ranking use the same construct.
   This is a code change, and it would change which designs pass — so it is a **new versioned artifact**,
   not an edit to `3370909a…`.
3. Before any native pilot, record `site1_35`'s correct-construct clash status as a known property.

## Hashes

| artifact | sha256 |
|---|---|
| `structural_recheck.json` | see `SHA256SUMS` in this directory |
| ranked file (unchanged) | `3370909acebba7300e017f468107ff8f8ebd57c3463b7186a6a9e51b0a9d0a6f` |
