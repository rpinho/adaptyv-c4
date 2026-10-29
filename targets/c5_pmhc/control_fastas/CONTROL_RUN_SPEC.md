# C5 positive-control run specification

Three AF2-Multimer predictions using the JDI TCR from 7OW6 as the binder.

## Expected outcome

| FASTA | expected | reason |
|---|---|---|
| control_tcr_vs_g12d | binding signal (ipTM/pTM) | KD 63 µM experimentally |
| control_tcr_vs_wt | weak/no signal | NB experimentally |
| control_tcr_vs_pepfree | weak/no signal | expected — TCR reads peptide |

## Pass criterion

ipTM(G12D) > ipTM(WT) AND ipTM(G12D) > ipTM(pepfree).
If this fails, the scoring pipeline cannot discriminate and must not be used.

## Compute note

Each prediction is a 5-chain, ~811 residue complex.
This is large for AF2-Multimer. Use at least an A40 (48 GB). Expected wall time: 30-90 min per prediction depending on MSA mode.

## MSA amortisation

The TCR alpha, TCR beta, heavy chain, and b2m MSAs are identical across all three runs.
Only the peptide changes (10 aa, single-sequence mode sufficient). Pre-compute the four
target-side MSAs once and reuse them.
