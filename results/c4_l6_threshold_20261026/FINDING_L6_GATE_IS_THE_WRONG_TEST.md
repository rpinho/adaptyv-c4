# FINDING: the L6 gate asks the wrong question. The repo documented the right threshold and the rule ignored it.

Written 2026-10-26 from the real clock. CPU only, **$0**, read-only. No GPU, no submission, no
organiser contact. No frozen file touched (4/4 verify).

I asked the four reviewers whether rule v4's `L6 < 0` gate has a calibration **offset** of order
0.2. Two answered that the offset is the wrong worry, and both were right. This records what they
found and what I verified.

## 1. The repo already specified the correct threshold, in writing

`scripts/ph_linkage.py` line 17, verbatim:

> "A binder that holds at pH 7.4 and lets go at pH 6.0 needs **L <= -2** (a hundredfold weaker
> binding at the lower pH)."

**Rule v4 gates at `L6 < 0`.** The boundary is misplaced by about **2 log units**, not 0.2.
Found by Fable 5; verified in place.

## 2. The same file contradicts itself — a FOURTH inverted path

| line | content | direction |
|---|---|---|
| 17 | `L <= -2` | correct |
| 128 | `switch_ok=bool(L >= 2.0)` | **inverted** |
| 164 | `"C4 gate needs L >= +2"` | **inverted** |
| 170-171 | sorts `reverse=True`, "most pH-selective first, C4 = L most positive" | **inverted** |

All four lines verified by direct read. The production path
(`c4_linkage_apical_recompute_20261025.py:124`) passes 6.0/7.4 explicitly and does not use the
inverted branch, so **the 342-design pool is not contaminated**. But any artifact ever produced
from `ph_linkage.py --csv` carries a reversed pass flag. This is the fourth inverted path found in
this repo, after `c4_linkage_audit_20261025.py`, `score_c4_acid.py`, and the v4 gate itself.

## 3. What `L6 < 0` actually admits (Astra's point, verified)

L = log10[Kd(7.4)/Kd(6.0)], so L is a **fold-change**, and the gate at zero admits *any* weakening
however small. Measured on the 342:

| weakening at pH 6.0 | L6 <= | n / 342 | + iPTM >= 0.60 |
|---|---:|---:|---:|
| 2x | -0.301 | 196 | 66 |
| 10x | -1.000 | 113 | 43 |
| **100x (the documented requirement)** | **-2.000** | **41 (12.0%)** | **18** |
| 1000x | -3.000 | 12 | 2 |

Pool median L6 = -0.508 = **3.22-fold**. Astra: *"1 nM becomes 3.22 nM. That is not 'no
detectable binding.'"* **35 designs pass the gate with under 2-fold weakening.**

## 4. Applied to the eight the reviewers unanimously chose to submit

Reproduced from source: `epitope_strategy.json` `footprint_human_auth` joined to
`native_context_gap.json` `per_position.rel_ctx4 >= 0.20`.
Chain: 342 -> acid-OFF 231 -> iPTM>=0.60 76 -> native area retention >=0.60 15 ->
>=10 natively-exposed contacts **8**.

| design | L6 | fold @ pH 6.0 | iPTM | retention | exposed contacts |
|---|---:|---:|---:|---:|---:|
| tfr1_c4_acid_site1_13 | -1.635 | 43.2x | 0.616 | 0.747 | 23 |
| tfr1_c4_acid_site1_147 | -1.609 | 40.6x | 0.943 | 0.668 | 10 |
| tfr1_c4_acid_site1_113 | -1.486 | 30.6x | 0.808 | 0.669 | 10 |
| tfr1_c4_acid_site2_24 | -0.826 | 6.7x | 0.782 | 0.686 | 15 |
| tfr1_c4_acid_site1_143 | -0.778 | 6.0x | 0.734 | 0.700 | 10 |
| tfr1_c4_acid_site2_128 | -0.555 | 3.6x | 0.640 | 0.634 | 15 |
| tfr1_c4_acid_site1_125 | -0.503 | 3.2x | 0.750 | 1.000 | 14 |
| tfr1_c4_acid_site2_55 | -0.236 | **1.7x** | 0.871 | 0.746 | 12 |

**None reaches the repo's own documented `L <= -2`.** Three reach 10x. One is under 2x.
Median 6.3x. Astra reached this from Kd arithmetic and Fable 5 from the documentation,
independently.

## 5. Zero is not a neutral boundary (Astra's requested statistic)

Astra rejected my first symmetry test as the wrong null and specified the right one: per design,
the half-sum of `L(pf+d, pf)` and `L(pf-d, pf)` — the expected L6 under **equiprobable shift
signs**, holding the observed shift magnitudes and pKa composition fixed. Measured on all 342:

| | mean L6 | n below zero |
|---|---:|---:|
| **observed** | **-0.506** | **231 / 342** |
| null: equiprobable shift signs | **+1.095** | **26 / 342** |
| null: single random sign draw | +1.079 | 86 / 342 |
| observed minus own null, per design | **-1.602** | 312 / 342 |

**The composition-matched null sits at +1.095, not 0.** Zero has no special status. Per Astra this
does **not** license subtracting the expectation; it establishes that the current boundary is
arbitrary.

Note the direction: designs sit ~1.6 units **below** their own null. Randomising shift direction
destroys the systematic burial effect, so this is a direct measurement of the mechanism in §6
rather than a refutation of it.

## 6. The deeper failure: the gate may not be selecting for pH switching at all

Fable 5, measuring a 66-design residue dump (the apical recompute did not persist per-residue
rows): HIS carries **-35.00** of a **-29.7** pool total; the other ~17 groups per design sit 2+ pH
units outside the 6.0-7.4 window and are mathematically inert. Median contributing histidines per
design: **3**.

> "`split_chains` deletes the partner's atoms from the *bound* coordinates, so the only thing it
> can see is the partner's electrostatic and desolvation effect on titratable groups in a frozen
> interface... burial destabilises the charged form of His and drops its pKa... A pool median of
> -0.508 is approximately three buried histidines doing nothing special. **The gate at L6 < 0 does
> not select for pH switches. It selects for 'has histidine near the interface,'** and it passes
> 231 of 342 because the design process put histidines near interfaces."

**VERIFIED on the full 342 pool** (`his_census_342.json`), and the claim holds with one
refinement:

| residue | n | total contribution | mean | median pKa_free |
|---|---:|---:|---:|---:|
| **HIS** | 2400 | **-363.8** | -0.1516 | **6.14** |
| LYS | 1461 | -130.9 | -0.0896 | 10.12 |
| N+ | 334 | -55.4 | -0.1659 | 7.98 |
| ARG | 712 | -0.0 | -0.0001 | 12.34 |
| TYR | 1594 | +0.3 | +0.0002 | 10.96 |
| ASP | 1564 | +87.6 | +0.0560 | 3.71 |
| GLU | 2470 | **+289.0** | +0.1170 | 4.61 |

| quantity | value |
|---|---:|
| pool total L6 | **-173.19** |
| pool total HIS contribution | **-363.77** (2.1x the net) |
| **median fraction of \|L6\| carried by His** | **0.998** |
| median contributing His per design | 6 of 31 groups |
| **designs where His alone flips the sign** | **207 / 342** |
| **corr(L6, number of contributing His)** | **-0.458** |

The refinement: Fable 5 called the non-His groups "mathematically inert", which is right for TYR
(+0.3 total) and ARG (-0.0) but **wrong for GLU/ASP/LYS**. GLU alone contributes +289.0 and acts as
a large partial counterweight to His's -363.8. The headline survives intact and is stronger than
the proxy suggested: **His carries 99.8% of the median design's \|L6\|, determines the sign outright
in 207 of 342 designs, and the more interface histidines a design has the more negative its L6
(r = -0.458).**

That closes the loop with the generation defect: both RFD4 configs set `C_SEQ_PROB_HIS: 0.15`.
A pool built to carry interface histidines sails through a gate that, measurably, scores interface
histidine count. The gate and the generator were selecting for the same thing.

## 7. The second gate is noise-compromised too

The only multi-seed data on disk (`results/boltz2_multiseed_results.json`, 5 designs, iPTM):

| design | single-seed | 5-seed mean | min | range | reliable |
|---|---:|---:|---:|---:|---|
| c4_v_185_4x | 0.929 | 0.860 | 0.698 | 0.239 | **false** |
| tfr1_e1_his015_148 | 0.924 | 0.805 | 0.727 | 0.188 | **false** |
| tfr1_e1_plain_98 | 0.923 | **0.617** | **0.308** | **0.559** | **false** |
| c4_v_58_4x | 0.916 | 0.740 | 0.561 | 0.358 | **false** |

Single-seed overstates the multi-seed mean in **4 of 4**, by up to 0.31. Rule v4's other gate is
`iPTM >= 0.60` on **single-seed** values, and two of the eight sit at 0.616 and 0.640 — inside
that overstatement band.

Fable 5 named the matching L6 job (recompute L6 across seeds; if within-design SD >= 0.5 the
111/231 split is model noise) as the highest-value free measurement. **It is not runnable:**
multi-seed structures were never saved for the 342, only these five iPTM summaries. Doing it
properly requires re-prediction, i.e. GPU, i.e. the open construct-scope decision.

## 8. Corrections to my own earlier claims

1. **I framed this as an offset problem. It is a magnitude problem.** Two reviewers said so
   independently and the repo's own line 17 had said so before either.
2. **I reported that the eight were not reproducible. They are.** My reconstruction used
   `frac_footprint_occluded` where the original used native relative SASA >= 0.20. The artifact
   was correct; my check was wrong. The reproduction is in §4 with the exact join.
3. **My P2 defence of the external references is void.** I argued they match the designs on
   contributing-group count (30/42/45 vs median 31). Fable 5: group count is ~85% inert residues;
   the comparison that matters is interface-His count and their pKa positions, where n = 3.
4. **My external sign audit is weaker than presented.** Astra: 3S9L carries iron in the N-lobe
   while the acid-retention biology concerns apo-transferrin, and the real cycle changes iron
   occupancy and conformation, whereas this calculation compares protonation on a fixed structure.
   That is closer to a state-matching error than the holo/apo caveat I recorded.

## 9. What is NOT claimed

- **No binding, affinity or pH response is measured anywhere here.** Every number is a PROPKA-based
  estimate on predicted complexes.
- **Nothing here says the designs are bad.** It says the gate that selected them does not test what
  the challenge asks.
- **Nothing here sets a correct threshold.** Astra: a calibrated gate needs measured Kd at both pH
  values across complexes spanning acid-OFF, neutral and acid-ON, construct- and state-matched.
  That is not achievable on CPU in six days.
- **§6 is verified on the full 342**; the mechanism is measured, but "generic burial" remains an
  interpretation of the measurement, not an independently demonstrated cause.
- The reviewers' recommendation to submit eight was made **before** §3-§4 existed. Whether it
  stands is theirs to revisit, not mine to adjust.
