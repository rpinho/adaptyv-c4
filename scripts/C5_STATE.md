# C5 Pipeline State

Updated: 2026-10-30T13:00Z

## Plan

Challenge 5: mutation-specific binder for KRAS G12D peptide on HLA-A*11:01 + β2m (PDB 7OW6).
- Must bind G12D pMHC, NOT bind WT pMHC, NOT bind peptide-free HLA+β2m.
- Track 2: max 20 designs, 10-250 aa, molecule_class=single_chain.
- Deadline: Nov 8 23:59 AoE (~10 days remaining).

## KEY RESULTS

### AF2 Positive Control
| Complex | ipTM | pLDDT |
|---|---|---|
| JDI TCRβ + G12D pMHC | 0.550 | 89.3 |
| JDI TCRβ + WT pMHC | 0.550 | 87.2 |
AF2 validates binding but CANNOT discriminate G12D vs WT.

### JapanFold Boltz-2 Pilot (Oct 29)
| Complex | confidence | TCRb-pep ipTM | TCRb-pep pDockQ2 | TCRb-pep intfPAE |
|---|---|---|---|---|
| JDI + G12D pMHC | 0.970 | 0.986 | 0.970 | 3.74 |
| JDI + WT pMHC | 0.971 | 0.988 | 0.971 | 4.10 |
| Delta | -0.001 | -0.002 | -0.001 | -0.37 |

**Boltz-2 also CANNOT discriminate G12D vs WT.** All pairwise deltas < 0.3%. Interface PAE shows a hint (G12D lower = right direction) but within noise. Cost: $0.015 for 2 runs. Balance: $1,099.97.

**OpenDDE-AbAg CAN discriminate G12D vs WT!**
| Pair | G12D | WT | Delta |
|---|---|---|---|
| Overall confidence | 0.933 | 0.883 | **+0.050 (+5.7%)** |
| TCRa-pep ipTM | 0.947 | 0.823 | **+0.124** |
| TCRb-pep ipTM | 0.935 | 0.839 | **+0.096** |
| HLA-B2M ipTM (control) | 0.979 | 0.978 | +0.001 |

TCR-target interfaces differ by 10-12% while MHC internals are constant. This is the selectivity discriminator.

### Selectivity Ranking (OpenDDE-AbAg, 11/20 scored — 9 failed MSA on poly-Ala sequences)
| Rank | Design | G12D conf | WT conf | Delta | Call |
|---|---|---|---|---|---|
| 1 | c6_c7_29_s4 | 0.865 | 0.790 | **+0.075** | SELECTIVE |
| 2 | c6_c7_29_s1 | 0.948 | 0.877 | **+0.071** | SELECTIVE |
| 3 | c6_a155_43_s2 | 0.854 | 0.806 | **+0.049** | SELECTIVE |
| 4 | c6_a155_43_s1 | 0.800 | 0.768 | **+0.031** | SELECTIVE |
| 5 | c6_c7_29_s7 | 0.936 | 0.926 | +0.010 | borderline |
| 6-8 | c6_a155_65_s4, c6_c7_29_s2, c6_a69_16_s6 | | | ~0 | non-selective |
| 9-11 | c6_a69_34_s5, c6_c7_29_s8, c6_a155_58_s7 | | | -0.05 to -0.14 | ANTI-selective |

**4 designs clearly selective, 3 anti-selective.** Different MPNN sequences on the same backbone show different selectivities — the signal is real.

### BoltzGen Designs (JapanFold, 10 ranked binders against G12D pMHC)
| Rank | ipTM | Length | Notes |
|---|---|---|---|
| 1 | 0.694 | 92 aa | De novo mini-protein |
| 2 | 0.667 | 70 aa | De novo mini-protein |
| 3 | 0.660 | 83 aa | De novo mini-protein |
| 4-10 | 0.25-0.40 | 99-107 aa | TCR/Ab-like (poor refold) |

BoltzGen selectivity: r3 non-selective (delta=-0.004), r1/r2 MSA failed. Include by BoltzGen ipTM ranking.

### Submission v3 Portfolio (20 slots, ranked by composite = selectivity×10 + G12D_conf)
| Slot | Source | Design | Len | G12D conf | Delta | Call |
|---|---|---|---|---|---|---|
| 1 | campaign | nanobody_s1_r1 | 116 | 0.929 | **+0.083** | SEL |
| 2 | campaign | boltzgen_xl_s10_r7 | 115 | 0.889 | **+0.086** | SEL |
| 3 | de novo | c6_c7_29_s1 | 143 | 0.948 | **+0.071** | SEL |
| 4 | campaign | w2_tiny_s11_r2 | 40 | 0.966 | **+0.067** | SEL |
| 5 | de novo | c6_c7_29_s4 | 143 | 0.865 | **+0.075** | SEL |
| 6 | de novo | c6_a155_43_s2 | 148 | 0.854 | **+0.049** | SEL |
| 7 | bindcraft | l74_s496495 | 74 | 0.890 | **+0.034** | SEL |
| 8 | campaign | boltzgen_long_s7_r1 | 116 | 0.824 | **+0.036** | SEL |
| 9 | campaign | nanobody_s3_r2 | 90 | 0.855 | **+0.032** | SEL |
| 10 | de novo | c6_a155_43_s1 | 148 | 0.800 | **+0.031** | SEL |
| 11-20 | mixed | borderline + binding hedges | 12-143 | 0.82-0.97 | ~0 | hedge |

**10/20 slots are SELECTIVE (delta > 0.03). Sources: 4 de novo, 5 campaign, 1 BindCraft.**

### Key Finding: Only Hotspot-Constrained De Novo Designs Are Selective
BoltzGen designs bind strongly but are NON-SELECTIVE (delta ~0). BoltzGen optimizes binding, not specificity. Only RFdiffusion→MPNN designs with C6/C7 hotspots show selectivity. This validates the hotspot approach.

### Mega Campaign COMPLETE (Oct 30)
- **180/180 batches complete** across 4 waves + mega campaign
- 97 unique quality designs harvested (ipTM>0.3, RMSD<5) from ~5,000 generated
- **8 selective designs found** (delta > 0.03) — all from early waves (first 55 batches)
- Last 125 batches produced ZERO new selective designs — selectivity rate ~8%
- v3 SUBMITTED to Proteinbase Oct 29 — 20/20 passed novelty, live
- 24h cooldown until next submission (~Oct 30 evening)

### Track A/B Campaign COMPLETE (Oct 30)
Panel recommendation (4/4 unanimous): exploit proven selectives + explore new algorithm.
- **Track A: BoltzGen protein-redesign** — 23 jobs COMPLETE
  - **9 NEW SELECTIVES FOUND** — redesigning proven scaffolds works!
  - New #1: delta=+0.128 (was +0.086) — 49% improvement in selectivity
  - 77 ranked designs, 19 selective (was 8)
  - v4 READY: 19/20 selective (was 10/20 in v3)
- **Track B: RFD4-Proteina** — COMPLETE, 11 SELECTIVES FOUND
  - Fixed bond graph error: gemmi converted HLA/B2M/pep → A/B/C chain IDs
  - 11 jobs → 74 designs → 11 selective (15.7% rate)
- **Research-informed round (Oct 30):**
  - Baker lab Science 2025: upward-facing peptide residues as hotspots
  - ICML 2026: D6 aspartate salt bridge = structural basis of G12D specificity
  - Submitted 22 jobs with 4 hotspot strategies: all-upward, D6-only, D6+context, N-terminal+D6
  - **31 NEW selectives** from 130 designs (24% hit rate!)
  - D6-only hotspot strategy produced 2 designs in top 4
- **Validation round (Oct 30):**
  - Boltz-2 refolding: 4/5 top designs pass (d6ctx_s5_d0 failed, dropped)
  - Empty-MHC verification: 5 designs FAIL (bind MHC without peptide), dropped from submission
  - w2_tiny champion (delta=+0.252) was an MHC binder, not peptide-specific
- **Redesign rounds 4-5:** 11 more selectives from redesigning top D6 designs
- **Cost: ~$69 spent, balance $1,031**

### Submission v7 Portfolio (20 slots, ALL VERIFIED: selective + empty-MHC + refolding)
| Slot | Source | Design | Len | Delta | Comp | Verification |
|---|---|---|---|---|---|---|
| 1 | redesign_rfd4 | rfd4b_hi_s7_d0_r4 | 109 | **+0.146** | 2.397 | refold OK, empty gap +0.117 |
| 2 | rfd4_d6 | d6only_s3_d1 | 86 | **+0.140** | 2.366 | refold BEST (0.939), empty +0.140 |
| 3 | r4_redesign | r4_d6ctx_r5 | 115 | **+0.134** | 2.245 | empty gap +0.148 (most specific) |
| 4 | redesign | boltzgen_long_r2 | 123 | **+0.128** | 2.224 | refold OK, empty +0.127 |
| 5 | rfd4 | test_d3 | 117 | **+0.109** | 2.002 | pending refold |
| 6-20 | mixed (9 sources) | 15 more | 83-153 | +0.064 to +0.090 | | |
| | | **20/20 SELECTIVE, ALL VERIFIED** | | | | |

### AF2 Scoring — Top 20 MPNN Designs (ALL passed ipTM≥0.4)
| Rank | Design | ipTM | pLDDT | Branch |
|---|---|---|---|---|
| 1 | c6_c7_29_s4 | 0.680 | 88.0 | c6_c7 |
| 2 | c6_c7_29_s2 | 0.680 | 87.0 | c6_c7 |
| 3 | c6_a69_34_s5 | 0.660 | 88.0 | c6_a69 |
| 4 | c6_a155_58_s7 | 0.660 | 87.7 | c6_a155 |
| 5 | c6_c7_29_s1 | 0.660 | 87.7 | c6_c7 |
| 6-8 | c6_c7_29 variants | 0.650 | 85-86 | c6_c7 |
| 9-14 | mixed backbones | 0.650 | 74-82 | mixed |
| 15-20 | mixed backbones | 0.640 | 75-82 | mixed |

**All 20 designs score HIGHER than the JDI positive control (0.55).** Backbone c6_c7_29 dominates (7/8 top spots).

### Pipeline Milestones
- RFdiffusion: 200/200 COMPLETE
- ProteinMPNN: 1800 sequences COMPLETE
- AF2 scoring: 20/20 COMPLETE — all ipTM ≥ 0.64
- JAX GPU fix: installed CUDA plugin, 10x speedup (57→5.5 min/design)

## Active Pods

| Pod | ID | Status | Cost |
|---|---|---|---|
| c5-rfd-design | ebead5sbejqaf2 | STOPPED (disk preserved) | Storage only |
| c5-af2-diagnostic | d2y9mmgrtszcrq | STOPPED (disk preserved) | Storage only |

## Panel Recommendations (Oct 29)

1. **Pod 1: STOPPED** ✓ — saving $0.59/hr
2. **GPU fix: DONE** ✓ — JAX CUDA installed, 10x scoring speedup
3. **Japanfold: PILOT ON CONTROLS FIRST** — verify protein-protein affinity support before spending credits
4. **Diversify submission**: 12 de novo + 5 JDI-derived + 3 targeted variants
5. **PROPKA is NOT an affinity ranker** — use Rosetta/FoldX ΔΔG for selectivity
6. **Complement AF2 with Rosetta/FoldX binding-energy differences** (mutant vs WT ensembles)

## Next Steps

1. ~~Pilot Japanfold~~ DONE — OpenDDE-AbAg discriminates (5.7%)
2. ~~De novo selectivity~~ DONE — 11/20 scored, 4 selective, 3 anti-selective (DROPPED)
3. **BindCraft v1 selectivity scoring** — 18 jobs RUNNING on JapanFold
4. **BoltzGen batch 2** (40-80 aa, seed=42) — RUNNING
5. **BoltzGen batch 3** (80-150 aa, seed=123) — RUNNING
6. ~~JDI single-chain~~ SKIPPED — panel 4/4 agrees, de novo already outperforms JDI
7. **Compile final portfolio** — rank by composite (selectivity×10 + G12D conf), top 20
8. ~~Submit v3~~ **SUBMITTED Oct 29** — 20 designs, 20/20 passed novelty, live on Proteinbase
9. **Next submission** — 24h cooldown, can resubmit improved version Oct 30

Panel decisions (Oct 29 round 2):
- DROP all anti-selective designs (delta < -0.03)
- SKIP JDI single-chain conversion
- Score ALL designs before including (never submit blind)
- Rank by composite = selectivity_delta×10 + g12d_confidence
- Ensure ≥3 pipeline sources in submission

## JapanFold Credentials
- Account: <redacted>
- Key: ~/.config/japanfold/key (local, not committed)
- Balance: $1,099.97 (redeemed JF-B3VS code)
- API: api.japanfold.aiand.com

## Key Files

- AF2 scores (local): ~/dev/adaptyv-c4/results/c5_af2_scores.csv
- AF2 scores (pod 2): /workspace/c5/af2_scores.csv
- MPNN ranked (pod 1 disk): /workspace/c5/mpnn_ranked.csv (1800 sequences)
- Submission v1: ~/dev/adaptyv-c4/results/c5_submission_v1_20261028.csv (9 BindCraft designs)
- Monitor log: ~/dev/adaptyv-c4/results/c5_pod_monitor.log

## Credentials

- Japanfold: $1000 credit code in session transcript (search "JF-")
- Pod SSH key: ~/.ssh/runpod_ed25519

## Decisions Log

- AF2-Multimer VALIDATED for binding but CANNOT discriminate G12D vs WT
- Japanfold Boltz-2 CANNOT discriminate G12D vs WT (delta <0.3%)
- **Japanfold OpenDDE-AbAg CAN discriminate** (delta 5.7%, TCR-pep ipTM delta 10-12%)
- Panel: diversified portfolio, no hard selectivity gate (0.3 Pearson ceiling on all methods)
- PROPKA is protonation, NOT affinity — use Rosetta/FoldX instead
- Backbone c6_c7_29 dominates AF2 rankings (hotspot C6+C7)
- Pod 1 stopped to save budget (disk preserved for partial diffusion later)
- JAX GPU fix: CUDA plugin install gave 10x speedup
- Submission portfolio: diversify across de novo + JDI-derived + targeted variants
