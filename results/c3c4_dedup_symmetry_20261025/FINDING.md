# The frozen "near-identical" rule had two defects — and they partially cancel (measured, $0)

**Status: evidence only.** Read-only re-measurement of the four 20-row artifacts. Nothing submitted,
nobody contacted, **no artifact modified** except one stale number in a README (see §5). Both advisors
were consulted before this was written; see §6 for who said what.

## 1. The rule as shipped

```
near-identical  iff  TM >= 0.65  AND  coverage >= 0.80
TM       = r.tm_norm_chain2            # normalised by CHAIN 2 only
coverage = min(non_gap_A/len_A, non_gap_B/len_B)
```

Two defects, both measured, both in `scripts/c3_round2_cluster_crosswalk_20261025.py:43-52`.

## 2. Defect 1 — the TM statistic is direction-dependent

`tm_norm_chain2` normalises by whichever structure is passed second. So `compare(a,b) != compare(b,a)`
when the chains differ in length, and dedup therefore depends on **selection order**.

Measured `|TM(a,b) − TM(b,a)|` across all 190 C4 pairs:

| chain-length gap | n pairs | mean abs TM delta | max |
|---|---:|---:|---:|
| 0–5 | 39 | 0.0075 | 0.0251 |
| 5–15 | 79 | 0.0359 | 0.0722 |
| 15–30 | 63 | 0.0824 | 0.1694 |
| 30+ | 9 | 0.1185 | 0.1421 |

**Pearson r = 0.940** between length gap and TM asymmetry. This is not noise; it is a systematic
length-dependent bias, and it is largest exactly where two designs differ most in size.

Concrete: `tfr1_c4_acid_site1_158` (L94) vs `tfr1_c4_acid_site2_1` (L70) gives
**TM = 0.7748 one way and 0.6054 the other** — the same pair straddles the 0.65 threshold purely on
argument order.

## 3. Defect 2 — the coverage clause is dead code

`coverage` counts non-gap residues **within each sequence separately**, including residues sitting
opposite a gap. That is sequence completeness, not *paired* alignment coverage. Consequence: it is
**exactly 1.0000 for all 760 pairs** across all four artifacts (verified: `min == 1.0` in every
artifact). The `coverage >= 0.80` half of the rule has **never rejected anything**.

True paired coverage on the pairs in question is far below 0.80 — and cannot exceed the length ratio:

| pair | lengths | as-shipped "coverage" | true paired coverage | max possible (min/max) |
|---|---|---:|---:|---:|
| filtered `p1_H8_52_64_910009_sd2` x `p4_H8_65_80_910127_sd1` | 52 / 76 | **1.0000** | 0.6711 | 0.6842 |
| robust `p1_H8_52_64_910015_sd1` x `p5_H4_65_80_910406_sd2` | 56 / 80 | **1.0000** | 0.6625 | 0.7000 |
| C4 `tfr1_c4_acid_site1_158` x `tfr1_c4_acid_site2_1` | 94 / 70 | **1.0000** | 0.734 | 0.7447 |

## 4. The four shipped artifacts: under the corrected rule, **their membership does not change**

*Scope: this section is about the four CSV artifacts only. Their **clustering** does change, and §4c is
where that matters.*

Corrected rule = **symmetric TM** (`max(TM_ab, TM_ba)`) **AND true paired coverage ≥ 0.80**, applied
to all 190 pairs of each shipped set — i.e. not in selection order, which removes the ordering bias too.

| artifact | as shipped | corrected | components (corrected) | document claim holds? |
|---|---:|---:|---:|---|
| C3 v8 | 1 pair | **1 pair** | **19** | yes (1 documented) |
| C3 filtered challenger | 0 pairs | **0 pairs** | **20** | **yes** |
| C3 all-five-robust | 0 pairs | **0 pairs** | **20** | **yes** |
| C4 direction-ranked | 0 pairs | **0 pairs** | **20** | **yes** |

My replay of the shipped selector reproduces the documented pair count for **all four** artifacts
(1/0/0/0), so no published component count is wrong.

**But this is membership only, and it is luck, not design — §4c shows the clustering changed a lot.** The two defects push in opposite
directions: the dead coverage clause makes the rule *too permissive*, while directional TM makes it
*permissive in one order and strict in the other*. The pairs that become violations under symmetric TM
(0.6541, 0.6879, 0.7748) all have true paired coverage **below 0.80**, so the coverage gate that was
supposed to catch exactly this class of length-mismatched containment was the one that was broken.
**These two bugs are not guaranteed to cancel in general.** Any future set with a long/short pair that
aligns well *and* pairs most of its residues would be admitted.

## 4b. Defect 1 alone already flags C4 — the selection order was hiding it

Worth separating, because it is the clearest demonstration of why order-dependence matters. Running the
**shipped** rule (chain-2 TM + the dead coverage clause) but evaluating it **order-free over all 190
pairs** instead of greedily in selection order:

| artifact | shipped rule, selection order | shipped rule, order-free | corrected rule |
|---|---:|---:|---:|
| C3 v8 | 1 pair / 19 | 1 pair / 19 | 1 pair / 19 |
| C3 filtered | 0 / 20 | **0 / 20** | 0 / 20 |
| C3 robust | 0 / 20 | **0 / 20** | 0 / 20 |
| C4 ranked | 0 / 20 | **1 pair / 19** | 0 / 20 |

So C4's published "0 near-identical pairs" is true **only in the exact order the selector happened to walk
the pool**. Evaluated order-free, the *same* rule finds the `site1_158` x `site2_1` pair (0.7748). The
claim happened to be saved by C4's selection order, not by the rule.

That pair then drops out again under the corrected rule because its true paired coverage (0.734) fails
the 0.80 gate — which is the coverage clause finally doing the job it was written to do, 24 Å-mismatched
residues later (94 vs 70).

## 4c. C3's clustering and the 08:00 verdict, re-run — with two of my own errors corrected

**First, two errors I made in the first version of this section, caught by Astra. Both changed numbers.**

1. **I compared on a truncated graph.** My first re-run used only the **153** edges that were *both*
   frozen and corrected, silently discarding the **56 corrected-only** edges. That is not the corrected
   rule; it is the intersection of two rules. Corrected figures below use **all 209**.
2. **My "no spurious links were added" check was tautological.** I intersected `keep` with `allf` where
   `keep ⊆ allf` by construction — it could not have found anything. In fact **56 corrected-only edges
   exist**, so the corrected graph is *not* a subset of the frozen one, and the "nothing was added"
   claim was simply false.

### The clustering result (holds)

Recomputed over all **19,900** pairs of the 200-key universe (`scripts/c3_recluster_corrected_rule_20261025.py`,
**126 s, $0**). Note the frozen matrix contained only **16,740 unique pairs**, so **3,160 pairs were
never measured before** — a corrected-only edge therefore cannot be attributed to symmetrisation alone.

| | frozen | corrected |
|---|---:|---:|
| near-identical edges | 284 | **209** |
| agree / frozen-only / corrected-only | — | **153 / 131 / 56** |
| connected components | 101 | **117** |
| **largest component** | **73** | **33** |

**131 of 284 frozen edges (46%) are rejected by the corrected predicate.** Examples:
`p2_H8_65_80_910115` x `p4_H8_52_64_910029` — sym TM **0.9575**, paired coverage **0.7571**;
`c3_glp1r_ecd_4` x `p3_H8_81_100_910206` — **0.7888 / 0.5000**. These are the length-mismatched
containment cases the dead coverage clause was meant to reject.

### The 08:00 verdict survives — but not at the magnitude I first reported

`scripts/c3_0800_comparison_corrected_rule_20261025.py`, same harness both columns, **full 209-edge set**:

| metric | frozen | corrected |
|---|---:|---:|
| v8 global components (PASS rows) | **9** | **9** |
| v8 own-matrix components (PASS) | 11 | 11 |
| depth own-universe components (all) | 53 | 64 |
| depth own-universe components (PASS) | 36 | 46 |
| **depth global components (PASS)** | **25** | **42** |
| depth PASS backbones novel vs v8 | 46 | 55 |

The harness **reproduces the published 25 vs 9 exactly** under the frozen rule — that is what validates
it. Under the corrected rule it is **42 vs 9**, not the **45 vs 10** I first published (truncated graph).

**Conclusion, scoped exactly:** the depth arm still wins rule 1 under both edge sets, and the corrected
rule moves the result in the depth arm's favour (25→42 while v8 stays 9). **I do not quote this as a
multiple** ("4.5×" was in my draft and is withdrawn). Astra's two standing caveats are adopted: adding
edges cannot inflate an absolute component count, so the *direction* is robust to graph bounds, but a
99-backbone pool compared with a 20-row artifact measures **qualifying units, not portfolio efficiency**.

### What the 73→33 change does and does not mean

Astra's correction, adopted: **components are single-linkage**, so a 73-member component can be entirely
legitimate (A matches B, B matches C, A need not match C). Size and fragmentation do **not** establish
that the members are distinct folds. The defensible statement is only:

> The frozen implementation connected 73 keys through threshold-passing paths, including edges the
> corrected coverage requirement rejects; under the corrected graph the largest component contains 33.

Also: 33 is the *size of the largest resulting component*, **not** the number of components it split
into — and that component includes one key outside the original 73. My earlier "splits it into 33" was
wrong phrasing.

### The validation gap Astra raised — now closed, and it was real
Astra: the "corrected" code took TM from both directions but **coverage from the forward alignment only**,
so the *full* predicate might still be directional. **It was.** Measured on all **19,900** pairs
(`scripts/dedup_predicate_symmetry_proof_20261025.py`):

| test | result |
|---|---|
| TM term directional | **19,457 / 19,900 pairs** |
| paired-column count differs by direction | **17** |
| coverage differs by direction | **17** |
| **final verdict differs by direction** | **67 pairs** |

So `max(TM)` alone did **not** make the predicate symmetric. Example: `A_26102500` x `c3_glp1r_ecd_4`
scores (TM 0.6850, cov 0.8182) one way and (TM 0.6050, cov 0.8182) the other — the *same pair* passes in
one order and fails in the other.

**A genuinely symmetric v2 was then written and tested** (`scripts/dedup_rule_v2_symmetric_20261025.py`):
`TM = max(tm_ab, tm_ba)`, `coverage = max(paired_cols_ab, paired_cols_ba) / max(len_a, len_b)` — both
terms are order-free, so the predicate is symmetric **by construction**. Result:

| | v1 (directional) | **v2 (symmetric)** |
|---|---:|---:|
| edges | 209 | **209** |
| edges changed | — | **0** |
| components | 117 | **117** |
| largest component | 33 | **33** |

**The 67 directional pairs are all non-linking**, so the published edge set, component count and largest
component are **unchanged** by the final symmetrisation. The numbers in this document stand; what changes
is that they now rest on a predicate that is symmetric by construction rather than one that happened to
agree on the edges that mattered.

## 5. What I changed, and what I did not

**Changed (one number, in prose):** `results/c4_ranked_20261025/README.md` said *"closest TM 0.646"*.
That figure is **stale** — it is from the superseded pre-gate run (`1d9aeaf2…`), not from the shipped
file. Under the corrected rule the closest **coverage-eligible** pair in the shipped C4 file is
**TM 0.6379**. The claim "0 near-identical pairs, 20 components" is unchanged and now verified.

**Not changed:** no CSV, no manifest, no component count, no selection. The artifacts are byte-unchanged
(hashes re-verified below). Correcting the *rule* would change selection and is not a $0 edit.

## 6. Advisors — both consulted independently, both converged

**Astra** (gpt-6-astra, full log `astra.raw.log`) diagnosed the coverage defect **before I had confirmed
it** — she flagged that `coverage = 1.00` was impossible for 52/76, 56/80 and 70/94, whose upper bounds
are 0.684 / 0.700 / 0.745, and predicted "the three new violations and 19-component conclusions would
not follow." **Her prediction was exactly right**, and it is why §4 concludes no change.
She also: called `max` "a defensible symmetric screen, not uniquely the correct definition"; warned
against claiming v8 unaffected "beyond TM symmetry"; and warned that component labels permit transitive
chains.

**Kimi** (`kimi.raw.log`, 911 tok in / 18517 out) independently reached the same coverage verdict in her
§E, and — importantly — ruled that on **dedup alone the two challengers are now tied**, so *"dedup cannot
be cited as an advantage for all-five-robust."* She also noted the filtered pair at **0.6541** is
genuinely marginal and that rounding-then-comparing is unsafe.

**Effect on the C3 ranking:** none on the verdict. `BEST = all-five-robust` was already justified
**solely** by 8/8 new seats clash-clean across five models, not by dedup. What changes is that the
ranking README's dedup row (`0 pairs` for both) is now stated as *"verified under the corrected symmetric
rule"* rather than *"as produced by the selector"* — same number, stronger footing.

## 7. The fix — a regression validator now exists (`scripts/dedup_rule_validation_20261025.py`)

Written and passing. It encodes both defects as executable invariants and computes all three rule modes
(`legacy` / `sym_tm` / `fixed`) over any artifact, so the old numbers stay reproducible for provenance.

- `--selftest` **PASS**: asserts the known asymmetric pair really is asymmetric, that symmetric TM flags
  it, that a 94-vs-70 pair cannot reach 0.80 paired coverage, and that the coverage bound holds.
- `--all` reproduces every number in this document from scratch in ~7 s at $0. Output saved to
  `validator_run.json`. Per-artifact maxima of `|TM_ab - TM_ba|`: **0.2216 v8, 0.2265 filtered, 0.2024
  robust, 0.1694 C4** — so the asymmetry ceiling is over 0.2 TM units, i.e. a third of the threshold.
- It also reports, per artifact, how many pairs had `legacy_cov == 1.0` **despite unequal chain lengths**
  (184–186 of 190). That is the dead-clause defect, now visible rather than silent.

Remaining recommendations for Ricardo (still his call — they change selection):

1. Treat the frozen rule as **deprecated for selection**. Reselection should use symmetric TM and true
   paired coverage, order-free.
2. Do **not** reselect the current four artifacts — and for C4 this is now **measured, not assumed**.
   I re-ran the full C4 selection pipeline (`iptm>=0.60`, `L6>0`, `L6(5.5/7.4)>0`, 2.0 Å rank-001 clash
   gate, greedy dedup, take 20) with the rule swapped to **symmetric TM + true paired coverage**. Result:
   **the same 20 designs, in the same order, with zero additions and zero drops.** The corrected rule
   rejects exactly one candidate (`tfr1_c4_acid_site2_99`, sym TM 0.7369, paired coverage 0.8529 — the
   same design the shipped rule already rejected) and accepts `site1_158` x `site2_1`, which it should.
   So **C4's membership is invariant to the rule fix**; only its prose needed the order-conditional
   qualifier. My earlier phrasing ("reselection would change membership") was an assumption I had not
   tested — it is **wrong for C4**, and is corrected here.

## 8. Artifacts and hashes (re-verified after this analysis)

All four CSVs byte-unchanged: v8 `631ed0bf…` · filtered `1450bcc7…` · robust `643490da…` ·
C4 ranked `3370909a…`. Raw measurements: `measure.json`, `rule_variants.json`, `corrected_closest.json`,
`length_gap.json`, `measure.log`, `rule_variants.log`, `corrected_closest.log`, `length_gap.log`,
`astra.raw.log`, `kimi.raw.log`.
