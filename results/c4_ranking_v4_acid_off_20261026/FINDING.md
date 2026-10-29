# C4 rule v4 (ACID-OFF): the ranking, and the sign fix is not worth anything measurable

Written 2026-10-26T14:35:28Z from the real clock (`date -u`). $0, CPU only. No GPU, no Modal, no RunPod,
no pod, no submission, no upload, no organiser contact. C3 not touched: it is closed and
submitted. v8/v7/v6/`submission_v9.csv` verified unchanged by hash before and after
(`manifest_frozen.json`, all four prefixes match).

**Headline: verdict (b), statistically indistinguishable — twice over.** The
direction-correct ranking fills all 20 seats, is nothing like the incumbent
(0/20 overlap), and is **better than the pool it rejects on nothing after
Bonferroni correction**. The gate was inverted and fixing it moves the seats
completely; it moves the *evidence* not at all.

---

## 1. The ranking

Rule: `results/decisions/20261026-c4-selection-rule-v4.md`, sha
`4c2ebc9a4179cfbbd83380dc3a0bb6a7fed0bbac5c550269126b52a7d189723d`.
Freeze taken **before** any scoring: `manifest_frozen.json`, sha
`b2d3bbf7431246ea182d25841039ea4c9d75215b35832fe1126e1c56c9a7d394`.

| quantity | value |
|---|---:|
| pool (records with a recomputed apical L6) | 342 |
| iPTM >= 0.60 alone | 116 |
| L6 < 0 (ACID-OFF) alone | 231 |
| **eligible (both, plus length/shape gates)** | **49** |
| **seats filled** | **20** |
| **shortfall** | **0** |

### Eligibility clauses verified, not assumed

The rule has two numeric gates and several identity gates. The numeric ones are above; the
identity one that matters most is "exact human P02786 189-383 target identity", because it
is what makes the whole apical construct defensible at all. Measured on the 342 complexes
actually scored (`target_identity_check.json`):

| check | result |
|---|---|
| receptor chain (B) sequence == P02786 189-383, extracted from Uniprot | **342/342 exact** |
| receptor chain length | 195 in 342/342 |
| negative control: P02786 1-195 (same length, wrong window) | does not match |
| negative control: P02786 189-343 and 260-414 (wrong lengths) | do not match |

So every scored complex really is the 195-residue apical domain and not a look-alike crop.
Note the scope this does **not** cover: it confirms what was *measured*, and the separate
finding `results/c4_native_context_20261026/FINDING_WRONG_DOMAIN.md` establishes that the
designs were *generated* against a different domain. A correct measurement of the wrong
construct is still the wrong construct, and that is why no file here is submittable.

The rule document quotes 76 for `iPTM>=0.60 AND L6<0`. That is the pre-clash number. My
49 = 76 minus **27 designs that fail the mandatory eligibility clause** "zero interchain
heavy-atom pairs <2.0 A in every available ORIGINAL human apical complex". No deferred
comparisons blocked any seat, so the shortfall is 0 and no row was dropped for missing
data — the shortfall rule is documented but never fired.

| # | design | S | L6 | iPTM | D | C | N |
|---:|---|---:|---:|---:|---:|---:|---:|
| 1 | `tfr1_c4_acid_site1_67` | 42.98 | -0.033 | 0.791 | 1.000 | 0.866 | 1.0 |
| 2 | `tfr1_c4_acid_site2_74` | 32.96 | -0.352 | 0.742 | 0.667 | 0.641 | 1.0 |
| 3 | `tfr1_c4_acid_site1_125` | 27.84 | -0.503 | 0.750 | 1.000 | 0.522 | 0.0 |
| 4 | `tfr1_c4_acid_site2_66` | 27.70 | -0.628 | 0.780 | 0.789 | 0.794 | 0.0 |
| 5 | `tfr1_c4_acid_site1_7` | 27.68 | -0.578 | 0.821 | 0.600 | 0.712 | 0.5 |
| 6 | `tfr1_c4_acid_site1_177` | 27.00 | -0.160 | 0.769 | 0.621 | 0.639 | 0.5 |
| 7 | `tfr1_c4_acid_site2_169` | 26.51 | -0.336 | 0.752 | 0.588 | 0.650 | 0.5 |
| 8 | `tfr1_c4_acid_site1_0` | 25.74 | -1.865 | 0.878 | 0.694 | 0.791 | 0.0 |
| 9 | `tfr1_c4_acid_site2_92` | 24.34 | -0.577 | 0.760 | 0.537 | 0.574 | 0.5 |
| 10 | `tfr1_c4_acid_site2_177` | 22.67 | -1.154 | 0.785 | 0.478 | 0.540 | 0.5 |
| 11 | `tfr1_c4_acid_site2_24` | 22.20 | -0.826 | 0.782 | 0.600 | 0.680 | 0.0 |
| 12 | `tfr1_c4_acid_site1_147` | 21.49 | -1.609 | 0.942 | 0.360 | 0.953 | 0.0 |
| 13 | `tfr1_c4_acid_site1_49` | 21.29 | -0.066 | 0.848 | 0.343 | 0.962 | 0.0 |
| 14 | `tfr1_c4_acid_site1_149` | 20.90 | -0.163 | 0.753 | 0.531 | 0.685 | 0.0 |
| 15 | `tfr1_c4_acid_site1_161` | 20.83 | -2.324 | 0.684 | 0.526 | 0.354 | 0.5 |
| 16 | `tfr1_c4_acid_site2_181` | 20.41 | -0.743 | 0.624 | 0.541 | 0.307 | 0.5 |
| 17 | `tfr1_c4_acid_site2_100` | 20.06 | -0.425 | 0.858 | 0.438 | 0.754 | 0.0 |
| 18 | `tfr1_c4_acid_site1_113` | 19.35 | -1.486 | 0.808 | 0.471 | 0.663 | 0.0 |
| 19 | `tfr1_c4_acid_site2_128` | 18.48 | -0.555 | 0.640 | 0.595 | 0.439 | 0.0 |
| 20 | `tfr1_c4_acid_site2_25` | 17.57 | -1.858 | 0.693 | 0.513 | 0.487 | 0.0 |

Max pairwise TM among the 20 seats is **0.8016**; the pair at
0.8016 (`site1_67` vs `site2_25`) has coverage 0.674, below the 0.80 half of the
near-identical predicate, so it is legitimately admissible. Every other pair is below
0.67 TM.

### 55 of the 100 weight points do nothing

`S = 30P + 25M + 20D + 15C + 10N`, measured across the seats:

| component | weight | value on the seats | discriminating? |
|---|---:|---|---|
| P (pH) | 30 | 0.0 for all 20 | **no — inert** |
| M (mouse) | 25 | 0.0 for all 20 | **no — inert** |
| D (diversity) | 20 | 0.343 – 1.000 | yes, 18 distinct |
| C (confidence) | 15 | 0.307 – 0.962 | yes, 20 distinct |
| N (novelty) | 10 | 0.0 / 0.5 / 1.0 | yes, 3 values |

The rule's two heaviest components — 55% of the score — are structurally silent on this
pool, for reasons that are about **missing evidence, not about the designs**:

- **P = 0** because `P=0.5` requires "a separately implemented **controlled** pH
  estimator". `results/c4_acid_off_scored.json` is a genuine second implementation
  (His-distance/burial switch score, nothing in common with the PROPKA linkage) but it
  records **no controls**, so it cannot grant 0.5. A working control run on that scorer
  is the single cheapest thing that would un-inert 30 points.
- **M = 0** because `M=0` is prescribed for *missing* mouse evidence, and the mouse audit
  covers **the 20 incumbent acid-ON seats and 0 of the 49 acid-OFF eligible designs**
  (verified: intersection = 0). Nothing here says these designs are mouse-incompatible;
  it says nobody looked.

So S is effectively `20D + 15C + 10N`. That is disclosed rather than patched, because
relaxing P or M to make the ranking look better is exactly the post-result tuning the
rule forbids.

### N is a screening flag, not a property of the design

Of the 49 eligible designs, **11 are in the 35-design novelty/provenance set and 38 are
not**. Of the 20 seats, 9 were screened and **11 score N=0 simply because nobody screened
them**. No claim that one seat is more novel than another is supported by this artifact.

---

## 2. The falsifiable test

Same 13 metrics, same two-sided permutation-on-medians, same 20,000 permutations as
`results/c4_independent_check_20261026/` — deliberately, so the two runs are comparable.
Join by `sha256(sequence)`, never by name.

### Controls — all four fired

| control | requirement | measured | fired |
|---|---|---|---|
| SHUFFLE | random relabelling fires at the nominal rate | 6.0% of 200 relabellings at p<0.05 (band 2.1–8.9%) | yes |
| PLANTED | perfect synthetic separation recovered | p = 5.0e-05 = the floor | yes |
| NOISE | an unrelated N(0,1) column is silent | p = 0.252 | yes |
| SEATS | every seat really is L6 < 0 | 20/20 | yes |

### Positive control on my own implementation

Before trusting a null I reproduced the **incumbent** result with this lane's own code:
`incumbent_positive_control.json`. **13/13 metrics match** the reference p-values and
directions, and the "better on none, worse on two" claim reproduces exactly — the two
significant metrics are `core_charged_frac` and `frac_acidic_DE`. The test is not broken.

### Result A — seats vs **all 293** v4-ineligible designs

| | value |
|---|---|
| Bonferroni-significant metrics (13) | **0** |
| directionally better / worse | 11 / 2 |
| smallest raw p | **0.204** (`frac_apolar`) |
| verdict | **(b) indistinguishable** |

### Result B — seats vs the **29** eligible-but-unseated designs (head-to-head inside the gate)

| | value |
|---|---|
| Bonferroni-significant metrics (13) | **0** |
| directionally better / worse | 7 / 6 |
| smallest raw p | 0.0349 (`frac_his`) — survives nothing at 13 tests |
| verdict | **(b) indistinguishable** |

### Robustness — the null is not a median artifact

Re-tested with difference of means, Cliff's delta and AUC (`robustness_v4.json`,
controls fired for all three statistics):

| comparison | min p over 13 metrics x 3 stats | after Bonferroni (39) | conclusion |
|---|---|---|---|
| A: vs all ineligible | 0.041 | 1.00 | null |
| B: vs unseated eligible | 0.005 | 0.195 | null |
| D: acid-OFF seats vs acid-ON seats | 0.013 | 0.488 | null |

Largest effects anywhere are Cliff's delta 0.19 (A) and 0.40 (D, `iface_apolar_frac`) —
"small" to "moderate", on one metric out of thirteen, with nothing surviving correction.

### The direct question: is the corrected set different from the one it replaces?

Seat overlap with the incumbent acid-ON 20 is **0/20** — the sets are disjoint. On the
incumbent's own reporting format (`component_audit_v4.json`):

> acid-OFF seats vs acid-ON seats: better on **1**, worse on **0** at raw p<0.05
> (`frac_acidic_DE`, p=0.036), **0** after Bonferroni.

One metric at raw p just under 0.05 out of thirteen is what chance produces. A completely
different 20 molecules, chosen by the opposite gate, are statistically
indistinguishable on every independent axis measured.

---

## 3. Verdict

**(b) Statistically indistinguishable — and that is the finding.**

The honest reading, stated so it cannot be softened later:

- The sign **was** wrong. That is arithmetic and both advisors confirmed it.
- Fixing the sign **does** replace the portfolio wholesale (0/20 overlap), so this is not
  a cosmetic change.
- Fixing the sign **buys no measurable quality**. The corrected seats beat the rejected
  pool on nothing that survives multiple testing, on either of two comparison groups, on
  four statistics. They are also not worse — which is why this is (b) and not (c).
- Therefore: **the gate direction was a real defect, but the surrogate design pool is
  dead regardless of which way the gate points.** That is the materially stronger
  argument the brief anticipated, and it is what the data says.

What this does **not** establish, and must not be read as establishing:

- That the acid-OFF designs are good. **Nothing measured here is a binding, affinity or
  pH-response measurement.** L6 is a PROPKA estimate on predicted complexes.
- That the ranking is meaningful biology. 55% of its score is inert; 11 of 20 seats carry
  N=0 purely for missing screening; mouse and pH evidence are absent, not negative.
- That v4 should be submitted. **No C4 file on disk is submittable.** v4 is gated
  correctly now but its own rule document records that **zero of 342 designs pass
  acid-OFF AND iPTM>=0.60 AND >=60% native area AND >=10 exposed contacts AND mouse
  identity >=0.90 together** — and separately, the evidence base still scores a domain
  the designs were not generated against
  (`results/c4_native_context_20261026/FINDING_WRONG_DOMAIN.md`).

The construct and file choice remain Ricardo's.

---

## 4. Corrections made inside this cycle (kept, not hidden)

1. **Frozen-file guard named nonexistent paths.** v1 of the manifest pointed at
   `results/submission_v8_c3_20261025.csv` etc. All four guard files exist under other
   names (`results/c3_submission_v8_20261025.csv`, `..._v7_20261025.csv`,
   `..._v6_20261023.csv`). Caught before any scoring ran. v1 preserved annotated as
   `manifest_frozen_v1_SUPERSEDED.json`; the freeze was regenerated and all four prefixes
   now verify.
2. **The dedup predicate was directional.** v1 took TM from the forward alignment only, so
   `site1_67 -> site2_25` scored 0.8016 and the reverse 0.5743 — seat order would have
   decided the verdict, and the ranking seated both while its own pairwise report flagged
   the violation. This is the identical defect Astra caught in the C3 lane
   (`scripts/dedup_rule_v2_symmetric_20261025.py`). Replaced with the repo's frozen
   symmetric predicate: `TM = max(tm_ab, tm_ba)`, `coverage = max(paired_cols_ab,
   paired_cols_ba)/max(len_a,len_b)`. Exactly one seat changed (`site1_149` in,
   `site1_86` out). v1 preserved as `v4_ranking_acid_off_DIRECTIONAL_DEDUP_SUPERSEDED.json`.
3. **Coverage term approximated paired columns.** v2 of that fix used
   `min(non-gap residues per sequence)` instead of paired columns. No seat changed; the
   intermediate is preserved as `v4_ranking_acid_off_v2_NONGAPCOV_SUPERSEDED.json`.
4. **The robustness script's permutation null was centered on the observed value.**
   `|stat_perm - obs| >= |obs|` instead of `|stat_perm - center| >= |obs - center|`. For
   Cliff's delta on 20 vs 293 (null spread ~0.03) that drove all 13 metrics to the p floor
   and would have reported a dramatic false positive in every comparison. Caught by
   noticing 13/13 significant with trivial effect sizes, discarded before being reported,
   fixed with explicit null centers (0 for mean/Cliff, 0.5 for AUC) plus noise and planted
   controls per statistic. The headline test never had this bug — its shuffle control
   calibrated it at 6.0% vs the nominal 5%.
5. **`incumbent_positive_control.py` first labelled "worse on two" as a direction count**
   rather than a raw-p<0.05 count. Recomputed like-for-like; the reference claim
   reproduces exactly.

## 5. Artifacts

| file | sha256 |
|---|---|
| `manifest_frozen.json` (freeze, pre-scoring) | `b2d3bbf7431246ea182d25841039ea4c9d75215b35832fe1126e1c56c9a7d394` |
| `v4_ranking_acid_off.json` (20 seats) | `75236da7cc1ae106c5649cc310a347e8564c4fbd525cfc7f4779abba6d613b3f` |
| `falsification_v4.json` | `1ae16d5e5d4178b61f1ab3231b41b140e1527e824cd4bcf57969a944d74246a2` |
| `robustness_v4.json` | `fae8d5eb22458df2ff9367bbe2436f491edc8965a6e911017de0b146090d9c2d` |
| `component_audit_v4.json` | `e7565a1f3d7e8af3df7f0b51aa632ba126a6fa84ee3bceea53a77a277c09f117` |
| `sensitivity_v4.json` | `5142d43dec66bd27168fab16c4185278fe602ebfdd0589f6b80228c89d8fa06a` |
| `incumbent_positive_control.json` | `a6467e7f838c05b63d5e362033d1792ed649b4776e60be8c05bacfc7dd031859` |
| `target_identity_check.json` | `611d9dfd5e6affecb9f12b87baf2f0d18ce52464577656ccdb6d2402576e4496` |

Scripts: `scripts/c4_v4_freeze_20261026.py`, `c4_v4_score_20261026.py`,
`c4_v4_falsification_20261026.py`, `c4_v4_incumbent_control_20261026.py`,
`c4_v4_robustness_20261026.py`, `c4_v4_component_audit_20261026.py`,
`c4_v4_sensitivity_20261026.py`.

### Reproduce

```sh
cd ~/dev/adaptyv-c4
python3 scripts/c4_v4_score_20261026.py          # ~70 s, writes nothing if frozen inputs moved
python3 scripts/c4_v4_falsification_20261026.py  # ~36 s
python3 scripts/c4_v4_incumbent_control_20261026.py
python3 scripts/c4_v4_robustness_20261026.py     # ~115 s
python3 scripts/c4_v4_component_audit_20261026.py
python3 scripts/c4_v4_sensitivity_20261026.py
```

## 6. If this is worth another $0 cycle

The one measurement that would un-inert 30 of the 100 score points is a **controlled run
of the existing acid-OFF switch scorer** (`scripts/score_acid_off.py`) on the same
positive/negative controls the PROPKA linkage already uses. That is CPU-only, needs no new
model, and would let P reach 0.5 where the direction is right — turning the largest
inert weight into a discriminating one. It is a proposal, not a result; I have not run it.
