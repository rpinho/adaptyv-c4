# C4 draft candidate file — 13 designs, unranked — DRAFT, NOT SUBMITTED

Built 2026-10-26T19:47:39Z from the real clock. This is a **draft candidate file for Ricardo's
own review**, built on the unanimous instruction of four reviewers (Fable 5, Opus 5.5, Astra,
Kimi; Fable 5.1 blocked by content classifier). **It has not been submitted, uploaded, or shown
to anyone outside this repo. Ricardo submits by hand, if and when he chooses.**

File: `results/c4_draft_13_unranked_20261026.csv`
sha256: `d673bb52bfe536ad9d333afe48b4e141e44c391272c4a08255ea9963483f14a7`

## Selection

Rule v5 (`results/decisions/20261026-c4-selection-rule-v5.md`) removed the pH-linkage quantity
L6 from eligibility and ranking. Rebuilt **without any L6 criterion**:

    342 designs -> iPTM >= 0.60 (surrogate construct): 116
                -> native footprint-area retained >= 0.60: 24
                -> >= 10 natively-exposed footprint contacts: 13

All 13 are included. No padding to 20 (unanimous reviewer instruction: free submission is not a
reason to fill the cap with filler). 13 unique sequences; 0 overlap with
`c3_submission_v8/v7/v6_20261025` or `submission_v9.csv`; every binder sequence was extracted from
its complex PDB and verified against its recorded `sequence_sha256` — **13 of 13 match**.

**Order is random** (seed 20261026), per the unanimous reviewer instruction that these 13 carry
no defensible ranking. If the submission form requires an order, this file's order should be
read as administrative only, not as a preference.

## MANDATORY disclosure, if this is ever submitted

All five reviewers converged on this content; use it verbatim or expand it, do not shorten it:

> Every selection metric reported for these designs (iPTM, native footprint-area retention,
> natively-exposed contact count) was computed on a 195-residue apical surrogate construct that
> the designs were **not** generated against — `results/c4_native_context_20261026/
> FINDING_WRONG_DOMAIN.md` establishes from the generator's own hotspot-atom selectors that they
> were generated against a different receptor region (author ~651-741). The generator was biased
> toward the wrong pH direction at the sequence level (`C_SEQ_PROB_HIS: 0.15` aimed at
> carboxylate hotspots). No valid pH-dependence estimator exists in this pipeline: the scorer that
> was used (L6) was withdrawn from selection by rule v5 after being shown to be noise-dominated
> near its own gate and uncalibrated at its zero point. No design has any measured binding, KD, or
> pH response. Mouse cross-reactivity support is sequence conservation only, not measured binding.

## What this is NOT

- **Not a claim that these 13 are good binders.** They are the designs that pass three
  surrogate-construct thresholds, nothing more.
- **Not validated against the official 640-aa target.** The construct-scope decision
  (`results/decision_requests/20261026-c4-construct-scope.md`) remains open; no design here has
  been scored against the real receptor.
- **Not a ranking.** Explicitly unordered.
- **Not submitted.** This document and the CSV beside it are a draft for review.

## Reviewer objections on record, unresolved

- The 0.60/0.60/10 thresholds were set while L6 carried the selection load and were never
  stress-tested as the sole criteria (Fable 5).
- The 5 designs admitted only after L6's removal were never independently reviewed; "uninformative"
  is not "equivalent to the designs it rejected" (Opus 5.5).
- Only 1 of 13 contacts the anchor residue His707; matching rule v5 is not validation (Kimi).

## Reproduce

```sh
cd ~/dev/adaptyv-c4
cat results/decisions/20261026-what-ships-outcome.md
cat results/decisions/20261026-c4-selection-rule-v5.md
```
