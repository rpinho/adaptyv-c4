# C3 sprint 8h — 2026-10-25 — deliverable index

**READ THIS FIRST — this index was stale until 20:11 ET.** It previously recommended
`c3_candidate_v9_B-plddt` (`704988ca`) from the 12:10 cycle and was never updated as the sprint
progressed. **That is no longer the recommendation.** Everything below the "Superseded" line is
retained as history and must not be acted on.

$0 spent on GPU. **SUBMITTED 2026-10-26 07:16 EDT** — candidate D (`ba5e3789`) to Challenge 3; see `SUBMISSION_RECORD_20261026.md`. This index line originally read "nothing submitted" and was true when written. v8 byte-identical throughout.

## CURRENT RECOMMENDATION (20:11 ET)

**`c3_candidate_D_swap18_20261025.csv`** — sha **`ba5e3789c0dfdf8c…`** — 18 rows, validator PASS
(`--slots 18`). Full reasoning: `RECOMMENDATION_D_FINAL_20261025.md`.

Decision path: `704988ca` (12:10) → … → `5313ce40` core18 → **`ba5e3789` D**.
Panel record: Kimi + Astra voted D on the real numbers; Opus 5.5 + Fable 5.1 voted D on a
paraphrased prompt and did **not** see D's figures. Treat it as 2 evidence votes + 2 analogy votes.

## Every candidate file, measured the same way

All figures recomputed from raw per-model ColabFold score JSONs, joined on `sha256(sequence)`.
`coverage` = share of rows that have a 5-model run on disk. **A row without one is unmeasured, not good.**

| id | sha | rows | coverage | mean min-5 | worst | on 0.85 floor | v8 rows retained |
|---|---|---:|---:|---:|---:|---:|---:|
| **D** | **`ba5e3789`** | 18 | **100%** | 0.8967 | 0.86 | 0 | 0 |
| core18 | `5313ce40` | 18 | 100% | 0.8961 | 0.85 | 1 | 0 |
| core17 | `8237096d` | 17 | 100% | **0.8988** | 0.86 | 0 | 0 |
| superset20 | `9fd48b44` | 20 | 100% | 0.8925 | 0.85 | 1 | 0 |
| B-plddt | `704988ca` | 20 | 40% | 0.8975* | 0.86* | 0* | 12 |
| v9_A | `b44a7854` | 20 | 35% | 0.9114* | 0.90* | 0* | 13 |
| v9_B | `643490da` | 20 | 40% | 0.9100* | 0.90* | 0* | 12 |
| v8 (frozen fallback) | `631ed0bf` | 20 | 0% | – | – | – | 20 |

\* measured only over the rows that have 5-model data — the missing rows are inherited v8 designs,
which have **no** 5-model run at all. These means are not comparable to D's and should not be quoted
as if they were. This is why the v9 challengers look better on paper than they are supported.

**Open side-by-side:** `8237096d` (core17) is D minus its weakest row, so it has a higher mean. It is
**not** recommended because 17 < 18 seats with no evidence the 18th hurts, and D is the artifact that
has been independently verified end-to-end. Recorded so the choice is visible, not hidden.

## Superseded (history — do not act on)

### Superseded recommendation (12:10 ET)

**`c3_candidate_v9_B-plddt_20261025.csv`** — sha `704988caef6ee56f…` — 12 v8 survivors + 8
replacements chosen under Astra's stated rule (**binder pLDDT, then identifier**) with pLDDT measured
on the **binder chain only**.

Cycle 2 recommended variant B (`643490da…`). B is **unchanged and still on disk**; it is the fallback
choice if the tie-break correction is rejected. The move to B-plddt is caused by fixing a measurement
error, not by a change of preference.

| | v8 (frozen) | B (cycle 2) | **B-plddt** |
|---|---:|---:|---:|
| rows | 20 | 20 | 20 |
| rows that clash (AF2m) | **7** | 0 | 0 |
| rows below 0.85 iPTM | 2 | 0 | 0 |
| iPTM median / min | 0.885 / 0.23 | 0.92 / 0.91 | 0.915 / 0.89 |
| mean binder pLDDT, raw | – | 96.84 | 97.42 |
| mean binder pLDDT, length-adjusted | – | +0.38 | **+0.44** |
| mean admitted binder length | – | 68.8 aa | 83.9 aa |
| admitted pod mix (NOT families — see caveat) | mixed | p1,p2,p4,p5 | p1,p2,p3,p4,p5 |
| **epitope clusters (19 scored)** | **2** | 1 | 1 |
| **novelty: total strict ref hits** | **944** | 7 | **7** |
| **novelty: max identity (strict)** | **0.662** | 0.351 | **0.351** |
| admitted rows screened for novelty | 20/20 | 12/12 | **12/12** (4 newly screened) |

### Three corrections adopted at 12:45 ET (from worker B's independent audits)
1. **"Families p1–p5" is not a backbone-family axis.** The prefix is a generation *pod* label, and it is
   almost perfectly confounded with the target window: `p1`→52–64 only, `p2`→65–80 only, `p3`→81–100
   only; `p4` is the balanced H8 pod and `p5` the only H4 pod (verified centrally over all 235 rows of
   `sequences.csv`). B-plddt's entire move from B is a **target-window shift toward 81–100**, not a
   backbone-family gain. Genuine backbone diversity is fine and separately verified: **12 admitted rows
   = 12 distinct backbone parents**, none reused.
2. **"944 → 7" is chiefly row deletion, not redesign.** 944 is dominated by `ecd_12` (641/0.662) and
   `ecd_29` (296/0.438) — rows the challenger removes. 12 of 20 rows are unchanged sequences. The
   direction is right; the 135× framing is not.
3. **Novelty had a real coverage gap, now closed.** The earlier screen covered 36 designs and included
   only 8 of the 12 admitted rows — so 4 admitted rows (**3 of them in this recommended file**) had
   never been screened. Worker B re-ran the identical strict pipeline with a positive control
   (`ecd_12` reproduces its per-database counts exactly): **all 12 admitted rows are strict-clean, 0
   hits in all four databases, 0 hits in the antibody-specific `ta`/`plabdab` DBs even at e ≤ 100.**
   Restricted to alignments covering ≥70% of the query, no admitted row exceeds **0.312** identity.

### The epitope caveat, stated properly
v8 is **already 18/19 in a single epitope cluster**; its only distinct epitope is
`c3_glp1r_ecd_29`, the 0.23-iPTM row. Replacing all 8 flagged rows therefore **does lose v8's sole
distinct cluster**. It does not collapse a diverse portfolio — v8 was concentrated to begin with — and
retaining a 0.23-iPTM design solely to count a second cluster is not supported by the evidence. No
claim of epitope diversity is made for B or B-plddt. Both advisors agree on this framing.

Variant A (7 replacements, keeps `ecd_29`) and A-bc (6 replacements, keeps a 0.31-iPTM 5-clash row)
are retained for contrast. **A-bc is not a recommendation.**

## Evidence chain
1. `native_vs_af2m_clash.py` → `native_vs_af2m.json` — the headline finding. Self-check reproduces
   all 8 published AF2m numbers exactly, so the measurement is trustworthy *before* it is
   interpreted. **5 of 6 `c3_bc_*` clash in BindCraft's own fold**, not just in AF2m.
2. `build_variants.py` → `challenger_variants_provenance.json` — the exact selection rule and the
   funnel (83 pool designs → 30 clash-clean on all five models → 25 iPTM ≥ 0.87 → 25 refold PASS →
   **19** admissible after the frozen diversity rule).
3. `verify_candidates.py` → `verifier_report.json` — V1–V8 structural verification.
3b. `scripts/c3_epitope_baseline_v8_20261025.py` → `epitope_baseline_v8.json` — scores **frozen v8**
   under the identical epitope rule; self-check reproduces variant B's published numbers exactly.
3c. `scripts/c3_tiebreak_corrected_20261025.py` → `tiebreak_corrected.json` — fixes the pLDDT scope
   bug (whole-complex vs binder-only) and re-applies Astra's rule over the full 19-design pool; its
   self-check reproduces the shipped 8 under the old rule.
3d. `scripts/c3_build_b_plddt_20261025.py` — builds the recommended CSV; deterministic (verified by
   rebuilding to an identical sha).
4. `acceptance_test.py` → `acceptance_report.json` — **Astra's pre-committed six conditions**,
   written before any candidate was finalized.
5. `provenance_table.py` → `PROVENANCE.md`, `provenance_table.csv`, `source_hashes.json` — per row.
6. Target identity — verified at the sequence level: both screens used the identical target (see Comparability).
7. `novelty_compare.py` → `novelty_compare.json` — re-groups already-measured rows only (**no new
   search**), so it under-covers the admitted set. Read it as a pointer, not the headline; the
   authoritative novelty evidence is item 7b. Residual exposure does fall 944 → 7 and max identity
   0.662 → 0.351, but this is **chiefly because the two most exposed rows were dropped**
   (`ecd_12` 641/0.662, `ecd_29` 296/0.438), not because the design method improved.
7b. **`novelty/FINDING.md`** — the audit that closes that gap. `novelty_compare.py` covered 36 designs;
   four of the eight admitted rows had **never been screened**. Re-run on the identical four databases
   and parameters: **all 12 admitted rows are clean (0 hits) in the strict regime**, and the control
   reproduces. Three corrections to how this file should be described: (i) "944 → 7" is dominated by
   **row removal**, not by more-novel replacements; (ii) "backbone families p1…p5" are **pods**
   confounded with the target window and an H8/H4 tag — the real backbone count is **12 distinct
   parents**; (iii) the 8 replacements are a **single-generation-method, helix-rich block**, so the
   file is more homogeneous than v8 in three specific senses. Novelty direction is unchanged
   (the challenger is **more** novel, not less); the wording needs fixing.

## Comparability (the objection that could have killed this)
Verified, not assumed: the depth arm and the v8 screen used a **byte-identical 111-residue target**
(sha `5c2bae2d…` in both `af2m_input.csv` and `c3_screen_input_20261024.csv`, matching the depth-arm
manifest), the same predictor family (ColabFold alphafold2_multimer_v3), 5 models, 3 recycles.

## What this does NOT establish
- Nothing here measures affinity or agonism. Clean clash ≠ binding.
- BindCraft's own fold is **a computational output**, so "clashes in both predictions" is a risk
  signal, not an experimental failure. No "nearly certain non-binder" language is used anywhere.
- The challenger is better on **audited computational criteria only**.
- **Nothing here establishes epitope diversity.** Both B and B-plddt are single-cluster. Family
  coverage (p3) is a diversification argument, not a demonstrated performance gain.
- pLDDT differences of ~0.5 points are a **deterministic tie-break**, not evidence of better binding.
- **The pLDDT edge is mostly a length artifact.** Binder pLDDT rises ~1.01 points per 30 residues
  (r = 0.477). B-plddt's admitted binders average 83.9 aa vs B's 68.8, so of the raw +0.58 pLDDT gap
  only about +0.06 survives length adjustment (~90% is length). B-plddt is recommended because it
  applies the pre-committed rule, **not** because it is demonstrably better. Sources:
  `plddt_length_confound.json`, independently reproduced from `sequences.csv` lengths.

## Decisions
`results/decisions/20261025-1145-c3-native-clash-and-challenger.md` (full record),
`20261025-1142-c3-sprint-astra.raw.md` (Astra verbatim),
`20261025-1150-c3-sprint-kimi.raw.md` (Kimi verbatim).

Both advisors independently ruled: replace the flagged rows, spend $0, keep v8 as baseline, require
mutual non-redundancy. Kimi overrode my initial preference to keep the 0.23-iPTM row; I accepted.

## Reproduce
```
cd ~/dev/adaptyv-c4
python3 results/c3_sprint_8h_20261025/native_vs_af2m_clash.py   # must print PASS self-check
python3 results/c3_sprint_8h_20261025/build_variants.py
python3 results/c3_sprint_8h_20261025/verify_candidates.py      # exit 0
python3 results/c3_sprint_8h_20261025/acceptance_test.py
python3 results/c3_sprint_8h_20261025/provenance_table.py
```

## Same-standard comparison (the honest way to state the clash claim)
Worker A correctly flagged that v8 was screened on **rank-001 only** while the pool was screened on
**all five models**. That is a real asymmetry in the *screening rule*, so the portfolio comparison is
stated both ways:

| standard | v8 | challenger B |
|---|---:|---:|
| rank-001 clash rows (v8's own standard) | **7** | **0** |
| rows with unknown rank-001 status | 0 | 0 |
| clash-clean on all five models | **not measurable** (only rank-001 archived) | 20/20 |

**The defensible claim:** at v8's own screening standard — rank-001, same predictor, same target —
the challenger has 0 clashing rows against v8's 7, with no unknowns. The all-five-model screening of
the pool is a *stricter* standard applied only to the replacements, so it cannot manufacture a pass.

**The residual gap that stays open (by decision, not oversight):** the 12 v8 survivors cannot be
re-screened under the all-five-model rule, because only their rank-001 model was ever archived.
Closing it needs ~12 designs × 4 extra models of AF2m inference — i.e. GPU. Both advisors ruled $0
for this sprint, so it is left open and *named*. It affects the **baseline**, not the challenger.
