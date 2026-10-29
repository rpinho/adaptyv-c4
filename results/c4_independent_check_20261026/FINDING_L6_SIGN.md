# FINDING: the v4 C4 ranking gates the pH switch in the direction `fa98723` identified as wrong

Written 2026-10-26, from the real UTC clock (`date -u` at write time: see `cycle_receipt.json`).
CPU only. No GPU, no submission, no organiser contact. No ranking was built. No frozen file touched.

> **AMENDMENT, same day, after reviewer feedback and a root-cause finding.**
> The gate inversion below **stands**, on the scorer's definition plus its three firing sign
> controls. It does not depend on the histidine correlation. Three corrections:
> - **Retracted:** "Histidine is the *only* residue whose pKa titrates between pH 7.4 and 6.0."
>   That is wrong. Buried glutamate pKa values span roughly 5.2–9.4 experimentally, and PROPKA
>   assigns environment-specific pKa values, so Asp/Glu can titrate in this window. The argument
>   that bulk acidic content "cannot supply the switch" is withdrawn.
> - **Downgraded:** the His/apolar/acidic correlations are **not independent corroboration**. L6 is
>   computed from the same structures and ionisable groups, so the association is partly circular;
>   scaffold family, burial and design selection could also generate it. They are consistent with
>   the reading, not evidence for it.
> - **Root cause found:** the acid-ON direction was baked into *generation*, not just gating. Both
>   RFD4 configs set `C_SEQ_PROB_HIS: 0.15` on the binder and aim every hotspot at a carboxylate —
>   a deliberate His–carboxylate pairing, exactly what `fa98723` called the wrong direction. And the
>   designs target the helical domain, not the apical domain, so the L6 values here were computed on
>   a non-target receptor. See `results/c4_native_context_20261026/FINDING_WRONG_DOMAIN.md`.

## Summary

The v4 ranking selects designs with **apical L6 > 0**. Under the scorer's own documented
convention that means *tighter binding at pH 6.0 than at pH 7.4* — bind in acid, release at
neutral. Challenge 4 asks for the opposite: **measurable KD at pH 7.4 and no detectable binding
at pH 6.0**. Ricardo committed that correction himself on 2026-10-24 (`fa98723`, "CRITICAL FIX:
Correct pH direction to acid-OFF"). The L6 scoring path written the following day reintroduced
the inverted direction, and v4 is built on it.

This is defect #1 from the handoff reappearing in a different code path after being fixed.

## The contradiction, with sources

| source | written | says C4 wants |
|---|---|---|
| `fa98723` commit message (author: Ricardo Pinho) | 2026-10-24 10:08 | **bind pH 7.4, NO binding at pH 6.0** ("the opposite of what we were scoring") |
| `C4_C5_ORCHESTRATION_BRIEF_20261026.md` §2, objectives verified from the live page | 2026-10-26 08:50 | **measurable KD at 7.4, no detectable binding at 6.0** |
| `scripts/score_acid_off.py` scoring logic (created by `fa98723`) | 2026-10-24 | His-driven acid-OFF: reward His near target Arg/Lys and buried interface His; penalise His near target Asp/Glu |
| `scripts/c4_linkage_audit_20261025.py` docstring | 2026-10-25 06:04 | "For C4's conditional-binding objective we need **L > 0** (tighter binding at the lower pH)" |
| `scripts/score_c4_acid.py` | earlier | "Positive L = net proton uptake on binding = acid-preferential (**correct for C4**)", at pH 5.5/7.4 |
| `HANDOFF_TO_CLAUDE_20261026.md` §4 defect 1 prose | 2026-10-26 08:30 | C4 is "bind 5.5/6.0, release 7.4" — **contradicts the brief and the commit** |

Four sources say acid-OFF (L < 0). Two scripts and one line of handoff prose say acid-ON (L > 0).
The two that say acid-ON are the ones v4 inherits.

## The L6 convention is internally consistent — the gate direction is what is wrong

`scripts/c4_linkage_audit_20261025.py` computes the Wyman linkage

    L = log10 Ka(pH_low) - log10 Ka(pH_high),  pH_low = 6.0, pH_high = 7.4

so **L > 0 means tighter binding at pH 6.0**. Its hard-coded sign controls all fire correctly
and are recorded in `results/c4_linkage_apical_20261025/linkage_apical.json`:

| control | L | expected | fired |
|---|---:|---|---|
| null (identical bound/free pKa) | 0.0000 | exactly 0 | yes |
| acid-ON (pKa raised 1.0 on binding) | +0.5918 | > 0 | yes |
| acid-OFF (pKa depressed 1.0 on binding) | −0.4539 | < 0 | yes |

So the scorer is not broken and its sign is not ambiguous. The defect is that v4's rule
(`"selection_rule": "iPTM>=0.60 AND apical L6>0, ranked by apical L6 desc then iPTM desc then name, take 20"`)
keeps the **acid-ON** half and ranks by how strongly acid-ON a design is.

A separate, already-documented label inversion in the same path is **not** a defect: the scorer is
called with binder label "B" and target label "A" while the files have chain A = binder and
chain B = receptor. The recompute script documents this at
`scripts/c4_linkage_apical_recompute_20261025.py:45-50` and the linkage sums over titratable
groups on both chains, so the labels only name temp files. Verified, not a finding.

## Measured consequence

Pool = 342 apical complexes with recomputed L6 (`linkage_apical.json`, sha
`f36cc0780430a00b16b9ca98b033c3d349e9361654a37d76bdb7536089579737`).

| quantity | value |
|---|---:|
| designs with L6 > 0 (acid-ON) | **111 / 342** |
| designs with L6 < 0 (acid-OFF) | **231 / 342** |
| designs with L6 exactly 0 | 0 |
| v4 seats with L6 > 0 | **20 / 20** |
| v4 seats with L6 < 0 | **0 / 20** |
| eligible under v4's rule (iPTM ≥ 0.60 **and** L6 > 0) | 40 |
| eligible in the acid-OFF direction (iPTM ≥ 0.60 **and** L6 < 0) | **76** |
| **overlap between v4's 20 seats and the acid-OFF eligible set** | **0** |

If acid-OFF is the correct direction, all 20 v4 seats are drawn from the wrong 111 of 342, and
the direction-correct eligible pool is nearly twice as large (76) and entirely disjoint from them.

## Independent corroboration: the histidine signal

This was found from the opposite end, by a signal that entered none of v4's selection. Across all
342 designs, L6 correlates with sequence composition (Spearman, 20,000-permutation p; the floor is
p = 5e-5 = 1/20001):

| correlation | rho | p |
|---|---:|---:|
| L6 vs **histidine fraction** | **−0.3824** | 5e-5 |
| L6 vs interface apolar fraction | −0.3611 | 5e-5 |
| L6 vs acidic (D+E) fraction | +0.2437 | 5e-5 |

All three survive Bonferroni correction for the 26 tests run (alpha = 0.00192).

Histidine is the **only** residue whose pKa (~6.0–6.5) titrates between pH 7.4 and pH 6.0, and it
is the residue `fa98723`'s acid-OFF scoring is built on. Asp/Glu (pKa 3.65/4.25) are essentially
fully deprotonated at both pH values and cannot supply a 7.4→6.0 switch. So ranking on high L6
selects **away from** the one residue that can physically deliver the required response, and
**toward** bulk acidic content that cannot. Histidine content by stratum:

| stratum | median His fraction | n |
|---|---:|---:|
| L6 < 0 (acid-OFF half) | **0.1053** | 231 |
| L6 > 0 (acid-ON half) | 0.0778 | 111 |
| v4's 20 seats | 0.0767 | 20 |
| natural human TfR1 apical 189–383 (reference) | 0.0205 | 1 |

## The repo's own acid-OFF scorer agrees

`results/c4_acid_off_scored.json` (343 designs, created by `fa98723`) scores the His-driven
acid-OFF mechanism directly. Joined to the v4 seats:

| | median | min | max |
|---|---:|---:|---:|
| v4's 20 seats, acid-OFF switch score | **0.5** | **−4.0** | 10.0 |
| whole pool, acid-OFF switch score | 0.5 | — | 22.0 (p90 = 6.4) |

By the direction-correct scorer the v4 seats sit exactly at the pool median, and at least one seat
scores **−4.0** — actively wrong-direction. The pool contains designs up to 22.0 that v4's gate
excluded.

Note also that frozen `results/submission_v9.csv` — the file `fa98723` built with acid-OFF
scoring — is recorded in the handoff as "mean L6 −0.919, 5/20 acid-ON". Under acid-OFF semantics a
negative mean L6 is the **correct** sign, and "5/20 acid-ON" means 15/20 are direction-correct.
Part of the stated case against v9 was an inverted reading of its sign. This does not make v9 good
(mean iPTM 0.418 and 10/20 passing the 2.0 Å gate are independently weak) but the sign element of
that rejection does not stand.

## What is NOT claimed

- **No binding, KD, or pH response is measured here.** L6 is a PROPKA-derived computational
  estimate, and the handoff records a prior "PROPKA calibration failure" (`726d9b5`). A design with
  the correct L6 sign is not thereby a pH-selective binder.
- **I did not read the official challenge page.** The acid-OFF direction rests on Ricardo's own
  `fa98723` commit and on the brief's statement that the objectives were verified from the live
  page on 2026-10-26. One look at the page, or one line from Ricardo, settles it definitively.
- **No replacement ranking was built.** The pre-registered rule requires the eligibility pool and
  model inventory to be frozen by sequence hash before any ranking, and v3 of that rule is itself
  written around mouse evidence, not around a corrected pH direction. Quantification only.
- The 231 acid-OFF designs are **not** thereby good. They are unscreened on clash, iPTM, mouse
  compatibility, foldability and novelty in that direction.

## Consequence for the cycle

1. **No C4 file on disk is submittable.** v4 is gated in the direction the challenge does not
   want; v9 has the right sign but weak confidence and a half-failing clash screen.
2. **Any future ranking must fix the direction first**, then freeze eligibility by sequence hash,
   then rank. The selection-rule document (v3) does not currently specify the L6 sign the gate
   should use; it says only "controlled human apical L6 > 0 (direction gate only)". **That line
   encodes the inverted direction and needs an amendment with its own reason, real timestamp and
   new SHA before it is used.**
3. The mouse audit's conclusions are unaffected — it measured geometry on fixed poses and made no
   pH claim — but it characterised the **acid-ON** 20. If the direction flips, the 20 molecules
   whose mouse compatibility matters are a different 20, and that work would need redoing on the
   corrected set.

## Reproduce

```sh
cd ~/dev/adaptyv-c4
date -u '+%Y-%m-%dT%H:%M:%SZ'
# the sign convention and its controls
sed -n '1,25p' scripts/c4_linkage_audit_20261025.py
python3 -c "import json;d=json.load(open('results/c4_linkage_apical_20261025/linkage_apical.json'));print(d['sign_controls'])"
# v4's rule, verbatim
python3 -c "import json;d=json.load(open('results/c4_ranked_v4_apical_L6_20261025/c4_acid_ranked_v4_apical_L6_20261025.json'));print(d['selection_rule'])"
# the committed correction
git show fa98723 --stat
# the counts in this document
python3 scripts/c4_independent_analysis_20261026.py --out /tmp/c4_reanalysis.json
```
