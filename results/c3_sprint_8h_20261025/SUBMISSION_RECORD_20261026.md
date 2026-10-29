# C3 SUBMISSION RECORD — 2026-10-26

**Status: SUBMITTED.** Challenge 3 (GLP1R agonist), Anthropic x Adaptyv Protein Design Competition.

## Submitted artifact

| field | value |
|---|---|
| file | `results/c3_sprint_8h_20261025/c3_candidate_D_swap18_20261025.csv` |
| sha256 | `ba5e3789c0dfdf8cf44e087f0d281b37335fe71fd370d0b64d11ea09d6d1f4d2` |
| size | 2,080 bytes |
| rows | 18 designs + header; columns `name,sequence,molecule_class` |
| line endings | CRLF (validator PASS, `--slots 18`) |
| backups | `~/adaptyv_c3_backup_D_swap18_ba5e3789.csv`; `results/c3_sprint_8h_20261025/backup/` |

Methodology submitted verbatim from
`results/c3_sprint_8h_20261025/independent/METHODOLOGY_D_swap18_20261025.txt`
(sha256 `6da47973814530c4aae581fe4a5c41edb5c1b159990e87a4598ee954e12d49db`, 2,669 bytes).

## Portal result

- Banners: **Submission successful**. Portal said the submission to the competition was successful.
- Confirmation: <https://proteinbase.com/collections/anthropic-adaptyv-protein-design-competition-submission-1-1bZrjd_jHr>
- Collection record shows **Submitted to** Anthropic x Adaptyv Protein Design Competition, **18 proteins**,
  all private to Ricardo Pinho.
- Record UUID `01a177ec-c245-0eb2-10a9-5d2f10ad7984`; the embedded 48-bit timestamp decodes to
  **2026-10-26 07:16:07 EDT (11:16:07 UTC)**. This is derived from the record ID, not a captured clock.
- Deadline was 2026-10-25 23:59 AoE = 2026-10-26 07:59 EDT. Submission preceded it by about 43 minutes.

## Portal-side verification at submission

| check | observed |
|---|---|
| Challenge selected | `SUBMITTING TO CHALLENGE 3` |
| Parsed / included | 18 / 18; 0 needing fixes; 0 invalid |
| Novelty check | 18 / 18 `Passed` |
| Methodology textarea | populated with the submitted text |
| Terms checkbox | checked (`#accept-competition-terms` = `on`, `aria-checked=true`) |
| Submit control | enabled, then returned the success page |

## Frozen artifacts — unchanged

`shasum -a 256 -c results/c3_sprint_8h_20261025/frozen_hashes.txt` returns **OK** for all four:

- v8 `631ed0bf6b10e0f8fbf64bd5d9730e294da15f0436af20f52a42bafc92229e12`
- v7 `1d99f843b9d5af6cedd74e5ccfbf6cfd260f26e5461b2bffb0a4ee2204349c55`
- v6 `8e9b2d741586e8b360fbed401b4de77317d579e183409cf4928de42661096891`
- submission_v9 `31d98a4133c558e1a3a9f582a607383bca8d7aad1b8bd820e684cb8fa231eadd`

## Why this candidate

Panel ratified D on post-audit evidence (Kimi receipt `e328d8c1`); the pre-registered decision rule
(`2ac0f66f`) never fired. The portal itself enforces the de novo rule in its page copy, and D is 18/18
computationally generated `p_` designs versus v8's 8/20 BindCraft. Full comparison of every candidate
under identical measurement: `README.md` in this directory.

## Standing scientific limits (as stated in the submitted methodology)

- Agonism is **unproven**; no design here tests activation or blockade.
- Predicted groove contact is a binding hypothesis, not a demonstration.
- All 18 designs share one RFD4 generation lineage, so their errors are correlated.
- Experimental binding, affinity, and function remain unmeasured.

Nothing in this record asserts a wet-lab result.
