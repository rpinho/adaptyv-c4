# FINDING: C3 has NOT been scored, and no feedback can arrive before C4 closes

Checked 2026-10-26T16:11:29Z from the real clock. Read-only browsing of pages Ricardo is already
signed in to. **No submission, no upload, no form, no organiser contact.** Requested directly by
Ricardo. $0.

## The question

Ricardo: "when we submitted C3 this morning, the format of the website seemed like it was oriented
towards giving feedback, like giving scores. Can you actually check and see if they scored our C3
submission? Because if they do then there is value in submitting once per day for C4."

## Answer: no scores exist yet, and none can exist inside the C4 window

### 1. His own C3 submission carries no scores

`https://proteinbase.com/collections/anthropic-adaptyv-protein-design-competition-submission-1-1bZrjd_jHr`
("Submission 1", 18 designs, badge `DESIGNATED`). Every scoring field is empty for all 18:

| column | value for all 18 |
|---|---|
| Target | `--` |
| ipSAE | `--` |
| pLDDT | `--` |
| Expressed | `--` |
| Mol. Weight | `--` |

About panel, verbatim: **"18 proteins designed"**, **"0 proteins validated"**,
Binding Strength "Unknown 100%", Expressed "Unknown 100%", Classification "Unknown 100%".

**Ricardo's read of the UI was correct** — the site IS score-oriented. The sort control even
defaults to "Boltz2 ipSAE". The columns are simply unpopulated: they fill from wet-lab validation,
not from an upload-time scorer. There is no automatic score returned on submission.

### 2. The challenge timeline puts C3 at stage 2 of 5

From the C3 page, verbatim stage labels: `Submissions open` -> **`NOW / Submissions closed`** ->
`Selection in progress` -> `Wet-lab validation in progress` -> `Results published`.
Selection has not finished. Nothing has been synthesised.

### 3. FAQ 9, verbatim

> "Once designs have been selected for a problem and DNA has been ordered, we will release a
> Collection on Proteinbase containing the designs that were selected. We estimate this Collection
> will be available **one week after the final competition challenge**."

The final challenge is C5, closing Sunday 2026-11-08. So the selection list is expected
**~2026-11-15** — two weeks after C4 closes.

### 4. FAQ 12, verbatim

> "We hope to complete experimental validation for each problem approximately one month after the
> submission window closes. **The first set of data will be released in December 2026.**"

### 5. FAQ 6, verbatim (deadline confirmation)

> "Deadlines: Challenge 1: Wednesday, October 7; Challenge 2: Sunday, October 18; Challenge 3:
> Sunday, October 25; **Challenge 4: Sunday, November 1**; Challenge 5: Sunday, November 8."

23:59 AoE (UTC-12) on Sunday 2026-11-01 = **2026-11-02 11:59 UTC**, which matches the
orchestration brief exactly. The challenge card's "Oct 19 - Nov 1" is the AoE date, not a
discrepancy.

## Consequence for the "submit once per day" idea

The premise does not hold. There is **no feedback channel that returns anything during a
submission window**:

| information | earliest available | vs C4 close (2026-11-02) |
|---|---|---|
| per-submission automated score | **does not exist** | — |
| which designs were selected | ~2026-11-15 (FAQ 9) | 13 days too late |
| experimental measurements | December 2026 (FAQ 12) | ~1 month too late |

So repeated submission yields **zero** information gain. Two further points of fact, not opinion:

- The design cap is **per challenge, not per day**: "Up to 40 for Track 1, up to 20 for Tracks 2
  and 3." Daily submission does not enlarge it.
- The collection carries a single `DESIGNATED` badge. What repeated submission does to that
  designation was **not tested**, and must not be tested by experiment — doing so would mean
  uploading, which is prohibited and is Ricardo's to do by hand regardless.

## What is NOT claimed

- **This does not say whether to submit anything for C4.** That is Ricardo's decision with the
  four decision models (Fable, Opus 5.5, Kimi, Astra). This document only removes one argument
  from the table: "submit daily to harvest feedback" has no feedback to harvest.
- Pages were read once each on 2026-10-26. Adaptyv states the FAQ is updated during the
  competition, so these quotes are a snapshot, not a permanent fact.
- C5 was observed incidentally as still "Not yet announced" (Oct 26 - Nov 8, locked). This was not
  a deliberate second C5 poll; the day's single permitted check was spent at 12:46:00Z and returned
  the same text.
