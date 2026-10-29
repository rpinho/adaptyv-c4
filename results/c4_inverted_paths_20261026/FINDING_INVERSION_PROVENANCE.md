# FINDING: the inversion has a single origin commit, and the "critical fix" never reverted it

Written 2026-10-26 from the real clock. CPU only, $0, read-only. No GPU, no submission, no
organiser contact. Frozen files 4/4 verify.

## Origin

`git blame` on the contradictory lines points at one commit.

| commit | date | author | line 128 of `scripts/ph_linkage.py` |
|---|---|---|---|
| `4280bc7` | 2026-10-20 | Ricardo Pinho | `L <= -2` |
| **`75c834b`** | **2026-10-23** | Ricardo Pinho | **`L >= 2.0`** |

`75c834b` is titled *"Fix C4 scoring direction: positive = proton uptake = bind at pH 5.5"*.
It changed two files and 12 lines in `ph_linkage.py`, flipping the pass flag, the printed banner
and the sort order to the acid-ON direction. It did not touch the module docstring, which is why
line 17 still carries the original `L <= -2` convention. **The file has held both conventions
since 2026-10-23.**

## The correction that did not correct

`fa98723` (2026-10-24), *"CRITICAL FIX: Correct pH direction to acid-OFF (bind pH 7.4, release
pH 6.0)"*, **does not touch `scripts/ph_linkage.py` or `scripts/complex_propka_switch.py`.**
Verified from its file list: it adds `results/c4_acid_off_scored.json`,
`results/submission_v9.csv`, `scripts/build_submission_v9_acid_off.py`,
`scripts/fetch_and_score_all.py` and refold JSONs.

So the Oct 24 fix built a **new parallel scoring path** rather than reverting the Oct 23
inversion. That is the mechanism behind the multiple conflicting conventions in this repository:
each correction added a path instead of repairing the one that was wrong, and the inverted ones
were left in place.

## The second file from the same commit is still inverted

`75c834b` also changed `scripts/complex_propka_switch.py`:

```
+    # C4: positive = proton UPTAKE on binding = correct (bind at pH 5.5)
+    result["ddg_switch_interface"] = round(-total_ddg, 3)
```

Current state, verified today: line 196 carries that comment and line 197 still negates the
value. `git log --follow` shows no commit has touched the file since `75c834b`. It sorts
descending by `ddg_switch_interface` (line 266), so it ranks the acid-ON direction first.

**No committed number in the current C4 lane comes from either file** — the production path is
`c4_linkage_apical_recompute_20261025.py`, which passes 6.0/7.4 explicitly. The hazard is
prospective: anyone running either file fresh gets output in the direction the challenge does not
want, with a comment asserting it is correct.

## Prior records of the same problem

- `results/REVIEW_SUMMARY.md:292` records `ph_linkage.py:128` as `L <= -2`, i.e. the pre-Oct-23
  state. That line is now stale.
- `results/REVIEW_SUMMARY.md:322` states *"FIX THE OBJECTIVE — DONE. switch_sites.py,
  ph_linkage.py, complex_propka_switch.py all corrected."* Measured today, two of those three
  still hold the acid-ON direction. **That claim does not hold.**

## A caller exists, so quarantine-as-specified does not apply cleanly

Opus 5.5 and Fable 5 both recommended moving the file rather than editing it, conditioned on
nothing importing it. `tests/test_install_switch.py:92` does
`from ph_linkage import score_one`, and patches `ph_linkage.run_propka`. Moving the module
breaks that test. This is recorded rather than worked around; the containment decision returns to
the reviewers with this fact attached.

## What is NOT claimed

- **No number currently in the C4 lane is contaminated.** This is provenance and prospective
  hazard only.
- **This is not a judgement about who was right on direction.** The Oct 23 commit was made when
  the objective was understood differently; the point is that the later correction did not reach
  these two files.
- **No file was moved, edited or deleted.** Containment is a decision for the four reviewers, and
  the caller above changes the facts they answered on.
