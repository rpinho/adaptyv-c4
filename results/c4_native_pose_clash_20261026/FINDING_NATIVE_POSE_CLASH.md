# FINDING: the native-pose clash test is INCONCLUSIVE — and the "0 of 343" it produced is a frame artifact, not a catastrophe

Written 2026-10-26 from the real clock. CPU only, **$0**, read-only inputs. No GPU, no submission,
no organiser contact. No frozen file touched (4/4 verify). All five controls fired.

Mandated by **Astra**, PART 3 of the 2026-10-26 reviewer round, and run on the **full 343** with no
shortlist fixed, as Opus 5.5, Fable 5.1 and Astra unanimously specified (process Q2 = A).

> Astra: "Rigidly superpose each shortlisted complex onto the correctly numbered native homodimer
> and measure nonbonded binder–receptor heavy-atom overlaps... **Flip B to A if all eight have
> extensive, unavoidable backbone interpenetration after alignment verification.**"

## Headline: the flip condition is NOT met

The raw output looks like it fires the condition — **0 of 343 designs clash-free**, median 382
heavy-atom overlaps below 2.0 Å, median 111 backbone–backbone. Taken at face value that would
flip the reviewers' unanimous "submit 8" to "submit nothing".

**It should not be taken at face value.** The repo's own rule — *a zero needs a control that
fires* — was applied, and the control shows the zero is produced by the reference frame, not by
the designs.

## The confound, measured

| binder clashes < 2.0 Å, by receptor region | median over 343 |
|---|---:|
| chain A, **the band they actually bind (651–741)** | **6** |
| chain A, outside the band | **311** |
| partner subunit B | 16 |
| glycans | 0 |

And the decisive control:

| | median |
|---|---:|
| binder heavy atoms clashing with chain A outside the band | **311** |
| **the crop's OWN non-band atoms clashing with the same region** | **316** |

The crop's own atoms interpenetrate the native receptor to the same degree as the binder's. The
binder is docked relative to the whole crop, so it inherits the crop's misplacement. The 311 is
therefore a property of **Boltz-2's distorted crop geometry**, not of the binder's pose.

Astra's condition requires interpenetration **after alignment verification**. Alignment
verification *fails* outside the band. So the condition is not met, and this check is
**inconclusive** for every region except the band.

Against the band itself, where alignment is verified at 1.21 Å: median **6** clashing heavy atoms,
p90 **46**, and **50 of 343 designs are entirely clash-free**. That is not a demonstration of
steric feasibility either — but it is the opposite of unavoidable interpenetration.

## The by-product worth keeping: defect #2, quantified

The first run refused to produce any result because its global superposition control failed. That
refusal was correct and the number behind it is a finding in its own right:

| fit region | residues | median CA RMSD |
|---|---:|---:|
| whole 154-residue crop | 154 | **17.46 Å** |
| apical fragment 260–330 | 34 | 8.62 Å |
| helical 532–750 | 120 | 15.19 Å |
| **real design band 651–741** | 64 | **1.21 Å** |

Boltz-2 reproduced the band the designs bind at near-experimental accuracy and **misplaced the
stitched fragments relative to one another by ~17 Å**. This is the first quantification of the
"six nonnative junctions" defect that `HANDOFF_TO_CLAUDE_20261026.md` §4 defect 2 described
qualitatively. Consequence: any measurement depending on the crop's **global** geometry is
unreliable, while band-local measurements are sound.

## Controls — 5/5 fired

| control | requirement | measured | fired |
|---|---|---|---|
| MAPPING_LEN | crop order is 154 residues | 154, author 260–750 | yes |
| **HIS707_REGRESSION** | mapping reproduces the independently recorded His707 contacts | **19/343, exactly as recorded** | yes |
| SUPERPOSE_RMSD | band-local fit is tight | median 1.21 Å (global 17.46 Å reported alongside) | yes |
| POS_DISPLACED | a binder translated +500 Å has zero contacts | 0 | yes |
| POS_NATIVE_PARTNER | the real partner subunit registers contacts against chain A | 220 atoms | yes |

`HIS707_REGRESSION` is the load-bearing one: it reproduces a count measured by a different script
in a different cycle, confirming the file-order → author-numbering mapping independently.

## What changed in method, and why

The fit frame was moved from the whole crop to author 651–741 **after** the global fit failed its
control, not to obtain a nicer answer. The reason is measured (17.46 Å vs 1.21 Å), the global RMSD
is still computed and reported per design, and the first run's refusal is preserved in this
document. Had the frame been changed silently, the 0/343 would have been reported as a result.

## What is NOT claimed

- **Not a binding, affinity or pH measurement.** Steric overlap on rigid coordinates only.
- **Band-local non-clashing is not feasibility.** 50 designs clash-free against the band says their
  pose does not interpenetrate *that* surface; the rest of the receptor remains untested because
  the frame cannot test it.
- **The receptor is treated as rigid** and modelled from one crystal form; real receptors flex.
- **This does not rescue any design.** It removes a reason to reject them, and leaves the
  reviewers' position where it stood.
- **A correct version of this test needs a correctly built target**, i.e. the designs re-folded
  against a properly gapped crop or the full native assembly. That requires GPU and is a decision
  that sits in `results/decision_requests/20261026-c4-construct-scope.md`.
