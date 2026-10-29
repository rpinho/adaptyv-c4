# The designs contact BOTH domains — the "no overlap at a single residue" claim is wrong

Written 2026-10-26, real clock from `date -u` at write time. CPU only, **$0 GPU**. No submission,
no organiser contact, no frozen file touched.

## What this corrects

The handoff headline of the previous cycle states:

> The 342 designs were GENERATED against TfR1 author residues ~651-741 (the HELICAL domain).
> Every evaluation from v3 onward scored them against author 189-383 (the APICAL domain).
> **The two do not overlap at a single residue.**

The first two sentences stand. **The third does not.** It is a statement about two residue *ranges*
being disjoint, which is true, but it was used to imply the designs make no apical contact, which
is false and measurable.

## The measurement

The design complexes renumber the receptor 1..154 in file order. Mapping each position back to the
author numbering of `targets/E1_his707.pdb` chain A (`pos i -> order[i-1]`) makes the two ranges
comparable. Position 1 is author 260; position 34 is author 330; position 35 is author **532**.

That mapping is the whole point: **positions 1-34 of the construct ARE apical residues 260-330.**
The 201-residue jump 330 -> 532 sits inside the construct, so a single contact index can land in
either domain. The two domains are disjoint as *ranges* and simultaneously present in the *same
file*.

Across all 343 design complexes at 4.5 A, contacts mapped to author numbering:

| band | total contacts | designs with >=1 |
|---|---:|---:|
| apical 189-383 | 2,915 | **294 / 343** |
| helical 621-750 | **8,568** | **342 / 343** |
| cysteine-rich 384-620 | 2,075 | 281 / 343 |

So 86% of designs make at least one apical contact, and 100% make at least one helical contact.
The helical band carries roughly **3x** the contact mass of the apical band.

## The sharper, defensible claim

"The designs were aimed at the wrong domain" survives, but the honest form is a **distribution**,
not an exclusion:

- Dominant contact mass is helical (8,568 vs 2,915), consistent with generation against the
  651-741 groove.
- But apical contact is not absent; it is present in 294/343 designs, via the construct's own
  260-330 segment.
- **Only 19 of 343 designs contact author HIS 707 at all** — the residue the campaign is named
  after and the anchor of the intended His707 epitope. 272/343 contact at least one of the named
  hotspot residues (684, 685, 699, 701, 702, 707, 709, plus 651/676/721/741).

The last line is the strongest single number in this file. 5.5% of designs touch the residue the
whole campaign is about. That is a far better-grounded objection than "no overlap at a single
residue", because it is per-design, per-residue, and it does not depend on a range-disjointness
argument that the data contradicts.

## Why this matters for the decision

The previous headline would have supported "the pool is entirely irrelevant to the stated target."
The measurement supports something narrower and still serious: **the pool is mostly pointed at a
different surface, and almost none of it touches the intended anchor residue.** Designs that do
bind would mostly bind the helical face. Whether any of them could still compete with transferrin
or cross the BBB is not answered here, and is not answered by any file on disk.

## Caveat, stated plainly

Contact is not binding. This measures proximity in a predicted complex, nothing more. A design
with zero His707 contact is not thereby non-functional, and a design with His707 contact is not
thereby functional. The counts describe where the models put the interface, not what the molecules
do.

## Artifacts

- `hotspot_contact.json` — per-design booleans: contacts His707, contacts any named hotspot,
  n_contacted, and the full sorted list of contacted author residues.
- `provenance_join.json` — the 343/343 receptor identity join this file builds on.
