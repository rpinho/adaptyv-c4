# Provenance join by sequence hash: all 343 C4 designs were built against `E1_his707` chain A

Written 2026-10-26T14:22:46Z (real clock, `date -u`). CPU only. **$0 GPU.** No submission, no
organiser contact, no frozen file touched. This is a *confirmation by hash* of a conclusion the
earlier finding had reached by config-name matching.

## Why this file exists

`results/c4_native_context_20261026/FINDING_WRONG_DOMAIN.md` concluded the C4 designs target the
wrong domain. Its stated basis was indirect: an exact config-name match, binder-length agreement
(60-100 configured vs 65-100 observed), and the His707 target files. The prior session flagged this
gap itself, verbatim: it identified the generator *"not from a per-design provenance join by
sequence hash to generator output."*

That join had not been run. It is now.

## Method

For every `results/c4_rebuild_20261024/remote_pdbs/*_complex.pdb`:

1. Parse all `ATOM` records, group by chain, order by residue number.
2. Select the longer chain as receptor and the shorter as binder (verified: receptor is 154 in all
   343 files, binder 60-100).
3. Compare the receptor sequence **byte-exactly** to chain A of `targets/E1_his707.pdb`.
4. Join key is `sha256(sequence)`. Design names are recorded but never used as the join key.

## Result

| quantity | value |
|---|---:|
| design complexes inspected | 343 |
| receptor chain == `E1_his707` chain A, exact | **343 / 343** |
| mismatches | **0** |
| unparsable | 0 |
| distinct binder sequences (by sha256) | 343 |
| binder length range | 60-100 |

Target identity: `targets/E1_his707.pdb` chain A, **154 residues, 260-750**, sha256
`1653daa4ab8cc1b36e111b7d91765c6e24c6bbede73cb163eef2dcf33775130e`.

So the claim is now a measurement, not an inference: **every design in the C4 pool was docked to
exactly one 154-residue construct**, and that construct is identical in all 343 files.

## What that construct actually is

It is a **stitched, multi-fragment crop**, not a clean domain. Chain A retains native numbering but
carries **6 nonnative junctions**:

| junction | residues skipped |
|---|---:|
| 268 -> 300 | 31 |
| 315 -> 322 | 6 |
| 330 -> 532 | **201** |
| 551 -> 616 | 64 |
| 630 -> 633 | 2 |
| 644 -> 678 | 33 |

Residue-level assignment against the apical domain (189-383):

| band | residues in construct | share |
|---|---:|---:|
| 260-383 (apical) | 34 | **22.1%** |
| 384-620 (cysteine-rich / helical) | 25 | 16.2% |
| 621-750 (helical, C-terminal) | 95 | 61.7% |
| **outside apical 189-383** | **120** | **77.9%** |

**78% of the modelled target is not the apical domain the campaign intended to hit.** The single
largest continuity block is 621-750, which is the helical domain.

## Two things this does *not* claim

- **It does not show the designs are worthless.** A binder raised against a real surface of the real
  protein may still bind. What it shows is that the designs were never selected against the domain
  the campaign named, so no claim about apical binding, TfR1 competition, or BBB transport follows
  from this pool. The measured abilities were measured on a different surface.
- **It does not by itself make any file submittable or unsubmittable.** That judgement rests on the
  combination of this, the inverted pH gate (`FINDING_L6_SIGN.md`, rule v4), and the native-context
  gap, not on this receipt alone.

A separate, unresolved point worth naming: 120 of 154 residues lie outside apical 189-383, but the
construct's own provenance note in `TARGET_PREP.md` describes a deliberate spatially-contiguous
crop. Whether that crop was intended as the design surface, or was a target-prep artifact that
silently became the design surface, is **not settled by this join** and I am not settling it here.
The join establishes *what was docked to*. It does not establish *what was intended*.

## Artifacts

- `provenance_join.json`, sha256 `74e2b8a72556d84b8747839c44dd9817c3c6c60dae2a13f202182c7d976c4e36`
  — 343 per-design records: design file, binder sha256, binder length, receptor sha256, match flag.

## Controls

- Positive: the receptor of a known `E1_his707`-derived complex matches exactly -> fires (343/343).
- Negative: comparing against `E2a_apical_202_213` and `E2b_apical_340_376` chains produced no
  exact match, so the equality is not trivially satisfied by any TfR1 crop.
- Distinctness: 343 distinct binder sha256 values, so the match is not an artifact of one file
  counted 343 times.
