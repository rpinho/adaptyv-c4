# Join gap closed: the two construct sets are the same 342 molecules

Written 2026-10-26T15:00Z+ from the real clock (`date -u` at write time). CPU only, **$0 GPU**.
No submission, no organiser contact, no frozen file touched. This closes the priority-1 item
recorded in the 14:56Z checkpoint.

## The gap that was open

Two construct sets existed and both self-verified, but it had **not** been proven they contained
the same molecules:

| set | receptor | n | used for |
|---|---|---:|---|
| `results/c4_rebuild_20261024/remote_pdbs/` | 154 res (`E1_his707` chain A) | 343 | my chemical-contact and hotspot findings |
| `results/c4_apical_pdb_audit_20261025/all343/` | 195 res (P02786 189-383) | 342 | Worker A's ranking and falsification |

The risk was concrete: if the binder sequences differed, the contact findings and the ranking
findings would describe different molecules and could not be cited together.

## The measurement

Every complex parsed, chains ordered by length, the shorter chain taken as binder, and binders
compared by `sha256(sequence)` — never by name.

| quantity | value |
|---|---:|
| rebuild set files | 343 |
| apical audit set files | 342 |
| distinct binder sha256, rebuild | 343 |
| distinct binder sha256, apical audit | 342 |
| **intersection by sha256** | **342** |
| only in rebuild | **1** |
| only in apical audit | 0 |
| shared design names with identical sequence | 342 / 342 |

**Conclusion: 342 of 343 binders are byte-identical across the two sets.** The contact findings and
the ranking findings describe the same 342 molecules. They differ only in which receptor construct
each was docked to — which is precisely the defect already documented in
`results/c4_apical_pdb_audit_20261025/FINDING.md` (the shipped ranking took iPTM/L6 from the
195-residue apical run while its clash screen read the 154-residue files).

## The one extra file, fully reconciled

`tfr1_c4_acid_site1_41` (binder len 100, sha `9c0a6ad3...`) exists only in the 154-receptor set.
It is not a silent drop. Traced end to end:

1. `results/c4_apical_consolidated_results.json` holds **343** records.
2. Exactly one carries an error: `tfr1_c4_acid_site1_41` -> `"error": "timeout"`.
3. 343 records - 1 timeout = **342** scored, which equals `linkage_apical.json`'s `n = 342`.
4. The retry log confirms an attempt was made:
   `Retrying ONLY tfr1_c4_acid_site1_41 (100 res) vs contiguous apical 195-mer; 342 valid designs untouched`.

So `all343` is misnamed: it contains **342** files because one apical refold timed out. The
exclusion is documented in the results data, not invented. **No design was dropped for an
unrecorded reason.**

## Latent landmine found and defused (no result affected)

`modal_boltz2_apical.py:191` sets `cached["skipped"] = True` on the **cache-reuse** path. So
`skipped: True` means *"reused from cache"*, not *"excluded"*. Of the 343 consolidated records,
**341 carry `skipped: True` while holding a valid iPTM**.

Any downstream consumer that treated `skipped` as an exclusion filter would have discarded 341 of
342 valid designs. Verified that this did **not** happen: the ranking's 342 valid names equal
`linkage_apical.json`'s n of 342, and the scored count was reproduced independently by Worker A.
**No reported number is contaminated.** But the field name is inverted from its meaning and should
not be trusted by future scripts without reading the cache path.

## What this does and does not change

- **Changes:** nothing numeric. This is a provenance repair, not a new measurement.
- **Does not change:** the ranking verdict (b), the contact-mass distribution, the His707 finding,
  or the conclusion that no C4 file on disk is submittable.
- **Does now permit:** citing the contact findings and the ranking findings together as describing
  the same 342 molecules.

## Artifacts

`join_gap.json` — per-set counts, receptor-length histograms, sha intersection, the single extra
design, and its reconciliation chain.
