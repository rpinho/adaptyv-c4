# Conditional antibody novelty evidence

CPU ANARCII2.0.8 annotation of all returned-hit-associated heavy references completed. This is a supplementary antibody-scale calculation; the organizer structural typer and VH/VHH distinction remain unmeasured.

- Inputs:50 complete-flanked H annotations among400 saved designs. The other350 have no CDRH3 result here.
- Restored109,342 unique heavy reference sequences from all raw paired/unpaired source associations;107,440 canonical references annotated,1,902 noncanonical sequences explicitly excluded.105,529 references have usable complete H loops.
- Unit Levenshtein identities were calculated with max-length, query-length and reference-length normalization; normalized indel similarity was tested separately. Boundaries70% and95% are inclusive named constants.22/50 measured-subset levels change with normalization.
- Primary max-length arithmetic over measured reference subsets:11 level2,34 level3,5 level4. These are not predictions under missing annotations.
- Allowing missing-reference completions and all four normalizers:49/50 scalar calls withheld; one conditional local4. That one has only two returned heavy references (one distinct usable loop). Finite full-chain retrieval is a major limitation; this is not a novelty-gate pass.
- All200 exported best-reference identities independently recomputed and source associations verified. Publisher self-loop, unequal-length and boundary controls passed.

Reference annotation took1,904.43s; conditional scoring11.30s, on local CPU. Raw per-reference annotations and associations remain in the cache SQLite file named by POOL.json; raw model output hashes and an immutable annotation checkpoint are retained.

Not established: organizer antibody class, organizer edit-distance normalization, complete PLAbDab CDR retrieval, organizer reference snapshots/global-identity policy or patent corpus. No binding, pH behavior, assay success, entry ranking or qualification is measured. See PREDICTIONS.json for per-design evidence and gaps, EVIDENCE.csv for compact rows, and VALIDATION.json for checks.

Source rule: [Adaptyv novelty](https://www.adaptyvbio.com/blog/novelty). Same public PLAbDab source snapshot as the sequence run; no replacement database was used.

Recorded own-lane spend: Runpod$56.589988 + Modal estimate$10.718947 =~$67.31. New GPU/cloud compute$0.00; separate generation/account-wide billing excluded.
