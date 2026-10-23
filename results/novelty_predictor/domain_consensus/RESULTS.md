# CPU domain consensus and broader retrieval: current provisional novelty evidence

The current report is **parent_context_run/predictions.json** (and its CSV).
For the 73 retired designs it predicts **20 local level 2, 40 level 3, 7 level 4,
6 coordinate-unavailable**. **Official novelty levels established: 0.** No design
is promoted or qualified for submission. Five of the seven local4 calls arise
from limited domain coverage despite finding matching domains; they are not
demonstrated new folds.

## What ran

All three domain tools now run locally on this Mac, with four numerical threads,
CPU explicitly selected, zero cloud/GPU/MPS use. The complete census analyzed
68 unchanged coordinate models: 67 retired designs and the native TfR1 positive
control. Every full sequence, backbone atom set and original coordinate source
was verified. No missing atoms, folds or sequences were repaired.

| CPU stage | Measured execution seconds |
|---|---:|
| Chainsaw, 68 inputs | 10.939 |
| Merizo, 68 inputs | 15.991 |
| UniDoc, 68 inputs | 3.838 |
| Consensus calculation | 1.396 |
| Initial domain search → full TM-align → exact export | 187.743 |
| Additional parent-context TM-align + exact export | 58.211 |

These are measured stage times, excluding acquisition/install/debugging time.
The original TED source is pinned to f07e24ef4a19a722f2c3c3001cbd03d477ba3f31;
all 208 published file objects verify unchanged. Merizo's three official weight
shards are pinned and checked against their Git objects. Chainsaw uses its shipped
weights. UniDoc/STRIDE were compiled natively from published sources. The newer
UniDoc 20250514 source archive is explicitly identified; the historical URL
returned 404. Runtime dependencies, native build and source receipts are retained.

Accepted high/medium consensus means agreement from three/two methods, using
the published 70% overlap and 25-residue-domain/5-residue-fragment postfilters.
Merizo's non-domain filtering feeds UniDoc; those methods are not fully independent.
There are **70 accepted domains**: 66 inputs have one, two have two; four inputs
have a discontinuous domain. No input is entirely unassigned. Low one-method
domains are excluded and retained in the raw census for inspection.

The novelty predictor now represents a discontinuous domain as its actual ordered
residue segments. It preserves gaps and source CA coordinates, normalizes TM by
the number of domain residues, and counts each covered protein residue once.
All 70 actual Foldseek database sequences match the declared residue sets;
all 70 query-coordinate sets preserve the source at PDB decimal precision.
The unchanged historical 1,562,678-chain structural DB, two discovery modes,
50,000 candidate cap/mode, full TM-align and exact TM export are used.

## Why the additional candidate pass was needed

The first domain-only search returned no hits for eight short domains. A local
point estimate based on those no-hits is weak retrieval evidence, not measured
TM zero. A prespecified follow-up applied to **every domain**, regardless of
its score: all returned reference targets from that design's frozen whole-chain
search were additionally aligned against the domain with actual full TM-align.

All **26,998 requested query/target pairs** were scored. They span 15,633 distinct
parent-reference targets. Original domain alignments were retained, with repeated
targets deduplicated by best query-normalized TM. The no-hit domain count fell
**8 → 0**; five local calls decreased (two 4→3 and three 3→2). Original best
scores were never lost. No score/coverage threshold, coordinate, domain, reference
snapshot or sequence evidence was changed in this follow-up.

| Retired local level | Whole-chain assumption, with public patents | Initial domain-only search | Current domains + parent-context candidates |
|---|---:|---:|---:|
| 2 | 16 | 17 | 20 |
| 3 | 47 | 41 | 40 |
| 4 | 4 | 9 | 7 |
| Unavailable coordinates | 6 | 6 | 6 |

Ten current calls differ from the whole-chain baseline. The same fixed domain
policy is used for every input; no more-favorable policy was selected per design.
The native control stays local1, with domain2–194/195 (98.974% coverage), best
query-domain TM **0.9927**. The separate patent-source control remains sequence-only
and has no numeric local level.

Of the seven local4 designs, five have accepted-domain coverage only **46.03–59.38%**
while their returned domain TMs are **0.5842–0.9514**. They stay below the published
whole-query structural-coverage threshold. That is a coverage-policy result, not
an absence of known structural elements. The other two have **94.74%** domain
coverage and best domain TMs **0.4807 / 0.4602** in this finite candidate pool.
Those are also provisional results, not proof that no closer structure exists.

## Software controls and remaining gaps

The original pinned TED consensus code shifts synthetic one-based ranges by one
residue. An explicit isolated one-line offset correction passes six consensus
fixtures; the unmodified code fails them. Five checks exercise the actual published
filters. Our initial bare-TER formatting error was also corrected, preserving all
68 ATOM blocks byte-for-byte. INDEXING_CORRECTION.md and failed-run receipts retain
the diagnoses. No biological or novelty threshold was loosened.

All **21 focused predictor tests pass**, including discontinuous coverage,
overlap/order rejection, query extraction, TM denominator and exact-export
behavior. Reclassification of all 68 frozen whole-chain cases with the new code
preserves every legacy best hit, class and local level. All named boundaries
remain unchanged. PARENT_CONTEXT_VALIDATION.json independently checks current
raw hits, controls, invariants and transitions; the comparison CSV contains every
coordinate-bearing input, not just changed examples.

Missing evidence persists: organizer ESMFold/ESMFold2-Fast structures and exact
domain/search versions, its patent corpus/snapshots, antibody typing/CDR policy,
exhaustive structural retrieval and actual portal novelty feedback. The previously
measured standalone TM-align disagreement (maximum0.00225 versus preregistered
0.001 tolerance) remains unresolved; this work does not establish a global error
bound. Public patent coverage is now available, explicitly partial versus the
unknown organizer corpus. No significant sequence hit is not measured zero identity.

A read-only local check of the other worktree found only the explicitly labeled
20-row do-not-submit format test, all matching existing retired analysis sequences.
No new sequence/coordinate design package was found in that inspected scope;
cloud generation state was not queried. Generation and submission remain in
their agreed lanes. The 6W3H gate-control result remains GATE VINDICATED.

New GPU/cloud CPU spend **$0.00**. Recorded Runpod **$56.589988** posted + our Modal
**$10.718947** estimated = **$67.308935 (~$67.31)**, excluding the separate RFD4
lane, storage/builds and consultations. No account, portal/submission, allocation,
generator, exporter/ranker or main-branch action occurred.
