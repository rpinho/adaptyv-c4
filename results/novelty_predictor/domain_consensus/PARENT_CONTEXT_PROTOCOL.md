# Whole-chain context supplies additional domain-alignment candidates

The first complete domain search returned no structural hits for several short
domains. It nevertheless emitted provisional local levels from finite discovery.
This is a measured retrieval gap, not evidence of TM-score zero or a novel fold.
The initial domain results remain frozen as novelty_run/ and VALIDATION.json.

Before inspecting the follow-up scores, apply one fixed augmentation to **every
one of the 70 domains**, independent of current level: take all returned reference
targets from the same design's frozen whole-chain search, and run actual full
TM-align on those references against each domain. Retain all original domain
alignments and combine by maximum query-domain-normalized TM per query/target.
Do not use only the best/top10 parent hits, select favorable candidates, or change
coordinates, domains, reference snapshot, scoring flags, class boundaries or
coverage thresholds. This can increase detected similarity; it cannot constitute
new evidence of novelty when no match was found.

Verify identical full protein sequences and coordinate sources between parent
and domain reports, complete raw-hit fingerprints, reference DB fingerprints,
and numeric query/reference key mappings. Use a real Foldseek prefilter DB and
the same full TM-align/exact export. Never reimplement the aligner. No extra
sequence or structure model inference is needed; existing domain query DB is
reused. Four CPU threads, finite per-stage timeouts, no GPU/cloud resources.

Report all per-domain score changes, no-hit counts, local level transitions and
remaining coverage gaps. Require original best domain scores to be retained and
native control to stay local1. A known matching domain covering too little of
the whole query may still produce a local4 under the published coverage rule;
call out that mechanism rather than describing it as a demonstrated new fold.

This remains a finite candidate set, not exhaustive alignment to all PDB chains
and not the organizer's unpublished discovery policy. All 73 designs remain
retired. Official levels remain unestablished. No automatic second augmentation.
Recorded Runpod $56.589988 + our Modal $10.718947 =~$67.31; new GPU/cloud CPU $0.
