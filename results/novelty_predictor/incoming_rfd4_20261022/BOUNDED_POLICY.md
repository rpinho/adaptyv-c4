# CPU continuation policy, before bounded alignment results

The five sequence searches and both full-reference structural discovery modes
completed. Their merged whole-chain candidate DB contains16,451,768 alignment
records. The original exhaustive TM-align stage reached its1200-second CPU
execution limit without completing. Its failed run is retained in run/; no
complete structural result or numeric novelty level is claimed from it.

Completed work is reused only after exact input/query, reference, binary/version,
flags and native DB fingerprint checks. CACHE_CHECKPOINT.json and
structure/DISCOVERY_CHECKPOINT.json cover completed discovery only. Partial
TM-align output is never reused as complete. The original predictor source is
preserved in Git at cbb72a1. A sequence checkpoint is now written before a future
structural stage, and large transient hit collections are released before alignment.

For this continuation, all400 inputs remain included. Use the same public sources,
coordinates, reference snapshot, two discovery modes,50000-candidate caps,
full TM-align (fast0), query-domain normalization, exact score export and fixed
classification boundaries. The new positive-evidence mode bounds **alignment
coverage**, explicitly reducing what can be concluded. It is not a completed
exhaustive comparison or a claim of identical best numerical TM-scores.

Two full TM-align passes seek positive matches at the existing0.8 and0.5 query
thresholds, stopping each query after1 accepted or64 rejected comparisons per
pass. Exact exported scores are still evaluated by the existing domain coverage
rules. No hit is unestablished; it is never interpreted as measured zero or a
new fold. Every unsampled comparison can still reveal a stronger known match.

A partial match set gives a lower bound on structural resemblance. Local point
levels are emitted only when positive evidence already establishes at least
moderate resemblance AND every completion of the unperformed comparisons yields
the same general-scale level. For example, low sequence/moderate structure leaves
levels2-or3 unresolved; moderate sequence/moderate structure fixes local2, and
high sequence/moderate structure fixes local1. High structural resemblance fixes
local1-or2 according to the observed sequence class. Otherwise numeric levels
are null and the unfinished-alignment interval is reported. This mode emits no
coverage-only local3/4, even where the fixed segmentation mathematically limits
the coverage. It cannot certify an eligible candidate from missing evidence.

The same policy applies to the whole-chain baseline, accepted consensus domains
and parent-context augmentation. Parent follow-up here adds every **returned
bounded positive** whole-chain reference, not every reference from a completed
exhaustive whole-chain alignment. That is a further explicit completeness gap.

Antibody typing/CDR-policy caveats remain separate (ANTIBODY_SCOPE.md). Exact
organizer levels:0. No functional, assay or submission qualification is made.
New GPU/cloud compute:$0.00; recorded Runpod$56.589988 +own Modal$10.718947
=~$67.31, other lane/account-wide spend excluded.
