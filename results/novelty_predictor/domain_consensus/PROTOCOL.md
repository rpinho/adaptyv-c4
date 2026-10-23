# CPU feasibility and explicit local TED-consensus comparison

This extends the CPU novelty lane. No generator, GPU/MPS, cloud CPU, account,
submission, exporter/ranker or other lane's file is used or changed.

The public novelty description borrows three-predictor consensus from TED but
does not publish exact program versions or settings. The local implementation
below is therefore an explicitly versioned approximation, never an organizer
score. Do not claim that identifying TED's programs identifies the organizer's
exact implementation.

Use original TED consensus source pinned at
psipred/ted-tools f07e24ef4a19a722f2c3c3001cbd03d477ba3f31. Its methods are
Merizo, Chainsaw and UniDoc. Merizo removes non-domain residues before UniDoc;
do not treat those two as entirely independent. Pin Merizo's three official
weight shards to psipred/Merizo 41d12fb84e6e8fdb586c2c859d12161dc7bb5bfd,
verify their Git blob identities, and record SHA256. Chainsaw weights ship in
the pinned TED tree. The historical UniDoc URL returned 404; use the explicit
20250514 publisher source archive identified by the maintained UCL pipeline,
and record this version difference. Compile native CPU UniDoc/STRIDE from
published sources rather than running bundled Linux binaries on this Mac.
Do not execute TED's setup.sh, which also installs unrelated notebook packages.

Before domain-based novelty claims, establish working CPU inference and validate
residue indexing/filters with source-level fixtures. Require full sequence and
coordinate correspondence to frozen inputs. Models receive complete available
backbone atoms, not the CA-only search files. Never synthesize missing backbone
atoms, repair folds, select a favorable domain prediction, or alter novelty
thresholds. Keep all raw predictor outputs, failures and unsupported intervals.

If preflight succeeds, analyze every available complex-coordinate input from the
frozen exact-score report (67 retired + native positive control), in the frozen
input order. Six unavailable retired inputs stay unavailable. Every method runs
on CPU explicitly, sequentially with at most four numerical threads and finite
timeouts. Report all outputs, not only agreeing predictions. Preserve the
published TED 70% domain-overlap rule and its 25-residue domain/5-residue fragment
filters. Discontinuous domains require faithful residue-set representation;
do not silently split or flatten them to fit the current interval schema.

Only a complete, checked three-method run may be described as local TED-style
consensus. Partial methods or failed indexing preflight remain a documented gap;
do not turn them into a guessed official level. The whole-chain novelty results
remain frozen as a separate comparison. Any later structural search using new
domains must preserve the same reference DB, exact scoring and level boundaries.

Recorded spend remains Runpod $56.589988 posted + our Modal $10.718947 estimated
= approximately $67.31. New GPU/cloud CPU cost: $0.00. All work is on this Mac.

Primary sources:
- https://www.adaptyvbio.com/blog/novelty
- https://github.com/psipred/ted-tools
- https://github.com/psipred/Merizo
- https://github.com/UCLOrengoGroup/domain-annotation-pipeline
- https://yanglab.qd.sdu.edu.cn/UniDoc/download/UniDoc_20250514.tgz
