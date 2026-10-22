# Independent TM-align check exposed an approximate score-export default

The remaining provisional local4 was checked against all26 reference matches
returned by the existing search, plus the known native-domain control and a
self-alignment. No new search, inference, relaxation or design was used here.
All query scores are normalized by the complete query length,95 residues for
the retired binder and195 for the native positive control.

| Measurement for the retired binder | Best query-normalized TM |
|---|---:|
|Earlier Foldseek exported score (approximate default)|0.47600|
|Standalone TM-align, default free alignment|0.42860|
|Standalone-I, fixed saved correspondence|0.44770|
|Standalone-i, saved correspondence as starting seed|0.42860|

All26 matches stayed below0.5 in both standalone modes. The initial maximum
score difference exceeded the predeclared0.01 criterion, so the software
repeatability check was implementation-sensitive, not vindicated. The self
control scored1.00000; the native control scored0.99284.

## Concrete scoring defect and correction

The final Foldseek convertalis step recomputes qtmscore/ttmscore. Its
--exact-tmscore default is0, approximate. A preceding tmalign stage with fast
mode disabled does not override that export default. The prior predictor
therefore exported approximate scores while describing them as full TM-align
scores. This was a defect in our pipeline.

The predictor now explicitly exports --exact-tmscore1. Primary evidence is
preserved in SCORE_EXPORT_SOURCE.json, foldseek_convertalis_help.txt and the exact-score
stage logs. At pinned Foldseek commit463739e0014a1549a527de589102cde98f802f37,
structureconvertalis.cpp constructs TMaligner using par.exactTMscore and computes
the requested qTM at export; TMaligner.cpp dispatches between approximate and
exact score calculations. Source is cached outside the repository.

The controlled exact export changes no candidate/alignment records. It reduces
the maximum absolute difference against standalone fixed correspondence from
0.05369 to0.00225 across26 matches.
The predeclared exact-score agreement tolerance was0.001: **it still
fails**. No tolerance was relaxed. The residual cause is unestablished; no claim
of numerical equivalence between implementations is made. None of those26
fixed-correspondence comparisons crosses0.5. The full-pool correction and new
local level counts are documented in ../exact_score/RESULTS.md.

## Reference identity and reproducibility

The first batch convert2pdb export produced filename collisions, including a
258-CA file where129 were expected. The guard stopped the run before scoring.
References were then exported one numeric DB key at a time; all27 actual header
identifiers and CA counts matched the saved search records. Full query/reference
sequences were also verified against every saved aligned substring before
padding only unaligned terminal residues with gaps. No author/label chain-letter
mapping was guessed and no coordinates were repaired or optimized. See
LOOKUP_CHECK.json, REFERENCE_MANIFEST.json and rejected export logs.

Standalone source: https://github.com/pylelab/USalign
Pinned commit: 23f14ffcf8647f2dc1c0e01043db3a7f3a661c28
Reported program version: TM-align20240303. Compiler, command, source/binary
SHA256 and local installation path are in INSTALLATION.json. The source was
compiled locally with-O3, without-fast or user-defined normalization/scaling.
Binary input fingerprints, all27 free outputs,54 fixed/seeded outputs, input
alignment FASTAs and commands are preserved. Run scripts are included. This is
a comparison of implementations sharing algorithm ancestry, not an independent
biological experiment or an organizer validation.

Exact organizer levels remain unestablished. The historical reference snapshot,
finite discovery, missing patent corpus, supplied coordinates instead of organizer
monomer predictions, whole-chain domain assumption, unknown organizer settings,
and antibody typing gaps persist. The retired pool remains retired.

All work was local CPU, at most4 parallel pair checks,30seconds per pair. No
GPU/cloud CPU allocation, generator/exporter edit, account or submission occurred.
New GPU spend$0.00; recorded Runpod$56.589988 +our Modal estimate$10.718947 =~$67.31,
excluding the separate RFD4 lane, storage/builds and consultations.
