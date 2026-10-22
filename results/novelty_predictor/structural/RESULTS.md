# CPU structural novelty screen completed; levels are provisional

Superseded scoring notice: this report used approximate TM export. Current exact
score results are in ../exact_score/RESULTS.md.
Corrected retired local counts:16 level2,47 level3,4 level4,6 unestablished.
Prior counts and score values remain below as historical evidence.

Current CPU follow-up supersedes the earlier search counts: see ../sensitivity/RESULTS.md.
The predictor now combines both discovery modes; current retired local counts are
29 level2,37 level3,1 level4,6 unestablished. All levels remain provisional.

The existing sequence predictor now runs **Foldseek candidate discovery followed
by real TM-align**, in the same `scripts/novelty_level_predictor.py`. The finite
database, supplied coordinates and declared domains are explicit in every result.
No organizer score is claimed. `novelty_audit.py` is unchanged.

## Measured results

The existing 73 retired sequences plus the native 195-aa control were retained.
Their earlier four-database sequence evidence was reused only after verifying
the input CSV hash, names, order and per-sequence hashes. All 73 retain their
previous no-significant-sequence-hit result; this is not measured zero identity.

| Retired designs | Structural evidence | Local provisional level |
|---:|---|---:|
| 26 | High similarity | 2 |
| 37 | Moderate similarity | 3 |
| 3 | Returned matches, all below moderate | 4 |
| 1 | No structural candidate match detected | 4, especially weak evidence |
| 6 | Exact coordinate models unavailable locally | Unestablished |

Thus **67/73** retired designs received local predictions: **26 level 2, 37 level 3,
4 level 4**. Sequence class is provisionally treated as low when no significant
sequence hit was detected. Every numeric level has `prediction_confidence=low`,
`level_status=provisional_local_prediction`, `organizer_level_established=false`,
the missing sequence database list and explicit possible levels. The JSON and
CSV both carry the caveat; these numbers must not be read as portal scores.

The native-domain positive control has **TM-score 0.9946**, normalized by its full
195-aa query domain, against **7zqs-assembly1_B**. Its aligned domain sequence is
identical, structural high coverage is **1.0**, and its local level is **1**, as
expected with the previously measured 100% sequence identity. The control is not
a candidate. The target-normalized score is only 0.3045 because that target chain
has 639 residues: using it would produce the wrong classification.

Across all 68 coordinate-bearing queries there were **12,715** TM-align rows;
per-query counts ranged from **0 to 1,383**. Among retired designs with returned
matches, best query-normalized scores ranged **0.4668–0.9699**. The three
below-moderate maxima were:

| Retired design suffix | Best query-domain TM | Best returned target |
|---|---:|---|
| `l96_9003044e98f7dfc2_candidate6` | 0.4668 | 7np7-assembly1_D6 |
| `l96_9003044e98f7dfc2_candidate10` | 0.4945 | 8fit-assembly2_C |
| `l95_c6a5afd6447f1fd0_candidate7` | 0.4760 | 8gzu-assembly1_69 |

`l95_c6a5afd6447f1fd0_candidate4` returned no candidate match. It has no measured
maximum TM-score; its score is null, not numerical zero. Nothing here proves a
new fold. The maximum nearest either threshold was 0.4945 near 0.5 and 0.8031 near
0.8; no best-hit classification depends on a printed score exactly at a threshold.

## Database actually obtained

Installed the official macOS universal CPU Foldseek binary, version/commit
**463739e0014a1549a527de589102cde98f802f37**. Installation/archive fingerprints are
in `INSTALLATION.json`. No GPU, Modal, Runpod or cloud CPU was used.

Downloaded the public **Foldseek PDB100** archive on **2026-10-22**, using
`foldseek databases PDB`. We searched **all 1,562,678 member chains in `pdb_seq`**,
not just its 324,204 cluster representatives. The searched sequence/3Di/CA/header/
index files total **3,725,663,586 bytes**; archive size **2,326,827,389 bytes**.
The archive's MD5 matches the publisher's version file, and all search inputs
have SHA256 fingerprints in `DATABASE.json`.

Dates must be distinguished: publisher version file **`PDB_DATE=250101`
(2025-01-01)**; database file build timestamps **2025-07-25**; HTTP archive
Last-Modified **2025-12-19**. The download date does not make this a 2026 PDB
snapshot. This historical database cannot establish novelty against later folds.

Source: https://foldseek.steineggerlab.workers.dev/pdb100.tar.gz

Version: https://foldseek.steineggerlab.workers.dev/pdb100.version

Documentation: https://github.com/steineggerlab/foldseek

## What ran and how coverage was calculated

The pipeline verifies each coordinate fingerprint, extracts the explicitly named
chain, and requires its complete CA sequence to equal the query sequence. Only
that chain's coordinates are written to Foldseek query files, excluding receptor
chains. Minimal BindCraft CIFs are parsed without requiring a nonexistent
occupancy column. A first attempt with a wrong chain manifest was rejected by
the sequence check before search; the corrected manifests explicitly select C
for the retired binders and native author A189–383 for the control.

The 67 short binders are each treated as **one whole-chain domain by a declared
assumption**, not by an implemented domain predictor. Each control/query domain
is its own Foldseek entry. Intervals are validated as disjoint, preventing double
counting when a protein has several domains or several hits for the same domain.
Missing-coordinate queries do not receive a fabricated model or structural class.

Foldseek's 3Di+AA search uses sensitivity 9.5, permissive E-value ≤10, no coverage
or identity filter and at most 10,000 prefilter candidates per query. Every
returned candidate is then realigned with **Foldseek's built-in TM-align**,
`--tmalign-fast 0`, no TM cutoff, full query-domain normalization and query-TM hit
ordering. This is a distinct TM-align stage, not simply the local alignment's
TM-like rescoring. The raw output includes both normalizations, aligned sequences,
query/target lengths and aligned positions. Commands and four logs are saved.
The built-in implementation may differ from an organizer standalone TM-align
version; neither the organizer version nor search settings were published.

For each domain, its highest query-normalized TM determines whether the domain
qualifies at 0.5 and 0.8. The lengths of disjoint qualifying domains are divided
by the full query sequence length. High requires ≥70% at TM≥0.8; moderate requires
>70% at TM≥0.5. This is domain coverage, not a weighted mean TM-score and not a sum
over redundant hits. A whole-chain domain meeting either threshold consequently
contributes coverage 1.0 under this local domain policy; it does not mean every
residue has a close atomic match. The published nine-cell sequence/structure
truth table then produces the local level. Equality is tested explicitly.

All stages used **4 CPU threads**, CPU search `--gpu 0`, and a 4G search split
memory limit. Coordinates already existed; no new inference or design was run.

## Published definition that remains unreproduced

The [published ProteinTyper description](https://www.adaptyvbio.com/blog/novelty)
predicts structures with ESMFold/ESMFold2-Fast and uses a consensus of three domain
predictors, followed by Foldseek and TM-align. This screen does not reproduce:

- Those sequence-derived structures: existing complex-predicted binder coordinates
  are used. A different monomer prediction can change the structural result.
- The unnamed three-predictor consensus: whole-chain minibinder domains are an
  explicit approximation. Arbitrary verified disjoint domain intervals are supported.
- The organizer structural database/provider/snapshot, which the article does not
  identify. Its exact corpus is unobtainable from the published specification;
  the public historical PDB database is a documented replacement, not a replica.
- The unidentified patent sequence collection already missing from the sequence
  half, and unpublished sequence/structural search parameters.
- Antibody detection and the separate antibody/CDR novelty scale.
- Structural evidence for six missing exact coordinate models (names/reasons are
  in `STRUCTURES.json`). No CPU structure inference was substituted.

**Exact organizer levels established: 0.** For the retired sequences, missing
sequence data and organizer structure/domain policy differences still permit
levels 1–4 in the broad uncertainty field. `possible_levels_if_structure_policy_accepted`
shows narrower bounds conditional on accepting this screen's structural result.
A finite prefilter/search and historical database can miss similar structures;
in particular the four local level-4 calls are provisional, not evidence of a
guaranteed novelty gate pass. Existing retired/admission decisions are unchanged.

## Reproduce

Use a fresh output directory and the local downloaded database:

```sh
.venv/bin/python scripts/novelty_level_predictor.py \
  results/novelty_predictor/analysis_inputs.csv \
  --sequence-report results/novelty_predictor/retired_pool_run/predictions.json \
  --structures results/novelty_predictor/structural/STRUCTURES.json \
  --structure-database results/novelty_predictor/structural/DATABASE.json \
  --foldseek /Users/ricardo/.local/opt/foldseek/foldseek/bin/foldseek \
  --out results/my_fresh_cpu_structural_screen
```

Alternatively provide `--databases` instead of `--sequence-report` to run both
halves on new sequences with corresponding coordinate/domain manifests. No tool
automatically downloads, guesses a chain, or infers a missing structure.

Fourteen unit tests pass: earlier sequence/threshold tests plus multi-domain
coverage, query-versus-target normalization, incorrect query-length rejection and
real native-CIF extraction. The real search/control and independent report
invariants are checked in `VALIDATION.json`. Raw alignments, query domains,
predictions, hashes and commands are saved in `run/`.
