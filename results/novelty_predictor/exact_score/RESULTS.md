# Exact TM-score export corrected; earlier local counts superseded

The predictor had a score-export defect: Foldseek convertalis recomputes TM-scores and defaults to approximate scoring. Running full TM-align earlier did not change that default. The exporter now explicitly sets --exact-tmscore1. This correction changes16 retired-pool local calls; it does not create or qualify any new design.

Current results are in run/predictions.json and run/predictions.csv. Older structural/sensitivity local counts used approximate export and are superseded. All current levels remain provisional; exact numerical TM scoring does not establish an organizer novelty score.

## What ran and measured results

A cache-only exact export reused the frozen query DB, reference DB and full TM-align result. Every candidate, aligned sequence, correspondence and other output field was identical to the earlier export; only the score fields changed. A separate fresh execution reran the full default two-mode discovery -> merged candidates -> full TM-align -> exact export pipeline. Its per-query/target exact scores exactly matched the cache-only rescore. Neither sequence evidence, coordinates, domain intervals, reference snapshot nor classification thresholds changed.

The fresh full pipeline completed in265.026 seconds,4 CPU threads,4G search-memory limit,300seconds per stage. It evaluated68 coordinate-bearing queries,42718 raw alignment rows and24653 distinct query/target matches. Four sequence databases were already searched; those fingerprint-verified sequence results were reused. All74 inputs are accounted for:73 retired designs plus one native positive control.

| Retired-pool local level | Earlier approximate export | Corrected exact export |
|---|---:|---:|
| 2 | 29 | 16 |
| 3 | 37 | 47 |
| 4 | 1 | 4 |
| Unestablished: coordinates missing | 6 | 6 |

The positive control remains local1, query-normalized TM0.9928 against7zqs-assembly1_B. It has100% sequence identity over its195-aa native query. It is not a candidate.

| Retired design whose local call changed | Approximate level / TM | Exact level / TM |
|---|---|---|
| c4_tfr1_groove_human_forced_l78_cc1057879d35a8ec_candidate2 | 2 / 0.8031 | 3 / 0.7537 |
| c4_tfr1_groove_human_forced_l98_ac2cbaba03c60f1c_candidate6 | 2 / 0.8101 | 3 / 0.7637 |
| c4_tfr1_groove_human_forced_l98_ac2cbaba03c60f1c_candidate8 | 2 / 0.8058 | 3 / 0.7556 |
| c4_tfr1_groove_human_forced_l84_fa5e4e94654ff9c1_candidate1 | 2 / 0.8232 | 3 / 0.7685 |
| c4_tfr1_groove_human_forced_l76_06cf7c724b17b60d_candidate1 | 2 / 0.8366 | 3 / 0.7782 |
| c4_tfr1_groove_human_forced_l76_06cf7c724b17b60d_candidate4 | 2 / 0.8429 | 3 / 0.7931 |
| c4_tfr1_groove_human_forced_l56_b9b1bd89e75cae01_candidate3 | 3 / 0.5560 | 4 / 0.4949 |
| c4_tfr1_groove_human_forced_l60_63126a071c041270_candidate1 | 2 / 0.8234 | 3 / 0.7615 |
| c4_tfr1_groove_human_forced_l60_63126a071c041270_candidate2 | 2 / 0.8334 | 3 / 0.7814 |
| c4_tfr1_groove_human_forced_l60_63126a071c041270_candidate3 | 2 / 0.8405 | 3 / 0.7891 |
| c4_tfr1_groove_human_forced_l95_2aeadb4dd8ac935e_candidate1 | 2 / 0.8395 | 3 / 0.7978 |
| c4_tfr1_groove_human_forced_l70_54a5b53444de2b49_candidate8 | 2 / 0.8352 | 3 / 0.7747 |
| c4_tfr1_groove_human_forced_l58_9f94b53776050aad_candidate4 | 2 / 0.8196 | 3 / 0.7681 |
| c4_tfr1_groove_human_forced_l58_9f94b53776050aad_candidate5 | 2 / 0.8227 | 3 / 0.7724 |
| c4_tfr1_groove_human_forced_l96_9003044e98f7dfc2_candidate10 | 3 / 0.5384 | 4 / 0.4854 |
| c4_tfr1_groove_human_forced_l95_c6a5afd6447f1fd0_candidate5 | 3 / 0.5273 | 4 / 0.4800 |

## Independent implementation check: result and remaining failure

An official standalone TM-align20240303 binary was compiled from pinned pylelab/USalign source. The original remaining-local4 was checked against all26 returned references, with query self and native-domain positive controls. All query/reference lengths, identities, coordinate fingerprints and saved aligned substrings were verified. Each numeric reference key was exported separately after a batch filename export was rejected by a CA-count guard before scoring.

The original approximate best query score was0.4760. Standalone default alignment gave0.4286; holding the saved correspondence fixed gave0.4477. Correct exact export also gives0.4477 for that winning match. Its broad category remains below moderate for this finite match set.

Across all26 fixed-correspondence comparisons, the maximum absolute difference decreased from0.05369 to0.00225. The predeclared agreement tolerance was0.001: **that numeric agreement check still fails**. It was not loosened. The residual cause is unestablished; no global numerical equivalence or error bound is claimed. Both implementations kept all26 matches below0.5. Full evidence and limitations are in ../independent_tmalign/RESULTS.md.

## Regression and data checks

- 16 focused unit tests pass. A behavioral regression exercises the complete structural command path against a test-only tool stub whose approximate/default score and exact score straddle a classification threshold. It checks the resulting level, not just a command string. The same regression deliberately fails with export flag0 and passes with1. These synthetic fixture scores are not biological measurements.
- The actual fresh export command contains --exact-tmscore1 and its report records exact_tmscore_export=1. All fresh exact scores equal the isolated cached scores; every old aligned-record field other than the two TM values is unchanged.
- Input CSV, sequence digests, coordinate/domain manifest, reference DB fingerprints and raw exact output SHA256 match the saved receipts. The class/threshold rule functions and native-context gate remain unchanged.
- Validation code and VALIDATION.json preserve these checks. The cross-implementation0.001 agreement failure is separate from the passing software regression and fresh/cache consistency checks.

## What remains unestablished

Exact organizer novelty levels established:0. Sequence no-hit remains provisional low similarity for local point estimates, not measured zero identity. The separate patent sequence corpus is unidentified/unavailable; SwissProt/PDB/antibody public snapshots differ from unknown organizer versions. Structural reference search is the historical1,562,678-chain public PDB database. Existing complex-predicted coordinates and explicit whole-chain domain assumptions are used; organizer monomer prediction and the unnamed three-domain-predictor consensus are not reproduced. Finite discovery may miss similar folds. Antibody typing/CDR scoring and six missing coordinate models remain gaps. Nothing reinstates any native-context-rejected design or qualifies the four local4 calls for submission.

Boundary choices remain unchanged and named in the script: sequence30% is low and70% moderate; TM0.5/0.8 inclusive; moderate structural coverage requires>70%, high coverage>=70%.

## Primary source and reproducibility

Pinned Foldseek CLI help reports that --exact-tmscore defaults0. The score-export source constructs TMaligner using par.exactTMscore and recomputes qTM normalized by res.qLen. Its helper dispatches to different approximate/exact routines. Source URL/commit/SHA receipts are in ../independent_tmalign/SCORE_EXPORT_SOURCE.json.

[Pinned score exporter](https://github.com/steineggerlab/foldseek/blob/463739e0014a1549a527de589102cde98f802f37/src/strucclustutils/structureconvertalis.cpp), [pinned score dispatcher](https://github.com/steineggerlab/foldseek/blob/463739e0014a1549a527de589102cde98f802f37/src/commons/TMaligner.cpp), [published novelty definition](https://www.adaptyvbio.com/blog/novelty).

Use the same CLI/manifests as the previous structural report with a fresh output directory; exact export is now mandatory by default. The fresh invocation, commands, logs, query domains, raw alignments and all local predictions are in run/. The cache-only command is in CACHED_EXPORT_COMMAND.json; its raw diagnostic file is in local cache with SHA256 in CACHED_EXACT_VALIDATION.json. validate.py checks the same inputs and scores without rerunning discovery.

No generator, GPU/Modal/Runpod service, exporter/ranker, existing novelty_audit.py, account, submission, main or KB-main was changed. New GPU/cloud CPU spend$0.00. Recorded Runpod$56.589988 posted +our Modal estimate$10.718947 =~$67.31, excluding separate RFD4 lane/storage/builds/consultations.
