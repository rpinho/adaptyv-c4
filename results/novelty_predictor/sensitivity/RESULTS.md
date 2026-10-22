# CPU follow-up: candidate discovery sensitivity

The earlier single-discovery structural calls were too optimistic for some retired designs. A focused eight-query test found two of four local level4 calls drop to3, and one local3 drop to2. This is a measured change between local methods, not an organizer false-positive rate.

The predictor now defaults to the union of Foldseek kmer0 and all-reference ungapped1 discovery, with50000 candidates per mode/query. mergedbs combines candidates for the same query database, then the existing full TM-align stage runs. Repeated target hits retain their highest query-domain TM and do not count as additional evidence. The old coordinates, domain intervals, thresholds and historical1,562,678-chain public PDB snapshot are unchanged.

Full follow-up:68 coordinate-bearing queries in214.201 CPU seconds; 24653 unique query/target matches from42718 raw alignments. Every previous best query-normalized TM was retained or improved. All74 input sequences are accounted for; the six missing exact coordinate models remain unestablished.

| Local level, retired designs | Before | After |
|---|---:|---:|
| 2 | 26 | 29 |
| 3 | 37 | 37 |
| 4 | 4 | 1 |
| Unestablished | 6 | 6 |

6 retired local calls were downgraded. Positive control remains local1 with query-normalized TM0.9946.

| Changed retired design | Before TM / level | After TM / level |
|---|---|---|
| c4_tfr1_groove_human_forced_l98_ac2cbaba03c60f1c_candidate8 | 0.7871 / 3 | 0.8058 / 2 |
| c4_tfr1_groove_human_forced_l84_fa5e4e94654ff9c1_candidate1 | 0.7869 / 3 | 0.8232 / 2 |
| c4_tfr1_groove_human_forced_l58_9f94b53776050aad_candidate4 | 0.7910 / 3 | 0.8196 / 2 |
| c4_tfr1_groove_human_forced_l96_9003044e98f7dfc2_candidate6 | 0.4668 / 4 | 0.6250 / 3 |
| c4_tfr1_groove_human_forced_l96_9003044e98f7dfc2_candidate10 | 0.4945 / 4 | 0.5384 / 3 |
| c4_tfr1_groove_human_forced_l95_c6a5afd6447f1fd0_candidate4 | no hit / 4 | 0.5546 / 3 |

No significant sequence hit is still provisionally treated as low sequence similarity. Exact organizer levels established:0. Missing patent sequences, older public snapshots, supplied complex-predicted coordinates instead of organizer monomer predictions, assumed whole-chain domains instead of the unnamed three-predictor consensus, and finite search sensitivity remain unresolved. The remaining local4 calls are not guaranteed passes. Nothing reinstates the retired pool or changes the separately vindicated native-context clash gate.

A broad exploratory CA-spacing diagnostic found0 flags in68 structures (adjacent distances2.5–4.5 Å). This is not a peptide/sidechain/folding or functional validator. GEOMETRY.json preserves every measurement.

Selection and CPU limits were recorded before each follow-up in PLAN.json and FULL_POOL_PLAN.json. The eight-query comparison is in COMPARISON.json; the full comparison, predictions, raw alignments, coordinate fingerprints, parameters and logs are preserved. Fifteen focused unit tests pass, including duplicate-target handling, domain coverage, normalization, sequence units and threshold equality.

Use the same predictor/coordinate manifests as the previous structural report. Default discovery is now union, cap50000. Explicit --structure-prefilter-mode 0 or 1 and --structure-max-seqs remain available for sensitivity checks; higher sequence/structural similarity can lower predicted novelty. Use fresh output directories. No new inference, generator, exporter or account/submission work occurred.

All stages ran locally,4 CPU threads,4G search memory,300seconds per stage. No GPU, Modal, Runpod or cloud CPU service was touched. New GPU spend$0.00. Recorded Runpod$56.589988 +our Modal estimate$10.718947 =~$67.31; separate RFD4 lane, storage/builds/consultations excluded.

The combined run additionally found that candidate10 on the96-residue family reached
TM0.5384 in the combined run, above both its original0.4945 and the focused
ungapped-only0.4967. The sole remaining local4 has best detected TM0.4760,
not an absent hit. It remains provisional under the unchanged database/domain gaps.
