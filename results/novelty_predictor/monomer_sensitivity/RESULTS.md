# Saved monomer versus complex coordinates: four paired local calls stable

CPU-only comparison completed with the production novelty predictor. All four
sequence-matched monomer predictions give **provisional local level 3**, as do
their frozen exact-score complex predictions: **0 of 4 level changes**. This
establishes stability only in this small, related sample. No official novelty
score, new qualified design, or binding result is established.

## What ran

The preregistered inventory examined 51 existing monomer CIF files, with zero
parse errors. Exactly four full-chain sequences matched the 73 retired analysis
designs. All four are variants of the same 78-residue parent. The other 69 retired
inputs had no exact saved-monomer match. No sequence was trimmed or repaired, no
coordinate model was selected using its score, and no new structure was predicted.
The six previously missing complex-coordinate cases remain unresolved.

The unchanged production predictor searched these four monomers and the unchanged
native positive control. It reused fingerprint-verified sequence evidence for all
74 analysis rows, marking the other 69 rows unavailable **for this comparison**.
It used the same historical 1,562,678-chain PDB structural database, whole-chain
domain assumption, union of discovery modes 0/1, 50,000-hit cap per mode, full
TM-align, and explicit exact TM-score export. Four CPU threads, 4G search memory,
300 seconds per stage; measured structural runtime **76.086 seconds**. The output
contains 8,437 raw alignment rows and 4,585 distinct query/target pairs.

## Measurements

Names below share the prefix
`c4_tfr1_groove_human_forced_l78_cc1057879d35a8ec_`.
RMSD fits every sequence-corresponding CA, monomer onto complex, with the actual
`run_context_ensemble.fit(mobile, fixed)`. TM values are each search's best match,
normalized by the complete 78-residue query; they are not monomer-to-complex TM.

| Design | Paired CA RMSD (Å) | Complex best qTM | Monomer best qTM | Complex local level | Monomer local level |
|---|---:|---:|---:|---:|---:|
| candidate1 | 0.409 | 0.6967 | 0.7017 | 3 | 3 |
| candidate2 | 0.480 | 0.7537 | 0.6961 | 3 | 3 |
| candidate3 | 0.256 | 0.6780 | 0.6691 | 3 | 3 |
| candidate5 | 0.345 | 0.7071 | 0.7100 | 3 | 3 |

The native positive control remains local level 1, qTM **0.9928**, with identical
structural evidence. Sequence evidence, database, domains, settings and production
predictor source are unchanged. Coordinate and raw-output fingerprints verify;
all 74 input rows are accounted for. `COMPARISON.json` contains exact RMSDs,
source hashes, best-reference identities and checks; `COMPARISON.csv` is the
paired table. The raw production outputs are under `run/`.

## What this cannot establish

Four related models with paired CA RMSD below 0.5 Å cannot establish robustness
across different folds, uncertain monomers, or scores near a class boundary.
These saved monomer predictions are not the organizer's ESMFold/ESMFold2-Fast
models. The patent sequence corpus, organizer snapshots/settings, and domain
consensus remain missing. Finite structural discovery and the previously measured
standalone TM-align residual disagreement remain unresolved. No retired design
was promoted, and no favorable coordinate form was selected for a submission.

The authoritative full-pool exact-score analysis remains
`../exact_score/run/predictions.json`: 16 local level 2, 47 level 3, 4 level 4,
6 coordinate-unavailable retired inputs. This paired experiment does not replace
that analysis or turn its local estimates into official novelty levels.

## Reproduction and cost

The frozen input inventory and structure manifest were produced by `prepare.py`
after `PREREGISTRATION.md` was written. The production command is:

```sh
.venv/bin/python scripts/novelty_level_predictor.py \
  results/novelty_predictor/analysis_inputs.csv \
  --sequence-report results/novelty_predictor/retired_pool_run/predictions.json \
  --structures results/novelty_predictor/monomer_sensitivity/STRUCTURES.json \
  --structure-database results/novelty_predictor/structural/DATABASE.json \
  --foldseek /Users/ricardo/.local/opt/foldseek/foldseek/bin/foldseek \
  --threads 4 --timeout-structure 300 --out /tmp/novelty-monomer-fresh-run
```

Use a fresh output path. `summarize.py` verifies the frozen run and writes the
paired report without rerunning searches. Runtime binary DBs may be moved to
local cache with their fingerprints retained in `RUNTIME_CACHE.json`; raw hits,
query PDBs, logs and JSON evidence remain in the repository.

New GPU/cloud CPU spend **$0.00**. Recorded Runpod **$56.589988** posted + our Modal
**$10.718947** estimated = **$67.308935 (~$67.31)**. This excludes the separate
RFD4 lane, storage/builds and consultations. No cloud service was started, no
submission/account action occurred, and no main branch was written.
