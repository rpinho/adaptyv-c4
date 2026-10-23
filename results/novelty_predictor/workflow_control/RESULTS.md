# Generic CPU workflow: real PDB end-to-end control passed

`scripts/run_cpu_novelty.py` now accepts an arbitrary `name,sequence` CSV and an
optional coordinate manifest. It runs the demonstrated sequence → whole-chain
context → three CPU domain tools → domain search → parent-context alignment
workflow without editing the retired-pool benchmark inputs. Usage and scope:
`../USAGE.md`.

The actual integration test supplied a monomer PDB, renamed the native195-aa
control to `GENERIC_NATIVE_CONTROL`, and deliberately omitted chain and source
hash. The workflow acquired/verified the fingerprint, identified chainA by the
unique full sequence, converted PDB for backbone processing, verified780 backbone
atoms, ran all3 CPU domain methods, searched the consensus domain, and scored
the additional whole-chain-context references. All6 workflow stages completed.

Measured result: **100% best sequence identity; local novelty level1; best domain
TM0.9927**, with the 2–194 consensus domain covering193/195 residues. JSON and CSV
agree. Official organizer level remains unestablished. This is a native software
control, not a challenge design or submission. `VALIDATION.json` binds actual
inputs, source versions, outputs and checks. The final report path is in
`run/STATUS.json`.

**23 focused tests pass.** Input-integrity checks demonstrate unique-chain
identification and rejection of ambiguous chains/changed expected fingerprints.
The structural command-path regression now also checks that reused numeric levels
clear to null when coordinates or accepted domains become unavailable. This
corrects a stale-value bug; it changes no named novelty threshold. Earlier domain
scoring evidence uses the predictor recorded at commitb32bdba, and its source
fingerprint is verified against that Git version rather than falsely treating it
as the newer code. Its measured counts remain20/40/7/6.

Absent coordinates give a sequence-only report with structural levels withheld.
No accepted consensus domains also withhold numeric structural novelty. The tool
does not infer protein structures from sequences on CPU, does not guarantee the
organizer score, and does not create a submission. It uses the explicit public
database snapshots and pinned local domain runtime; no missing source is replaced.

New GPU/cloud CPU spend **$0.00**. Recorded Runpod **$56.589988** + our Modal
estimate **$10.718947** = **~$67.31**, excluding the separate RFD4 lane, storage,
builds and consultations. No generation, account, submission, cloud allocation,
exporter/ranker or main-branch action occurred.
