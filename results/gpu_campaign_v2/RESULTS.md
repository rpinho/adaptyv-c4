# Completed GPU pilot: no accepted scaffold

Three trajectories, zero accepted designs. One reached annealing but missed interface confidence (i_pTM0.31); two failed initial binder confidence (pLDDT0.44 and0.53). No pH switch, mouse cross-reactivity or measured affinity was established. Three trials do not establish that the epitope is unusable.

The POD was stopped automatically and independently verified EXITED. Combined runtime compute estimate across both GPU attempts: $1.7952. Fresh posted billing for the migrated pod: GPU $1.736618, pod disk $0.009028, total $1.745646; billing can lag, and network-volume storage is separate. No workers or GPU are left running by this completed attempt.

The CPU continuation started before the supervisor wrote its initial status file and exited with FileNotFoundError. Fixed with a missing/partial-record retry, atomic queue state writes and explicit follow-up failure reporting. Added two regressions; all17 tests pass. Replayed the queue after verified stop: no accepted scaffolds, so no pKa models were invented or scored. Original failure log preserved in continuation.log.

Next decision is in ../gpu_campaign_v3/DECISION.md: test native-target optimization before changing epitope or weakening filters. The original $3.58 total compute cap remains binding.
