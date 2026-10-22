# Astra second read: priority change

Recorded from Astra's consultation messages2026-10-22; its longer review may
follow separately. No GPU allocation or additional paid consultation occurred.

- Stop same-recipe seeds. No current rejected structure is an obvious native
  binder being withheld only by a confidence threshold: saved native poses are
  roughly46-61A from the requested patch and have abnormal peptide geometry.
- Put a tiny hard-sequence refolding diagnostic first. Predeclare native v9
  trajectory487e943cff8955d6(69aa) and replicate3/a207b66be0c44216(67aa).
  The second was selected for least adjacent-CA compression among the five
  independent-seed saved structures, before seeing any refold outcome.
- Withhold all designed binder coordinates. Fold each sequence alone first;
  predict a full-native-target complex only if a plausible monomer emerges.
  This tests hard-sequence folding versus coupled soft-sequence search and
  creates validation needed by a later domain-search experiment.
- Held-out AF2 monomer checkpoints1_ptm/2_ptm differ from the multimer design
  pool, but are not a different architecture. Do not overstate independence.
- Monomer-first is resource triage for stable minibinders, not a biological
  impossibility proof; induced folding remains possible.
- Make6W3H optional/time-boxed geometry calibration. Its Fc architecture and
  rearranged receptor fragment cannot validate our minibinder generator and
  should not delay candidate work. Domain generation remains a conditional
  next experiment, not a multi-day prerequisite pipeline project.

Root updated NEXT_MOVES.md accordingly. No submission recommendation or
experimental result is created by this review. Runtime estimates unchanged:
Modal~$9.30; Runpod$35.6652includinghistory; consultations$1.1922999;
setup/storage unreconciled. GPUs remain stopped.
