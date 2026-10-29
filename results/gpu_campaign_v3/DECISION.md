# Paired native-target control after the unsuccessful pilot

**Provider/budget update:** the host-capacity queue below is cancelled. Ricardo
explicitly authorized paying more and using his existing Modal account. The
same prepared experiment is migrating to Modal; the old $3.58 cap is superseded
for that two-hour pilot. Fresh prices and live paths:
`../modal_campaign_v3/README.md`. The text below records the original plan.

The completed pilot accepted 0/3 trajectories: one reached annealing with i_pTM0.31, two failed initial binder confidence (pLDDT0.44/0.53). Three tries cannot establish an epitope failure. The old forced-targeting run replaced371/639 target residues with lysine during optimization.

Test that representation before changing the epitope or relaxing acceptance gates: keep the complete native monomer, the same patch, binder-length distribution, campaign seed, contact losses/weights, accessible-termini objective and all finite acceptance filters. Disable lysine replacement. This is a controlled computational hypothesis, not a prediction that native sequence will improve binding. Hotspot coverage >=0.5 remains mandatory. Two trajectories maximum. Capture trajectory structures for diagnosis. Opposite-subunit restoration, native sequence verification and clash checks still precede downstream pKa scoring.

Total compute cap remains $3.58 across all attempts; deduct the prior $1.7952 runtime estimate. Rate must be checked before resume. No account, entry, main or KB writes. Stop the POD on completion, failure or cap.
