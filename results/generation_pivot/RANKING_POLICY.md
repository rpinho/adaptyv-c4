# Provisional ranking policy for the measured generation arms

Registered before looking at the new scored-arm outputs. Reuse and adapt the
canonical ranker; do not create a second ranker or submission emitter.

Join the canonical sequence export to inline metric rows by the exact design
ID within each run, retaining the arm and source hashes. Reject duplicate IDs
or conflicting sequences rather than picking one. Preserve missing/error
metric rows explicitly; do not turn them into measured zeros.

The first pass is **geometry triage**, not a prediction of the competition
ranking. Order fully measured rows lexicographically by descending
`hotspot_obedience.hotspot_obedience`, ascending
`clash_metrics.frac_atoms_clashing_1p5A`, then ascending
`bond_length_metrics.bond_length.frac_unreasonable`. Stable source/name ties;
retain the existing declared family-diversity policy and expose its settings.
Keep incomplete rows in the evidence table with the missing fields, rather
than claiming their score is zero or that they failed an official gate.

These three measures have different meanings; no invented weighted affinity
or pH score will be computed. Do not relabel hotspot obedience as the old
Hotspot_Contact_Fraction without establishing metric equivalence. No iPTM,
pLDDT, pAE, interface ddG or developability values are supplied by `geom`.
Neither total histidine count nor its global prior establishes a pH switch.

Apply the real restored-context audit to the actual shortlist, retaining all
audit values and unresolved fits. Distinguish computational triage, CSV-format
validity, organizer novelty acceptance and experimental binding. Missing
experimental confirmation is not itself an official submission prohibition.
Ricardo alone decides and submits the final entry; no automatic upload.

No new generation, sequence modification, ranking output or qualification is
claimed by this policy document. Two metric arms are approved after the
validated novelty report commit; see GPU_APPROVAL.json. Recorded
Runpod$56.589988; own prior total about$67.31; new GPU/cloud compute$0.00.
