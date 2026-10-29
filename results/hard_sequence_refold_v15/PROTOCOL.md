# Fixed two-sequence refolding diagnostic

Ricardo instructed implementation and continued work after Astra's second read.
manifest.json pins both source artifacts and hard sequences before any refold.
The GPU receives only sequence strings, a native receptor template and metadata;
no saved binder coordinates or templates. CPU preparation verifies exact inputs,
hard one-hot representation, SEQUENCE-only flags and all-zero positions/masks.

First predict BOTH sequences with BOTH held-out model_1_ptm/model_2_ptm AF2
checkpoints,3recycles,dropoutfalse,large_complex_initializationfalse. Reuse pinned
predict_validation_ensemble without stage_filters (no checkpoint early exit).
All4 monomer jobs precede any large complex prediction. Monomer output is chainA.

Advance only if BOTH models have binder pLDDT>=.70, credible backbone geometry,
no severe nonadjacent backbone overlap and mutual internal-fold CA RMSD<=3.5A.
Predeclared geometry ranges(A): N-CA1.2–1.75,CA-C1.2–1.8,C-O1.05–1.5,
peptideC-N1.1–1.7,adjacentCA2.7–4.2; nonadjacent backbone atoms>=1.5.
These are compute triage, not experimental validation or impossibility proof.

Each advancing sequence receives2 native639-residue monomer-target complex
predictions. This known-feasible memory scope is explicitly NOT dimer-aware
refolding. After actual GPU release, restore/check native11QC subunitB and NAG
clearance using existing canonicalize(). Complex output A=target/B=binder.
Binder-only pLDDT is reported, plusiPTM/pTM, geometry and hotspot distance.
Mutual internal-fold agreement and target-frame binder pose agreement are named
separately. No automatic candidate qualification, pH scoring or submission.

At most8 prediction jobs (up to32 internal network evaluations with3recycles).
Worker30min/provider35min/localmonitor40min; one single-use H100/A10080 allocation,
zero retries/mincontainers0, actual resource-zero verification before CPU work.
Existing authorization applies. H100GPU$3.9492/hr; requestedCPU/RAM-inclusive
upper$4.649328/hr.45min planning reserve$3.486996 plus setup/storage. Prior
Modal runtime$9.2964416754; Runpod$35.6652inclhistory; consultation$1.1922999.

Then follow the observed result: plausible new native complex→robustness/evidence;
monomer-only success→interface-generation pivot; monomer failure→no rescue sweep,
one bounded domain-generation hypothesis and assess a documented fallback.
No same-recipe seed extensions. A stopped campaign releases the actual GPU.
