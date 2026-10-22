# Supplied-coordinate sensitivity: saved monomers versus complexes

Registered before searching any monomer in this comparison. CPU only; no new
fold prediction, generation, GPU allocation, or submission.

Question: does replacing a saved complex-predicted binder with an already saved
monomer prediction change its provisional local novelty level? This tests one
coordinate-policy gap. It cannot validate the organizer's ESMFold/ESMFold2-Fast
policy or official novelty score.

Inventory every existing CIF under results/ with "monomer" in its path, excluding
this novelty work. Read first-model CA sequences with the production reader.
Include every exact full-chain sequence match to analysis_inputs.csv, regardless
of its model confidence or prior novelty result. Never repair or trim a sequence
to make it match. Report nonmatches and parsing failures. If multiple distinct
monomers match one input, stop rather than choose based on their scores.

Keep all 74 analysis rows in the input to reuse the fingerprint-verified sequence
report. Supply the matched monomers and the unchanged native positive control;
mark other inputs unavailable for this comparison. Compare each matched design
against its frozen exact_score/run/predictions.json result. Use the same local
PDB database and explicit whole-chain domain assumption. Use default union
prefilter modes 0 and 1, 50,000 candidate cap per mode, full TM-align, exact score
export, four CPU threads, 4G search memory, and 300 seconds per stage. No changes
to thresholds, source code, search settings, or conformer selection after results.

Report paired sequence-ordered CA RMSD after rigid fitting, each best query TM,
local level, and changed-level count. A changed level is evidence of coordinate
sensitivity in this paired sample; zero changes establishes only stability in
this sample, not general coordinate independence. Keep missing coordinates and
unmatched sequences explicit. Do not promote any retired design or select a
submission on these scores. The native positive control must remain local level1;
otherwise the run is invalid for comparison and requires diagnosis.

Recorded spend: Runpod $56.589988 posted + our Modal $10.718947 estimate =
approximately $67.31, excluding the separate RFD4 lane, storage/builds and
consultations. New GPU/cloud CPU spend: $0.00.
