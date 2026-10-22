# CPU novelty predictor: current evidence and limits

Current five-source sequence + exact structural comparison:
[patent_extension/RESULTS.md](patent_extension/RESULTS.md). Public NCBI pataa was
obtained and searched; it is not verified as the organizer's patent corpus.
Retired local counts remain16 level2,47 level3,4 level4,6 coordinate-unavailable.
Exact organizer levels established:0. Updated sequence-source manifest:
patent_extension/DATABASES.json. Earlier manifests are preserved as history.

The exact-export correction and residual numerical disagreement are documented
in exact_score/RESULTS.md. Four paired saved-monomer calls remain local3 in
monomer_sensitivity/RESULTS.md. Domain-consensus CPU feasibility work is separate
under domain_consensus/; no unmeasured consensus is assumed in these counts.

Everything below describes earlier sequence-only/approximate-score checkpoints.
Its missing-source statements and old counts are historical, not current results.

## Task1: existing6W3H real-gate result reviewed first

GATE VINDICATED for deposited FcB/TfRC: unchanged audit target CA RMSDs0.963/0.950 Å,
zero restored-context backbone pairs<2 Å, omitted-context minima10.717/10.737 Å.
Whole-native-dimer backbone minima2.884/2.899 Å include the original interface.
Alternate FcA/TfRD is INCONCLUSIVE: target fit3.709 Å, irrespective of zero clashes.
All mapped CAs, both copies and both native subunit orientations were retained.
The independent annotation reconciliation is157 native-mapped positions +9 explicitly
annotated, unresolved expression-tag/linker positions=166; no mapping segment is
missing. All332 per-copy construct positions are in the saved mapping table.

What ran: the actual run_context_ensemble.fit(mobile,fixed), both fit directions,
and unchanged audit_full_context.main() with an isolated ROOT fixture. Before/after
hashes match. This was completed in commitfd85fa1 before the new lane-boundary
message; this task reviewed those existing artifacts rather than rerunning it.
The original66 stay retired; the originalone with fit>2 Å remains an inconclusive
context transfer. What was not established: broad sensitivity/specificity from
one experimental system, experimental binder performance, or minibinder admission
of the >150-residue Fc. Complete control evidence: ../real_gate_control/RESULTS.md.

## Task2: what ran and what it measured

New file: ../../scripts/novelty_level_predictor.py. Neither novelty_audit.py nor
the exporter/ranker, RFD4 files, Modal apps or GPU campaigns were changed.
Installed CPU MMseqs2 version18-8cc5c locally. Public downloads were normalized
with source URLs, HTTP Last-Modified dates, SHA256 fingerprints and explicit
per-family status in DATABASES.json. No local reference set was substituted.
Identical sequences were deduplicated within each database only; ambiguous
reference letters were retained, reference case normalization counted, and PDB
nonprotein polymers excluded by their declared mol:protein marker.

| Named family | Distinct sequences searched | HTTP Last-Modified | Actual obtained scope |
|---|---:|---|---|
| SwissProt |487476|2026-09-03|Full public canonical FASTA snapshot,575748 entries before deduplication|
| PDB |188017|2026-10-17|Full public SEQRES protein snapshot,1109692 protein chains before deduplication|
| Therapeutic antibodies |2211|2025-03-18|Thera-SAbDab public CSV: all four heavy/light and bispecific sequence columns|
| PLAbDab |372631|2024-08-28|Both paired heavy/light and unpaired GBSeq_sequence public downloads|
| Separate patent sequences |not obtained|unknown|Article does not identify provider, corpus or release; PLAbDab is not a replacement|

These are downloaded snapshots, not verified organizer snapshots. In particular,
Thera-SAbDab/PLAbDab public bulk downloads are older than today; recent entries may
be absent. This gap persists even though their downloads succeeded. The article's
phrase therapeutic-antibody database is interpreted as Thera-SAbDab and explicitly
labelled as that interpretation; the organizer has not identified its exact source.
No account or access-controlled resource was used.

Input was the actual73 retired project sequences plus one exact native195-aa
TfR1 sequence as a positive control; none are recommended submission candidates.
Analysis CSV has only name,sequence. INPUT_SCOPE.json preserves source fingerprint.
The four database searches all completed in20.358 seconds total at4 CPU threads,
4G search memory limit, --gpu0, sensitivity7.5, real alignment identity mode3,
E-value<=0.001, maximum10000 prefilter hits/query, no identity or coverage filter.
Identity is computed from integer nident/alnlen, avoiding rounded fident values
at the decision boundaries; alignment backtraces are enabled. The initial search
before this precision change is retained separately. The best identity among
qualified returned hits is selected, rather than identity
of the highest bit-score hit. Raw hit tables, commands, logs, alignments, query/
target coverage and per-hit E-values are saved in retired_pool_run/.

Results:

- 73/73 retired sequences: no significant hit under these settings; best detected
 identity is null, not numerical zero. Zero final novelty levels established.
 Missing patent and structural evidence leave possible general levels1–4.
- Positive control:100% identity,195/195 aligned residues, query coverage1.0,
 PDB11qc_A positions68–262, E-value5.458e-127. It returns sequence class high,
 possible levels1–2 and an explicit risk of failing the>=3 gate. No exact level
 is claimed without structure. SwissProt32 and PDB34 qualifying alignments were
 returned for this control; antibody database searches returned zero.

No significant sequence hit does not prove identity<=30% against all sequences,
prove a new fold, or qualify a design. Source MMseqs search parameters, versions,
coverage filters and database snapshots are unpublished, so our search policy is
an approximation. The10000 prefilter cap can miss a high-identity hit that scores
poorly; cap saturation is flagged. The native positive control establishes real
search wiring, not recall near30% for short de novo proteins. No portal predictions
were compared with actual upload feedback.

## Boundary policy and rule evaluator

Named constants in the new script fix equality explicitly:

| Quantity | Equality policy |
|---|---|
|Sequence30%|low; moderate requires>30%|
|Sequence70%|moderate; high requires>70%|
|TM-score0.5|qualifies for moderate-domain coverage|
|TM-score0.8|qualifies for high-domain coverage|
|Moderate structural coverage70%|does not qualify; requires>70%|
|High structural coverage70%|qualifies; requires>=70%|

The general-scale rule evaluator implements all nine sequence/structure combinations,
with high structure treated as at least moderate for level1. Its structural helper
sums disjoint query-domain coverage at each TM cutoff; it refuses overlapping domains
and scores normalized by anything except query-domain length. This helper is only
a rule evaluator exercised with synthetic boundary fixtures, not measured structure
search evidence. The main CLI intentionally does not accept unverified structural
category assertions. It writes predicted_level=null when the evidence permits
multiple levels, with possible_levels and conditional levels for each structural
class. A missing database can raise sequence similarity; it cannot lower a real hit.

Ten unit tests pass, including four sequence equality cases, all nine level
combinations, six TM/coverage edge cases, overlap/normalization rejection,
missing-evidence bounds, identity-unit parsing, exact integer identity counts versus rounded output,
highest-identity versus highest-score selection, and strict query input validation.
The real four-database run plus positive-control assertions also passed.

## Structural half: documented gap, not a prediction

Not run: sequence-to-structure prediction, the consensus of three domain predictors,
Foldseek structural search, or query-domain-normalized TM-align. CSV sequences alone
do not provide coordinates/domains. I did not achieve that full pipeline within
this CPU-only task; its prediction/runtime assets and organizer reference snapshot
are not installed or identified, and the article does not name the three predictors.
Foldseek and TM-align can themselves run on CPU when coordinates and domains exist;
this is not a claim that those searches intrinsically require a GPU. No GPU was
requested or started. General scale only: antibody typing and the separate CDR
novelty scale are also unimplemented.

Consequently this deliverable is a working, measured sequence half and a tested
rule/bounds evaluator. It is not yet a complete ProteinTyper replica or a dependable
portal gate predictor. No design has been assigned an invented level3 or4.

## Use

CPU invocation on a new analysis CSV (use a fresh output directory):

    python3 scripts/novelty_level_predictor.py designs.csv --databases results/novelty_predictor/DATABASES.json --out results/my_novelty_screen

The manifest references local cached FASTA files and verifies their hashes before
searching. On another machine fetch the same sources or provide explicitly named,
fingerprinted replacements in your own manifest. Download/normalization scripts
and raw-source fingerprints are saved here; live URLs may later change. Missing
families must remain explicit. The main predictor never downloads, substitutes
or uploads data automatically. Output: predictions.json and predictions.csv plus
raw search evidence. Empty predicted_level means unestablished, not a passing score.

Rules: https://www.adaptyvbio.com/blog/novelty
MMseqs documentation: https://github.com/soedinglab/MMseqs2/wiki
PDB FASTA documentation: https://www.rcsb.org/docs/programmatic-access/file-download-services

## Exporters, lane boundaries and spend

Read-only comparison is in EXPORTER_COMPARISON.md (three lines). Their emitter
lacks overwrite protection and source/output SHA256 fingerprints. No missing
portal-format capability found; our emitter is superseded and work on it stops.
The single canonical exporter/ranker and all generation belong to the other lane.

No GPU work or service allocation occurred. New GPU spend$0.00. Recorded cumulative
Runpod$56.589988 posted +our Modal$10.718947 conservative estimate =~$67.31;
other RFD4 lane, storage, CPU builds and consultations excluded. Both tasks are
CPU-only; no cloud CPU compute was rented. Everything remains on the Codex branch,
with no main/KB-main changes, account creation or submission.
