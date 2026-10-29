# Method: proton-coupled two-state design, with the numbers checked

## Oct23 independent-review response: no new optimization executed

The new request called for sequence-specific filtering, ranking and generation
for targeted blood-brain-barrier delivery. Those phases were not executed under
the assistance boundary stated in this session. No new GPU/cloud allocation was
started. Existing sequences, ranks, coordinates, thresholds and CSV are unchanged.
Read-only canonical format/uniqueness checks passed for the existing20-row CSV;
no new submission packet was produced and no entry was submitted. See
outbox/adaptyv-c4-independent-review-status.html and
results/independent_review_response_20261023/STATUS.json for completed work,
review-source fingerprints, scope corrections and unexecuted phases.
New GPU/cloud compute spend$0. Recorded Runpod$56.589988; own tracked prior
estimate$73.311827 plus unreconciled Modal host charges, excluding the other lane
and account-wide billing. Current listed deployed apps have0 tasks; prior stopped
generation-app release receipts are retained.

## Oct23 exact shortlist audit and CPU His-scan complete

The exact final shortlist20 was independently rerun through the unchanged
6W3H-validated omitted-context gate and existing full-native staging helper.
20/20 pass both;0 omitted-context backbone pairs below2A; target CA fit
0.2291–0.3467A. Every numeric gate row matches the earlier80-pose subset.
Minimum omitted-protein backbone clearance20.970A; minimum full-native-target
backbone clearance3.013A. Gate and staging source hashes remain unchanged.

The actual upstream TMolMutationDdgMetric completed all10 preregistered
contacting-His sites in the top5 designs, on local CPU:30 mutation deltas and
40 scored WT/H-to-ALA/H-to-GLN/H-to-HIS control legs. Complete observed native
receptor dimer1278 residues retained. Every leg retained expected residue
identities; maximum measured backbone displacement0.0016735A, below the declared
0.002A audit tolerance. Scored coordinate snapshots are preserved by fingerprint
in the local cache. CPU runtime is pinned:TMol0.2.0, torch2.14.1, AtomWorks3.0.4;
installed metric/helper bytes match the upstream source. Two instrumentation
preflight failures are preserved separately and are not biological failures.

Both His-removal deltas are positive at3/10 sites, belonging to original
shortlist ranks1 and3. Other sites give mixed or negative signs. Maximum absolute
H-to-HIS self-substitution delta0.35712 score units; do not assume a zero control.
These numbers are relaxed SITE/receptor interaction-energy changes in TMol
beta2016 score units, not calibrated kcal/mol, rigorous binding ddG, pKa, KD or
pH selectivity. Neither pH nor protonation equilibrium was varied. No automatic
functional pass threshold or new ranking was introduced.

The unchanged canonical build_submission.py emitted the20-row CSV for Ricardo's
review, preserving every original name, sequence and rank and attaching actual
scan evidence. Other15 designs are explicitly not_scanned, with missing energy
fields left blank. This establishes format readiness, not organizer novelty,
independent refolding, expression/developability, mouse binding or biological
qualification. Twin-Strep exclusion and upload replacement/accumulation remain
for Ricardo to resolve. No account, upload or entry was created.

Read outbox/adaptyv-c4-his-scan-review.html; CSV and complete numeric evidence:
results/rfd4_his_scan_20261023/{submission_for_ricardo_review.csv,RESULTS.json}.
The finite CPU controller and its awake keeper both exited. No GPU/cloud task
was started; new GPU/cloud spend$0.00. Recorded Runpod$56.589988; own tracked
prior estimate$73.311827 plus unreconciled Modal host charges, excluding the
other lane and account-wide billing. Prior first-pair stopped/zero-task proof
is retained; that app aged out of the current recently-stopped listing. The
higher pair remains explicitly listed stopped with0 tasks. See RELEASE_CHECK.json.

## Oct23 generation and CPU review complete; all owned GPU apps released

Four approved inline-geom arms completed200 outputs each,800 total and798
unique amino-acid sequences. All three ranking geometry fields are finite
for800/800 outputs; target sequence violation count is0. Designs with at least
three His: plain7/200, prior0.10 33/200, prior0.20 135/200, prior0.30 183/200.
Mean His count rises0.66/1.28/6.325/24.97. This measures conditioning response,
not improved binding, pH selectivity or a benchmark pass rate.

Canonical per-arm geometry ranking retained20 rows each for an80-pose pool.
All80 passed the UNCHANGED omitted-context gate and full-native staging;
target CA fit RMSD range0.2127–0.3467A. The registered adaptive contact-hypothesis
policy selected20 with global sequence-family cap3, retaining all80 annotations
and every final-ranker deferral. All20 have descriptive His-N/basic-receptor
proximity within4.5A; this is not an energy, protonation or affinity calculation.
Seven selected sequences contain at least30% His, a visible composition concern,
not a validated developability threshold. Many geometry ties remain.
The canonical emitter is unchanged; its format-valid CSV, workflow draft,
full evidence and limitations are linked in outbox/adaptyv-c4-generation-review.html.
This review set is not a final recommendation to consume20 competition slots.
No independent native-complex refolding, interface ddG/His scan, pH-binding assay,
mouse-binding validation or fresh complete novelty screen of these800 outputs
has been performed. Organizer acceptance and novelty remain unmeasured.
Ricardo alone decides and uploads any actual submission; none was uploaded.

The first pair finished generation but blocked in sampler cleanup; raw failure
and completed-output recovery are preserved separately. The canonical sampler
now has bounded shell-free cleanup. The higher pair returned successfully and
automatically stopped its app. Its raw controller status incorrectly treated
CLI “already stopped” as failed release; independent positive stopped/zero-task
state is saved in RELEASE_VERIFIED.json. The controller now accepts that positive
state even when the CLI returns nonzero, without accepting empty-account evidence.
Both OWNED generation apps are verified stopped with zero tasks. No paid successor
was launched. Four shutdown-control tests pass; the CPU runtime switch has numeric
parity on overlapping audited poses. See results/rfd4_final_review_20261023/VALIDATION.json.

Recorded Runpod$56.589988; prior own Modal estimate$10.718947. Four-arm GPU-time
upper estimate after release$6.002892: first pair$3.129680, higher pair$2.873211.
Prior tracked total plus that GPU estimate$73.311827, PLUS unreconciled host
CPU/RAM/build/storage charges. These are own-lane estimates, exclude Claude's
separate lane and are not an account-wide invoice. The final review used local
CPU only. New GPU spend for the CPU review$0.00.

## Oct23 measured completion; novelty work stops here

Recovered cached exact-score export completed1469.78s; all3,523,855 requested
pairs returned and400-row final validation passed. Local general-scale counts
73 level1/92 level2/212 level3/12 level4/11 unknown. Preserved143 earlier
points unchanged. Ten no-candidate domains affect7 unresolved inputs;4 other
inputs have no accepted domain. Official score/type/corpus and experimental
success remain unmeasured. See adaptive_completion/RESULTS.md and receipts.
Commit this completed report and pivot to the two explicitly approved
canonical inline-geom generation arms. No further novelty campaign.
Recorded Runpod$56.589988; own prior total about$67.31; new cloud compute$0.00.

## Oct23 export recovery and generation priority

Full native alignment completed; its exact-score export exceeded1200s.
run_completion.py now supports guarded --resume-export with7200s export
timeout, preserves the failed attempt and fingerprints the completed native
cache. Export-only retry is active;31 focused tests pass. Do not repeat the
alignment or declare resolved levels before validation. After committing
this report, stop expanding novelty work and follow the existing generation
handoff. Read-only geometry preflight and bounded invocation are prepared
under results/generation_pivot/. GPU yes-after-quote remains required and
approved for the two named arms after the validated report commit; see
generation_pivot/GPU_APPROVAL.json. Runpod$56.589988 +own Modal estimate$10.718947 =~$67.31; new cloud
compute$0.00, other-lane/account-wide billing excluded.

## Current CPU classification: measured bounds, then remaining full pairs

Actual400-input bounded screen completed:73 local1/70 local2/257 unknown,
zero organizer scores or qualified designs. Preserve those143 positive-bound
point records unchanged. Full TM-align/exact export is now running locally
on the domain plus whole-parent retrieval union:3,523,855 unique pairs.
Ten unresolved domains have no candidates and retain unknown evidence.
Eight-hour finite CPU limit; no GPU. Plan, native-source audit,96-pair exact
seed check and automatic validated finalizer:adaptive_completion/.30 software
checks pass, plus a real400-input simulated-empty-export negative control.

Separate antibody reference analysis annotated107,440 canonical heavies,
excluding1,902 noncanonical sequences.50 query H loops have both IMGT flanks.
22/50 subset levels change with edit-distance normalization;49/50 scalar
calls remain withheld under missing-reference completions and normalization
sensitivity. One conditional4 has only two retrieved references. This does
not establish organizer typing or a novelty-gate pass. Evidence/gaps:
results/novelty_predictor/antibody_cdr_reference/. No sequence generation,
functional improvement, ranking or qualification measured. Recorded Runpod
$56.589988 +own Modal estimate$10.718947 =~$67.31; new GPU/cloud compute$0.00,
separate generation/account billing excluded.

## Actual incoming400:CPU novelty evidence, partial alignment continuation

Both saved RFD4 arms are now local and all400 sequences match verified chains.
Five sequence searches completed. The16,451,768-record structural discovery
exceeded the original exhaustive alignment's1200-second stage budget. Its
failed run is retained, not described as a complete result. A separate bounded
positive-evidence continuation reuses verified completed work and withholds
conclusions that unfinished comparisons could change. Policy, active status
and receipts:results/novelty_predictor/incoming_rfd4_20261022/.

CPU ANARCII annotation and controls completed separately;85 assigned variable
domains among400, not official structural types or antibody novelty levels.
See results/novelty_predictor/antibody_scope_audit/RESULTS.md. No generator,
mutation, account or submission action. Recorded Runpod56.589988 +own Modal
estimate10.718947 =~67.31; new GPU/cloud compute0.00, separate lane excluded.

## Arbitrary-input CPU novelty workflow verified

scripts/run_cpu_novelty.py now accepts CSV+saved PDB/CIF coordinates and performs
the demonstrated sequence, whole-chain-context,3-method domain and expanded-domain
alignment workflow. Renamed native PDB control passes end to end: local1,
qTM.9927,193/195 domain coverage,780 verified backbone atoms.23 tests pass;
missing coordinates/consensus clear numeric levels rather than retaining stale
whole-chain values. No structure inference, GPU/cloud resource or organizer-score
claim. See results/novelty_predictor/USAGE.md and workflow_control/RESULTS.md.
Retired-pool measured counts below remain unchanged; no design promoted.

## Current measured CPU domain policy (Oct22)

See results/novelty_predictor/domain_consensus/RESULTS.md. All3 domain tools ran
on all68 available models;70 domains include4 discontinuous inputs. Explicit
index correction/source checks and21 predictor tests pass. Whole-chain code
compatibility checked on all68 frozen cases. Additional whole-parent targets
were aligned against every domain,26,998 pairs, removing8 retrieval no-hits.
Current provisional retired counts20 local2/40 local3/7 local4/6 missing; native
control remains1/qTM.9927.5 local4 calls arise from insufficient whole-query
domain coverage despite matching domains; do not equate them with new folds.
All thresholds/reference/sequence evidence/coordinates remain fixed. No official
score or qualified design established. Earlier counts below are frozen alternatives.

## Public patent evidence and CPU domain preflight (Oct22)

NCBI pataa public collection now obtained and searched explicitly:3,982,722
records/787,058,342 residues,0 qualifying retired-design hits,110/110 identity
source control. It is not verified as the organizer patent corpus. Full68-domain
exact rerun is identical to frozen evidence; provisional counts16/47/4/6 unchanged.
Current report/manifest: results/novelty_predictor/patent_extension/RESULTS.md
and DATABASES.json. Initial partial-coordinate preparation error preserved;
corrected setup inherits and checks all68 frozen coordinate domains.

TED-style CPU domain preflight passed all3 tools on native control, consensus
2-194 of195. Published consensus offset bug reproduced and explicitly corrected
in an isolated compatibility copy,6 consensus/5 real-filter checks pass; no IoU,
size or novelty threshold changed. Full census is still in progress; no new
domain-based novelty level claimed. See domain_consensus/INDEXING_CORRECTION.md.

## Saved monomer coordinate sensitivity — CPU comparison (Oct22)

results/novelty_predictor/monomer_sensitivity/RESULTS.md preregisters and measures
an unchanged-predictor comparison using every exact saved-monomer sequence match.
51 CIFs examined;4 matched variants of one78aa parent. All4 remain local3 versus
their complex coordinates,0 changed levels, paired CA RMSD.256–.480A; native
control remains local1/qTM.9928. This is sample stability only, not organizer
monomer-policy validation. Full-pool exact-score analysis remains authoritative;
no favorable conformer selection, new prediction or design promotion occurred.

## Current CPU novelty checkpoint — exact score-export defect corrected (Oct22)

Generation, RFD4 files, Modal apps and GPU campaigns belong to
claude/submission-and-control. Do not edit/duplicate its build_submission.py,
rank_candidates.py, modal_rfd4.py or configs/rfd4/. Exporter comparison done;
our emitter is superseded. NO GPU without Ricardo yes and quoted hourly rate.

6W3H remains GATE VINDICATED for compatiblecopy, alternatecopy INCONCLUSIVE.
Its complete mapping/unchanged gate proof is in results/real_gate_control/.

A standalone TM-align audit exposed our score-export defect: convertalis defaults
approximate even after full TM-align. Production now forces --exact-tmscore1.
Do NOT use earlier approximate-score structural/sensitivity counts for ranking.
Current evidence: results/novelty_predictor/exact_score/RESULTS.md and
run/predictions.json. Cached score-only rescore preserves every alignment record;
a fresh full pipeline265.026s gives exactly the same24653 query/target scores.
Current67 retired local calls:16 level2,47 level3,4 level4;6 missing coordinates.
Positivecontrol local1,exact queryTM.9928.16 local calls changed because of the
export correction. This creates/qualifies no candidate; retired pool stays retired.

16 focused unit tests pass, including behavioral exact/default export regression;
negative control fails with approximate0 and passes exact1. Separate cross-tool
fixed-alignment agreement check still FAILS: max difference.00225 vs prereg.001,
down from.05369. Winner's exact.4477 matches standalone fixed.4477, but numerical
equivalence across all26 pairs is not established. See independent_tmalign/.
No thresholds or agreement tolerance were relaxed.

All levels remain provisional: missing patent corpus, older public snapshots,
supplied coordinates, whole-chain domains vs unnamed consensus, finite discovery
and organizer-policy/antibody gaps persist. No exact organizer level established.
Existing novelty_audit.py, exporters/ranker and all generation files unchanged.

New GPU/cloud CPU spend0.00. Recorded Runpod56.589988 posted +our Modal10.718947
estimate =~67.31, excluding separate RFD4 lane/storage/builds/consultations.
No account, submission, main or KB-main write. No automatic GPU successor.

## Latest instruction — generator lane STOPPED; CPU real-gate control COMPLETE (Oct22)

Ricardo explicitly superseded this generator lane with the separate approved
RFD4-Proteina work on claude/submission-and-control. Do not restart RFdiffusion,
MPNN, AF2 successors or any GPU from earlier plans below. No generator duplication.
All our Modal functions zero runners/inputs/backlog; all3 Runpod pods EXITED.
Release occurred before replying. No successor is authorized in this lane.

New CPU-only task completed: results/real_gate_control/RESULTS.md. Independent
P02786 annotation accounting:157 mapped residues plus9 unresolved tag/linker
residues=166. Actual fit() used; unchanged audit_full_context.main() exercised
in isolated fixtures; all source/helper gate hashes unchanged. GATE VINDICATED
for FcB/TfRC, CA RMSD.963/.950, zero omitted-context BB pairs<2A, min10.717/10.737A;
full-native BB min2.884/2.899A. Alternate FcA/TfRD is INCONCLUSIVE, fit3.709A.
All copies/orientations retained, no outlier trimming or threshold adjustment.
The66 remain retired; originalone>2A fit is an inconclusive context transfer.

Our lane Modal10.718947 estimated +Runpod56.589988 posted =~67.31 total;
new control GPU0.00. Other RFD4 lane/consultations/builds/storage excluded.
V18 did complete50 real diffusion steps,0 generator-accepted backbones; its
sparse-contact rejection and exploratory-only coordinates are historical.
No MPNN sequence generated. Earlier pending/next-run plans below are superseded.
Portal: single_chain; once/24h; novelty>=3/4; closesNov2 11:59UTC.
Ricardo handles every account and entry. No main or KB-main writes.

## Current checkpoint — sequence-only refolds completed; native-domain experiment prepared

V15 completed eight coordinate-free fixed-sequence prediction jobs; actual GPU
release independently verified. Both free folds pass both checkpoint geometry
checks (pLDDT .796-.935; internal model agreement .98/1.28 A). Native receptor
complexes remain off-patch (52.38-70.50 A), iPTM .148-.197; three restored-context
clashes, zero credible TfR1 binders. This separates plausible sequence folding
from failed interface generation. No pH or mouse-binding claim.

Astra's outcome and domain review are in results/astra_review/. Exactly one
native contiguous A189-383 generation experiment is prepared: only the target
representation changes vs native v9; same seed, lengths, losses and gates.
Installed-image CPU frames, native identity, exact hotspots and author-numbered
CIF->PDB output roundtrip pass; 80 local checks pass. A direct CIF reader uses
label IDs and failed the first CPU build; the corrected test uses the actual
CIF->PDB author-ID downstream path. No GPU was allocated during this correction.
Protocol: results/apical_context/DOMAIN_PROBE_DECISION.md. Launch/live evidence
belongs in results/native_domain_v16/. The post-release extraction prepares one
predeclared hard sequence even if its raw search pose fails; full-native refolding
still requires its own manifest-specific CPU image preflight. Nothing is admitted
or submitted automatically. No same-recipe seed successor is scheduled.

Modal completed runtime estimate $9.5172782; Runpod runtime/history estimate
$35.6652; consultation usage $1.1922999, setup/storage unreconciled. H100 GPU
$3.9492/hr, requested CPU/RAM upper $4.649328/hr combined. Continued spending is
already explicitly authorized. Runpod remains stopped. Earlier live states below
are historical; use current supervisor files for actual allocation/execution.


**Current qualification (2026-10-21):** the old 73-residue receptor crop omitted
protein that all 65 well-aligned source candidates overlap at backbone level.
Their linkage results diagnose that crop; they do not validate binding to the
full receptor. See `results/full_context_audit/RESULTS.md`. Independently
prepared free-binder sidechains also remove the strongest proposed H46 free-pKa
advantage in the 72-model scout: `results/free_networks_v3/RESULTS.md`. These
findings, and fresh full-context generation, supersede advancing the old numerical
winners directly to submission. The equations below remain model definitions.

## The metric

`scripts/ph_linkage.py` scores a complex by Wyman linkage. For every titratable group with
free-state pKa `pf` and complex pKa `pb`:

    L = log10 Ka(pH_low) - log10 Ka(pH_high)
      = SUM_i [ t(pb,pf,pH_low) - t(pb,pf,pH_high) ],   t(pb,pf,pH) = log10((1+10^(pb-pH))/(1+10^(pf-pH)))

L is in log10 units of affinity. **L <= -2 is our internal 100-fold model target,
not an official challenge cutoff or a calibrated detector limit.** The lab defines
measurable versus undetectable binding. PROPKA values are computed on complex,
isolated binder and isolated target and matched by chain and residue number.
An extracted chain retains the bound conformation; it does not establish the
independently relaxed free-state ensemble. See `WIN_PLAN.md` for the current strategy.

## Calibration (this corrects the first version of the plan)

The plan originally claimed "2-3 interface histidines with pKa depressed 1-1.5 units gives about
2 log units". **That is wrong, and the first real scorer output is what exposed it.** On a PDL1
smoke design, a binder histidine with pKa 6.27 free -> 5.37 bound (a clean 0.9-unit depression)
contributed only **-0.34** log units, not ~1.

`scripts/linkage_budget.py` computes why. Contribution per group, pH 6.0 vs 7.4:

| pKa_free | depression 1.0 | 2.0 | 3.0 | 5.0 |
|---|---|---|---|---|
| 6.0 | -0.24 | -0.28 | -0.28 | -0.28 |
| 6.6 | -0.49 | -0.62 | -0.63 | -0.63 |
| 7.0 | -0.61 | -0.86 | -0.89 | -0.90 |
| 7.4 | -0.61 | -1.02 | -1.11 | -1.12 |
| 8.2 | -0.32 | -0.95 | -1.28 | -1.34 |

Two things follow, and both change the design instructions:

1. **A histidine saturates at about -1.1 to -1.3 log units**, and only when its *free-state* pKa is
   around 7.4-8.2 **and** the depression on binding is 2 units or more. A histidine sitting at the
   textbook 6.0-6.5 with a 1-unit shift buys about -0.25, so a design would need five or six of
   them. The real target is therefore **2 to 4 interface histidines with free-state pKa 7.0-8.0
   whose pKa drops by at least 2 units on binding**, which is a much sharper instruction than
   "install some histidines".
2. **The optimum is not free.** A histidine with free pKa 7.4 is already half protonated at pH 7.4,
   so if protonation disfavours binding, affinity at the operating pH is reduced roughly twofold.
   That is an acceptable price for a hundredfold switch, but it has to be stated, and it means the
   pH 7.4 binding energy has to be strong enough to absorb it.

## What that implies for the mechanism families

Raising a histidine's *free* pKa and depressing it hard on *binding* are two separate asks, and
they point at a specific construction:

- Put a carboxylate (Asp/Glu) next to the switch histidine **on the binder's own surface**, not
  across the interface. In the free state that pairing raises the histidine's pKa into the useful
  7.0-8.0 window; on binding the pair is desolvated and its geometry broken, which depresses the
  pKa sharply. This is a refinement of mechanism family 1 and it is what the numbers above argue
  for, rather than the across-the-interface dyad the plan first described.
- Mechanism family 3 (burying an unpaired histidine against an apolar interface) gives the large
  depression but tends to leave the free-state pKa near 6.5, so it should be paired with family 1
  rather than used alone.
- Mechanism family 6 (the target's own His707 or His603) is unaffected by this: the target supplies
  a histidine whose free pKa is whatever the receptor makes it, and the binder's job is to depress
  it. Worth measuring the free-state pKa of His707 and His603 on the apo receptor before committing.

## The two-estimator rule still holds

PROPKA is empirical. The second-estimator requirement remains, but the previous identification of
PROPKA-driven PDB2PQR/APBS as an independent pKa estimator was incorrect. PDB2PQR
can reuse PROPKA's titration assignment; APBS then evaluates electrostatics for
those charges. A genuinely independent titration calculation must be configured
and checked before claiming independent agreement. PDB2PQR is available locally;
an APBS executable was not found on PATH in the strategy review. See the
[PDB2PQR documentation](https://pdb2pqr.readthedocs.io/en/stable/using/index.html).

`scripts/install_switch.py` supplies a neutral-state structural filter: negative
reference/variant InterfaceAnalyzer dG and ddG <= 3 REU. It is neither a measured
pH 7.4 affinity nor a second pKa estimator.

## The acid budget: the term the first two versions of this method missed

**Measured 2026-10-20, on real BindCraft2 output, not assumed.** The pH layer had never been
executed end to end. Running it on the two ranked Challenge 3 designs (`results/c3_ranked/`, the
first real designs this project produced) gave the result that changes the selection rule.

Design `c3_glp1r_ecd_groove_forced_l99_e26affcf9b8831e4_seq0`:

| group | pKa free | pKa bound | dpKa | contribution to L |
|---|---|---|---|---|
| **B33 HIS** (designed) | 6.43 | 3.76 | **-2.67** | **-0.521** |
| **B26 HIS** (designed) | 6.20 | 1.87 | **-4.33** | **-0.386** |
| A29 N+ | 7.96 | 7.76 | -0.20 | -0.049 |
| B34 GLU | 5.05 | 5.89 | +0.84 | +0.192 |
| **A127 GLU** (target) | 4.83 | 6.49 | +1.66 | **+0.534** |
| **A68 GLU** (target) | 4.92 | 7.91 | **+2.99** | **+1.255** |
| | | | | **total L = +1.20** |

Two buried histidines delivered exactly the mechanism this method is built on - a 2.7 and a 4.3
unit depression, worth **-0.91** between them - and the design still came out with the **wrong
sign**, because three carboxylates on the *target* had their pKa pushed **up** on binding and
contributed **+1.91**.

The physics is symmetric and I had only used half of it. Burying a neutral histidine against a
desolvated interface depresses its pKa, which releases a proton on binding and makes binding weaker
at low pH. Burying a *charged* carboxylate does the mirror image: desolvation and neighbouring
acids raise its pKa, so it **accepts** a proton on binding and makes binding **stronger** at low
pH. A glutamate whose pKa rises from 4.9 to 7.9 is a proton acceptor sitting right in the
pH 6.0-7.4 window, and a single one of those outweighed both designed histidines.

Three consequences, all of them acted on:

1. **The design metric is a budget, not a count of histidines.** `scripts/ph_linkage.py` now writes
   a per-residue CSV (`--csv` also emits `*_residues.csv`) and prints both ends of the sum, because
   the positive tail is where a design is lost. Never read the His rows alone.
2. **Epitope choice is the only lever on the target's acids** - the receptor cannot be mutated.
   `scripts/epitope_charge.py` censuses every titratable group a binder would bury within 10 A of
   the hotspots. Measured on the candidate set:

   | epitope | patch | Asp/Glu | His | Lys/Arg | formal net | acid:His |
   |---|---|---|---|---|---|---|
   | **E1 His707 (human, live campaign)** | 49 | **2** | **2** | 6 | **+4** | **1.0** |
   | E1 His710 (mouse ortholog) | 67 | 4 | 2 | 6 | +2 | 2.0 |
   | E2a apical 202-213 | 53 | 5 | **0** | 7 | +2 | 5.0 |
   | E2b apical 340-376 | 100 | **10** | 1 | 11 | +1 | 10.0 |
   | E4 His603 | 33 | 5 | 1 | 3 | **-2** | 5.0 |

   The epitope the campaign is running on is the best of the set on this axis, and it is the best by
   a wide margin: two acids against two target histidines, in a patch carrying six basic residues.
   A base-rich patch is also what depresses a buried histidine's pKa, so the same property helps
   twice. This is an **independent confirmation of a choice made for structural reasons** (dimer
   occlusion and 7/7 human/mouse epitope identity), arrived at from different physics.
3. **Both apical sites are now dead for a second, independent reason.** They were already
   compromised by 0.59 human/mouse identity across the patch. E2a carries **no histidine at all**,
   so there is no lever to install against, and E2b buries ten carboxylates. E4 His603 is
   acid-dominated and formally negative. Nothing in the portfolio should be spent on them.

### What sets the sign, measured rather than reasoned

The first draft of this section said to declare the interface acids as `coldspots`. That is wrong
on its own terms - A714 GLU **is one of the chosen hotspots** - and the per-residue table says
something sharper. Sorting all 33 groups by contribution:

| group | dpKa | contribution | why |
|---|---|---|---|
| A88 TYR | **+3.84** | **+0.000** | pKa 13.2 -> 17.0, never near the window |
| A67 ASP | **-1.57** | **-0.000** | pKa 2.95 -> 1.38, never near the window |
| A138 GLU | -0.66 | -0.009 | acid, pKa went DOWN, helps |
| A68 GLU | +2.99 | **+1.255** | pKa 4.92 -> **7.91**, lands inside the window |

**Size of the shift is almost irrelevant. What counts is how much of the pH 6.0-7.4 window the
titration crosses.** A 3.8-unit shift outside the window is worth exactly nothing; a 1.7-unit shift
into it costs half a log unit. Acids are not the enemy - acids whose pKa is pushed *into* the window
are, and an acid whose pKa goes *down* (salt-bridged to a binder cation) contributes the right sign.
The rule is therefore **charge-compensate every buried carboxylate**, not exclude it.

### The ceiling on the target's own histidines (this kills mechanism family 6 alone)

propka3 on the apo trimmed human target (`results/apo/E1_his707_tight.pka`) gives the free-state
pKa values the plan said to measure next:

| group | apo pKa | model | ceiling at total desolvation (pKa -> 2.0) |
|---|---|---|---|
| HIS 707 A | **6.09** | 6.50 | -0.328 |
| HIS 699 A | **5.81** | 6.50 | -0.205 |
| HIS 684 A | 6.35 | 6.50 | -0.473 |
| GLU 694 A | 4.60 | 4.50 | (+1.18 if it rises as A68 GLU did) |
| GLU 714 A | 4.55 | 4.50 | (+1.16 if it rises as A68 GLU did) |
| GLU 533 A | 4.55 | 4.50 | (+1.16 if it rises as A68 GLU did) |
| **all three His together** | | | **-1.01** |

The gate needs **L <= -2.00** and the target can supply at most **-1.01**, at a physical limit no
binder reaches. **Mechanism family 6 cannot clear the gate on its own** - the receptor's histidines
are a useful bonus of roughly -0.3 each, not a strategy. The plan must not rest on them.

Worse, a *single* interface glutamate rising as far as A68 GLU actually rose in a real design
(+2.99) contributes **+1.18**, which cancels one whole well-designed histidine. The downside term is
the same size as the upside term.

### Where the leverage actually is: free-state pKa, not dpKa

The same arithmetic, applied to a binder histidine:

| binder His free pKa | driven to | contribution each | two of them |
|---|---|---|---|
| 6.0 | 2.0 | -0.28 | -0.56 |
| 7.0 | 3.5 | -0.90 | -1.79 |
| **7.4** | **3.5** | **-1.12** | **-2.23 (clears the gate)** |
| 7.8 | 3.5 | -1.26 | -2.52 |

This explains the measured design completely. Its histidines were depressed by 2.67 and 4.33 units -
far more than the 1.5 units needed to saturate - and returned only -0.52 and -0.39, because they
*started* at 6.43 and 6.20. **The design objective is not a large dpKa. It is a high free-state pKa
that collapses on binding.** Everything past about 1.5-2.0 units of depression is wasted effort.

So the binder specification is two-state and explicit:

- a histidine whose **free-state pKa is 7.4 or above**, which means siting it where the *unbound
  binder* presents carboxylate neighbours and solvent exposure;
- **buried against the target's cations when bound** - and the E1 His707 patch carries six Lys/Arg
  (net formal +4), which is exactly the environment that collapses a histidine's pKa. The property
  that makes this epitope acid-poor is the same property that makes it a good switch site;
- **two such sites minimum**, three if any interface carboxylate cannot be charge-compensated;
- A699 and A707 HIS counted as a bonus of about -0.5 together, never as load-bearing.

**Honest statement of what this is.** PROPKA is empirical and single-structure, and these pKa shifts
are computed on predicted complexes, not crystal structures. The budget is a *ranking* tool for
selecting among designs and epitopes, not a prediction of the assayed pH ratio. That is why the
two-estimator rule above is not optional, and why nothing is claimed about a design until the
independent APBS arm agrees on the sign.

## Status of the layer (2026-10-20)

- `scripts/cif2pdb.py` **new and needed**: BindCraft2 writes an `_atom_site` loop with no
  `occupancy` column, which Biopython's `MMCIFParser` rejects with a bare `KeyError`. The converter
  reads the loop by its own header order, never by fixed columns, so it survives a writer change.
- `scripts/ph_linkage.py` **executed on real designs for the first time**, propka3 3.5.1 in the
  local CPU venv (`.venv/`, Python 3.9). The pH layer deliberately runs on CPU here rather than on
  the rented GPUs, which stay on trajectories.
- `scripts/epitope_charge.py` **new**, see above.
- `scripts/switch_sites.py` **new**. `install_switch.py` ranked candidate positions by *interface
  burial*, which optimises only the bound half of the two-state problem and is precisely the error
  the measurement exposed. This ranks them on both halves from coordinates alone - exposure plus
  carboxylate neighbours in the free binder against burial plus target cations in the complex - so
  the expensive Rosetta pass runs on a short list. On the only two positions where an independent
  propka answer exists it agrees in order and within about 20% in magnitude (B33 HIS: estimate
  -0.62, propka -0.52; B26 HIS: estimate -0.50, propka -0.39). Two points is consistency, not
  validation, and it is a pre-screen, never a pKa claim.
  It also rated four positions in that design *above* the two histidines BindCraft2 placed, which
  is the use case: pick switch sites deliberately instead of accepting whichever histidines the
  generator happened to leave at the interface.
- `scripts/install_switch.py` remains **untested on real designs** - no Challenge 4 design has been
  accepted yet - and must be re-pointed at `switch_sites.py`'s ranking instead of its own burial
  sort before it is used. It is next, the moment a design lands.
- Best achievable pair on the two existing designs scores about -1.4 and -1.65 estimated, short of
  the gate, which is consistent with the budget above: two sites are the floor, not the plan.

## Where the binder actually sits, measured (2026-10-21)

The layer above was built and calibrated before any design had been observed docking. Thirty-seven
refolded complexes now exist, and they move one assumption from reasoned to measured.

**The epitope is the hydrophobic groove, not the His707 patch.** 37/37 designs contact
LEU709/LEU713/LEU716/PHE685/PHE701/TRP702/PHE731/LEU735/THR739/ILE742 and 0/37 contact
SER704/GLY705/SER706/HIS707/ALA711. HIS707 sits 4.0 A from LEU709, at the rim.

**What this changes for the proton budget, and it is an improvement.** The method already concluded
that the leverage is the *binder's* free-state pKa, not the target's histidines (ceiling -1.01 for
all three, against a -2.00 gate). A binder histidine buried against a hydrophobic groove loses more
solvation than one pressed onto a polar patch, so the bound-state pKa drops further and each
histidine is worth more than the -1.115 estimated against the polar patch. The groove is also
base-rich (LYS717, LYS720, ARG719, LYS534, ARG732 against GLU714, GLU728, GLU533), and the measured
rule from the previous section is that bases help and acids fight: this epitope's charge census is
on the right side of the ledger, which the apical sites were not.

**What must now be re-measured rather than inherited.** Every number in the sections above was
computed against designs docked on the polar patch, or against the apo target. The per-residue
contributions, the acid budget and the `switch_sites.py` rankings all have to be recomputed on
groove-docked designs before any of them is quoted again. The *rules* carry over - sign is set by
how much of pH 6.0-7.4 a titration crosses, acids can go either way and must be charge-compensated
rather than excluded - but the magnitudes do not.

**Honest status.** The groove's own effect on a designed histidine's bound pKa is a prediction from
the desolvation term, not yet a measurement. It is testable the moment designs clear the filters:
run `ph_linkage.py` on groove-docked complexes and compare the per-residue contributions against
the -0.52 and -0.39 measured on the polar patch. If the groove does not deepen them, the
improvement claim in this section is wrong and should be struck.

## The groove prediction was wrong, and the reason is the rule I failed to apply (2026-10-21)

The section above predicted that burying a designed histidine against a hydrophobic groove would
deepen its contribution beyond the -0.52 and -0.39 measured on the polar patch, and stated that if
it did not, the claim gets struck. It did not. **Struck.**

Three groove-docked designs, same backbone, same histidine position (C60):

| design | binder His C60 free pKa | bound pKa | dpKa | contribution | total L |
|---|---|---|---|---|---|
| cand3 | 6.69 | 6.15 | -0.54 | **-0.335** | -0.643 |
| cand2 | 5.52 | 5.00 | -0.52 | **-0.079** | -0.452 |
| cand5 | 5.42 | 4.95 | -0.47 | **-0.061** | -0.399 |

Against -0.52 and -0.39 on the polar patch, the groove made the binder histidine **worse**.

**Why, and it is the method's own measured rule turned back on me.** This layer already established
that the leverage is the **free-state** pKa, not dpKa. Desolvation does not act only on the bound
state: a histidine positioned in or beside a hydrophobic groove is *already partly buried in the
free monomer*, so the same apolar environment that depresses the bound pKa depresses the free pKa
too. The free pKa lands at 5.42-6.69 instead of near 7.4, the titration therefore barely crosses the
pH 6.0-7.4 window, and dpKa collapses to -0.47 to -0.54. The two effects cancel, and on net the
groove loses. I predicted a one-sided benefit from a two-sided mechanism.

**What did work, for a different and already-documented reason.** Total L moved from **+1.20 on the
polar patch to -0.64 here** - a 1.84-log improvement that also **flips the sign to the required
direction**. That gain is the **acid budget**, not desolvation: the polar patch's three glutamates
contributed +1.91 and dominated everything, and the base-rich groove simply does not have them. So
the epitope change was right, and right for the reason the acid-budget section gave, not the reason
the groove section guessed.

**Where His707 actually ended up.** It is contacted by none of the designs, yet contributes
-0.103 to -0.121 in all three, because at 4.0 A from LEU709 the binder perturbs it from the rim.
HIS699 contributes more (-0.120 to -0.225). Small, consistent, and in the right direction - but
together the target's histidines still sit inside the -1.01 ceiling measured on the apo structure.

**The gate is still missed and the route to it is now specified by data.** L = -0.64 against a
-2.00 gate leaves about 1.4 log units. The requirement is explicit: a designed histidine must be
**solvent-exposed in the free binder** (free pKa near 7.0-7.4) and buried **only** on binding.
C60 fails the first half. This is exactly the quantity `switch_sites.py` was written to score - free
-state exposure and acidic neighbours for a high free pKa, against burial and target cations for a
low bound pKa - and `install_switch.py` has never been wired to use it, still running its own
burial sort. That wiring is now the highest-value remaining step, and these three accepted
backbones are the first real substrate for it.

**Caveat on the arithmetic.** The target is split across chains A (532-743) and B (686-693) and this
run scored chain A only, so chain B's eight residues are absent from the sum. They are not at the
interface the designs use, but the figures above are chain-A-complete rather than target-complete.

## The free-state pKa rule, now confirmed prospectively on 66 designs (2026-10-21)

The campaign completed on its own at its configured 25 final designs (177 trajectories, 14%
acceptance, against 0 of 99 on the two mis-aimed epitopes). Sixty-six passing complexes were
scored. This is the first time the layer has had a population rather than a handful, and the
population settles the central rule.

| quantity | value |
|---|---|
| best L | **-1.175** |
| median L | -0.387 |
| designs with L < 0 (correct direction) | **62 / 66** |
| designs clearing the -2.00 gate | 0 / 66 |
| binder histidines scored | 68 |
| their free pKa | 3.54 - 7.16, mean **6.17** |
| of those, free pKa >= 7.0 | **3** |

**The rule predicts which histidines matter, before looking at their contributions.** The three
deepest single-histidine contributions in the whole set are **-0.822, -0.659 and -0.658**, and all
three are the *same position* - C24 on the l58_9f94b53776050aad backbone - at free pKa **7.14, 7.16
and 7.14**. They are also precisely the three histidines in the set with free pKa >= 7.0. The other
65 histidines, averaging free pKa 6.17, produce the mediocre bulk. Sorting by free-state pKa and
sorting by contribution give the same answer, on 68 independent cases. The rule was inferred from
two residues on 10-20; it now holds on a population, and it is what `switch_sites.py` scores.

**And the design carrying that histidine is the best design overall** (L = -1.175), which is the
useful form of the result: the single free-state pKa of one residue, computable before any linkage
sum, ranked the whole campaign.

**Two levers, and the data separates them.** C24 goes 7.14 -> 5.67. The earlier estimate for a
histidine making the full excursion (7.4 -> 3.5) was -1.115 each. We are getting -0.822 because the
free pKa is 7.14 rather than 7.4 *and* because the bound state only reaches 5.67 rather than 3.5 -
the burial is incomplete. So the gap to the gate decomposes into: raise the free pKa toward 7.4
(position the histidine in a solvent-exposed, acid-flanked site in the free monomer) and deepen the
bound-state burial (place it where the groove and the target's cations desolvate it). Those are the
two terms `switch_sites.py` already scores separately, and nothing in the design loop is using
either yet.

**Honest arithmetic on the gate.** At -0.822 per well-placed histidine, the -2.00 gate needs about
2.4 of them, or two that make a fuller excursion. Two histidines at the full 7.4 -> 3.5 excursion
would give -2.23 and clear it. That is a design target, not an achievement: no design in this set
has two histidines above free pKa 7.0, and whether two can be placed on one backbone without
spoiling the interface is untested.

**What this does not show.** All 66 designs come from one campaign on one epitope, and the 25 ranked
designs derive from far fewer distinct backbones, so these are not 66 independent samples. Chain B
(686-693) is excluded from every sum, as noted above. And no orthogonal estimator has been run:
every figure here is PROPKA3.

## Ranked two-site installation, now executed on all 25 ranked C4 backbones (2026-10-21)

`scripts/install_switch.py` now uses `switch_sites.score_design` directly. There is
no second burial sort. The CPU campaign tested the three best estimated pairs
from eight eligible sites on each ranked backbone, in best historical-L order.
One best-scored candidate represents each backbone; its sequence can differ from
ranked seq0, and both identities are recorded. Existing His sites participate in
HH and QQ twins rather than being dropped as no-op edits.

Real validation exposed two problems the unexecuted installer had concealed:

1. **Reference preparation mattered more than the mutation.** The best raw source
   gave InterfaceAnalyzer dG +804.75 REU when only separated partners were globally
   repacked. Common whole-complex repacking and side-chain minimisation, with the
   backbone fixed, gave -44.66. Local mutation-site packing alone compared a strained
   bound input with globally repacked free partners. The installer now prepares one
   common complex, then applies identical local protocols to WT, HH and QQ. Both
   reference and variant dG must be negative, in addition to ddG <= 3 REU. These
   remain neutral-state structural surrogates, not pH 7.4 affinity measurements.
2. **An omitted His is not a zero contribution.** PROPKA suppresses the C1 ring in
   two variants by coupling it to N+. Those variants retain energy-filter status
   but are explicitly **unscorable** for complete two-site linkage. Partial sums
   are preserved only for diagnosis and excluded from ranking. A regression test
   prevents an apparently passing partial sum from entering selection.

The completed screen uses **AB as target, C as binder**, and all matched groups
(`min_shift=0`). Historical A-only scores order the substrates; comparisons below
use matched prepared AB references, not historical campaign numbers.

| Quantity | Result |
|---|---:|
| Ranked backbones processed | 25/25 |
| HH pairs modelled, with matched QQ attempts | 75 |
| HH passing neutral-energy filters | 38 |
| HH with complete installed-site linkage | 36 |
| Complete-linkage HH with retained QQ control | 25 |
| HH energy-retained but unscorable | 2 |
| Complete HH with L <= -2 | **0/36** |
| Complete HH with both selected free pKa >= 7 | **0/36** |
| Best complete HH L | **-1.530** |

The best pair is **C:D4H/T46H** on `l63_59b0f47ecf134650`: its matched reference L
is -0.351, HH L -1.530, and QQ L -0.444. HH therefore improves L by **-1.179**
against its reference and **-1.086** against QQ. Rosetta ddG is **+1.74 REU**,
inside the 3-REU filter, with negative reference and mutant dG. Selected free pKa
values are **6.33 and 7.25**: one useful high-free-pKa site, still not two. This is
a computed screen result, not experimental binding or clearance of the pH gate.

The native-coordinate heuristic has a limited real-data check: on 27 native
interface histidines across 17 selected backbone IDs, estimated versus PROPKA
contribution Spearman correlation is 0.567 and mean absolute error 0.153 log units.
Free-pKa estimate mean absolute error is 0.711. This does not validate a mutation
pKa prediction. Rotamers also matter: native C24 on the original best source moved
from free pKa 7.14 to **6.71 in its first prepared reference**, despite unchanged
backbone coordinates. Native-position ranking is a pre-screen, not evidence that
mutating its top two sites creates two elevated free pKa values.

Every saved source/output checksum, mutation identity, target sequence, energy
filter and complete linkage sum was audited. **181 PDBs across all 25 backbones
passed; maximum backbone atom displacement was 0.000 A.** Seven regression tests
passed. Diagnostic runs are retained and clearly marked invalidated/superseded.

Limits: three pairs per backbone, not all 28 shortlist combinations; related
backbones; original trimmed receptor with artificial fragment boundaries; ranking
on source rotamers before preparation; one pKa estimator. PROPKA warns about
backbone groups at A683 and A728. Full-context confirmation, APBS, mouse validation
and orthogonal refolding remain outstanding. No account or entry was created or
submitted. Both rented pods stayed stopped; **new GPU spend $0.00**, historical
last recorded pod spend $33.87, persistent storage not freshly reconciled.

Full evidence: `results/switch_install/campaign/RESULTS.md`, `variants.csv`,
`residues.csv`, `integrity.json`, per-design manifests and PDBs. Protocol diagnostics
and reproduction commands are in `results/switch_install/VALIDATION.md` and
`runtime_setup.txt`.

## Coordinate initialization needs an explicit relative-placement preflight

The v4 attempt exposed a separate initialization problem, not evidence for the
cause of v3 rejection: enabling upstream initial_atom_pos retained native target
coordinates while unresolved binder atoms stayed at origin,188.816A from the
hotspot centroid. That attempt was cancelled and its GPU released. No comparison
of accepted-model quality is inferred from the aborted attempt.

The replacement target is the same entire639-residue monomer, rigidly translated
so origin lies15A outward from the hotspot centroid, with7.955A nearest-atom
clearance against the complete native dimer. Atom identities, residue numbers,
sequence and all CA pair distances are preserved to PDB precision. This changes
the relative starting placement of the unresolved binder, not receptor shape.
It remains an optimization hypothesis, not a physical binding measurement.
Opposite-subunit restoration/native-sequence and clash checks remain mandatory
before any accepted structure receives CPU linkage scoring. Exact transform and
integrity measurements: results/modal_campaign_v5/human_11qc_local_frame.audit.json.
