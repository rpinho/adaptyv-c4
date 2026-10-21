# Switch installation — CPU validation, 2026-10-21

Status: completed and audited, on branch `codex/install-two-state-switch`. No KB changes.

Runpod: neither pod started or otherwise modified. Incremental GPU spend $0.00;
last recorded historical pod spend $33.87 (handoff). Persistent-volume storage is
not included in this session's GPU figure; no fresh billing API read is available.

## Implementation and substrate

`install_switch.py` now calls `switch_sites.score_design` directly and preserves
its order and exclusion flags. Positive/zero heuristic contributions and sites
with no burial gain are excluded. It models two sites together, with matched HH
and QQ twins, including native His/Gln sites rather than dropping no-op mutations.
Default effort is one pair from eight eligible sites; `--max-pairs` can expand it.
The pair sum orders work only; it does not predict total L or account for coupling.

`run_switch_campaign.py` selects one best historical-L candidate per backbone
from the 66 scored complexes. All 25 ranked backbone IDs are covered exactly once,
ordered by their best L. The ranked seq0 and best candidate can differ; both IDs
are recorded. This uses the best scored sequence/complex for each ranked backbone
rather than silently treating a candidate score as a seq0 score.

C4 chains are **AB target, C binder**. PDB chain/number identifiers are mapped via
`pdb2pose`, never treated as pose indices. Historical L used chain A only; the new
run includes both A and B. Compare twins with their **relaxed AB reference**, not
with the historical chain-A score.

The reference, HH and QQ receive the same sequential local mutation/repacking
and fixed-backbone side-chain minimisation, starting from the same input complex.
InterfaceAnalyzer uses the same scorefunction for all three; ddG > 3 REU rejects
a twin. This is a neutral-state structural surrogate, not measured pH 7.4 affinity.
A retained switch without a retained QQ control is not a complete matched pair.

## Evidence written so far

- `planning/campaign.json` and 25 per-backbone `plan.json` files: complete CPU-only
  ranking and pair provenance. No measured mutant pKa or affinity is claimed.
- Five regression checks passed under Python 3.9: real C4 scorer order/native His,
  missing target fragment rejection, empty selection, PDB-to-pose identity mapping,
  and exact 25-backbone coverage in best-L order.
- `ph_linkage.py` now fails a score if either free-state PROPKA run fails, rather
  than silently building a partial linkage sum. This is required for mutant QA.

## Runtime and references

A separate `.venv-rosetta` uses Python 3.11 and the native macOS arm64 quarterly
PyRosetta build, leaving the existing Python 3.9 scoring environment intact.
PyRosetta is CPU software; a rented GPU is unnecessary for installation/PROPKA.

Official runtime and API references consulted:
- [PyRosetta downloads](https://www.pyrosetta.org/downloads): quarterly builds and pip mirror.
- [Apple Silicon](https://www.pyrosetta.org/downloads/m1-macs): native arm64 builds.
- [InterfaceAnalyzer API](https://graylab.jhu.edu/PyRosetta.documentation/pyrosetta.rosetta.protocols.analysis.html).
- [BindCraft PyRosetta interface scoring](https://github.com/martinpacesa/BindCraft/blob/main/functions/pyrosetta_utils.py): compatibility with the newer DockingPartners setter.

PROPKA remains a single empirical estimator on predicted, locally relaxed
structures. APBS and orthogonal refolding remain outstanding. Rosetta energy and
PROPKA pKa agreement do not constitute two independent pKa estimators.

## Unmodified-source control (complete)

`source_full_target.json` rescored all 25 selected source complexes with the existing
PROPKA 3.5.1 environment, target **AB**, and `min_shift=0.0` (all matched groups).
This succeeded on all 25. The first three original-source L values are -1.085,
-1.013 and -0.991, versus the historical A-only/truncated sums -1.175, -1.104 and
-1.081. Inclusion of B affects free target A pKa as well as adding B terms, so
simply adding the B residue contributions does not correct a historical A-only sum.
The new campaign retains every matched group's contribution, even below the old
0.1-pKa shift cutoff, so installed-site free pKa remains available when a shift is
small. These unmodified values are separate from the relaxed references used for ddG.

The native C24 example has heuristic free/bound pKa 7.43/6.38 and estimated
contribution -0.638, versus PROPKA 7.14/5.67 and -0.822. It ranks first on that
backbone; this is a useful site-order consistency check, not calibration of the
heuristic on new mutations. Mutation outcomes still require real modelling.

## Native-site ranking check (complete)

`native_his_calibration.json` pairs the heuristic with PROPKA for 27 native
interface histidines on 17 of the selected backbones. Spearman rank correlation
between estimated and PROPKA contributions is **0.567**; mean absolute contribution
error is **0.153 log10 units**, and mean absolute free-pKa error is **0.711**.
See `native_his_calibration_summary.json`. These are selected, related native
structures, not independent validation samples or mutation predictions. The
measured free-pKa rule does not by itself validate this coordinate heuristic.

## First real Rosetta run exposed a reference-preparation error

The first unprepared smoke run is preserved under `smoke_unprepared/` and explicitly
invalidated there. Its reference dG was **+626.91 REU**; a relative ddG threshold
alone was therefore insufficient. Diagnostics on the unmodified best-L input gave:

| Bound-state preparation | InterfaceAnalyzer dG (REU) |
|---|---:|
| Input coordinates, separated partners repacked | +804.75 |
| InterfaceAnalyzer also packs input | +8.93 |
| Shared global repack plus side-chain minimisation, backbone fixed | -44.66 |

The bound input had unrelaxed rotamers, while `pack_separated=True` globally
repacked free partners. Local mutation-site repacking did not remove input/core
strain elsewhere. Fix: prepare the entire complex once, preserve that prepared
PDB, then apply the same local protocol to the reference and both twins from that
common prepared pose. PROPKA scores the actual saved coordinates. Merely enabling
InterfaceAnalyzer input packing on a clone would score coordinates different from
the saved variants, so that shortcut is not used.

Both reference and variant must now have **negative dG**, as well as ddG <= 3 REU,
to retain a twin. These are structural surrogates, not measured affinity. This
preparation sensitivity is itself a limitation on interpreting absolute energies.
Ranking uses the original source coordinates; actual mutant pKa is recomputed
after side-chain preparation and installation. No backbone coordinates are
optimised. The final screen covers the **three highest-ranked pairs per backbone**,
not an exhaustive search over the eight-site shortlist.

## Installed-site completeness is now required

The first full pass exposed two retained HH variants whose **C1 histidine ring was
absent from PROPKA's summary in both states** (backbones `l78_cc1057879d35a8ec`
and `l60_63126a071c041270`). This was not a regex/parser error. The installed
PROPKA 3.5.1 API shows the ring coupled to the N-terminal N+ group and suppressed
from the independent-group summary. See `terminal_his_diagnostic/coupling.json`
and the raw `.pka` files. Internal suppressed-group pKa is diagnostic only; it
cannot be inserted into an independent-site Wyman sum by hand.

`score_model` now checks for every installed histidine in the residue terms. A
missing site preserves `partial_linkage` for diagnosis, marks
`linkage_status=incomplete_selected_histidines`, and removes the total L and its
reference delta from selection. Structural retention is a separate fact. Such a
variant is **unscorable**, not a pH-gate failure or a passing two-site design.
The old coverage pass is preserved and marked superseded in `coverage_diagnostic/`.
A fresh full campaign applies the guard from the start. The seventh regression
test proves even an apparently passing partial L cannot enter complete scoring.

The standalone `ph_linkage.py` summary alone is not evidence that both installed
rings were scored: use the installation manifests and `verify_switch_outputs.py`.

## Trimmed-target limitation

The receptor remains the original trimmed input, including three fragments inside
chain A and the eight-residue B fragment. PROPKA reports backbone-group atom-count
warnings at fragment boundaries A683 and A728. No target sequence or backbone is
changed to hide those warnings. The fixed-structure calculations remain ranking
estimates on a trimmed predicted complex; full-context/orthogonal confirmation
is outstanding. Neither PROPKA nor Rosetta scores justify a binding claim.

## Final checkpoint

The audited campaign completed 25/25 backbones, 75 HH/QQ pair attempts. There
are 38 energy-retained HH variants, 36 complete linkage scores, 25 complete scores
with retained QQ controls, and two explicitly unscorable terminal-His variants.
No complete score clears L <= -2; no complete HH pair has both selected free
pKa >= 7. Best L is -1.530. See `campaign/RESULTS.md`.

Seven regression tests passed. The output audit passed on **181 PDBs across
25 backbones**, with unchanged target sequences, exact intended mutations,
verified source/output/code SHA256, valid energy gates and linkage sums. Maximum
backbone atom displacement is **0.000 A**. `campaign/integrity.json` is the audit.
The local CPU campaign is stopped at its configured count. Both Runpod pods
were left stopped. New GPU spend **$0.00**; historical anchor **$33.87**.
