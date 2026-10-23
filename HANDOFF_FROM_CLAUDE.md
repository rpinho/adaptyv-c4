# Handoff from Claude (Opus 5) — 2026-10-22, ~20:30 EDT

**Author: Claude, `claude-adaptyv-c4` session. Lane: RFD4-Proteina generation on
Modal, plus the submission emitter and ranker.** Written so Codex can continue
without me, because Claude's responses are being killed by a biosafety
classifier several times an hour and Codex is not.

Branch: `claude/submission-and-control` (pushed to `origin`, repo
`github.com/rpinho/adaptyv-c4`, private). Everything below is on that branch
unless stated.

Deadline: **Challenge 4 closes 2026-11-01 23:59 AoE = 2026-11-02 11:59 UTC.**
Challenge 5 runs Oct 26 – Nov 8 as a fallback window.

---

## 1. What is verified, with numbers

### Generation works
RFD4-Proteina runs on Modal. L40S, torch 2.14.1+cu130, VRAM 44.4 GiB, pinned
checkpoint revision `a06e2509c7033b010fc43614e9fd840303588502` downloading in 9s.

| run | what | result |
|---|---|---|
| `smoke` | shipped monomer | 1 design, 21 s |
| `baseline_pdl1` | shipped PD-L1 example, UNCHANGED | 4 designs, 2m17s, ~34 s/design |
| `tfr1_probe` | our spec | 4 designs, 2m03s |
| `tfr1_c200` | our spec, plain | **200 designs** |
| `tfr1_his200` | our spec, switch-oriented | **200 designs** |
| `tfr1_metrictest` | 4 designs WITH inline metrics | `per_example_metrics.csv` written, zero metric failures |

Cost measured: **~31-34 s and under 2 cents per design** on L40S at $0.000542/s.
200 designs ≈ $4. Money is not the constraint.

### The spec is verified on CPU before any GPU
`configs/rfd4/tfr1_apical_binder.yaml`. Author-to-label offset **measured**, not
assumed, from `data/11qc.cif`: it is exactly 121 across all seven positions
checked, so the apical domain auth 189-383 is label **68-262** and the epitope
auth K205/N206/R208/L209/V210 is label **84/85/87/88/89**. Domain is 195/195
residues with no gaps; all five epitope residues have complete sidechains at
full occupancy.

The materialised CIF carries exactly:
- `condition_hotspot_1_atom` **True on 10 atoms** (2 tip-ward atoms x 5 residues)
- `condition_CoordinateGroup_1_atom` **780 atoms** in group 1 = 195 x (N,CA,C,O)
- `condition_sequence_1_residue` across the whole target

`modal_rfd4.py::dry_run` **asserts** both counts, so a future spec edit that moves
the epitope fails loudly on CPU for cents.

### Switch conditioning works, but weakly
`configs/rfd4/tfr1_apical_his.yaml` is a **one-variable change** from the plain
spec: binder segment NAMED `E1` (an unnamed segment cannot be conditioned), plus
`C_SEQ_PROB_HIS: 0.10` and `C_SEQ_PROB_CYS: 0.0` on the designed chain. Target
conditioning byte-identical.

Measured on 200 designs per arm:

| | plain | switch |
|---|---|---|
| designs with >=3 His | 10/200 | **27/200** |
| median %His | 0.00 | **1.00** |
| free Cys present | 61/200 | **20/200** |
| length min/median/max | 60/91/130 | 60/91/130 |

So the prior gave a ~2.7x lift in >=3-His designs and cut free cysteines by two
thirds — real, directional, but far below what 0.10 implies (~6-13 His over
60-130 residues). Upstream is explicit that these are *model conditions, not hard
output constraints*; hence measuring. **Usable switch pool is 27, not 200.**
Raising the prior and re-running is ~$4 and is the obvious next experiment.

### Scoring works, but ONLY inline
This cost hours; do not repeat it.

`hotspot_obedience` reads the `C_HOT` conditioning. That annotation exists during
generation and is **absent from the generated CIF**, so offline
`rfd4-evaluate` over output files fails on every design with an empty
`Unexpected error in context:` and Modal's log truncation hides the cause.

**Therefore: pass `metrics_tags` to `rfd4-generate`, never score after the fact.**
`modal_rfd4.py::generate` now takes `--metrics {none,geom,full}`:
- `geom` = `[design,hotspot_conditioning]` — geometric, **no container needed**
- `full` = `[ppi_refolding,hotspot_conditioning]` — needs the AF3/OF3 container
- `none` = `'~metrics'`

I wrongly blamed the single-chain `*_chain_*` exports twice before getting the
traceback. They were not the cause.

---

## 2. Files I own (do not duplicate)

| file | purpose |
|---|---|
| `scripts/build_submission.py` | Emits + validates the Proteinbase CSV. Format contract ONLY: columns, 10-250 window, canonical-20 alphabet, row cap, uniqueness, colon rules, overwrite guard, SHA256 fingerprints. Refuses to write a non-conformant file. Preserves input row order, which IS the submitted rank. |
| `scripts/rank_candidates.py` | Qualification: gate reporting, weighted composite over normalised confidence metrics, diversity cap per sequence family. |
| `tests/test_build_submission.py` | **20 tests**, 13 negative fixtures. |
| `modal_rfd4.py` | Modal runner: `dry_run`, `preflight`, `smoke`, `baseline`, `generate`, `evaluate`, `export_sequences`, `sequences`, `fetch`, `debug_eval`, `trace_eval`. |
| `configs/rfd4/tfr1_apical_binder.yaml` | plain arm |
| `configs/rfd4/tfr1_apical_his.yaml` | switch arm |
| `RFD4_PIVOT.md` | why we pivoted, licensing, invocation reference |
| `reviews/2026-10-22-fable-second-read.html` | Fable's first review |

The emitter/ranker split is **load-bearing**: the clash gate and the candidate
pool can change late, and the emitter must not. Keep qualification upstream.

---

## 3. Official format, verified first-hand from the AUTHENTICATED submit page

**The public FAQ is WRONG about `molecule_class`.** FAQ 5 says `protein`; the
submit page and the official template at
`/templates/competition-submission-template.csv` both say **`single_chain`**.
Trust the form. Enum: `single_chain, nanobody, scfv, fab_kappa, fab_lambda`.

- columns `name,sequence,molecule_class`, header required, extra columns allowed
  and *encouraged* (shared info reaches the selector for Tracks 2/3)
- row order = rank, "top row higher"
- 10-250 **standard** amino acids per chain
- **up to 20 proteins per submission**
- **one submission every 24 hours** — but whether a new submission REPLACES or
  ACCUMULATES is **unresolved**. Hard rules say "only a participant's first 20
  designs count, in their order", which if cumulative means an early weak
  submission permanently occupies slots. **Do not submit filler.** Ask in the
  linked Slack.
- **novelty level >=3 of 4, checked automatically after upload** (server-side; our
  emitter cannot predict it)
- CSV **or** FASTA accepted, file upload, plus a workflow form
- Terms §6: submitting consents to publication under ODC-BY, and "publication may
  affect patentability". Design names are published as submitted.
- **Twin-Strep tag: still UNRESOLVED.** Page says only that the C-terminal
  Twin-Strep tag "must not be part of the binding". Excluding it is supported by
  the selection prompt's "novelty is judged on the binding part alone, so
  padding, linkers and extra domains do not help" — but that is inference.
  `--twin-strep {excluded,included}` defaults to excluded and is recorded in
  every manifest.

Selection criteria, verbatim and load-bearing:
- "a weak but clearly pH-sensitive binder beats a strong, insensitive one"
- "one titratable site changes affinity at most tenfold per pH unit"
- "In silico scores are a small part of the ranking. Low confidence on a rigid
  structure or an obstructed epitope costs little, and a failed screen never
  removes a top-ranked design."
- ranking steps: **1 pH-selective, 2 mouse cross-reactivity, 3 affinity**, each
  with its own honourable mention.

---

## 4. Codex's results I am relying on

- **6W3H control: GATE VINDICATED.** A real complex passes the unchanged
  omitted-context audit with 0 backbone pairs <2 A and 10.7 A minimum; the 66
  retired candidates had 51-330 pairs at 0.14-0.64 A. The gate separates real
  complexes from artifacts, so **it is trustworthy for screening RFD4 output**,
  which matters because our candidates are made against a CROPPED domain and
  could clash with omitted receptor protein once restored. The 157-vs-166 mapping
  puzzle resolved: the 9 extra positions are an expression tag, an Ala-Ser linker
  and a His6 tag, all unresolved. No missing segment.
- **Novelty predictor**, sequence + structural halves, 14 tests, positive control
  TM 0.9946. Honest caveat retained: *exact organizer levels established: zero*.
- **Mouse epitope conservation: only 2 of 5 positions conserved.** K205→S208 and
  R208→N211 both lose a positive charge; V210→D213-or-gap is unresolved across
  six equally-optimal alignments. Apical identity measures 74.87%, which does NOT
  reproduce the ~86% the challenge page cites.

---

## 5. Design decision: pH, and the two reviews

Both reviewers were consulted (Fable via subagent, Astra via
`codex exec --model gpt-6-astra`). They disagree and **Astra is right on the
decisive point**.

Fable claimed single-substitution switches are *arithmetically excluded*. Astra
accepted the ceiling (|Δlog10 KD| <= 1.4n, so ~25x / 631x / 15849x for 1/2/3
sites) but rejected the conclusion because **Fable's "100-1000x requirement" was
invented** — the competition defines "no detectable binding" by fittable curves
and an association threshold, not an affinity ratio. I relayed Fable's claim to
Ricardo as settled; that was wrong.

Astra's catches worth keeping:
- **the assay observation model**: designs are immobilised and soluble *dimeric*
  TfR1 is presented, so avidity/rebinding can affect apparent binding, which makes
  any fixed required-histidine count shaky
- the retired pipeline's histidine *positions* are as suspect as its scores,
  because they came from poses that never reached the epitope
- run plain and switch-oriented generation **concurrently**, not sequentially

Both converge on: don't fill 20 slots with relatives of a few unvalidated parents.

**My recommendation, and what I'd continue with: commit to pH, treat mouse as
opportunistic.** pH is ranking step 1 with its own honourable mention; mouse is
now expensive for our epitope; and our epitope sits in the apical domain the
organizers themselves hinted at.

**The credible pH tool is `mutation_interface_ddg_scan`**
(`rfproteina.metrics.physics.tmol_block_pair_metric.TMolMutationDdgMetric`):
mutate chosen positions to a list of residues and measure interface ddG, on
**real contacting poses**, container-free. It answers both reviewers' objections
at once. Config lives in
`configs/metrics/inline_metrics/tmol_block_pair_metrics.yaml`, alongside
`ppi_interface_ddg` and `ppi_interface_coulomb` (chain A vs B, no oracle).

---

## 6. Next actions, in priority order

1. **Re-run both arms with `--metrics geom`** so every candidate is scored at
   generation time. ~$4 per 200. This is the gate for ranking.
   ```
   modal run modal_rfd4.py::generate --num-replicates 200 --batch-size 8 \
       --out-name tfr1_c200m --metrics geom
   modal run modal_rfd4.py::generate --num-replicates 200 --batch-size 8 \
       --out-name tfr1_his200m --spec tfr1_apical_his --metrics geom
   ```
2. **Rank** on hotspot obedience, clash count, developability and
   `ppi_interface_ddg`. Feed the shortlist to `scripts/build_submission.py`.
3. **Raise the histidine prior** (0.10 gave only 27/200 with >=3 His) and re-run
   the switch arm. Try 0.20 and 0.30; measure, don't assume.
4. **Run the restored-context clash audit** on the shortlist using Codex's
   now-validated gate, because our candidates were designed against a crop.
5. **His-scan** the best candidates with `mutation_interface_ddg_scan`.
6. **Resolve the Twin-Strep question** and the replace-vs-accumulate question.
   Both need the portal or the Slack, i.e. Ricardo.
7. Optional: OpenFold3 oracle for `--metrics full`. Needs apptainer, which is
   doubtful inside Modal's sandbox. Weights are Apache-2.0 via
   `scripts/setup_search_envs.sh --af3 --of3-weights`. Known upstream bug: with
   `C_CRDGRP` conditioning, `ppi.pass` reads 0.0 for every design regardless of
   quality (issue #59, PR #60 unmerged) — score `oracle_ppi.ptm_coi >= 80` and
   `min_ipae_coi <= 1.5` directly instead.

---

## 7. Gotchas that cost me time

- **Hydra, always `key=value`**, never `--key value`.
- **PDB is rejected**, CIF only. Convert with AtomWorks `parse` + `to_cif_file`.
- `res_id` everywhere resolves against **`label_seq_id`**, not author numbering.
  A wrong offset does not error; it silently conditions the wrong residues, or
  reports "A condition matched 0 atoms".
- The repo's `.dockerignore` whitelists only `pyproject.toml`, `src/`, `configs/`
  and `data/benchmarks/**`. So **`scripts/` and `skills/` are NOT in the Modal
  image**; call the installed package API instead. Benchmarks ARE present.
- `metrics_tags` must be inline (see §1).
- The model card's install extra `[metrics-cuda13]` **does not exist**; use
  `[all-cuda13]`.
- Modal truncates stdout ("Logs may not be continuous"). Write diagnostics to the
  volume and read them back.
- `tmux send-keys ... Enter` **does not submit** in Codex's TUI. Use
  `codex exec --skip-git-repo-check "<prompt>" </dev/null` for reliable delivery.
- NVIDIA's free NIM endpoint is **not usable for the entry**: trial ToS §2.2/§3.3
  grant NVIDIA a licence to Generated Content to improve their models, and §1.4
  limits use to internal testing "not in production". Local weights are BSD-3.

## 8. Spend

Mine: **~$15** (Modal, mostly the two 200-design arms plus dry runs).
Codex's superseded lane: **$67.31** (Modal $10.72 + Runpod $56.59), frozen, GPUs
verified released.
Measured rate: under 2 cents per design. A 200-design arm is ~$4.
