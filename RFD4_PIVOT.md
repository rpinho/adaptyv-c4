# Generation pivot to RFD4-Proteina — 2026-10-22

Approved by Ricardo. This replaces the BindCraft-style gradient pipeline as the
primary generator. The old pipeline is not deleted; it is no longer the plan.

## Why

Our acquisition failure is a pathology of gradient descent through AF2, not a
tuning problem. Across all five v14 replicates (250 gradient steps) the binder
never came within 10 A of a hotspot (median 54-57 A) while `interface_contacts`
sat saturated at its 20 A cutoff, contributing *zero* gradient toward the
epitope, and the backbone collapsed to 0.16-0.79 A adjacent-CA spacing. A
diffusion model conditioned on hotspots places the binder at the epitope by
construction and does not have this failure mode.

Authors' reported benchmark: **10.6% clustered pass rate for binder design vs
5.9% for BindCraft, at ~42x throughput.** BindCraft is precisely the stack that
produced our 0-for-5. Wet-lab support reported: 12/18 calmodulin binders
<1 uM; 158/188 sub-uM TrkA; 112 pM YAP1 and 527 pM Tau disordered-peptide
binders. These figures are from the authors via launch coverage and are NOT yet
confirmed against a primary paper — treat as motivation, not evidence.

## Licensing — resolved, and the reason we run locally

| Route | License / terms | Verdict |
| --- | --- | --- |
| **Local weights** | **BSD-3-Clause.** Repo `LICENSE` and the Hugging Face `LICENSE` are byte-identical (1,533 bytes), Copyright (c) 2026 Institute for Protein Design, University of Washington. No non-commercial clause, no gating. | **Use this** |
| NVIDIA NIM free endpoint | NVIDIA API Trial ToS: §2.2 grants NVIDIA "a worldwide license with respect to User Content and Generated Content, to use, host, store, reproduce, modify, create derivative works"; §3.3 to use it "to improve NVIDIA products and services, including AI models"; §1.4 limits use to "internal testing and evaluation purposes, not in production". | Dry runs only, never the entry |
| Hugging Face Spaces | 10 generations/hour hard cap; inputs "may be retained for monitoring and improvement". | Sanity checks only |
| subseq.bio | Third party, unaffiliated with IPD or NVIDIA. Permissive terms (you own outputs, commercial use allowed), ~$1.99/GPU-hr. But no biosafety screening is mentioned, no company identity is published, and the site shows a "DEV MODE" badge. | Avoid for the entry |

NVIDIA's model page cites a "NVIDIA Community Model License" in one place and
BSD-3 in another; first-party file reads say BSD-3. The Community license
attaches to NVIDIA's hosted deployment, not to the downloaded weights.

## Sources

- Model: `github.com/RosettaCommons/RFD4-Proteina` — public, BSD-3-Clause, 331
  stars, last pushed 2026-10-22. Ships `.project-root`, `pyproject.toml`,
  `configs/inference.yaml`, `setup.sh`, `Makefile`, and its own `CLAUDE.md`,
  `AGENTS.md` and `skills/`.
- Weights: `huggingface.co/nvidia/RFD4-Proteina`
- Agent skills: `github.com/NVIDIA-BioNeMo/bionemo-agent-toolkit`, under
  `open-models-skills/rfd4-proteina/` (nine skills) and
  `nim-skills/rfd4-proteina-nim/`. Skill text CC-BY-4.0, helpers BSD-3-Clause.

The nine open-model skills: `setup`, `annotated-cif`, `generate`, `search`,
`evaluate`, `autoresearch`, `commitment-curves`, `campaign`, `package`.
Recommended order: **setup -> annotated-cif -> generate -> optional search ->
evaluate**.

## How it runs

Hydra CLI. Arguments are `key=value`, **never** `--key value`. Run from the
repository root (the directory containing `.project-root`).

```bash
# install
cd "$RFPROTEINA_ROOT" && make setup && source .venv/bin/activate
# env only, no weights:  bash setup.sh --skip-weights
# checkpoints land in DATA_PATH (default ./data_path); export DATA_PATH first

# readiness
bash <skill>/rfd4-proteina-setup/scripts/preflight.sh --repo-root "$RFPROTEINA_ROOT"
# must print PREFLIGHT PASSED before anything else

# smoke test
rfd4-generate inputs=data/benchmarks/monomer \
    subset_to_examples=monomer-short num_replicates=1 out_dir=./runs/smoke

# protein binder, the shape we need
rfd4-generate \
    inputs=data/benchmarks/ppi/ppi_all_atom/pdl1.cif.zst \
    num_replicates=8 batch_size=8 out_dir=./runs/ppi \
    metrics_tags='[ppi_refolding,hotspot_conditioning]'
```

Parameters that matter: `out_dir`, `num_replicates` (8-80; 1 is a smoke test),
`batch_size` (start 8, but confirm the longest input fits at `batch_size=1`
first), `metrics_tags`. Extras: `seed=<int>`, `skip_existing=true` to resume,
`save_chain_of_interest=true` **required** if we want foldseek clustering of
designed binders — which we do, for the diversity axis.

Diversity knobs (backbone and sidechain/sequence temperature):
`generation.model.bb_ca.simulation_step_params.sc_scale_noise` (default 0.3) and
`generation.model.local_latents.simulation_step_params.sc_scale_noise` (default
0.1). Roughly (0.1,0.1) = high designability, (0.5,0.5) = UniProt-like,
(0.9,0.9) = maximum diversity. `recycle_strategy=thinking_high` gives the 1.5x
recycle-depth variant. Default 1x schedule is 60 model evaluations / 185 trunk
passes.

## Hotspot conditioning — how the epitope is specified

Hotspots are **not** a CLI string. They are the `C_HOT` annotation on target
atoms inside an annotated CIF (AtomWorks `StandardAnnotation` columns). So the
target input must be authored, not just converted:

1. Convert the TfR1 target to CIF. **PDB is rejected** by the local package and
   by all three hosted routes.
2. Author `C_HOT` on the epitope atoms, plus coordinate/sequence conditioning on
   the target, via the `annotated-cif` skill's Python recipes
   (`references/python-authoring.md`, `annotators-and-selections.md`).
3. Validate: `scripts/inspect_cif.sh <file>` then
   `scripts/check_annotated_cif.sh <file> --strict`.
4. Pass it with `inputs=<file>.cif`.

The designed binder chain is the "chain of interest" and is identified
automatically from the conditioning annotations.

Our case: chain A residues **189-383** (the 195-residue apical domain) inside a
639-residue receptor numbered A122-760, hotspots **A205, A206, A208, A209,
A210**. Whether the rest of the receptor can be present-but-untargeted is the
open authoring question — Fable's domain-crop argument composes naturally here,
since conditioning is per-atom rather than per-file.

## Unknowns before a campaign

- **VRAM and runtime per design: NOT STATED** in any doc read. Must be measured.
- Checkpoint download size: not stated. `configs/inference.yaml` names the pair.
- Whether `*_refolding` metrics are worth it: every `*_refolding` tag runs the
  **AF3 oracle**, which needs an apptainer/docker container. Plain generation
  needs no container. Do not add refolding to the first smoke test.
- Whether the local path performs or requires any biosafety screening, and
  whether it runs fully offline.

## Plan

1. Provision one NVIDIA GPU host (CUDA 12/13, Python 3.12). Quote the hourly
   price before creating it.
2. Clone, `make setup`, download checkpoints, `preflight.sh` to PASSED.
3. Monomer smoke test. Then the shipped `ppi_all_atom` binder example
   **unchanged**, to establish a working baseline and measure VRAM and runtime
   per design on a known-good input.
4. Only then author our own TfR1 annotated CIF and run a bounded campaign.
5. Keep `scripts/build_submission.py` and `scripts/rank_candidates.py` as-is:
   they take sequences and are generator-agnostic by design, which is exactly
   why the emitter was kept free of qualification logic.

Step 3 before step 4 is deliberate. Running the shipped example first separates
"the model works here" from "our input is authored correctly" — the two failure
modes that would otherwise be indistinguishable.
