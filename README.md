# adaptyv-c4

Track 3 entry for **Challenge 4** of the
[Anthropic × Adaptyv Protein Design Competition](https://proteinbase.com/competitions/anthropic-adaptyv-2026).

**Target:** human transferrin receptor 1 (TfR1, UniProt P02786) — a pH-sensitive
conditional binder for blood–brain barrier transcytosis.

**Design criteria** (from the challenge specification):
- Bind at pH 7.4, no detectable binding at pH 6.0
- Mouse TfR1 (Q62351) cross-reactivity preferred
- Single-chain, 10–250 amino acids

## Method

Proton-coupled two-state binder design using RFdiffusion-Proteina
with histidine-enriched sequence conditioning, scored by Wyman linkage
analysis (PROPKA3 pKa shifts across pH). Candidates are screened through
a multi-stage pipeline: geometric metrics at generation time, structural
novelty via Foldseek/TM-align against PDB + patent databases, restored-context
clash audit against the full receptor, and mouse epitope conservation analysis.

See [METHOD.md](METHOD.md) for the full technical description with calibration
tables and scoring details.

## Running

```bash
python -m venv .venv && source .venv/bin/activate
pip install -e ".[all-cuda13]"  # or .[all-cpu] for local work

# CPU dry run (validates the spec, costs cents)
python modal_rfd4.py dry_run

# Generate 200 designs with inline geometric metrics
modal run modal_rfd4.py::generate --num-replicates 200 --batch-size 8 \
    --out-name tfr1_run --metrics geom

# CPU novelty screen
python scripts/run_cpu_novelty.py input.csv --databases databases.json \
    --structures structures.json
```

Requires a [Modal](https://modal.com) account for GPU generation and a local
install of [Foldseek](https://github.com/steineggerlab/foldseek) for structural
novelty screening.

## Disclaimer

This is an independent entry. It is not affiliated with, endorsed by, or
sponsored by Adaptyv Bio, Anthropic, or any competition organiser. All designs
are original work by [Ricardo Pinho](https://github.com/rpinho).

## License

Code in this repository is provided for reference. Designed sequences submitted
to the competition are governed by the
[competition terms](https://proteinbase.com/competitions/anthropic-adaptyv-2026).
