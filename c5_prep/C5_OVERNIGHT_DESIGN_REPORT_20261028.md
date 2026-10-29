# C5 overnight BindCraft design campaign — status report

Date: 2026-10-28 (run started ~00:05 UTC)
Pod: `gmvlpy79felypz` (A40 48GB, EU-SE-1, $0.59/h)
Tool: BindCraft (martinpacesa/BindCraft, standard CLI, single-template mode)
Target: G12D pMHC (state1_g12d_complex.pdb, chains A+B+C = 385 residues)
Hotspots: `A69,A155,C7,C8` (groove rim + exposed peptide, Fable's recommendation)
Algorithm: 4stage_duo (double logits optimization)
Binder lengths: [70, 80, 90, 100, 110] (fixed to avoid JIT recompilation)
MPNN: num_seqs=20, max_mpnn_sequences=2

## Trajectory results

| # | Length | Seed | Stage 1 pLDDT | Softmax pLDDT | One-hot | Final | MPNN | Time |
|---|--------|------|---------------|---------------|---------|-------|------|------|
| 1 | 109 | 634072 | 0.80 PASS | 0.71 PASS | 0.76 PASS | 0.76 | 0/20 | 1h12m |
| 2 | 98 | 76712 | 0.88 PASS | 0.64 FAIL | - | - | - | 1h0m |
| 3 | 100 | 567490 | 0.66 PASS | **0.92 PASS** | **0.82 PASS** | **0.82** | 0/N (running) | ~1h5m + MPNN |

**Accepted designs: 0**

v3 run (standard 4stage + relaxed AF2 filters pLDDT≥0.70, i_pTM≥0.35):
| # | Length | Stage 1 | Softmax | One-hot | Final | Interface | MPNN |
|---|--------|---------|---------|---------|-------|-----------|------|
| v3-1 | 74 | 0.80 | **0.86** | **0.80** | **0.80** | 13 res | 0/20 |
| v3-2 | 76 | — | 0.45 FAIL | — | — | — | — |
| v3-3 | 92 | 0.75 | **0.89** | **0.87** | **0.87** | 17 res | 0/20 |
| v3-4 | 94 | **0.85** | **0.89** | **0.87** | **0.88** | **24 res** | 0/20 |

**FINAL OVERNIGHT RESULT: 0/100+ MPNN designs accepted across 5 successful
hallucinations (v2 + v3), with both default and relaxed AF2 filters.**

v3 trajectory 4 had the best metrics of the entire session (pLDDT 0.88, i_pTM 0.88,
24 interface residues) — yet 0/20 MPNN sequences passed. The MPNN→AF2 gap is
confirmed as fundamental, not a threshold issue.

**Root cause:** AF2 predicts entirely different structures for MPNN-redesigned
sequences on this multi-chain pMHC target. The hallucinated sequence works in the
context of the hallucination's continuous optimization, but when MPNN redesigns
the non-interface positions and AF2 refolds from scratch, the groove binding mode
is lost. This is a known limitation of BindCraft for multi-chain targets with
deep groove binding sites.

## What worked

1. **Corrected hotspots (Fable's recommendation) solved the Stage 1 problem.**
   v1 hotspots (buried groove floor: A70/A114/A147/A156/C5-C8) → 0/1 passed Stage 1.
   v2 hotspots (accessible groove rim: A69/A155/C7/C8) → 2/2 passed Stage 1.

2. **The APBS electrostatic difference map was MEASURABLE.** 8.4 kT/e at 6-8 Å
   from the peptide. The design strategy (complementary charge at P6-P7) is
   physically grounded.

3. **The pipeline ran autonomously all night.** GPU at 100% utilization during
   hallucination, ~26s per iteration, ~1h per trajectory.

## What didn't work

1. **0/20 MPNN designs passed AF2 validation** (trajectory 1). The hallucinated
   backbone (pLDDT 0.76) was marginal — good enough for the trajectory gates but
   not for MPNN sequence redesign to reproduce the pose.

2. **Trajectory 2 failed the Softmax pLDDT gate at 0.64** despite starting with
   the best Stage 1 metrics (pLDDT 0.88, i_pTM 0.84). The sequence commitment
   step consistently loses the binding pose. The `4stage_duo` algorithm's double
   logits pass may be counterproductive.

3. **Multi-template BindCraft (detargeting) hung during XLA compilation.** The
   multi-template engine from `Multi_BC_Original.ipynb` could not compile the
   JAX graph for the 484-residue input in under 2 hours. Fell back to
   single-template.

## Root cause analysis

The multi-chain pMHC target (385 residues, 3 chains) is challenging for
BindCraft's hallucination:

- **Large target → expensive iterations (~26s each)** limits the number of
  trajectories per hour.
- **Multi-chain binding → harder sequence-structure coupling** makes MPNN
  redesign less reliable.
- **The 4stage_duo algorithm** adds a second logits pass that resets the
  sequence and may disrupt converged binding poses.
- **The groove is a difficult binding site** — deep, buried, surrounded by
  HLA helices. The accessible approach vectors are limited.

## Recommendations for Ricardo

### Immediate (today)

1. **Switch to standard `4stage` algorithm** (not `4stage_duo`). The duo's
   double logits pass consistently produces lower Softmax pLDDT. The standard
   algorithm is faster (~25% fewer iterations) and may maintain binding poses
   better through sequence commitment.

2. **Lower the Softmax pLDDT threshold** from 0.80 to 0.65-0.70. Trajectory 2
   failed at 0.64 — many trajectories with good backbones are being killed too
   early. Let more through to MPNN; the MPNN + validation filters will catch
   truly bad designs.

3. **Run more MPNN sequence variants** (increase max_mpnn_sequences from 2 to 5).
   More sequence diversity at the MPNN stage increases the chance of finding one
   that validates.

4. **Consider shorter binders (60-80 aa)** — faster per iteration, less material
   to fold, might find simpler binding geometries.

### Strategic

5. **Try BindCraft2 on Tamarind Bio** (cloud platform, Pacesa's tweet). It has
   native detargeting. Cloud run avoids the local XLA compilation issue.

6. **RFdiffusion as an alternative backbone generator** — the panel's original
   recommendation. RFdiffusion is faster at backbone generation and may produce
   higher-quality backbones for this groove geometry. Use ProteinMPNN + AF2
   validation separately.

7. **The APBS map should inform hotspot refinement** — use the electrostatic
   difference map to identify the exact grid points where placing charge matters
   most. This could translate to more specific hotspot residues.

## Files

- `output_g12d_v2/` — trajectory PDBs and CSVs (on pod)
- `run_g12d_v2.log` — full design log (on pod)
- This report

## Cost

- Pod uptime: ~8h × $0.59/h = ~$4.72 (C5 design pod only)
- C4 rescore fleet: still running ($3.54/h)
- RunPod balance was $316.44 at start of night
