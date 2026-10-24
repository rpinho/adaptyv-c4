# Fable Consultation: Boltz-2 ipTM Results Review
Date: 2026-10-23

## Early Results (3/20 at time of consultation)

| Design | Size | ipTM | pTM | Confidence |
|--------|------|------|-----|------------|
| his200m | 78 | 0.741 | 0.834 | 0.757 |
| his015_96 | 107 | 0.422 | 0.622 | 0.640 |
| his015_115 | 74 | **0.895** | 0.927 | 0.747 |

## Fable's Assessment

### Is ipTM=0.895 suspicious?
YES — exercise caution. De novo binders rarely exceed 0.7 against novel targets without:
- (a) Structural homology to known binder in MSA
- (b) Boltz-2 hallucinating a confident but incorrect interface

TfR1 is well-represented in training set (HFE, transferrin). Boltz-2 may over-interpolate.

**Action items:**
- BLAST his015_115 and his015_148 against PDB to check structural mimicry
- Run with --diffusion_samples 5 — if ipTM variance >0.15, score is noise
- A true 0.895 on 74-res de novo binder would be publication-grade

### ipTM vs Confidence Discrepancy
Moderate confidence (0.747) despite high ipTM (0.895) = **pLDDT dragging it down**
- Likely poor local confidence in loop regions or interface
- Check per-residue pLDDT at binding interface
- If interface residues <70 pLDDT, high ipTM is unreliable

### ipTM Threshold
Use **≥0.35** for PROPKA analysis
- Captures everything above shuffled-control noise (~0.235)
- Don't over-filter — need full distribution to calibrate
- Designs below 0.35 still worth keeping if good pH-switch metrics

### Compact Designs Scoring Higher
Expected — NOT concerning
- Smaller binders have fewer degrees of freedom
- Boltz-2 converges more reliably on compact designs
- Larger designs may need more diffusion samples

## Bottom Line
his015_115 is lead candidate pending BLAST check and multi-seed validation.
Analyze everything ≥0.35 with complex PROPKA.
