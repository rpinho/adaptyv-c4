# Fable Consultation: PROPKA + Spatial Analysis Review
Date: 2026-10-23

## Critical Finding: ΔΔG_monomer ≈ His Count

Regressing −ΔΔG_monomer on His count: **R² = 0.80, slope = 1.19 kcal/mol per His**.
The ranking is a histidine counter, not a pH-switch predictor.

Mean PROPKA shift from 6.5 model: -0.25 ± 0.51, 93% within |1.0|.
PROPKA accuracy on crystal structures is ~0.8-1.0 pKa RMSD; these are unrelaxed ESMFold.
Per-His pKa values are model value + noise.

## Monomer ΔΔG is the Wrong Quantity

pH-dependent binding depends on **ΔpKa = pKa_bound - pKa_free**, not on free pKa:
```
ΔG_bind(pH) = ΔG_ref - RT Σ ln[(1+10^(pKa_b,i-pH)) / (1+10^(pKa_f,i-pH))]
```
A His with pKa 6.5 free AND 6.5 bound contributes **zero** to the switch.
Our metric and actimo-labs' metric are DIFFERENT quantities on different structures.

## High-His Density is Mechanistically Wrong

1. Can't help by construction — only interface His matter
2. Charge coupling suppresses the switch (PROPKA under-counts at this density)
3. +6 to +8 net charge at pH 5.5 destabilizes fold and invites non-specific binding
4. ESMFold pLDDT at 95+ on 15.8% His is not reassuring (outside training distribution)

## c4_v_185_4x Has pI = 5.5

At isoelectric point at working pH = minimum solubility = fail on expression.
Drop it.

## Replacement: his015_118 > his015_176

his015_118 (66 res, 9 His, Δq=7.9, pI=6.6, 2 HH doublets) beats
his015_176 (121 res, 13 His, Δq=11.1, pI=6.8, 4 HH doublets).
Shorter, less supercharged, fewer coupling sites.

## His-D/E Analysis Too Loose

6Å Cα cutoff is not a salt bridge. Use 4.0Å heavy-atom Nδ1/Nε2→Oδ/Oε.
Prediction: 40-70% collapses to <15%.

## Action Plan (Fable's Recommendation)

1. Days 1-2: Boltz-2 complex prediction (ALREADY RUNNING)
2. Day 3: PROPKA on complex+monomer → true ΔΔG_switch at interface His
3. Day 3 (parallel): redo His-D/E at 4.0Å
4. Days 4-6: swap designs scoring 0 or negative on real metric
5. Days 7-8: charge-cap designs with q(5.5) > +6 via non-interface H→Q/N
6. Day 9: buffer

## Bottom Line

Pipeline is sound engineering pointed at the wrong observable.
Fix = complex structure prediction (Boltz-2) → real ΔΔG_switch.
This is fixable in 9 days with existing infrastructure.
