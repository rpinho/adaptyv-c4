# Fable Review — C4 Acid-Site Strategy (2026-10-23)

## Top-line

The mechanism is sound and correctly reasoned. Two steps in the funnel are actively working against us: the **≥2 His filter is on the chain, not the interface** (97% pass rate — it filtered nothing), and the **top-30 selection by *highest* free-state His pKa is pointing the wrong direction**. Both are fixable before Nov 1 and matter more than any threshold tuning.

## 1. Is L ≥ +2 the right gate?

The magnitude is defensible, the use as a hard gate is not.

Physical ceiling: for a single proton-linked site, the maximum achievable Δlog K across pH 7.4→5.5 is 1.9 log units. **L ≥ +2 is mathematically unreachable with one effective His.** You need ≥2 well-placed His, each near-ideal.

L is derived from PROPKA ΔpKa; PROPKA's His error is ~1 pKa unit, so σ(L) over 2-3 His is roughly 1-1.5 log units. Make L a **ranking** signal, and make the hard gates things that are robust to pKa error:
- Sign constraint: ΔpKa > 0 for ≥2 His
- Geometry: each contributing His within 4.0 Å of target Asp/Glu carboxylate in complex
- **Two-sided requirement**: gate absolute affinity at pH 5.5 too (L = +4 with Kd(5.5) = 10 µM fails)

## 2. ipTM vs L weighting

**ipTM as gate, L as rank within the gated set** — not a weighted sum.

L is computed on top of the predicted complex. If the pose is wrong, His–carboxylate geometry is wrong, so ΔpKa is wrong, so L is noise. L is only interpretable conditional on the pose being approximately right.

Suggested: hard gate ipTM ≥ ~0.75, rank by L within that. Reserve 3-5 of 20 slots for high-L / ipTM 0.6-0.75 as a mechanism hedge.

## 3. Free-state pKa clustering at 6.0-6.5 — CRITICAL FINDING

**Good that they cluster there. Our selection was INVERTED.**

For proton uptake on binding: want pKa_free **low** (His neutral at 7.4) and pKa_bound **high** (His protonated at 5.5). Top design at avg free pKa 7.05 is substantially protonated at pH 7.4 in the free state — the switch is largely spent before the target shows up.

260/340 designs sit inside one 0.5-wide bin around the model value. We ranked on 0.3-0.5 unit differences inside PROPKA's error. The top-30 was effectively random w.r.t. the mechanism, with a mild bias toward the wrong end.

**Better pre-screen**: band pKa_free to ~5.5-6.5, require His to be solvent-exposed and on the predicted interface face.

## 4. PROPKA on de novo designs

- His is PROPKA's worst residue class; benchmark RMSD ~0.8-1.0 units
- PROPKA parameterized on experimental natural-protein structures; RFD4 backbones have idealized geometry → systematic over-prediction of shifts
- Sensitive to rotamer and His ring flip/tautomer

Mitigations:
1. **Run FcRn/IgG positive control** (native pH switch in our direction) through the pipeline
2. Run PROPKA across multiple Boltz-2 samples per design
3. Count solvent-accessible His–(Asp/Glu) contacts under 4 Å at interface (geometric check)

## 5. MARGINAL designs?

No. With 343 PASS and 20 slots, no scarcity. Exception: check per-residue pLDDT at interface/His positions, not just mean.

## Additional risks

- **"Must NOT bind at 7.4" half is unscored** — nothing penalizes pH-independent affinity
- **Nonspecific polycation risk** — His-rich binder targeting acidic patch may stick to anything acidic
- **Interface-based His filter needed** — ≥2 His filter on whole chain did no work (97% pass)
- **Funnel is thin** — 400 backbones → 20 submissions for conditional mechanism
- **Own the pH-blindness in write-up** — no step predicts structure at 5.5 vs 7.4

## Actions taken

1. ✅ Launched batch 2 with corrected selection: free pKa 5.5-6.5, high pLDDT
2. ✅ Updated scoring to support ipTM gate + L ranking
3. ⬜ FcRn positive control (if time permits)
4. ⬜ Geometric His–carboxylate contact check after Boltz-2
