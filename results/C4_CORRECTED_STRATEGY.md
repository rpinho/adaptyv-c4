# C4 Corrected Strategy — Post Sign-Inversion Discovery

## Summary

On 2026-10-23, Fable consultation identified that our ENTIRE C4 scoring and generation
pipeline was optimized for C3's endosomal-release mechanism (bind at pH 7.4, release at
pH 5.5), not C4's conditional-binding objective (bind at pH 5.5, no binding at pH 7.4).

**All 20 v8 designs have net WRONG-direction pH-switch.** 1/20 (his015_96) has a single
interface His with weak correct-direction linkage (+0.293 kcal/mol). 0/20 reach the
actimo-labs benchmark of +2.84 kcal/mol.

## C4 Physics (Correct Direction)

**Objective**: Binder must take up protons on binding (Wyman linkage L ≥ +2).

**Mechanism**: Binder His placed near target Asp/Glu (carboxylates):
- Free state: His in solvent, pKa ~6.5 (normal)
- Bound state: His near target carboxylate → ionic salt bridge stabilizes protonated
  form → pKa RISES above 7.4 → His is protonated at 5.5 AND 7.4 when bound...
  
  Wait — this needs more careful analysis. Let me reconsider.

Actually the correct mechanism for a pH 5.5 binder is:
- At pH 5.5: His is protonated → forms salt bridge with target Asp/Glu → binding ON
- At pH 7.4: His is deprotonated → no salt bridge → binding OFF
- This works when His pKa in the complex ≈ 6.0-6.5 (between the two pH values)

But pKa_bound > pKa_free just means the bound state stabilizes protonation more.
The key for selectivity is that the His is PROTONATED at pH 5.5 (always true for
reasonable pKa) and DEPROTONATED at pH 7.4 (true when pKa < 7.4).

The Wyman linkage captures this: if binding requires protons (Δν > 0, proton uptake),
then binding is more favorable at lower pH. L > 0 means stronger binding at pH_low.

For the switch to work, we need His pKa_bound in the 5.5–7.4 window, and ideally
pKa_bound > pKa_free (because raised pKa means the salt bridge exists over more of
the pH range).

## TfR1 Apical Domain Acidic Epitope Analysis

TfR1 apical (154 res): 3D + 7E = 10 acidic, 7R + 7K + 5H = 19 basic.
The surface is NET BASIC — challenging for proton-uptake design.

### Acidic residues confirmed at interface (from Boltz-2 complexes):

| Residue | In his015_96 complex | In his015_115 complex | Notes |
|---------|---------------------|----------------------|-------|
| E36 | 7.1Å (interface) | 22.7Å (distal) | Best in 96, not in 115's docking site |
| D41 | 2.8Å (contact!) | 28.8Å (distal) | Direct contact in 96! |
| E98 | 5.0Å (interface) | 6.6Å (interface) | Interface in both |
| D17 | 14.3Å | 5.6Å (interface) | Interface in 115 |
| E132 | 29.8Å | 7.3Å (interface) | Interface in 115 |

### Two distinct docking sites:
1. **E36-D41-E98 cluster** (his015_96 site): has the best acidic contacts but also
   R84 (2.8Å), H88 (3.0Å), K97 (6.4Å) nearby
2. **D17-E98-E132 region** (his015_115 site): more dispersed, less basic contamination

## Regeneration Plan (9 days to Nov 1)

### Phase 1: Calibration (Day 1 — immediate)
- [x] Fix scoring scripts (complex_propka_switch.py, ph_linkage.py)
- [ ] Run actimo-labs C2 sequences through our pipeline to verify methodology
  - Their ΔΔG should come out NEGATIVE in our C4 convention (they bind at pH 7.4)
  - This validates our pipeline without needing positive C4 controls
- [ ] Scramble control: predict 5 scrambled binder sequences → calibrate ipTM noise

### Phase 2: Epitope Selection (Day 1-2)
- [ ] Download TfR1 crystal structure (PDB 1CX8 or 3S9L) for true surface accessibility
- [ ] Identify best acidic patch with BioPython SASA analysis
- [ ] Design RFdiffusion hotspots targeting D41+E36+E98 (site 1) or D17+E98+E132 (site 2)
- [ ] Consider both sites — diversify across submission

### Phase 3: Generation (Day 2-5)
- [ ] RFdiffusion: 200+ backbones targeting each acidic site
- [ ] ProteinMPNN: His-biased pool (sparse, 2-3 His per design)
- [ ] Selection: only accept His within 5Å of target D/E, >8Å from target R/K
- [ ] ESMFold refolding + quality gate

### Phase 4: Scoring (Day 5-7)
- [ ] Boltz-2 complex prediction (corrected pipeline)
- [ ] Complex PROPKA with CORRECTED sign
- [ ] Gate: L ≥ +2 (corrected direction)
- [ ] Multi-seed validation on top candidates

### Phase 5: Submission (Day 7-9)
- [ ] Replace wrong-direction designs with corrected ones
- [ ] Keep his015_96 as positive-direction reference
- [ ] Include 1-2 deliberate L ≈ 0 negative controls
- [ ] Final v9 submission CSV

## Key Design Constraints for Corrected Binders

1. **Sparse His**: 2-3 per design (actimo-labs pattern, confirmed by thermodynamics)
2. **His near target D/E**: within 5Å of Asp/Glu side chains
3. **His far from target R/K**: >8Å from Arg/Lys/His (avoids pKa depression)
4. **His >8Å apart**: prevents anti-cooperative coupling
5. **No carboxylate near binder His in free state**: maintains high free-state pKa
6. **pI < 8 at pH 7.4**: prevents non-specific electrostatic binding

## What We Keep From v8

- **his015_96**: only correct-direction design (L = +0.21, weak but right way)
- **Boltz-2 infrastructure**: Modal scripts, cache volumes all reusable
- **PROPKA pipeline**: fixed scripts ready to score new designs
- **Multi-seed validation**: methodology proven, apply to new designs
- **All ipTM data**: useful as binding signal even if noisy
