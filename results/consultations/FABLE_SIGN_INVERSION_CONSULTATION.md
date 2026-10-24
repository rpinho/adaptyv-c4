# Fable Consultation: ΔΔG_switch Sign Inversion (2026-10-23)

## Critical Finding: ENTIRE C4 PIPELINE IS SIGN-INVERTED

### The Problem

Challenge 4 requires: **bind at pH 5.5, no binding at pH 7.4** (proton UPTAKE on binding).

Our pipeline was designed for C3 (endosomal release: bind at pH 7.4, release at pH 5.5):
- `switch_sites.py`: selects positions for pKa DEPRESSION on binding (C3 objective)
- `ph_linkage.py`: gates at L ≤ -2 (weaker binding at low pH — C3)
- `complex_propka_switch.py`: ranks descending by ddg_switch (positive = pKa depression = C3)

For C4, the correct direction is:
- His pKa RAISED on binding (proton uptake — nearby target carboxylates)
- Gate: L ≥ +2 (stronger binding at low pH)
- Ranking: most NEGATIVE ddg_switch (or flip sign)

### Evidence

**his015_131** (our "best" at ΔΔG=+5.703):
- All 8 interface His have negative ΔpKa (pKa depressed on binding)
- H30: pKa drops from 6.82 → 4.22 (−2.60) — strongest contributor
- This is 15,000× selectivity FOR pH 7.4 binding (WRONG for C4)

**his015_96** (our "worst" at ΔΔG=−0.293):
- H43: pKa RISES from 6.85 → 7.17 (+0.32) — proton UPTAKE
- This is the ONLY design with net correct-direction linkage
- But only 2× selectivity — far below actimo-labs minimum of 100×

### Corrected Picture
- 1/20 designs has net correct direction (his015_96, corrected +0.293)
- 0/20 designs reach actimo-labs benchmark (+2.84 to +4.56)
- The generation objective in switch_sites.py is inverted for C4

### Root Cause
- `switch_sites.py` explicitly states: "a switch position must be elevated-pKa when free AND depressed when bound"
- This is the C3/endosomal-release mechanism, carried into C4 unchanged
- The epitope selection targets cation-rich patches on TfR1, which DEPRESS bound-state pKa
- For C4, need to target ACIDIC patches (Asp/Glu) that RAISE bound-state pKa

### Required Fixes (Fable recommendations)
1. Flip switch_sites.py: reward target carboxylates, penalize cations
2. Change gate in ph_linkage.py from L ≤ -2 to L ≥ +2
3. Fix ddg_switch ranking direction (or negate the score)
4. Re-target epitope to acidic patch on TfR1
5. Regenerate with sparse His (2-3) placed against target carboxylates
6. Calibrate: run actimo-labs sequences through corrected pipeline → should get L ≈ +2.1 to +3.3

### Additional Code Bugs Found
1. `complex_propka_switch.py:~170`: `interface_his = find_interface_his(...) or list(...)` — empty list is falsy, silently falls back to all His
2. `pka_b = binder_his_bound.get(resnum, 6.5)`: default 6.5 for missing His can manufacture spurious ΔpKa
3. `find_his707_contacts`: chain-B numbering runs 1-154, no His707 — mislabeled but harmless
4. pH window: switch_sites.py uses 6.0-7.4 but C4 is 5.5-7.4 (1.9 units, not 1.4)

### Fable's Recommended 9-Day Plan
1. Confirm sign (minutes) — DONE
2. Calibrate against actimo-labs sequences (hours)
3. Re-pick epitope to acidic TfR1 patch + regenerate (bulk of compute)
4. Scramble control: predict scrambled sequences to test ipTM bias hypothesis
5. Multi-seed the ΔΔG, not just ipTM
6. Keep one L ≈ 0 design as internal negative control

### Strategic Framing
Finding an inverted objective before submission is a BETTER Track 3 result than
a clean ranking built on it. The sign error is genuinely easy to make — both
directions are "pH switches," the C3 machinery was correct for C3.
