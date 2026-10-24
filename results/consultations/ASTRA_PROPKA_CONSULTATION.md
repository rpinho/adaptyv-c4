# Astra Consultation: PROPKA + Strategy Review
Date: 2026-10-23

## Critical Finding: No Evidence of Binding

"You have no evidence that any design binds TfR1."
ESMFold pLDDT = monomer folds. It says nothing about binding.
De novo binder hit rates are single-digit % for hard targets.
20-slot portfolio with no complex evidence = lottery.

## ΔΔG_monomer is Not a Switch Metric

- Our -19 kcal/mol = pH-dependence of FOLDING stability, not binding
- actimo-labs' 2.84-4.56 = actual binding free energy difference between pH values
- Attribution of gap: ~70% from His not at interface, ~20% from His count, ~10% from monomer vs complex
- "Drop -19 kcal/mol from your narrative. It reads as misunderstanding."

## Charge Swing is Not a Binding Switch

- +1.5 charge at pH 7.4 vs +13 at pH 5.5 = global electrostatic change, not specific binding
- Nonspecific stickiness at pH 5.5 will look like "binding" but binds everything negative
- If counter-screen shows biphasic kinetics, we lose

## Mechanistic Inversion Warning

"Your biggest mechanistic inversion: His-biased binders tune the binder's own charge
rather than exploiting the target's switch."

At pH 5.5, His707 is protonated (+). Binder should present D/E or polar acceptor opposite it.
A binder His at that contact is REPELLED by protonated target His.

## Boltz-2 Priorities (in order)

1. Does it bind? (ipTM, interface PAE, multi-seed pose consistency)
2. Is epitope at His707? (atoms within 4-5Å of target His707)
3. Chemistry at His707 contact (D/E opposite protonated His707, not His)
4. Complex PROPKA (ΔpKa per interface His)
5. Two-state energy (Rosetta with neutral vs protonated His)
6. Portfolio selection (epitope diversity, 2-4 interface His, charge limits)

## BindCraft: Keep as Binding Hedge

- Best published wet-lab hit rates among our methods
- Consider adding 2-4 His at interface positions → His scan, not redesign
- Keep 2-3 such hedges in the 20

## Biggest Blind Spots

1. No binding evidence for any design
2. Goodharting PROPKA (selecting for His count)
3. Developability (His-rich → IMAC/Ni-NTA interference?)
4. Mechanism assumption (is His707 really the dominant lever?)
5. No negative control or specificity check

## Revised Timeline

- Days 1-2: Boltz on WIDE set (all 20), multi-seed
- Day 3: epitope + His707 contact filter
- Days 3-4: complex PROPKA + two-state interface energy
- Day 5: BindCraft His scan
- Days 6-7: portfolio selection (epitope diversity, 2-4 interface His, charge limits)
- Remaining: developability checks, writeup
