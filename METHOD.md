# Method: proton-coupled two-state design, with the numbers checked

## The metric

`scripts/ph_linkage.py` scores a complex by Wyman linkage. For every titratable group with
free-state pKa `pf` and complex pKa `pb`:

    L = log10 Ka(pH_low) - log10 Ka(pH_high)
      = SUM_i [ t(pb,pf,pH_low) - t(pb,pf,pH_high) ],   t(pb,pf,pH) = log10((1+10^(pb-pH))/(1+10^(pf-pH)))

L is in log10 units of affinity, so the challenge's first criterion — measurable at pH 7.4,
undetectable at pH 6.0 — needs roughly **L <= -2**. pKa values come from PROPKA3 run three times
(complex, free binder, free target) and matched by chain and residue number.

## Calibration (this corrects the first version of the plan)

The plan originally claimed "2-3 interface histidines with pKa depressed 1-1.5 units gives about
2 log units". **That is wrong, and the first real scorer output is what exposed it.** On a PDL1
smoke design, a binder histidine with pKa 6.27 free -> 5.37 bound (a clean 0.9-unit depression)
contributed only **-0.34** log units, not ~1.

`scripts/linkage_budget.py` computes why. Contribution per group, pH 6.0 vs 7.4:

| pKa_free | depression 1.0 | 2.0 | 3.0 | 5.0 |
|---|---|---|---|---|
| 6.0 | -0.24 | -0.28 | -0.28 | -0.28 |
| 6.6 | -0.49 | -0.62 | -0.63 | -0.63 |
| 7.0 | -0.61 | -0.86 | -0.89 | -0.90 |
| 7.4 | -0.61 | -1.02 | -1.11 | -1.12 |
| 8.2 | -0.32 | -0.95 | -1.28 | -1.34 |

Two things follow, and both change the design instructions:

1. **A histidine saturates at about -1.1 to -1.3 log units**, and only when its *free-state* pKa is
   around 7.4-8.2 **and** the depression on binding is 2 units or more. A histidine sitting at the
   textbook 6.0-6.5 with a 1-unit shift buys about -0.25, so a design would need five or six of
   them. The real target is therefore **2 to 4 interface histidines with free-state pKa 7.0-8.0
   whose pKa drops by at least 2 units on binding**, which is a much sharper instruction than
   "install some histidines".
2. **The optimum is not free.** A histidine with free pKa 7.4 is already half protonated at pH 7.4,
   so if protonation disfavours binding, affinity at the operating pH is reduced roughly twofold.
   That is an acceptable price for a hundredfold switch, but it has to be stated, and it means the
   pH 7.4 binding energy has to be strong enough to absorb it.

## What that implies for the mechanism families

Raising a histidine's *free* pKa and depressing it hard on *binding* are two separate asks, and
they point at a specific construction:

- Put a carboxylate (Asp/Glu) next to the switch histidine **on the binder's own surface**, not
  across the interface. In the free state that pairing raises the histidine's pKa into the useful
  7.0-8.0 window; on binding the pair is desolvated and its geometry broken, which depresses the
  pKa sharply. This is a refinement of mechanism family 1 and it is what the numbers above argue
  for, rather than the across-the-interface dyad the plan first described.
- Mechanism family 3 (burying an unpaired histidine against an apolar interface) gives the large
  depression but tends to leave the free-state pKa near 6.5, so it should be paired with family 1
  rather than used alone.
- Mechanism family 6 (the target's own His707 or His603) is unaffected by this: the target supplies
  a histidine whose free pKa is whatever the receptor makes it, and the binder's job is to depress
  it. Worth measuring the free-state pKa of His707 and His603 on the apo receptor before committing.

## The two-estimator rule still holds

PROPKA is empirical. The second, independent estimator (PDB2PQR titration states feeding APBS
Poisson-Boltzmann at both pH values) is installed in the `phscore` env and must agree before a
design is submitted. `scripts/install_switch.py` supplies the orthogonal *binding* check: Rosetta
InterfaceAnalyzer `dG_separated` on the neutral state, so a variant that wins the switch by
destroying pH 7.4 binding is discarded (`--max-ddg`, default 3.0 REU).
