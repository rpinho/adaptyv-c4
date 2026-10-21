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

## The acid budget: the term the first two versions of this method missed

**Measured 2026-10-20, on real BindCraft2 output, not assumed.** The pH layer had never been
executed end to end. Running it on the two ranked Challenge 3 designs (`results/c3_ranked/`, the
first real designs this project produced) gave the result that changes the selection rule.

Design `c3_glp1r_ecd_groove_forced_l99_e26affcf9b8831e4_seq0`:

| group | pKa free | pKa bound | dpKa | contribution to L |
|---|---|---|---|---|
| **B33 HIS** (designed) | 6.43 | 3.76 | **-2.67** | **-0.521** |
| **B26 HIS** (designed) | 6.20 | 1.87 | **-4.33** | **-0.386** |
| A29 N+ | 7.96 | 7.76 | -0.20 | -0.049 |
| B34 GLU | 5.05 | 5.89 | +0.84 | +0.192 |
| **A127 GLU** (target) | 4.83 | 6.49 | +1.66 | **+0.534** |
| **A68 GLU** (target) | 4.92 | 7.91 | **+2.99** | **+1.255** |
| | | | | **total L = +1.20** |

Two buried histidines delivered exactly the mechanism this method is built on - a 2.7 and a 4.3
unit depression, worth **-0.91** between them - and the design still came out with the **wrong
sign**, because three carboxylates on the *target* had their pKa pushed **up** on binding and
contributed **+1.91**.

The physics is symmetric and I had only used half of it. Burying a neutral histidine against a
desolvated interface depresses its pKa, which releases a proton on binding and makes binding weaker
at low pH. Burying a *charged* carboxylate does the mirror image: desolvation and neighbouring
acids raise its pKa, so it **accepts** a proton on binding and makes binding **stronger** at low
pH. A glutamate whose pKa rises from 4.9 to 7.9 is a proton acceptor sitting right in the
pH 6.0-7.4 window, and a single one of those outweighed both designed histidines.

Three consequences, all of them acted on:

1. **The design metric is a budget, not a count of histidines.** `scripts/ph_linkage.py` now writes
   a per-residue CSV (`--csv` also emits `*_residues.csv`) and prints both ends of the sum, because
   the positive tail is where a design is lost. Never read the His rows alone.
2. **Epitope choice is the only lever on the target's acids** - the receptor cannot be mutated.
   `scripts/epitope_charge.py` censuses every titratable group a binder would bury within 10 A of
   the hotspots. Measured on the candidate set:

   | epitope | patch | Asp/Glu | His | Lys/Arg | formal net | acid:His |
   |---|---|---|---|---|---|---|
   | **E1 His707 (human, live campaign)** | 49 | **2** | **2** | 6 | **+4** | **1.0** |
   | E1 His710 (mouse ortholog) | 67 | 4 | 2 | 6 | +2 | 2.0 |
   | E2a apical 202-213 | 53 | 5 | **0** | 7 | +2 | 5.0 |
   | E2b apical 340-376 | 100 | **10** | 1 | 11 | +1 | 10.0 |
   | E4 His603 | 33 | 5 | 1 | 3 | **-2** | 5.0 |

   The epitope the campaign is running on is the best of the set on this axis, and it is the best by
   a wide margin: two acids against two target histidines, in a patch carrying six basic residues.
   A base-rich patch is also what depresses a buried histidine's pKa, so the same property helps
   twice. This is an **independent confirmation of a choice made for structural reasons** (dimer
   occlusion and 7/7 human/mouse epitope identity), arrived at from different physics.
3. **Both apical sites are now dead for a second, independent reason.** They were already
   compromised by 0.59 human/mouse identity across the patch. E2a carries **no histidine at all**,
   so there is no lever to install against, and E2b buries ten carboxylates. E4 His603 is
   acid-dominated and formally negative. Nothing in the portfolio should be spent on them.

Two further things this measurement hands us, both free:

- **A694 and A714 GLU are the two acids to keep out of the interface** on the live epitope. They go
  in as `coldspots` on the next campaign generation rather than being discovered after the fact.
- **A699 and A707 HIS are target-supplied levers.** Mechanism family 6 does not need the binder to
  carry a histidine at all; a binder that depresses the receptor's own His707 wins the same
  negative contribution for free. Their free-state pKa on the apo receptor is the number to measure
  next, and it decides how much depression is available.

**Honest statement of what this is.** PROPKA is empirical and single-structure, and these pKa shifts
are computed on predicted complexes, not crystal structures. The budget is a *ranking* tool for
selecting among designs and epitopes, not a prediction of the assayed pH ratio. That is why the
two-estimator rule above is not optional, and why nothing is claimed about a design until the
independent APBS arm agrees on the sign.

## Status of the layer (2026-10-20)

- `scripts/cif2pdb.py` **new and needed**: BindCraft2 writes an `_atom_site` loop with no
  `occupancy` column, which Biopython's `MMCIFParser` rejects with a bare `KeyError`. The converter
  reads the loop by its own header order, never by fixed columns, so it survives a writer change.
- `scripts/ph_linkage.py` **executed on real designs for the first time**, propka3 3.5.1 in the
  local CPU venv (`.venv/`, Python 3.9). The pH layer deliberately runs on CPU here rather than on
  the rented GPUs, which stay on trajectories.
- `scripts/epitope_charge.py` **new**, see above.
- `scripts/install_switch.py` remains **untested on real designs** - no Challenge 4 design has been
  accepted yet. It is next, the moment one lands.
