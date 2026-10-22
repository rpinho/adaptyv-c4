# Method: proton-coupled two-state design, with the numbers checked

## Current checkpoint — sequence-only refolds completed; native-domain experiment prepared

V15 completed eight coordinate-free fixed-sequence prediction jobs; actual GPU
release independently verified. Both free folds pass both checkpoint geometry
checks (pLDDT .796-.935; internal model agreement .98/1.28 A). Native receptor
complexes remain off-patch (52.38-70.50 A), iPTM .148-.197; three restored-context
clashes, zero credible TfR1 binders. This separates plausible sequence folding
from failed interface generation. No pH or mouse-binding claim.

Astra's outcome and domain review are in results/astra_review/. Exactly one
native contiguous A189-383 generation experiment is prepared: only the target
representation changes vs native v9; same seed, lengths, losses and gates.
Installed-image CPU frames, native identity, exact hotspots and author-numbered
CIF->PDB output roundtrip pass; 80 local checks pass. A direct CIF reader uses
label IDs and failed the first CPU build; the corrected test uses the actual
CIF->PDB author-ID downstream path. No GPU was allocated during this correction.
Protocol: results/apical_context/DOMAIN_PROBE_DECISION.md. Launch/live evidence
belongs in results/native_domain_v16/. The post-release extraction prepares one
predeclared hard sequence even if its raw search pose fails; full-native refolding
still requires its own manifest-specific CPU image preflight. Nothing is admitted
or submitted automatically. No same-recipe seed successor is scheduled.

Modal completed runtime estimate $9.5172782; Runpod runtime/history estimate
$35.6652; consultation usage $1.1922999, setup/storage unreconciled. H100 GPU
$3.9492/hr, requested CPU/RAM upper $4.649328/hr combined. Continued spending is
already explicitly authorized. Runpod remains stopped. Earlier live states below
are historical; use current supervisor files for actual allocation/execution.


**Current qualification (2026-10-21):** the old 73-residue receptor crop omitted
protein that all 65 well-aligned source candidates overlap at backbone level.
Their linkage results diagnose that crop; they do not validate binding to the
full receptor. See `results/full_context_audit/RESULTS.md`. Independently
prepared free-binder sidechains also remove the strongest proposed H46 free-pKa
advantage in the 72-model scout: `results/free_networks_v3/RESULTS.md`. These
findings, and fresh full-context generation, supersede advancing the old numerical
winners directly to submission. The equations below remain model definitions.

## The metric

`scripts/ph_linkage.py` scores a complex by Wyman linkage. For every titratable group with
free-state pKa `pf` and complex pKa `pb`:

    L = log10 Ka(pH_low) - log10 Ka(pH_high)
      = SUM_i [ t(pb,pf,pH_low) - t(pb,pf,pH_high) ],   t(pb,pf,pH) = log10((1+10^(pb-pH))/(1+10^(pf-pH)))

L is in log10 units of affinity. **L <= -2 is our internal 100-fold model target,
not an official challenge cutoff or a calibrated detector limit.** The lab defines
measurable versus undetectable binding. PROPKA values are computed on complex,
isolated binder and isolated target and matched by chain and residue number.
An extracted chain retains the bound conformation; it does not establish the
independently relaxed free-state ensemble. See `WIN_PLAN.md` for the current strategy.

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

PROPKA is empirical. The second-estimator requirement remains, but the previous identification of
PROPKA-driven PDB2PQR/APBS as an independent pKa estimator was incorrect. PDB2PQR
can reuse PROPKA's titration assignment; APBS then evaluates electrostatics for
those charges. A genuinely independent titration calculation must be configured
and checked before claiming independent agreement. PDB2PQR is available locally;
an APBS executable was not found on PATH in the strategy review. See the
[PDB2PQR documentation](https://pdb2pqr.readthedocs.io/en/stable/using/index.html).

`scripts/install_switch.py` supplies a neutral-state structural filter: negative
reference/variant InterfaceAnalyzer dG and ddG <= 3 REU. It is neither a measured
pH 7.4 affinity nor a second pKa estimator.

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

### What sets the sign, measured rather than reasoned

The first draft of this section said to declare the interface acids as `coldspots`. That is wrong
on its own terms - A714 GLU **is one of the chosen hotspots** - and the per-residue table says
something sharper. Sorting all 33 groups by contribution:

| group | dpKa | contribution | why |
|---|---|---|---|
| A88 TYR | **+3.84** | **+0.000** | pKa 13.2 -> 17.0, never near the window |
| A67 ASP | **-1.57** | **-0.000** | pKa 2.95 -> 1.38, never near the window |
| A138 GLU | -0.66 | -0.009 | acid, pKa went DOWN, helps |
| A68 GLU | +2.99 | **+1.255** | pKa 4.92 -> **7.91**, lands inside the window |

**Size of the shift is almost irrelevant. What counts is how much of the pH 6.0-7.4 window the
titration crosses.** A 3.8-unit shift outside the window is worth exactly nothing; a 1.7-unit shift
into it costs half a log unit. Acids are not the enemy - acids whose pKa is pushed *into* the window
are, and an acid whose pKa goes *down* (salt-bridged to a binder cation) contributes the right sign.
The rule is therefore **charge-compensate every buried carboxylate**, not exclude it.

### The ceiling on the target's own histidines (this kills mechanism family 6 alone)

propka3 on the apo trimmed human target (`results/apo/E1_his707_tight.pka`) gives the free-state
pKa values the plan said to measure next:

| group | apo pKa | model | ceiling at total desolvation (pKa -> 2.0) |
|---|---|---|---|
| HIS 707 A | **6.09** | 6.50 | -0.328 |
| HIS 699 A | **5.81** | 6.50 | -0.205 |
| HIS 684 A | 6.35 | 6.50 | -0.473 |
| GLU 694 A | 4.60 | 4.50 | (+1.18 if it rises as A68 GLU did) |
| GLU 714 A | 4.55 | 4.50 | (+1.16 if it rises as A68 GLU did) |
| GLU 533 A | 4.55 | 4.50 | (+1.16 if it rises as A68 GLU did) |
| **all three His together** | | | **-1.01** |

The gate needs **L <= -2.00** and the target can supply at most **-1.01**, at a physical limit no
binder reaches. **Mechanism family 6 cannot clear the gate on its own** - the receptor's histidines
are a useful bonus of roughly -0.3 each, not a strategy. The plan must not rest on them.

Worse, a *single* interface glutamate rising as far as A68 GLU actually rose in a real design
(+2.99) contributes **+1.18**, which cancels one whole well-designed histidine. The downside term is
the same size as the upside term.

### Where the leverage actually is: free-state pKa, not dpKa

The same arithmetic, applied to a binder histidine:

| binder His free pKa | driven to | contribution each | two of them |
|---|---|---|---|
| 6.0 | 2.0 | -0.28 | -0.56 |
| 7.0 | 3.5 | -0.90 | -1.79 |
| **7.4** | **3.5** | **-1.12** | **-2.23 (clears the gate)** |
| 7.8 | 3.5 | -1.26 | -2.52 |

This explains the measured design completely. Its histidines were depressed by 2.67 and 4.33 units -
far more than the 1.5 units needed to saturate - and returned only -0.52 and -0.39, because they
*started* at 6.43 and 6.20. **The design objective is not a large dpKa. It is a high free-state pKa
that collapses on binding.** Everything past about 1.5-2.0 units of depression is wasted effort.

So the binder specification is two-state and explicit:

- a histidine whose **free-state pKa is 7.4 or above**, which means siting it where the *unbound
  binder* presents carboxylate neighbours and solvent exposure;
- **buried against the target's cations when bound** - and the E1 His707 patch carries six Lys/Arg
  (net formal +4), which is exactly the environment that collapses a histidine's pKa. The property
  that makes this epitope acid-poor is the same property that makes it a good switch site;
- **two such sites minimum**, three if any interface carboxylate cannot be charge-compensated;
- A699 and A707 HIS counted as a bonus of about -0.5 together, never as load-bearing.

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
- `scripts/switch_sites.py` **new**. `install_switch.py` ranked candidate positions by *interface
  burial*, which optimises only the bound half of the two-state problem and is precisely the error
  the measurement exposed. This ranks them on both halves from coordinates alone - exposure plus
  carboxylate neighbours in the free binder against burial plus target cations in the complex - so
  the expensive Rosetta pass runs on a short list. On the only two positions where an independent
  propka answer exists it agrees in order and within about 20% in magnitude (B33 HIS: estimate
  -0.62, propka -0.52; B26 HIS: estimate -0.50, propka -0.39). Two points is consistency, not
  validation, and it is a pre-screen, never a pKa claim.
  It also rated four positions in that design *above* the two histidines BindCraft2 placed, which
  is the use case: pick switch sites deliberately instead of accepting whichever histidines the
  generator happened to leave at the interface.
- `scripts/install_switch.py` remains **untested on real designs** - no Challenge 4 design has been
  accepted yet - and must be re-pointed at `switch_sites.py`'s ranking instead of its own burial
  sort before it is used. It is next, the moment a design lands.
- Best achievable pair on the two existing designs scores about -1.4 and -1.65 estimated, short of
  the gate, which is consistent with the budget above: two sites are the floor, not the plan.

## Where the binder actually sits, measured (2026-10-21)

The layer above was built and calibrated before any design had been observed docking. Thirty-seven
refolded complexes now exist, and they move one assumption from reasoned to measured.

**The epitope is the hydrophobic groove, not the His707 patch.** 37/37 designs contact
LEU709/LEU713/LEU716/PHE685/PHE701/TRP702/PHE731/LEU735/THR739/ILE742 and 0/37 contact
SER704/GLY705/SER706/HIS707/ALA711. HIS707 sits 4.0 A from LEU709, at the rim.

**What this changes for the proton budget, and it is an improvement.** The method already concluded
that the leverage is the *binder's* free-state pKa, not the target's histidines (ceiling -1.01 for
all three, against a -2.00 gate). A binder histidine buried against a hydrophobic groove loses more
solvation than one pressed onto a polar patch, so the bound-state pKa drops further and each
histidine is worth more than the -1.115 estimated against the polar patch. The groove is also
base-rich (LYS717, LYS720, ARG719, LYS534, ARG732 against GLU714, GLU728, GLU533), and the measured
rule from the previous section is that bases help and acids fight: this epitope's charge census is
on the right side of the ledger, which the apical sites were not.

**What must now be re-measured rather than inherited.** Every number in the sections above was
computed against designs docked on the polar patch, or against the apo target. The per-residue
contributions, the acid budget and the `switch_sites.py` rankings all have to be recomputed on
groove-docked designs before any of them is quoted again. The *rules* carry over - sign is set by
how much of pH 6.0-7.4 a titration crosses, acids can go either way and must be charge-compensated
rather than excluded - but the magnitudes do not.

**Honest status.** The groove's own effect on a designed histidine's bound pKa is a prediction from
the desolvation term, not yet a measurement. It is testable the moment designs clear the filters:
run `ph_linkage.py` on groove-docked complexes and compare the per-residue contributions against
the -0.52 and -0.39 measured on the polar patch. If the groove does not deepen them, the
improvement claim in this section is wrong and should be struck.

## The groove prediction was wrong, and the reason is the rule I failed to apply (2026-10-21)

The section above predicted that burying a designed histidine against a hydrophobic groove would
deepen its contribution beyond the -0.52 and -0.39 measured on the polar patch, and stated that if
it did not, the claim gets struck. It did not. **Struck.**

Three groove-docked designs, same backbone, same histidine position (C60):

| design | binder His C60 free pKa | bound pKa | dpKa | contribution | total L |
|---|---|---|---|---|---|
| cand3 | 6.69 | 6.15 | -0.54 | **-0.335** | -0.643 |
| cand2 | 5.52 | 5.00 | -0.52 | **-0.079** | -0.452 |
| cand5 | 5.42 | 4.95 | -0.47 | **-0.061** | -0.399 |

Against -0.52 and -0.39 on the polar patch, the groove made the binder histidine **worse**.

**Why, and it is the method's own measured rule turned back on me.** This layer already established
that the leverage is the **free-state** pKa, not dpKa. Desolvation does not act only on the bound
state: a histidine positioned in or beside a hydrophobic groove is *already partly buried in the
free monomer*, so the same apolar environment that depresses the bound pKa depresses the free pKa
too. The free pKa lands at 5.42-6.69 instead of near 7.4, the titration therefore barely crosses the
pH 6.0-7.4 window, and dpKa collapses to -0.47 to -0.54. The two effects cancel, and on net the
groove loses. I predicted a one-sided benefit from a two-sided mechanism.

**What did work, for a different and already-documented reason.** Total L moved from **+1.20 on the
polar patch to -0.64 here** - a 1.84-log improvement that also **flips the sign to the required
direction**. That gain is the **acid budget**, not desolvation: the polar patch's three glutamates
contributed +1.91 and dominated everything, and the base-rich groove simply does not have them. So
the epitope change was right, and right for the reason the acid-budget section gave, not the reason
the groove section guessed.

**Where His707 actually ended up.** It is contacted by none of the designs, yet contributes
-0.103 to -0.121 in all three, because at 4.0 A from LEU709 the binder perturbs it from the rim.
HIS699 contributes more (-0.120 to -0.225). Small, consistent, and in the right direction - but
together the target's histidines still sit inside the -1.01 ceiling measured on the apo structure.

**The gate is still missed and the route to it is now specified by data.** L = -0.64 against a
-2.00 gate leaves about 1.4 log units. The requirement is explicit: a designed histidine must be
**solvent-exposed in the free binder** (free pKa near 7.0-7.4) and buried **only** on binding.
C60 fails the first half. This is exactly the quantity `switch_sites.py` was written to score - free
-state exposure and acidic neighbours for a high free pKa, against burial and target cations for a
low bound pKa - and `install_switch.py` has never been wired to use it, still running its own
burial sort. That wiring is now the highest-value remaining step, and these three accepted
backbones are the first real substrate for it.

**Caveat on the arithmetic.** The target is split across chains A (532-743) and B (686-693) and this
run scored chain A only, so chain B's eight residues are absent from the sum. They are not at the
interface the designs use, but the figures above are chain-A-complete rather than target-complete.

## The free-state pKa rule, now confirmed prospectively on 66 designs (2026-10-21)

The campaign completed on its own at its configured 25 final designs (177 trajectories, 14%
acceptance, against 0 of 99 on the two mis-aimed epitopes). Sixty-six passing complexes were
scored. This is the first time the layer has had a population rather than a handful, and the
population settles the central rule.

| quantity | value |
|---|---|
| best L | **-1.175** |
| median L | -0.387 |
| designs with L < 0 (correct direction) | **62 / 66** |
| designs clearing the -2.00 gate | 0 / 66 |
| binder histidines scored | 68 |
| their free pKa | 3.54 - 7.16, mean **6.17** |
| of those, free pKa >= 7.0 | **3** |

**The rule predicts which histidines matter, before looking at their contributions.** The three
deepest single-histidine contributions in the whole set are **-0.822, -0.659 and -0.658**, and all
three are the *same position* - C24 on the l58_9f94b53776050aad backbone - at free pKa **7.14, 7.16
and 7.14**. They are also precisely the three histidines in the set with free pKa >= 7.0. The other
65 histidines, averaging free pKa 6.17, produce the mediocre bulk. Sorting by free-state pKa and
sorting by contribution give the same answer, on 68 independent cases. The rule was inferred from
two residues on 10-20; it now holds on a population, and it is what `switch_sites.py` scores.

**And the design carrying that histidine is the best design overall** (L = -1.175), which is the
useful form of the result: the single free-state pKa of one residue, computable before any linkage
sum, ranked the whole campaign.

**Two levers, and the data separates them.** C24 goes 7.14 -> 5.67. The earlier estimate for a
histidine making the full excursion (7.4 -> 3.5) was -1.115 each. We are getting -0.822 because the
free pKa is 7.14 rather than 7.4 *and* because the bound state only reaches 5.67 rather than 3.5 -
the burial is incomplete. So the gap to the gate decomposes into: raise the free pKa toward 7.4
(position the histidine in a solvent-exposed, acid-flanked site in the free monomer) and deepen the
bound-state burial (place it where the groove and the target's cations desolvate it). Those are the
two terms `switch_sites.py` already scores separately, and nothing in the design loop is using
either yet.

**Honest arithmetic on the gate.** At -0.822 per well-placed histidine, the -2.00 gate needs about
2.4 of them, or two that make a fuller excursion. Two histidines at the full 7.4 -> 3.5 excursion
would give -2.23 and clear it. That is a design target, not an achievement: no design in this set
has two histidines above free pKa 7.0, and whether two can be placed on one backbone without
spoiling the interface is untested.

**What this does not show.** All 66 designs come from one campaign on one epitope, and the 25 ranked
designs derive from far fewer distinct backbones, so these are not 66 independent samples. Chain B
(686-693) is excluded from every sum, as noted above. And no orthogonal estimator has been run:
every figure here is PROPKA3.

## Ranked two-site installation, now executed on all 25 ranked C4 backbones (2026-10-21)

`scripts/install_switch.py` now uses `switch_sites.score_design` directly. There is
no second burial sort. The CPU campaign tested the three best estimated pairs
from eight eligible sites on each ranked backbone, in best historical-L order.
One best-scored candidate represents each backbone; its sequence can differ from
ranked seq0, and both identities are recorded. Existing His sites participate in
HH and QQ twins rather than being dropped as no-op edits.

Real validation exposed two problems the unexecuted installer had concealed:

1. **Reference preparation mattered more than the mutation.** The best raw source
   gave InterfaceAnalyzer dG +804.75 REU when only separated partners were globally
   repacked. Common whole-complex repacking and side-chain minimisation, with the
   backbone fixed, gave -44.66. Local mutation-site packing alone compared a strained
   bound input with globally repacked free partners. The installer now prepares one
   common complex, then applies identical local protocols to WT, HH and QQ. Both
   reference and variant dG must be negative, in addition to ddG <= 3 REU. These
   remain neutral-state structural surrogates, not pH 7.4 affinity measurements.
2. **An omitted His is not a zero contribution.** PROPKA suppresses the C1 ring in
   two variants by coupling it to N+. Those variants retain energy-filter status
   but are explicitly **unscorable** for complete two-site linkage. Partial sums
   are preserved only for diagnosis and excluded from ranking. A regression test
   prevents an apparently passing partial sum from entering selection.

The completed screen uses **AB as target, C as binder**, and all matched groups
(`min_shift=0`). Historical A-only scores order the substrates; comparisons below
use matched prepared AB references, not historical campaign numbers.

| Quantity | Result |
|---|---:|
| Ranked backbones processed | 25/25 |
| HH pairs modelled, with matched QQ attempts | 75 |
| HH passing neutral-energy filters | 38 |
| HH with complete installed-site linkage | 36 |
| Complete-linkage HH with retained QQ control | 25 |
| HH energy-retained but unscorable | 2 |
| Complete HH with L <= -2 | **0/36** |
| Complete HH with both selected free pKa >= 7 | **0/36** |
| Best complete HH L | **-1.530** |

The best pair is **C:D4H/T46H** on `l63_59b0f47ecf134650`: its matched reference L
is -0.351, HH L -1.530, and QQ L -0.444. HH therefore improves L by **-1.179**
against its reference and **-1.086** against QQ. Rosetta ddG is **+1.74 REU**,
inside the 3-REU filter, with negative reference and mutant dG. Selected free pKa
values are **6.33 and 7.25**: one useful high-free-pKa site, still not two. This is
a computed screen result, not experimental binding or clearance of the pH gate.

The native-coordinate heuristic has a limited real-data check: on 27 native
interface histidines across 17 selected backbone IDs, estimated versus PROPKA
contribution Spearman correlation is 0.567 and mean absolute error 0.153 log units.
Free-pKa estimate mean absolute error is 0.711. This does not validate a mutation
pKa prediction. Rotamers also matter: native C24 on the original best source moved
from free pKa 7.14 to **6.71 in its first prepared reference**, despite unchanged
backbone coordinates. Native-position ranking is a pre-screen, not evidence that
mutating its top two sites creates two elevated free pKa values.

Every saved source/output checksum, mutation identity, target sequence, energy
filter and complete linkage sum was audited. **181 PDBs across all 25 backbones
passed; maximum backbone atom displacement was 0.000 A.** Seven regression tests
passed. Diagnostic runs are retained and clearly marked invalidated/superseded.

Limits: three pairs per backbone, not all 28 shortlist combinations; related
backbones; original trimmed receptor with artificial fragment boundaries; ranking
on source rotamers before preparation; one pKa estimator. PROPKA warns about
backbone groups at A683 and A728. Full-context confirmation, APBS, mouse validation
and orthogonal refolding remain outstanding. No account or entry was created or
submitted. Both rented pods stayed stopped; **new GPU spend $0.00**, historical
last recorded pod spend $33.87, persistent storage not freshly reconciled.

Full evidence: `results/switch_install/campaign/RESULTS.md`, `variants.csv`,
`residues.csv`, `integrity.json`, per-design manifests and PDBs. Protocol diagnostics
and reproduction commands are in `results/switch_install/VALIDATION.md` and
`runtime_setup.txt`.

## Coordinate initialization needs an explicit relative-placement preflight

The v4 attempt exposed a separate initialization problem, not evidence for the
cause of v3 rejection: enabling upstream initial_atom_pos retained native target
coordinates while unresolved binder atoms stayed at origin,188.816A from the
hotspot centroid. That attempt was cancelled and its GPU released. No comparison
of accepted-model quality is inferred from the aborted attempt.

The replacement target is the same entire639-residue monomer, rigidly translated
so origin lies15A outward from the hotspot centroid, with7.955A nearest-atom
clearance against the complete native dimer. Atom identities, residue numbers,
sequence and all CA pair distances are preserved to PDB precision. This changes
the relative starting placement of the unresolved binder, not receptor shape.
It remains an optimization hypothesis, not a physical binding measurement.
Opposite-subunit restoration/native-sequence and clash checks remain mandatory
before any accepted structure receives CPU linkage scoring. Exact transform and
integrity measurements: results/modal_campaign_v5/human_11qc_local_frame.audit.json.
