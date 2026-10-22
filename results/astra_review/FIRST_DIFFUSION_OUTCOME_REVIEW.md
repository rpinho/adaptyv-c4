# First diffusion outcome review — superseded

**Status: historical review only.** Ricardo canceled this generator lane before this review was finalized. The provisional exploratory MPNN/validation and subsequent batch recommendations sent to root are withdrawn from the active plan. No sequence sampling, refolding allocation or further generator work is authorized by this file. Root reported all project GPU resources released.

## Observations preserved

The fixed v18 backbone completed50 diffusion steps. `BACKBONE_AUDIT.json` records credible backbone geometry but rejects the preregistered coverage criterion: two of five hotspot Cα positions within8Å, fraction.4 versus required.5. Its `sequence_design_allowed=false` and zero qualification remain unchanged.

Separate `EXPLORATORY_NATIVE_CONTEXT.json` records a native195-residue fit of0.0887Å RMSD, complete restored receptor backbone clearance4.1436Å, observed-NAG clearance24.0858Å and no native-context contradiction under those checks. This exploratory assessment did not revise the primary result.

I independently inspected the restored pose onCPU. Only native subunitA Val210 lies within5Å of the poly-Gly backbone by any heavy-atom distance. Thus the geometry supplies limited local proximity, not a well-established interface. A poly-Gly backbone has no designed side-chain interface and no binding, pH or mouse evidence.

The official challenge permits binders outside the recommended epitope; our exact five-hotspot coverage criterion is an internal pilot filter. RFdiffusion's documented hotspot convention uses Cβ proximity and its published filtering uses sequence-designed structures and AF2 with target templates plus an initial guess. These are different from a five-position Cα coverage threshold or our coordinate-free complex test. Removing the initial guess imposes an additional docking task; it is not automatically a better calibrated biological rejection rule. [Official challenge](https://proteinbase.com/competitions/anthropic-adaptyv-2026/challenges/tfr1), [official RFdiffusion guidance](https://github.com/RosettaCommons/RFdiffusion#practical-considerations-for-binder-design).

Before cancellation, I judged a separately labeled exploratory follow-up scientifically interpretable if the original rejection remained intact. That provisional recommendation is now superseded, along with the proposed generator batch. No implementation or follow-up experiment was performed by this reviewer.

Active work has moved to root's separate CPU-only6W3H mapping/calibration task. This document is not a plan for that task and does not change admission or geometry gates.

Reported spending at cancellation: Modal runtime upper estimate$10.7189472, Runpod posted$56.589988, consultations$1.87865; additional costs as tracked by root. No new GPU spending from this review.
