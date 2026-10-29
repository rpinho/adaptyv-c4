# FINDING: L6 is noise-dominated near the gate, and v4's top seat is not acid-OFF

Written 2026-10-26 from the real clock. CPU only, $0, read-only inputs. No GPU, no submission,
no organiser contact. Frozen files 4/4 verify. All four controls fired.

Mandated by Opus 5.5 and Fable 5 in the 2026-10-26 reviewer round, both of whom called it the
highest-value free measurement available.

## Method

Perturb structures already on disk by amounts that should not change the answer, and re-score:

- **rigid** global rotation + translation. A control: pKa is frame-invariant.
- **jitter** isotropic Gaussian coordinate noise at 0.05 / 0.10 / 0.25 A RMS per atom, far
  below the coordinate uncertainty of a predicted complex.

Fable 5 proposed the fold-seed route. Multi-seed structures were never saved for the 342 pool,
so that route needs GPU. This is Opus 5.5's alternative and tests the same property.

## Pool result, 70 designs x 8 replicates

| jitter | median sd | p90 sd | max sd | designs flipping sign | sd exceeds abs(L6) |
|---|---:|---:|---:|---:|---:|
| 0.05 A | 0.061 | 0.211 | 0.612 | 17 / 70 | 14 / 70 |
| **0.10 A** | **0.158** | 0.347 | 0.748 | **19 / 70** | **23 / 70** |
| 0.25 A | 0.298 | 0.546 | 0.989 | 31 / 70 | 26 / 70 |

The write-precision floor measured **0.00013**, so these are real, not rounding artifacts.

Opus 5.5 predicted: *"If per-design sigma is even 0.3, the sign of 16 of 116 iPTM-eligible designs
is undetermined and 4 of 20 seats are coin flips."* Measured sigma reaches 0.3 at 0.25 A and
19 of 70 designs already flip at 0.10 A.

## The decision-relevant sets, 10 replicates at 0.10 A

### v4 acid-OFF seats

| design | L6 | sd | range | flips / 10 |
|---|---:|---:|---|---:|
| site1_161 | -2.324 | 0.267 | -2.685 .. -1.937 | 0 |
| site1_0 | -1.865 | 0.416 | -2.518 .. -1.503 | 0 |
| site2_25 | -1.858 | 0.082 | -1.943 .. -1.672 | 0 |
| site2_177 | -1.154 | 0.071 | -1.190 .. -0.987 | 0 |
| site2_181 | -0.743 | 0.129 | -0.863 .. -0.396 | 0 |
| site2_66 | -0.628 | 0.107 | -0.781 .. -0.454 | 0 |
| site1_7 | -0.578 | 0.252 | -0.706 .. +0.147 | 1 |
| site2_92 | -0.577 | 0.214 | -0.869 .. -0.129 | 0 |
| site2_100 | -0.425 | 0.284 | -1.132 .. -0.172 | 0 |
| site2_74 | -0.352 | 0.040 | -0.432 .. -0.297 | 0 |
| site2_169 | -0.336 | 0.093 | -0.410 .. -0.106 | 0 |
| site1_149 | -0.163 | 0.251 | -0.443 .. +0.436 | **8** |
| site1_177 | -0.160 | 0.159 | -0.241 .. +0.289 | 1 |
| site1_49 | -0.066 | 0.254 | -0.331 .. +0.519 | **6** |
| **site1_67 (SEAT #1)** | **-0.033** | **0.430** | **+0.045 .. +1.292** | **10** |

**v4's top-ranked seat is positive in every single perturbed replicate.** Its qualifying value of
-0.033 is indistinguishable from noise, and the perturbed range does not include it.
5 of 15 seats have undetermined sign; 3 of 15 have sd exceeding abs(L6).

### The eight the reviewers chose

| design | L6 | sd | range | flips / 10 |
|---|---:|---:|---|---:|
| site1_13 | -1.635 | 0.196 | -1.739 .. -1.229 | 0 |
| site1_147 | -1.609 | 0.191 | -1.618 .. -1.058 | 0 |
| site1_113 | -1.486 | 0.341 | -2.341 .. -1.434 | 0 |
| site2_24 | -0.826 | 0.101 | -0.953 .. -0.620 | 0 |
| site1_143 | -0.778 | 0.127 | -0.922 .. -0.488 | 0 |
| site2_128 | -0.555 | 0.388 | -0.529 .. +0.646 | **4** |
| site1_125 | -0.503 | 0.074 | -0.508 .. -0.289 | 0 |
| site2_55 | -0.236 | 0.189 | -0.498 .. +0.104 | 1 |

**2 of 8 undetermined, 0 of 8 with sd exceeding abs(L6).** The top five are stable under every
perturbation. The reviewers selected this set on native-accessibility and confidence grounds
without any noise information; it is independently more robust than the operative v4 ranking.

## Consequences

1. **The L6 < 0 gate is noise-dominated for designs near zero.** Any design with abs(L6) below
   roughly 0.3 has an undetermined direction at perturbations far smaller than the uncertainty in
   the structures themselves.
2. **v4 seat #1 should not be presented as acid-OFF.** This is a property of the ranking, not of
   the molecule: nothing here says the design is bad, only that the value which qualified it does
   not survive a perturbation that should not matter.
3. **It strengthens the case for the eight over the v4 seats**, on evidence the reviewers did not
   have when they chose.
4. Combined with the magnitude finding, the usable rule is a band on abs(L6) well away from zero,
   not a sign test. Setting that band remains a decision for the four and requires a documented v5.

## Controls, all fired

| control | result |
|---|---|
| NULL_ZERO_JITTER | zero-magnitude jitter reproduces L6 exactly |
| RIGID_INVARIANCE | rigid-body transform moves L6 by 0.00013, the write-precision floor |
| JITTER_MONOTONIC | median sd rises 0.061 -> 0.158 -> 0.298 with jitter |
| STRATIFIED_PORTION_REPRESENTATIVE | stratified subset reproduces the pool median within 0.30 |

Two controls initially failed and were corrected before any result was reported: RIGID_INVARIANCE
used an unphysical 1e-6 tolerance when PDB coordinates are written to 3 decimals, and the sample
control required the whole sample to match the pool median while deliberately oversampling the
gate region. Both corrections are in the script.

## What is NOT claimed

- **No binding, affinity or pH response is measured.** Coordinate perturbation plus PROPKA only.
- **Coordinate jitter is not the same as fold-seed variation**, which would likely be larger. This
  is a lower bound on the instability, not an estimate of it.
- **A stable L6 is not a correct L6.** Designs whose sign survives perturbation are still scored
  on the apical construct, which the designs were not generated against.
- **No ranking was built and no design was selected or rejected here.**
