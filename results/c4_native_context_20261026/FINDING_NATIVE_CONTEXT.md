# FINDING: the 195-aa apical surrogate is not a cheap approximation of the C4 target — it is a different surface

> **AMENDMENT, same day. The measurements below stand; the mechanism offered for them is WITHDRAWN.**
> The designs were generated against human TfR1 author residues ~651–741 (the helical domain, as a
> two-subunit groove around His707), not against the apical domain at all. See
> `FINDING_WRONG_DOMAIN.md` in this directory.
>
> - **ρ = +0.679 stands as a measurement.** Its explanation is not "truncation exposed hydrophobic
>   faces and the design method was drawn to them." It is that Boltz-2 was asked to co-fold a
>   helical-domain binder against an apical construct containing none of its intended epitope; with
>   no real epitope available the binder settles on that construct's artificial packing faces.
> - The hydrophobic-attraction story is withdrawn. Astra's review flagged it as "plausible, not
>   demonstrated" and named imposed hotspots as a candidate confound; the hotspots were in fact
>   imposed (`forced_targeting: true`, `min_hotspot_contact_final: 0.5`) on a different region.
> - The four-context SASA table, the 13.0% / 4.9% / 1.5% decomposition, the 22 lost positions and
>   the 36.62 Å membrane-proximity retirement are all unaffected — they describe the apical domain
>   itself and are independent of which designs bind where.
> - Per Astra, "fraction of footprint SASA retained" is a screen, not the ideal statistic: whole-residue
>   accessibility can survive on the wrong face, and water accessibility is not binder accessibility.
>   A better statistic is the fraction of actual binder-buried target area retained, with absolute
>   retained interface area and a check that surviving contacts form a coherent interface. Not computed.
> - The "0 of 342 survive" line means zero pass **my** thresholds, which are conventions and not
>   calibrated. It is not a demonstration of biological futility, and the nine designs at ≥80%
>   retention preclude any blanket impossibility claim.

Decision evidence for the one reserved choice: model the official 640-aa tag-free homodimer, or
continue on the 195-aa apical surrogate.

Written 2026-10-26 from the real UTC clock (see `native_context_gap.json:created_at_utc`).
CPU only, $0, read-only. No ranking, no submission, no organiser contact. All four controls fired.

## The headline

**The design campaign bound preferentially to apical surface that exists only because the apical
domain was cut out of its subunit.** Across the 195 apical positions, how often a position is
contacted by the 342-design pool correlates with how much of its area disappears in the native
homodimer at

    Spearman rho = +0.6789,  permutation p = 5e-5 (the 20,000-permutation floor),  n = 195

Sixteen of the twenty most-contacted positions are materially occluded (≥50% area loss) in the
native assembly. The median design retains only **44.5%** of its footprint surface area. **Nine of
342 designs** retain ≥80%.

This is not a ranking defect and no reranking fixes it. It is a target-definition defect affecting
every design on disk.

## Why it happens, mechanistically

Cutting residues 189–383 out of a 639-residue subunit exposes the faces that were packing against
the rest of the protein. Those faces are hydrophobic interior, so a design or docking method finds
them highly attractive. The most-contacted positions are exactly that: F298, F288, F316, L295,
I290, F313, F297, M283, F264 — an apolar patch that is **protein interior** in the real receptor.

This also explains two earlier observations that looked reassuring:

- The designed interfaces had respectable apolar burial (pool median interface apolar fraction
  0.653, against 0.546 for the real experimental TfR1–HFE interface). They are burying apolar
  surface because they are binding an exposed core face.
- The footprints looked highly mouse-conserved (median 84.5% identical). Buried residues are
  conserved *because* they are buried. Restricted to natively exposed positions, apical mouse
  identity drops to **56.2%**.

## Measured: the apical surface in four contexts

All four contexts are built from the same experimental coordinates (11QC chain A), so the only
variable is context, not conformation. 11QC is P02786 122–760 per subunit (639 aa); the official
target is 121–760 (640 aa). Residue 121 is absent locally and is not measured.

| context | apical SASA (Å²) | |
|---|---:|---|
| CTX1 apical 189–383 alone — **the surrogate** | 10076.0 | 100% |
| CTX2 full 639-aa monomer | 8763.6 | 87.0% |
| CTX3 homodimer | 8268.0 | 82.1% |
| CTX4 homodimer + 6 observed NAG | 8114.4 | **80.5%** |

Decomposition of the 19.5% lost:

| lost to | Å² | % of surrogate area |
|---|---:|---:|
| **the rest of its own subunit** | **1312.3** | **13.0%** |
| the partner subunit | 495.6 | 4.9% |
| observed glycans | 153.6 | 1.5% |

**The dominant risk is intra-subunit, not the dimer interface.** This corrects the framing in the
orchestration brief, which emphasised partner subunit, glycans, dimer interface and
membrane-proximal geometry. Measured:

| risk named in the brief | measured verdict |
|---|---|
| rest of own subunit (not named in the brief) | **largest effect: 13.0% of apical area, 2.6× the partner subunit** |
| partner subunit | real but secondary: 4.9% of area; 14 positions within 8 Å; closest approach 3.14 Å |
| glycan proximity | real but small: 1.5% of area; 20 positions within 8 Å; closest approach 1.44 Å |
| membrane-proximal geometry | **retired.** Closest apical atom to residue 122 is **36.62 Å**. The apical domain is not membrane-proximal in this construct. |

## Measured: exposure, position by position

| quantity | value |
|---|---:|
| apical positions | 195 |
| exposed (rel. SASA ≥ 0.20) in the surrogate | 102 |
| exposed in the native homodimer + glycans | 80 |
| **exposed in the surrogate but NOT in the native target** | **22** |
| materially occluded (≥50% area loss) | 32 |

The 22 positions the surrogate presents as available and the native target does not:

    190, 252, 263, 288, 290, 295, 298, 300, 304, 305, 306, 307,
    311, 312, 313, 314, 315, 316, 317, 319, 320, 322

Nine of the twelve most-contacted positions in the entire campaign are in that list. Worst cases:
**P314 loses 100.0%** of its area, G306 97.7%, S315 95.6%, I190 91.9%, G312 86.5%, T305 84.3%,
L295 82.5%, I290 82.2%, F298 80.3%.

## Measured: per design

| quantity | pool (n=342) | v4's 20 seats |
|---|---:|---:|
| footprint area retained in native context, median | **0.4451** | 0.5033 |
| p10 / p90 | 0.3391 / 0.6282 | — |
| min / max | 0.2230 / 1.0000 | 0.2603 / 0.7732 |
| fraction of footprint contacts materially occluded, median | 0.600 | 0.533 |
| designs with ≥1 occluded contact | **335 / 342** | 20 / 20 |
| designs with ≥25% of footprint occluded | **330 / 342** | — |
| retains ≥80% of footprint area | **9 / 342** | 0 / 20 |
| retains ≥60% of footprint area | 43 / 342 | 4 / 20 |
| natively-exposed contacts per design, median | **8** | — |
| designs with fewer than 5 natively-exposed contacts | 43 / 342 | — |

## The intersection with every objective applied at once

Quantification only. **This is not a ranking and must not be used as one** — the pre-registered
rule requires eligibility and model inventory frozen by sequence hash first, and its L6 line still
encodes the inverted pH direction.

| filter (cumulative) | n |
|---|---:|
| pool | 342 |
| acid-OFF (L6 < 0, the direction C4 asks for) | 231 |
| + iPTM ≥ 0.60 | 76 |
| + retains ≥ 60% of footprint area in native context | 15 |
| + ≥ 10 natively-exposed contacts | 8 |
| + exposed-footprint mouse identity ≥ 0.90 | **0** |
| *for contrast:* v4's actual rule (iPTM ≥ 0.60 and L6 > 0) | 40 |

**Zero designs on disk satisfy the challenge's three objectives simultaneously** once native
accessibility is required. Eight survive to the penultimate line — direction-correct, confident,
and binding mostly real surface — but none of those eight has a strongly mouse-conserved exposed
footprint.

## What this means for the reserved decision

The two options are not "correct but expensive" versus "cheap but approximate". Measured:

**(a) Continue on the 195-aa apical surrogate.** Everything on disk is here. But the surrogate
over-exposes 1312 Å² of intra-subunit packing surface, the campaign preferentially bound it
(rho = +0.68), and the median design loses 55% of its footprint in the real assembly. Further
ranking inside the surrogate re-ranks designs whose epitopes are substantially unavailable. The
cheapness is real; the evidence it produces about the actual target is weak and this is now
quantified rather than suspected.

**(b) Model the official 640-aa tag-free homodimer.** This is the only context in which the
available surface is the real one. It is the only path to evidence about the actual target, and it
would redirect design toward the 80 natively exposed positions instead of the 102 the surrogate
advertises.

A third option the measurement raises, which was not on the table:

**(c) Re-dock or re-design against the 639-aa monomer rather than the full dimer.** 13.0% of the
19.5% area correction comes from the design's own subunit and only 4.9% from the partner. A
single-subunit 639-aa context captures two thirds of the correction at a fraction of a dimer's
cost. Whether that is a worthwhile intermediate is a modelling-cost question, not one this
measurement settles.

**Timing note.** C4 closes 2026-11-02 11:59 UTC. Option (b) means new GPU and a new pipeline on
measured cost of $1.83–2.79 per design; a 20-design run is $37–56, and the only prior C4 GPU
attempt returned zero usable rows. Option (b) is unlikely to produce a validated 20-design file
inside the remaining window. That is an argument about what can be *finished*, not about which
context is correct.

## What is NOT claimed

- **No binding, affinity, KD or pH response is measured anywhere in this document.** This is
  solvent accessibility and interatomic distance on experimental coordinates.
- **The binder poses were not modelled in the native assembly.** This measures whether each
  design's *human epitope* remains available there. A design whose footprint is occluded might
  still bind elsewhere, or bind a shifted pose — untested either way. Occlusion is a reason the
  current pose cannot be retained, not proof the molecule cannot bind.
- **Conformational change is not modelled.** The apical domain is treated as rigid across contexts.
  Real receptors flex, and TfR1's apical domain in particular may move; 11QC is one crystal form.
- **Glycans are only the ones observed in 11QC** (3 NAG per subunit). The assay protein will carry
  fuller, heterogeneous glycans, so the 1.5% glycan figure is a floor, not an estimate.
- **Residue 121 is absent locally**, so the measured construct is 639 of the official 640 aa.
- The 2.0 Å clash screens, iPTM values and L6 values elsewhere in this repo were all computed in
  the surrogate context and are **not** corrected by this document.

## Reproduce

```sh
cd ~/dev/adaptyv-c4
date -u '+%Y-%m-%dT%H:%M:%SZ'
OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
  python3 scripts/c4_native_context_gap_20261026.py --out results/c4_native_reproduce_20261026
# exits non-zero unless all four controls fire:
#   CTX1_MATCHES_SEQ, GLYCAN_PRESENT, BURIAL_MONOTONIC, SELF_CONSISTENT
```

Controls, this run: 195/195 apical residues at 189–383; 6 NAG found; **zero** monotonicity
violations (accessibility never increased as context was added, across all 195 positions); zero
self-consistency violations on a full CTX4 recomputation.
