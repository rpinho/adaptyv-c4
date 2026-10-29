# FINDING: mouse was never the blocker — the real epitope is well conserved, and all three of its histidines are natively inaccessible

Written 2026-10-26T14:56:01Z from the real clock (`date -u`). CPU only, **$0**, read-only on every
input. No GPU, no Modal, no RunPod. No submission, no upload, no organiser contact. No ranking
built. No frozen file touched (4/4 verify). All **seven** controls fired.

This closes the gap `FINDING_WRONG_DOMAIN.md` left open in its own words:

> "Mouse conservation of the real epitope, author 651-741 across both subunits, has not been
> measured by anyone."

It is now measured, and it changes the strategic answer.

---

## Headline

**Two results that point in opposite directions, and the second is the binding one.**

1. **Mouse cross-reactivity is not the obstacle for this design pool.** The real epitope
   (author 651-741) is mouse-identical at **76/91 = 83.5%**, against **74.9%** for the apical
   domain. Restricted to positions actually presented in the native homodimer the gap widens
   sharply: **80.0% (32/40)** at the real epitope versus **56.2%** at apical. **All ten named
   generator hotspots are 100% mouse-identical**, His707 included. Mouse was a real problem for
   the apical surrogate. It is not a problem for the surface these designs were actually aimed at.

2. **The pH-switch geometry does not exist natively.** All three histidines in the band are
   natively inaccessible: **H699 is completely buried (0.00 A^2)**, H684 retains 21.8 A^2
   (rel. SASA 0.097), and **H707 — the residue the whole campaign is named after — presents only
   28.5 A^2 (rel. SASA 0.127) in every context, including the band alone.** The four hotspots that
   *are* natively exposed are **all carboxylates** (D659, E694, E714, E728) — precisely the
   acid-ON anchors `fa98723` identified as the wrong pH direction.

So the campaign's own anchor residue was barely solvent-exposed from the start, and the only
natively available anchors drive the switch the wrong way. That is upstream of the construct
choice and upstream of the gate direction: **no construct and no reranking reaches it.**

---

## 1. Mouse conservation, by band

Full-length global BLOSUM62 alignment of human P02786 (760 aa) against mouse Q62351 (763 aa),
open -10 / extend -0.5 — the same parameters as `results/apical_context/`. **6 equally-best
alignments**, score 3100.5; all six were enumerated and a position counts as conserved only if it
is identical in **all** of them.

| band | identical | n | identity | gapped in some optimum |
|---|---:|---:|---:|---:|
| **real design epitope 651-741** | **76** | 91 | **0.8352** | 0 |
| partner band 676-696 (chain B in both contigs) | 18 | 21 | 0.8571 | 0 |
| apical 189-383 (what everything was scored on) | 146 | 195 | 0.7487 | 3 |
| official target 121-760 | 499 | 640 | 0.7797 | 3 |

The 15 substituted positions in the real band: 661, 667, 672, 673, 674, 678, 680, 693, 700, 710,
713, 720, 721, 723, 724. None is a named hotspot.

### Conservation restricted to what the native assembly presents

The apical cycle's lesson was that burial inflates apparent conservation, so the same split is
applied here:

| | n | identical | identity |
|---|---:|---:|---:|
| natively exposed (rel. SASA >= 0.20 in dimer + glycan) | 40 | 32 | **0.8000** |
| natively buried | 51 | 44 | 0.8627 |
| whole band | 91 | 76 | 0.8352 |

The burial artifact is present but small (0.800 exposed vs 0.863 buried). At apical the same split
collapsed from 0.845 to **0.562**. The real epitope's conservation is genuine surface conservation,
not a burial illusion.

---

## 2. The named hotspots: all conserved, mostly inaccessible

Per-context SASA, A^2, from identical 11QC chain A coordinates:

| pos | aa | mouse | band alone | monomer | dimer | +NAG | rel. native | exposed? | role |
|---:|---|---|---:|---:|---:|---:|---:|---|---|
| 659 | D | **D** | 54.1 | 54.1 | 54.1 | 54.1 | 0.280 | **yes** | carboxylate |
| 682 | E | **E** | 105.0 | 2.3 | 2.3 | 2.3 | 0.010 | no | carboxylate |
| 685 | F | **F** | 62.8 | 1.2 | 1.2 | 1.2 | 0.005 | no | groove |
| 694 | E | **E** | 157.3 | 152.0 | 152.0 | 152.0 | 0.682 | **yes** | carboxylate |
| 701 | F | **F** | 91.5 | 20.0 | 20.0 | 20.0 | 0.083 | no | groove |
| 702 | W | **W** | 195.9 | 35.8 | 35.8 | 35.8 | 0.126 | no | groove |
| **707** | **H** | **H** | **28.5** | **28.5** | **28.5** | **28.5** | **0.127** | **no** | groove anchor |
| 709 | L | **L** | 65.1 | 13.3 | 13.3 | 13.3 | 0.066 | no | groove |
| 714 | E | **E** | 114.5 | 114.5 | 114.5 | 114.5 | 0.513 | **yes** | carboxylate |
| 728 | E | **E** | 92.0 | 77.8 | 77.8 | 77.8 | 0.349 | **yes** | carboxylate |

**10/10 mouse-identical. 4/10 natively exposed, and all four are carboxylates.**

Two structural facts fall straight out of this table:

- **The partner subunit occludes none of them.** `lost_to_partner_subunit` is **0.0 A^2 at all ten
  hotspots.** Both generator contigs included chain B 676-696 and the campaign was described as a
  two-subunit groove, but at the residues it actually named, the second subunit contributes
  nothing. The "composite two-chain site" is not supported at the hotspots.
- **H707 is not occluded by anything — it is intrinsically barely exposed.** 28.5 A^2 in the band
  alone, 28.5 A^2 in the full homodimer with glycans. The construct is irrelevant to it. The five
  other groove hotspots (682, 685, 701, 702, 709) lose 50-160 A^2 each, and every A^2 of that loss
  is to **their own subunit** — the crop they were generated against presented them as exposed
  (rel. SASA 0.26-0.69) when the real monomer does not.

### All three histidines in the band

| pos | band alone | native | rel. native | exposed? | mouse |
|---:|---:|---:|---:|---|---|
| 684 | 100.5 | 21.8 | 0.097 | no | H (conserved) |
| **699** | 21.6 | **0.0** | 0.000 | no | H (conserved) |
| **707** | 28.5 | 28.5 | 0.127 | no | H (conserved) |

C4 asks for binding at pH 7.4 and none at pH 6.0. Histidine is the obvious residue to build that
on, and it is what this campaign was named for. **At the real epitope, receptor-side histidine is
not available in any construct.**

---

## 3. Native accessibility of the band

Four nested contexts, identical coordinates, context the only variable:

| context | band SASA (A^2) | |
|---|---:|---|
| CTX1 band 651-741 alone | 7727.0 | 100% |
| CTX2 full 639-aa monomer | 5082.3 | 65.8% |
| CTX3 homodimer | 4087.6 | 52.9% |
| CTX4 homodimer + 6 NAG | 3985.9 | **51.6%** |

| lost to | A^2 | % of band |
|---|---:|---:|
| **the rest of its own subunit** | **2644.8** | **34.2%** |
| the partner subunit | 994.7 | 12.9% |
| observed glycans | 101.7 | 1.3% |

| quantity | value |
|---|---:|
| band residues resolved in 11QC chain A | 91 of 91 |
| exposed in the band alone | 74 |
| exposed in the native assembly | **40** |
| exposed in the crop but not natively | **34** |
| materially occluded (>=50% area loss) | 42 |

Same shape as the apical result: **intra-subunit burial dominates** (34.2%) over the partner
(12.9%) and glycans (1.3%). This is now measured at both epitopes, which strengthens the
construct-scope recommendation rather than merely repeating it — see §6.

---

## 4. Per design

Each design's own contact footprint, already mapped to author numbering in
`results/c4_provenance_join_20261026/hotspot_contact.json`, joined by design and scored on the
real band:

| quantity | all 343 | the v4 acid-OFF 20 seats |
|---|---:|---:|
| designs with >=1 real-band contact | 339 / 343 | 20 / 20 |
| real-band contacts, median | 17 | 11.5 |
| real-band contacts that are natively exposed, median | **5** | **3.5** |
| designs with >=1 natively-exposed real-band contact | 317 / 343 | 17 / 20 |
| mouse identity of real-band footprint, median | **0.857** | 0.875 |
| mouse identity of natively-exposed footprint, median | **0.750** | 0.833 |
| **designs contacting H707** | **19 / 343** | **0 / 20** |

The footprints are well conserved in mouse — median 85.7%, and 90th percentile 1.000. But the
median design makes only **5** contacts on natively-exposed real-band surface, and **not one of
the 20 direction-correct seats touches H707 at all.**

---

## 5. What this does and does not change

### It overturns the pessimistic mouse reading — for the right reason

Task (A) of the brief asked whether a single fixed scaffold family binding one human patch can be
mouse cross-reactive when the patch differs at 6 of 9 residues. The answer given earlier in this
cycle (apical: 56.2% identity among natively exposed positions) was measured on **a domain these
designs were not generated against**. At the surface they *were* generated against, the number is
**80.0%**, with 10/10 named hotspots and 3/3 histidines mouse-identical. **Mouse cross-reactivity
is not what blocks this pool.** The earlier conclusion was not wrong about apical; it was answering
the question about the wrong surface.

### It does not rescue the pool

The pool is still not submittable, and this finding adds a reason that is *harder* to fix than the
previous ones:

| defect | fixable by | status |
|---|---|---|
| inverted pH gate | re-gating (done, rule v4) | fixed, bought no measurable quality |
| wrong evaluation construct | re-running on 640-aa monomer (costs GPU) | open, Ricardo's call |
| **no natively accessible histidine at the real epitope** | **nothing — it is receptor geometry** | **not fixable** |
| acid-ON baked into generation (`C_SEQ_PROB_HIS: 0.15`, carboxylate hotspots) | a new campaign | not fixable in this pool |

### The constructive part: where a redesign would aim

The 32 positions in 651-741 that are **both natively exposed and mouse-identical** are a real,
measured asset, and they are the only defensible target list this band offers:

    651R 654S 655R 657T 658T 659D 662N 664E 665K 666T 668R 670V 676D 688P
    694E 695S 696P 704S 705G 706S 714E 717K 718L 719R 722N 725A 728E 729T
    730L 732R 733N 737L

Rich in D/E/K/R/N/S/T — a polar, charged, well-conserved patch. Note what is **absent**: no
histidine. Any pH switch built here must come from the **binder** side, not the receptor.

---

## 6. Effect on the open construct-scope decision

`results/decision_requests/20261026-c4-construct-scope.md` asks A (stay on the 195-aa surrogate),
B (640-aa monomer) or C (official 640-aa homodimer + glycans), and I recommended B. This
measurement sharpens that materially, in **both** directions:

- **B over C is now strongly supported, not merely suspected.** At the real epitope the partner
  subunit removes 994.7 A^2 band-wide and **0.0 A^2 at every one of the ten named hotspots**;
  glycans remove 1.3%. B captures the entire hotspot-level correction. C's extra cost buys a term
  measured at zero where it matters.
- **But the upside of either is now lower than I estimated when I recommended B.** I argued B
  "makes future measurements honest", which is still true. What this finding adds is that an honest
  measurement of *this pool* at *this epitope* has a known ceiling: no accessible receptor
  histidine, only carboxylate anchors available, 0/20 seats touching the anchor residue, and a
  median of 5 natively-exposed footprint contacts. Spending $37-56 to measure that more precisely
  is a worse buy than it looked an hour ago.

**I am not changing the recommendation on my own authority** — the construct choice is reserved for
Ricardo and the request is already filed. I am recording that the evidence behind it moved, and
that "none yet — keep doing $0 work" has become more defensible than when the request was written.

---

## 7. What is NOT claimed

- **Nothing here is a binding, affinity, KD or pH-response measurement.** It is solvent
  accessibility on experimental coordinates plus sequence identity from one alignment family.
- **"Natively inaccessible" is not "cannot participate".** H707's 28.5 A^2 is small, not zero, and
  rel. SASA 0.20 is a convention, not a threshold with physical meaning. A partially buried residue
  can still contribute to a binding site, and receptors flex — 11QC is one crystal form and the
  domain is treated as rigid here.
- **A pH switch does not require receptor histidine.** Per Astra's correction earlier in this
  cycle, buried glutamate pKa values span roughly 5.2-9.4 and PROPKA assigns environment-specific
  pKa, so Asp/Glu *can* titrate between 7.4 and 6.0; binder-side histidine is also available. The
  claim here is narrower and stands: the *receptor-side histidine* route the campaign was built on
  is not available at this epitope.
- **Mouse identity is not a binding probability.** 83.5% identity does not predict
  cross-reactivity, and the three substitutions inside the partner band (678, 680, 693) are
  unweighted by contact.
- **The alignment is one family of six optima.** Conservation is reported only where all six agree;
  no position in the real band is gapped in any optimum, so the band is unambiguous, but the
  apical band has 3 gapped positions.
- **Glycans are only the 6 observed in 11QC**, so the 1.3% figure is a floor.
- **No design was re-scored against the real target.** This measures whether the real epitope is
  *available* and *conserved*, not whether any binder engages it. That experiment is still the
  missing one.
- **Contact footprints come from complexes predicted against the apical construct.** They are used
  here only to locate which author-numbered receptor positions each binder reached in the file it
  was folded in; they are not evidence of a native pose.

---

## 8. Controls — all seven fired

| control | requirement | measured | fired |
|---|---|---|---|
| ALIGN_SELF | P02786 vs itself identical everywhere | 760/760 | yes |
| **NUMBERING_11QC** | 11QC chain A author numbering == P02786 numbering | **639 checked, 0 mismatched** | yes |
| HOTSPOT_IDENTITY | 707 HIS, 659 ASP, 682/694/714/728 GLU | 6/6 | yes |
| **APICAL_REGRESSION** | apical identity reproduces the prior cycle's 146/195 | **146** | yes |
| NEG_WINDOW | a 200-residue-shifted window differs | 0.8352 vs 0.8462 | yes |
| BURIAL_MONOTONIC | accessibility never increases as context is added | 0 violations / 91 positions | yes |
| GLYCAN_PRESENT | 6 NAG in 11QC | 6 | yes |

Two notes on control strength, so they are not read as stronger than they are:

- **NUMBERING_11QC is the load-bearing one.** Every band number in this document depends on author
  numbering being UniProt numbering, and it is confirmed residue-by-residue at 639/639 with zero
  mismatches. Combined with HOTSPOT_IDENTITY (6/6) this also re-derives the root-cause finding's
  numbering conclusion from a different direction.
- **NEG_WINDOW is weak.** It fires — the shifted window gives a different number — but only just
  (0.8352 vs 0.8462), because TfR1 is fairly uniformly conserved. It establishes that the mapping
  is position-sensitive; it does not establish much else, and it should not be cited as validating
  the band numbers. NUMBERING_11QC does that.

## 9. Artifacts

| file | sha256 |
|---|---|
| `real_epitope_mouse.json` | `619620b4b44560505dfd066961f3423e41023c21370e47ca4fae619502fd318a` |
| `real_epitope_per_position.csv` (91 rows) | `07e1d8c8e45bef003ab90ed423902b130e8695d4fb353956620fa6ca287ef1b2` |
| `real_epitope_per_design.csv` (343 rows) | `50039e31de1b891403be7a23aa9986678e545747ce2425d47aa8d8c89bb8fc76` |

Script: `scripts/c4_real_epitope_mouse_20261026.py`

### Reproduce

```sh
cd ~/dev/adaptyv-c4
date -u '+%Y-%m-%dT%H:%M:%SZ'
# controls only, writes nothing:
python3 scripts/c4_real_epitope_mouse_20261026.py --self-test --out /tmp/unused
# full run (~3 min); refuses to overwrite and exits non-zero unless all 7 controls fire:
OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
  python3 scripts/c4_real_epitope_mouse_20261026.py --out results/c4_real_epitope_reproduce_20261026
```
