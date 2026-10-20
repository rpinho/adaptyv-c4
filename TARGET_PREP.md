# Target prep and epitope choice (2026-10-20)

Project root for the Challenge 4 entry. Deliberately **outside the KB**: nothing here goes onto
`main`. The entry plan document lives on the `adaptyv-c4` branch of the KB worktree only.

Scripts: `scripts/align_tfr1.py`, `scripts/epitope_scan2.py`, `scripts/hotspots.py`,
`scripts/his_patches.py`. Results: `results/*.json`. Structures and UniProt records: `data/`.

---

## 1. Sequence facts (computed here)

Global Needleman-Wunsch, UniProt **P02786** (human, 760 aa) vs **Q62351** (mouse, 763 aa):

- 77.5% identity overall, 76.5% over the ectodomain 89-760, **78.0% over the assay construct 121-760**
- fully conserved 20-aa windows include 461-486 and 725-744
- least conserved: 352-374 (35-40%) and 568-589 (40%)
- native histidines in 121-760: 186, 234, 300, 302, 318, 401, 475, 515, 603, 684, 699, 707

## 2. Structural basis (upgraded; v1 used the wrong structures)

A direct RCSB search API query returned **24 human TfR1 entries**. That changed the structural
basis materially:

| Use | PDB | Method | What it is |
|---|---|---|---|
| **Design template** | **11QC** | cryo-EM **2.4 A** | TfR1 ectodomain homodimer, chains A+B, residues **122-760** = exactly the assay construct |
| Transferrin footprint | **3S9L** | X-ray 3.22 A | TfR1-transferrin complex |
| HFE footprint | **1DE4** | X-ray 2.8 A | HFE/beta2m-TfR1 complex |
| Designed-binder epitopes | **6WRV, 6WRW, 6WRX** | X-ray 2.47-3.07 A | **computationally designed de novo binders** (3DS18, 2DS25.5, 2DS25.1) bound to TfR1 |
| Viral epitope | **3KAS** | X-ray 2.4 A | Machupo virus GP1 on TfR1 |
| Peptide epitope | **9GH7** | X-ray 2.08 A | bicyclic peptide on TfR1 |
| Noted, not used | 6W3H, 6OKD, 6Y76, 8P0Z | - | Fc brain-delivery construct and apical-domain redesigns; their numbering is offset or they contain no full TfR1, so no footprint was taken |

v1 had used 1CX8 (3.2 A) as the template and 1SUV (7.5 A) for the transferrin footprint. Both are
superseded above. **11QC is also the structure the official challenge page cites.**

Footprints computed at a 5 A heavy-atom cutoff, in TfR1 numbering:

- **transferrin (3S9L), 37 residues:** 121,123,125,195,196,197,199,201,211,212,215,292,294,619,622,623,626,629,640,643,644,646,647,648,650,651,658,661,...
- **HFE (1DE4), 27 residues:** 318,527,615,616,619,620,622,623,626,629,630,638,640,641,643,644,646,648,650,651,653,654,655,657,658,661,662
- **designed binders (6WRV/6WRW/6WRX), 23-38 residues each:** all three converge on **197-215 + 294 + 340-348 + 371-376**
- **Machupo GP1 (3KAS):** 197-217 + 293,294 + 339-350 + 370,371 - the same apical cluster
- **bicyclic peptide (9GH7):** a separate site, 150-161 + 190-195 + 219 + 382,384 + 412,414 + 458 + 569

N-glycosylation (NAG-linked ASN): **251, 317, 727**. Tag vicinity: 29 residues within 15 A of
residue 760. Dimer interface is handled by computing SASA in the A+B dimer, not the monomer.

## 3. The strategic finding

**The site that is proven designable is not the site that is conserved.** Measured on the same
structure:

| Epitope set | n | human/mouse identity | mean rel SASA | hydrophobics | native His |
|---|---|---|---|---|---|
| Apical consensus of the three designed binders, minus transferrin/HFE/glycan/tag | 17 | **0.59** | 0.304 | 6 | - |
| Machupo-only apical residues | 7 | 0.71 | 0.277 | 1 | - |
| Bicyclic-peptide site | 21 | 0.71 | 0.243 | 8 | - |
| **His707 patch (688,690,695,704,705,706,711,733,737)** | 9-10 | **1.00** | 0.299 | 3-4 | **707** |

So the apical hotspot where every known designed binder, the virus and the engineered BBB shuttles
all bind is **human-specific** (0.59 identity) — which is consistent with the known biology that
human-TfR1 BBB shuttles do not cross-react with mouse TfR1. Criterion 2 of the challenge rewards
mouse cross-reactivity; criterion 3 rewards affinity and a functional epitope. **Those two pull in
opposite directions on this target**, and the ranking is lexicographic after the pH gate, so the
entry should not bet everything on one epitope.

## 4. Epitope portfolio for the 20 designs

| ID | Epitope | Residues (BindCraft hotspots, chain A) | Why |
|---|---|---|---|
| **E1** | **His707 patch** | 688,690,695,704,705,706,711,714,730,733,737 | **Primary.** 100% human/mouse identity (criterion 2), carries a **mouse-conserved native His707** so the target supplies half the proton linkage (mechanism family 6), 19.5 A from both the transferrin and HFE footprints, and **no overlap with any known binder site in the PDB** (novelty). |
| **E2** | Apical, de-risked | 202,203,204,205,206,207,208,209,210,213 and 340,343,344,348,371,374,376 | Proven designable (three independent designed binders plus a virus bind here) and outside the transferrin footprint; likely the functional/transcytosis epitope class. Costs criterion 2 (identity 0.59). Sub-cluster it: 202-213 and 340-376 are 27 A apart, too wide for one interface. |
| **E3** | His475 patch | 471,472,476,477 | Second conserved native-His option (patch identity 1.00, His475 conserved). Small and only 2 hydrophobics, so lower designability. |
| **E4** | His603 patch | 122,124,443,600,603,604,605,608 | Third native-His option, identity 0.88, 3 hydrophobics. |
| E5 | Bicyclic-peptide site | 150-161,190-195,219,382,384,412,414,458,569 | Fallback: ligandable by a real molecule, 8 hydrophobics, identity 0.71. |

Native-His screen (all 12 histidines in the construct, `results/his_patches.json`): **His707 is the
best option** (exposed 0.142, conserved, patch identity 1.00); His475 (1.00) and His603 (0.88)
follow; **His318 is unusable — it sits in the HFE footprint**; His186, 234, 302, 401, 684, 699 are
buried or poorly conserved.

## 5. What is still weak

1. Domain boundaries were never taken from the Lawrence 1999 paper; everything above is by residue
   number and footprint, which is sufficient for hotspot choice but not for domain claims.
2. E2's span (27 A) must be split into two hotspot sets before use.
3. The scan's ranking score is hand-weighted. It shortlists; it is not evidence.
4. 11QC is cryo-EM at 2.4 A: good, but side-chain placement at the periphery is still approximate.
5. **Compliance note for the write-up:** 6WRV/6WRW/6WRX are existing designed binders. Their
   *epitope* informs E2; their *scaffolds and sequences* must not. Plan a novelty audit that scores
   every final design against all known TfR1 binders (these three, the bicyclic peptide, HFE,
   transferrin, Machupo GP1) and reports the maximum identity, so the de novo/zero-shot claim is
   demonstrated rather than asserted.

---

## 6. Dimer-occlusion check: the portfolio changes again (2026-10-20, later)

Relative SASA computed in the dimer is not enough. A residue can look exposed and still be
unreachable, because the partner subunit sits over the approach vector. `scripts/trim_target.py`
reported the epitope-to-chain-B closest approach, and `scripts/check_dimer_occlusion.py` then
measured, per residue, the SASA gained when chain B is deleted (how much the partner covers) and
the distance to chain B. Results (`results/dimer_occlusion.json`):

| Epitope | closest approach to chain B | verdict |
|---|---|---|
| E1 as first defined (688-737) | **3.3 A** | **partly occluded** - 688 (0.33 of its surface), 733 (0.16) and 737 (0.19) only appear on the isolated monomer. Drop those three. |
| **E1' refined: 695,704,705,706,707,711,714 (+730)** | 7.2 A | **usable.** Zero occlusion on every member; His707 is 9.6 A from chain B with rel SASA 0.142 and no occlusion. |
| E3 (His475: 471,472,475,476,477) | 3.7-5.1 A | **dropped.** 472 loses half its surface to the partner (0.496) and His475 itself is 4.7 A from chain B. |
| **E4 (His603: 600,603,604,605,608)** | 8.7 A | **usable.** Zero occlusion on all five, all mouse-identical. |
| E2a (202-213) | 39.2 A | irrelevant - nowhere near the dimer interface. |
| E2b (340-376) | 30.4 A | irrelevant. |

**The second honest tension: the conserved native-His patches are polar.** E1' reads S,S,G,S,H,A,E
and E4 reads K,H,D,V,N — one hydrophobic residue each. AF2 is documented to design hydrophilic
interfaces worse than hydrophobic ones, and the RFdiffusion guidance is to prefer sites with three
or more hydrophobics and avoid charged polar sites. The apical E2 sites carry 6 hydrophobics and no
dimer problem, but only 0.59 identity. So:

- **designability** favours E2 (apical, human-specific)
- **the scoring criteria** favour E1'/E4 (conserved, native His, novel)

That maps cleanly onto the mechanism families instead of forcing a single choice:

| Mechanism families | Epitope | Rationale |
|---|---|---|
| 1-5 (the switch lives in the **binder**: His-Asp dyad, His-His pair, buried unpaired His, scaffold fold-switch, charge ladder) | **E2a (202-213)** and **E2b (340-376)** | These families only need a designable interface; the apical sites are the ones the generators can actually hit, and they are the proven-designable, likely-functional epitope class. |
| 6 (the switch is the **target's own** His) | **E1' (His707)** and **E4 (His603)** | These families require a native, mouse-conserved histidine at the interface. Expect a lower pass rate from the generators and budget more trajectories. |

Trimmed, generator-ready targets are in `targets/` (template 11QC, 12 A neighbourhood, contiguous
fragments, hydrogens stripped):
`E1_his707_clean.pdb`, `E2a_apical_202_213.pdb`, `E2b_apical_340_376.pdb`, `E4_his603.pdb`.
E1' keeps 21 chain-B residues in the trim, so **that run must pass both chains as the target** or
the design will grow into space the partner subunit occupies.

## 7. Baseline reproduction (the prior-art-first step)

Pod `yy6bbyxlp7vwi7` (A100 80GB PCIe, CA-MTL-3, $1.79/hr), volume `adaptyv-c4` (100 GB) at
`/workspace`. BindCraft installed from source with mamba (the conda solve stalled for 25 minutes and
was abandoned); env verified: jax 0.6.0 on `CudaDevice(id=0)`, colabdesign 1.1.3, PyRosetta
2026.29, biopython 1.88. The AF2 weight extraction "failed" only because `tar` could not chown on
the network volume; all 15 `.npz` files are present and the 5.6 GB tar was removed afterwards.

Smoke test running: the shipped PDL1 example, `default_4stage_multimer` settings with
`default_filters`, one accepted design requested. Stage 1 optimisation confirmed on GPU.

**Licensing note for the write-up:** PyRosetta prints that commercial use requires a licence. This
entry is non-commercial, and the competition terms leave commercialisation rights with the
designer — so if anything here were ever taken commercial, the PyRosetta licence would have to be
settled first. BindCraft itself carries the same condition.
