# C5 scoring strategy — four-model advisory panel

Date: 2026-10-27
Context: AF2-Multimer positive-control validation FAILED (ipTM G12D = WT = 0.84).
Question: which scoring strategy for the C5 design campaign?

## Panel

| Model | Source | Method |
|---|---|---|
| Claude Fable 5.1 | Local (Agent) | Measured the PDB structures, ran SASA/contact analysis |
| Claude Opus 5.5 | Local (Agent) | Reasoning from the data |
| OpenAI GPT-6 Astra | OpenRouter | Reasoning from the data |
| Moonshot Kimi K3.1 Pro | OpenRouter | Deep reasoning (879 prompt / 10180 completion tokens) |

---

## Fable 5.1 — the structural bombshell

**Recommendation: Option E (new) + reframed D, NOT option D as written.**

Fable actually measured the structures and found Option D is **geometrically impossible**:

- Asp6 side-chain SASA = **8.8 Å² (10.4% relative)** — ~90% buried in the HLA groove
- Asp6 OD2 SASA = **0.00 Å²** — carboxylate completely inaccessible
- Asp6 is already salt-bridged to HLA **Arg114** + contacts to Gln70, Gln156, Trp147
- Nearest TCR atom to Asp6 carboxylate = **7.5 Å** — the TCR doesn't touch it
- TCR contacts at P6 are all via TCRβ **His100** to **backbone** atoms only

**Key insight:** "Design for direct Asp contacts" is not an available strategy. The side
chain is buried. The only discrimination mechanism is the buried charge perturbing the
electrostatic field above the groove — exactly what Poole et al. describe.

**Proposed Option E: target-side electrostatic map ($0, CPU, day 1)**

Run APBS (Poisson-Boltzmann) on state1 vs state2, difference the electrostatic potential
at the solvent-accessible surface. If the Asp→Gly perturbation is measurable at the
surface, the difference map IS the design map — it tells you where to place complementary
charge. If it isn't measurable, no in-silico method will discriminate.

**Second insight: C4 pH machinery DOES transfer.** state.md says "C5 has no pH axis" but
the discrimination mechanism is protonation-coupled electrostatics through His100. propka
pKa shifts are exactly the right readout. Fix the state.md line.

**Design strategy:** Place the binder over the P6-P8 bulge and Arg114/Gln70/Gln156/Trp147
microenvironment, with a cation or His ~7-8 Å above P6 (where TCRβ His100 sits). RFD4
hotspot conditioning on those residues, not on Asp6.

**Risks:** The PB signal may be smaller than method error. State2 is an in-silico D→G
mutant with G12D-context rotamers — may need relaxation or a real WT structure.

---

## Opus 5.5

**Recommendation: D as primary + A as lightweight validation.**

- Design binders with DIRECT polar contacts to the Asp carboxylate
  (Arg/Lys salt bridges, H-bonds from Asn/Ser/Thr/Tyr)
- Use AF2 for fold quality only
- FoldX validation on JDI TCR as day-1 experiment — if it discriminates, use as
  supplementary filter; if not, rely on design strategy alone
- Manual inspection of all 20 designs in PyMOL (feasible at this scale)

**Rejects:** B (can't calibrate with one control), C (structure prediction isn't the
problem — electrostatics is)

**Missing option:** Rosetta design with explicit interaction-energy optimization for
Asp contact. Also: manual PyMOL inspection for 20 designs.

**Note:** Opus did NOT measure the structures and assumes direct Asp contacts are
possible. Fable's measurement shows they are not — Opus's recommendation needs revision.

---

## GPT-6 Astra

**Recommendation: D + bounded A. Reject B. Skip C.**

"Spend the remaining time designing a credible reason for WT not to bind, rather
than searching for a predictor that reports WT nonbinding."

Key points:
- ipTM is structural confidence, not affinity prediction
- Don't chase composite scores with arbitrary weights on uncalibrated quantities
- Calibrate any physics term on the control in ONE working day — if it fails, delete it
- Use gates and a shortlist, not a magic scalar:
  1. Basic viability
  2. Credible target interface
  3. Asp dependence mechanism
  4. Supporting energetic evidence (if calibrated)
  5. Portfolio diversity
- Missing option: explicit negative design against WT and HLA-driven binding
- Stress-test finalists with a small panel of other HLA-A*11:01 peptides

**Budget allocation:**
- $20 validation benchmarks
- $120 candidate screening
- $40 finalist checks
- $20 reserve

**Note:** Like Opus, Astra assumes direct Asp contacts are available. Fable's
measurement changes this — but Astra's gating framework and negative-design
emphasis remain valid.

---

## Kimi K3.1 Pro

**Recommendation: D + narrowly scoped A. Reject B and C.**

"Stop trying to make ipTM certify specificity."

Key points:
- Use **sequential gates**, not weighted composites:
  1. AF2 for structural plausibility (fold, interface PAE, model agreement)
  2. Hard mechanism filters (direct Asp contact, peptide contribution to BSA,
     no alternative HLA-only binding pose, loss of contacts in WT)
  3. Physics as veto + secondary ranking (FoldX broadly, Rosetta on final set)
- ipTM is pose quality, not affinity. A high-ipTM WT model doesn't mean WT binds.
- FoldX D→G ΔΔG is NOT automatically a binding-selectivity ΔΔG (contains
  complex-stability effects). Treat carefully.
- "Repack and minimally relax WT. Scoring a frozen structure after deleting the
  Asp side chain can exaggerate the penalty."

**Missing option: explicit multistate negative design during generation.**
Lightweight version: design on G12D → repack/rescore WT and HLA-only →
mutate/resequence failures → repeat once. Moves selectivity upstream without
new infrastructure.

**Critical risk #1: "Asp accessibility or approach geometry may be unfavorable.
Check this in 7OW6 immediately."** — This is exactly what Fable proved: the
carboxylate is 90% buried, SASA = 0.00 Å². Kimi identified the risk; Fable
measured it.

**Risk #7: "Overconcentration on one solution."** If the shared scaffold or
Asp-recognition motif fails, 20 sequence variants may all fail together.

12-day plan: Day 1 freeze problem + half-day FoldX test on JDI (stop if it
doesn't work). Days 2-5 generate + triage (thousands → ~100 AF2 → ~40
finalists). Days 6-8 counter-screen WT. Days 9-10 select diverse 20.
Days 11-12 QC + submission.

---

## Convergence (all four models)

1. **Drop ipTM as a selectivity metric.** It's structural confidence, not binding
   free energy. The FAIL is expected and fundamental, not fixable.

2. **Option B (composite scoring) is rejected.** Arbitrary weights on uncalibrated
   quantities with one control point = false precision.

3. **Option C (alternative predictors) is deprioritized.** AF2's insensitivity to
   point mutations is well-established and not an MSA-depth problem.

4. **Design for selectivity, don't score for it.** Engineer the specificity into the
   binder's mechanism rather than hoping a confidence score discovers it.

5. **Physics (A) only as bounded validation.** Give it one day on the control.
   If it works, use as supplementary evidence. If not, delete it.

6. **AF2 keeps its role** for fold quality (pLDDT), structural plausibility (pTM),
   and binding mode confidence (ipTM for G12D only — not differential).

## But Fable changes everything

Fable's structural measurement — that Asp6 is 90% buried and the carboxylate is
completely inaccessible — invalidates the design premise that Opus and Astra share.
You cannot design hydrogen bonds to something no protein can reach.

**The revised consensus is:**

1. **Run the APBS electrostatic difference map first (day 1, CPU, $0).** This is the
   experiment that determines whether any in-silico method can see the G12D→WT
   difference at the binder-accessible surface.

2. **If APBS shows a measurable surface perturbation:** design binders that complement
   the G12D-specific field (cation/His ~7-8 Å above P6, over the groove bulge).
   Use the APBS difference map as the design map.

3. **If APBS shows no surface perturbation:** specificity is computationally invisible.
   Accept this, design for binding quality + portfolio diversity, and let the
   experimental screen sort selectivity.

4. **C4 propka/His-scan transfers.** The pKa-shift readout matches the mechanism.
   Correct state.md.

5. **State2 needs relaxation.** It's an in-silico D→G mutant with G12D rotamers.

## Decision for Ricardo

**Immediate next step (no spend):** Run APBS on state1 vs state2. This is the
$0 experiment that sizes the rest of the campaign.

**Then:** the 12-day plan follows from that result.
