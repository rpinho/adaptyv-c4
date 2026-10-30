# C5 Research Ideas — Outside-Input Sweep 2026-10-30

## IDEA 1: Use upward-facing peptide residues as hotspots (Baker lab method)
- **Source:** Baker lab, Science Jul 2025 — [Design of high-specificity binders for peptide-MHC-I complexes](https://www.bakerlab.org/wp-content/uploads/2025/08/science.adv0185.pdf)
- **Insight:** RFdiffusion conditioned on upward-facing peptide residues (pointing out of the MHC groove) generates binders that contact the peptide, not the MHC. This is the structural basis for peptide specificity. "Only a small fraction of trajectories yield designs that contact the peptide while avoiding the MHC" — so high selectivity hit rate requires correct hotspot targeting.
- **Action:** Identify which residues in VVVGADGVGK are upward-facing in PDB 7OW6 (or 9UV8). Our current hotspots (C4-C8) may include groove-buried residues that drive MHC binding, not peptide specificity. Refine hotspots to upward-facing only.
- **Status:** OPEN

## IDEA 2: Exploit the aspartate salt bridge for G12D specificity
- **Source:** ICML 2026 — [LLM-Assisted vs Agentic De Novo Minibinder Design for KRAS G12D](https://icml.cc/virtual/2026/70812)
- **Insight:** The agent recovered "the aspartate-anchored salt bridge linked to G12D specificity." The D at position 6 (the G12D mutation itself) forms a salt bridge that WT glycine cannot. Binders that form hydrogen bonds or salt bridges specifically to D6 should be inherently selective.
- **Action:** Design binders with explicit D6 contact — narrow hotspots to C6 only (the mutant residue), or design with positively charged residues (Arg/Lys) at the interface to form salt bridges with the negatively charged aspartate.
- **Status:** OPEN

## IDEA 3: Use PDB 9UV8 (newer G12D/HLA-A*11:01 structure)
- **Source:** ICML 2026 paper uses PDB 9UV8 as target — [Lacuna listing](https://lacuna.tiptreesystems.com/work/llm-assisted-versus-agentic-approaches-to-de-novo-minibinder-design-for-a-kras/wrk_9f2db8b2e2620a4af5599c8dc92b0192)
- **Checked:** 9UV8 is 1.79 Å crystal structure BUT has the **9-mer** peptide (VVGADGVGK), NOT our 10-mer (VVVGADGVGK). Challenge target is 7OW6 (10-mer). Different peptide register = different D position in groove.
- **Action:** Keep 10-mer target per challenge rules. 9UV8 is useful as reference for understanding D exposure, not as a design target.
- **Status:** INFORMATIONAL — no design action needed

## IDEA 4: 9-mer vs 10-mer peptide register affects immunogenicity
- **Source:** Comms Biology 2026 — [Structure guided analysis of KRAS G12 mutants in HLA-A*11:01](https://www.nature.com/articles/s42003-025-09285-0)
- **Insight:** G12D has a "length-encoded immunogenic advantage" — the 9-mer (VVGADGVGK) and 10-mer (VVVGADGVGK) registers have different structures and immunogenicity profiles. G12D's buried aspartate destabilizes the complex.
- **Action:** Check whether our peptide (VVVGADGVGK, 10-mer) vs the 9-mer (VVGADGVGK) matters for design. The ICML paper used the 9-mer. Consider designing binders against both registers.
- **Status:** OPEN

## IDEA 5: Adaptyv binder-prescreen v0.5 known-binder benchmark
- **Source:** [adaptyvbio/binder-prescreen v0.5](https://github.com/adaptyvbio/binder-prescreen/releases/tag/v0.5)
- **Insight:** 51 known TCR binders for KRAS G12D from IEDB, VDJdb, PDB, patents. Can benchmark our designs against known binders. May also reveal structural motifs for G12D recognition.
- **Action:** Download and run prescreen on our top designs. Extract structural features from known binders.
- **Status:** OPEN

## IDEA 6: Ahn et al. TCR-like antibody workflow (KRAS G12D)
- **Source:** Mol Ther 2026 — [Discovery of TCR-like antibodies to KRAS G12D](https://pubmed.ncbi.nlm.nih.gov/42237541/); bioRxiv preprint [2025.08.31.673313](https://www.biorxiv.org/content/10.1101/2025.08.31.673313.full.pdf)
- **Insight:** In silico docking from existing TCR-pMHC structures (PDB 6ULK), then limited CDR sequence design. Negative selection with a different peptide/HLA tetramer confirmed specificity. On HLA-C*08:02 (not A*11:01), but same peptide concept.
- **Action:** Could use PDB 6ULK or similar TCR-pMHC structures as scaffolds for our HLA-A*11:01 target. Structural alignment approach.
- **Status:** LOW PRIORITY (different HLA)

## IDEA 7: NYU monobody structures (KRAS G12D selective)
- **Source:** NYU dataset — PDB [8EZG](https://www.rcsb.org/structure/8EZG), [8F0M](https://www.rcsb.org/structure/8F0M)
- **Insight:** Monobodies selective to KRAS G12D (not pMHC, but the protein itself). Different binding mode but validates that G12D-specific contacts exist.
- **Action:** LOW PRIORITY — these bind KRAS directly, not the pMHC complex.
- **Status:** CLOSED (wrong target type)

## IDEA 8: Binder-peptide iPTM as selectivity metric
- **Source:** Own analysis of OpenDDE-AbAg pair_chains_iptm data (2026-10-30)
- **Insight:** Binder→peptide iPTM delta (G12D minus WT) correlates strongly with overall selectivity. test_d3 (#2 design, delta_conf=+0.109) has binder-peptide delta=+0.279. All selectives have positive binder-peptide delta. Designs with high binder→HLA iPTM but low binder→peptide iPTM are NON-selective (they bind MHC, not peptide). This validates Baker lab's upward-facing hotspot approach.
- **Action:** Use binder→peptide iPTM delta as a SECONDARY selectivity filter alongside overall confidence delta. Designs with binder-peptide iPTM delta > 0 and binder-HLA iPTM delta ~ 0 are the strongest selectivity candidates.
- **Status:** VALIDATED — incorporated into analysis

## IDEA 9: ipSAE_min metric (from C4 session literature review)
- **Source:** Dunbrack lab, bioRxiv 2025.02.10.637595; Overath et al. meta-analysis bioRxiv 2025.08.14.670059
- **Insight:** ipSAE_min (iPTM filtered by PAE cutoff) is the best single predictor of wet-lab binding in a 3,766-binder meta-analysis. But PAE field is EMPTY in our OpenDDE-AbAg results — not computable from current data.
- **Action:** NOT ACTIONABLE with current scoring model. Would need AF3 or Boltz-2 with PAE output.
- **Status:** BLOCKED (no PAE in OpenDDE-AbAg output)

## IDEA 10: Design diversity > better scoring (Adaptyv EGFR retrospective)
- **Source:** Adaptyv EGFR retrospective, bioRxiv 2025.04.17.648362v2 (from C4 session)
- **Insight:** Diversity of design approach quadrupled Round 2 hit rate, mattered more than better scoring.
- **Action:** Already executing — we have 5 pipeline sources. Research-informed hotspot strategies add 4 more design substrategies. Current approach is diversity-first.
- **Status:** VALIDATED — no new action needed

## Priority ranking
1. **IDEA 1** (upward-facing hotspots) + **IDEA 2** (D6 salt bridge) — EXECUTED: 7 new RFD4 jobs submitted
2. **IDEA 8** (binder-peptide iPTM metric) — VALIDATED: use for secondary filtering
3. **IDEA 3** (PDB 9UV8) — CHECKED: 9-mer structure, informational only
4. **IDEA 4** (9-mer vs 10-mer) — challenge target is 10-mer per 7OW6
5. **IDEA 5** (prescreen benchmark) — TODO: validate against known binders
6. **IDEA 10** (design diversity) — VALIDATED: already executing multi-source strategy

## IDEA 11: D6 aspartate is BURIED, not exposed — salt bridge is with MHC, not available for binder
- **Source:** Claude Science structural analysis of 7OW6 (session 7cb19b28, 2026-10-30)
- **MEASURED:** Asp6 side chain SASA = 7.9 Å² in 7OW6 (vs Val6 SASA = 83.0 Å² in G12V/7STF). Asp6 carboxylate is captured by HLA Arg114 (2.95-3.04 Å polar contacts) and HLA Gln70 (3.08 Å). Nearest TCR atom to Asp6 side chain = 5.73 Å (too far for contact).
- **Implication:** The D6-only hotspot strategy CANNOT work by direct salt bridge to Asp6 — the carboxylate is occupied by the MHC. BUT our D6-only designs ARE selective (delta +0.140), so they achieve selectivity through an INDIRECT mechanism: the Asp6→Arg114 salt bridge changes the peptide conformation and groove shape (positions 4-5-7 SASA differ between G12D and G12V), which the binder reads.
- **Revised design logic:** Don't target D6 directly. Target the CONFORMATIONAL CONSEQUENCE of D6: the altered backbone and side chain positions of neighboring residues (positions 5, 7, 8) that differ between G12D and WT due to the D6-Arg114 interaction. The "upward-facing" hotspot set (C5, C7, C8) may be better than D6-focused.
- **Status:** CRITICAL — revises understanding of the selectivity mechanism
