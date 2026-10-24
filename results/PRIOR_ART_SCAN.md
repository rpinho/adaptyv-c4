# Prior Art Scan — Adaptyv × Anthropic Protein Design Competition (2026-10-23)

## Scope
Web search + Codex GitHub scrape. 35 public repos found, plus preprints, guides, and X posts.

## C3 (GLP1R Agonist) — Competitive Landscape

### Public participant repos
| Repo | Tools | Designs | Strategy | Key insight |
|------|-------|---------|----------|-------------|
| [Abdulrahman-Neamah/glp1r-agonist-binder-design](https://github.com/Abdulrahman-Neamah/glp1r-agonist-binder-design) | BoltzGen, Boltz, Chai, Claude | 20 submitted | **Full receptor (chains R,S,T)** groove binders | Targets full receptor, not just ECD — may achieve agonism through TMD-adjacent binding |
| [jianurr/glp1r-binder-design](https://github.com/jianurr/glp1r-binder-design) | BindCraft, AF2, ColabFold, Chai, Rosetta | 2 submitted | ECD groove (BindCraft, Relaxed filter) | Minimal entry, single trajectory |
| [CCNY-McCann-Koder](https://github.com/jimjmccann/CCNY-McCann-Koder) | Chai | 3 (C3 symmetric backbones) | Symmetric backbone + 12 candidates | Academic lab, also covers C1/C2 |

### Relevant preprints
- **De novo miniprotein agonists/antagonists for GPCRs** (Nature, May 2026): RFdiffusion-designed. GLP1R antagonist IC50 7.9 nM. No strong agonism achieved for GLP1R specifically. [bioRxiv 2025.03.23.644666v3]
- **EvoBind cyclic peptide agonists** (bioRxiv, June 2025): 18-mer cyclic peptide, EC50 32 nM at GLP1R (liraglutide control: 2.4 nM). Sequence-only design against GCGR gave cross-activity. Key: targeted activation surface, not just binding.
- **AI GLP-1 RAs** (Advanced Science, 2025): ~10,000 designs, 2 leads with half-lives ~3× semaglutide.

### Our position
- **Strength**: 20 designs, 3 methods (BindCraft + variants + RFD4), WY-enriched for groove stacking, all ESMFold-validated
- **Gap**: ECD-only targeting — Abdulrahman-Neamah targets the full receptor which could engage TMD for agonism. We can't pivot in 2 days.
- **Gap**: No Boltz-2 complex prediction (everyone else uses it)

## C4 (TfR1 Conditional Binder) — Competitive Landscape

### Public participant repos
| Repo | Tools | Designs | Strategy | Key insight |
|------|-------|---------|----------|-------------|
| [anushajain2/adaptyv-challenge4](https://github.com/anushajain2/adaptyv-challenge4) | BoltzGen, PROPKA, OpenMM, Modal | 8 draft (not submitted) | His pairs at interface periphery, PROPKA pKa shift scoring | Shift = pKa_bound − pKa_free. Filter: both shifts ≤ −0.5, free pKa 6-7, ratio ≥ 2. No validated switch. |

### pH-switch methodology from C1/C2 entries (transferable to C4)
| Repo | Tools | pH method | Key metric |
|------|-------|-----------|------------|
| [actimo-labs/tnf-a-comp2](https://github.com/actimo-labs/tnf-a-comp2) | RFdiffusion + ProteinMPNN + Boltz-2 + PROPKA3 | **Thermodynamic linkage**: ΔΔG = ΔG_bind(6.0) − ΔG_bind(7.4) from pKa shifts | ΔΔG 2.84-4.56 kcal/mol, His >8Å apart |
| [Vijayavallabh/anthropic-adaptyv-2026-tnf](https://github.com/Vijayavallabh/anthropic-adaptyv-2026-tnf) | RFdiffusion + ProteinMPNN + LigandMPNN + AF2 + Chai + Boltz + PROPKA + Foldseek | Physics layer on top of folding models | "Folding models are blind to pH" |
| [byrneemma78-hub/tnf-ph-switch](https://github.com/byrneemma78-hub/tnf-ph-switch) | BindCraft2 + Chai + Boltz + PROPKA | "His placement, not His count" | Heading is the thesis |
| [Laurencejaboro19/Adaptyv-Bio-x-Anthropic](https://github.com/Laurencejaboro19/Adaptyv-Bio-x-Anthropic-competition-) | BindCraft2 + ProteinMPNN + Chai + Boltz + PROPKA + Modal | PROPKA3 / Wyman-linkage | Computational-only |

### Relevant preprint
- **Computational design of pH-sensitive binders** (bioRxiv 2025.09.29.678932): Two mechanisms — interface destabilization (electrostatic repulsion) and monomer destabilization (buried charged networks). Notes: most work historically uses random His mutagenesis + display selection. **Could not read full text (rate-limited).**

### Our position
- **Strength**: 20 designs, His-enriched (8-16%), charge-optimized variants, pI-validated, mutation-clustering fix
- **Strength**: Only public C4 entry with this many designs and iteration depth
- **Gap**: No PROPKA pKa calculations or thermodynamic linkage estimates — everyone else uses PROPKA
- **Gap**: No Boltz-2 complex prediction for interface validation
- **Gap**: No Wyman-linkage ΔΔG estimate (we use charge-swing Δq as a proxy)

## Cross-cutting insights

### Tools used by competitors (frequency across 35 repos)
- **PROPKA**: ~10 repos — the de facto pH gate standard
- **Boltz-2**: ~12 repos — dominant complex predictor
- **BindCraft/BindCraft2**: ~12 repos — most popular generator
- **Chai-1**: ~10 repos — orthogonal structure check
- **BoltzGen**: ~6 repos — de novo generator
- **RFdiffusion/RFd3**: ~5 repos — diffusion backbone generator
- **ProteinMPNN**: ~8 repos — sequence design
- **ESMFold**: ~5 repos — single-sequence folder (our main tool)
- **Foldseek**: ~4 repos — novelty screening
- **PROPKA + Wyman linkage formula**: The consensus pH-switch metric

### Key lessons
1. **PROPKA is essential** — "folding models are blind to pH" (Vijayavallabh). Every serious pH-switch entry uses PROPKA for pKa shift estimation.
2. **Boltz-2 over ESMFold for interface validation** — ESMFold is single-sequence, can't predict complexes. Boltz-2/Chai-1 give ipTM and interface quality.
3. **Shuffled-sequence controls** calibrate scoring thresholds (eliminate Boltz-2 false positives).
4. **His placement > His count** — "why histidine placement, not histidine count" (byrneemma78).
5. **Error correction is a feature** — harishrajaram tracked 13 errors and published the corrections. The community values transparency.
6. **Most entries are unvalidated** — every repo disclaims experimental results. All evidence is computational. This is a level playing field.

## Priority improvements for our entries

### C4 (9 days remaining) — HIGH PRIORITY
1. **Add PROPKA analysis**: Run PROPKA on ESMFold PDB outputs to get pKa shifts for all His residues. Calculate thermodynamic linkage ΔΔG.
2. **Add Boltz-2 complex prediction**: Fold our designs against TfR1 to get ipTM/ipSAE interface scores.
3. **Add Wyman-linkage scoring**: Replace our Δq proxy with actual ΔΔG estimates.

### C3 (2 days remaining) — LIMITED TIME
1. **Add Boltz-2 complex prediction** for top 5 designs (if time permits)
2. Cannot pivot to full-receptor targeting in 2 days

## Sources
- [Adaptyv × Anthropic competition page](https://proteinbase.com/competitions/anthropic-adaptyv-2026)
- [Intuition Labs guide](https://intuitionlabs.ai/articles/anthropic-adaptyv-protein-design-competition-guide)
- [Anthropic announcement](https://x.com/AnthropicAI/status/2100701582744797347)
- [Adaptyv Bio on BindCraft2](https://x.com/adaptyvbio/status/2101983777593364823)
- [EGFR competition paper](https://www.biorxiv.org/content/10.1101/2025.04.17.648362.full.pdf)
- [EvoBind cyclic peptide preprint](https://www.biorxiv.org/content/10.1101/2025.06.06.658268v1.full)
- [LigandAI GLP1R news](https://runtimewire.com/article/ligandai-free-peptide-design-glp1r-challenge)
- [GEM × Adaptyv RBX1 competition](https://www.gembio.ai/gemadaptyv2026)
