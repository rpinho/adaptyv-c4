# Independent review of the Challenge 4 evidence

23 October 2026

**Assessment: this is a computational hypothesis set with documented geometric screening, not an established set of pH-selective binders. The evidence does not establish performance on any of the three stated ranking criteria.** That is an evidence limitation, not proof that every design will fail.

Scope: I reviewed the four requested files and checked CSV format, aggregate annotations, and its recorded SHA-256. I did not rerun structural calculations or independently validate the supplied clusters. This review addresses scientific validity and reporting. I cannot provide sequence-specific selection, reranking, mutations, or generation instructions to improve targeted blood-brain barrier transport. Accordingly, it does not supply the requested submission subset or redesigned campaign.

## 1. Overall credibility and submission readiness

There is meaningful computational quality control here. RESULTS.json reports that all 20 pass the unchanged omitted-context and native-staging checks, with no omitted-context backbone pairs below 2 Å. Restoring receptor context addresses a real class of modeling artifact. However, the reported 0.229–0.347 Å target-fit RMSD measures receptor alignment; it is not evidence that a binder adopts its proposed fold or binds that receptor.

All 20 sequences are unique, contain canonical amino-acid letters, and are 73–121 residues long. They therefore satisfy the stated sequence-length requirement and are represented as single-chain entries. This does not establish their physical oligomeric state or organizer eligibility.

| Claim | Evidence in these files |
|---|---|
| Compatibility with the specified backbone clash gates | Reported for 20/20 |
| Binding at pH 7.4 | Unmeasured |
| No detectable binding at pH 6.0 | Unmeasured |
| Mouse TfR1 binding | Unmeasured |
| Affinity | Unmeasured |
| Independent folding, expression, and developability | Unmeasured |
| Organizer novelty and acceptance | Unmeasured |
| BBB transcytosis | Not established |

I would not endorse a description of these as demonstrated conditional binders. They can accurately be described as unvalidated computational designs. A design competition may intentionally accept untested hypotheses; missing experimental validation alone is not an automatic disqualification. Whether to spend submission slots on them is a different question from whether their claimed mechanism has evidentiary support.

## 2. Histidine enrichment and the proposed pH mechanism

First, the prompt's composition grouping needs correction. Ranks 6 and 11 contain approximately **1.90% and 0.85% His**, respectively, despite coming from the 0.30-conditioning arm. Rank 16 contains **27.27% His**. Across the entire shortlist, His content ranges from **0.85% to 54.62%**, and seven sequences contain at least 30% His. Conditioning strength is not the realized composition of an individual output.

An approximately half-histidine sequence is a substantial composition bias. It warrants uncertainty about folding and solution behavior; it does not justify a categorical assertion that the protein cannot fold. Composition alone cannot establish a unique stable structure, soluble expression, aggregation, or their absence. The available evidence cannot answer “Will these fold?”

Histidine's ionization makes a pH response chemically plausible in general. It does not establish the direction, magnitude, or molecular origin of that response. Experimental work shows that histidine ionization depends on its protein environment and can differ between conformational states. Thus a pH-dependent signal could reflect altered conformational stability rather than the intended change in interface affinity. These are general principles, not predictions for these sequences. See [Hansen and Kay, measurement of histidine ionization in protein states](https://pubmed.ncbi.nlm.nih.gov/24733918/) and [Edgcomb and Murphy, environmental dependence of histidine pKa](https://pubmed.ncbi.nlm.nih.gov/12211010/).

The current scan does not test the proposed switch:

- It changes residue identity, not solution pH or protonation equilibrium. Removing a histidine is not equivalent to protonating that histidine.
- Its output is a relaxed site-to-receptor interaction score in TMol beta2016 units. It is neither a calibrated binding free energy nor a measured affinity.
- It covers 10 sites in five preselected designs. The other 15 are unscanned; this is missing evidence, not evidence of failure.
- Three sites have positive deltas for both substitutions. That observation does not establish pH selectivity. Conversely, negative mutation deltas do not prove pH insensitivity.
- The largest absolute self-substitution delta is 0.3571 score units. A nominally unchanged residue can therefore produce a nonzero result under this procedure. This is a control observation, not a universal error bar or significance threshold.

Nor does His–basic-group proximity establish a switch: the annotation is a geometric hypothesis, as the files themselves acknowledge. **Neither a high His percentage nor a favorable mutation score demonstrates binding at one pH and undetectable binding at the other.** “Undetectable” also depends on the assay's detection limit and conditions; it cannot be equated to an arbitrary modeled affinity ratio.

## 3. Diversity

The CSV contains **15 cluster labels across 20 unique sequences**. Cluster sizes are two groups of three, one group of two, and twelve singletons. Thus the supplied annotations do not support calling the entire shortlist merely a few sequence families. They also do not establish 20 independent structural hypotheses.

These are supplied labels, not independently verified families. The four files do not establish the clustering threshold, fold diversity, or interface diversity sufficiently to interpret every singleton as a distinct solution. Four conditioning arms share a generator, target region, and broad mechanistic premise; changing the conditioning prior does not itself produce an independent mechanism. Twelve entries come from the 0.30 arm, five from 0.20, two from 0.10, and one from plain generation.

The defensible conclusion is **some annotated sequence diversity, with structural and mechanistic diversity unresolved**. Sequence uniqueness also does not establish organizer novelty.

## 4. Red flags and potential embarrassment

The primary risk is overclaiming, rather than an identified format disqualification.

**Proxy measurements could be mistaken for biological results.** “Clash-audited,” “His-scanned,” and “format-ready” each describe completed operations. None establishes binding, a pH switch, or transport. Even demonstrated receptor binding and pH dependence would not, by themselves, demonstrate BBB transcytosis.

**METHOD.md mixes current results with historical modeling narratives.** Older sections use words such as “measured” and “confirmed” for computational pKa/linkage outputs and contain superseded conclusions about epitopes. Those results concern different design sets and cannot validate this shortlist. The current opening limitations are appropriately explicit, but the document as a whole is not a clean, self-contained methods account for these 20 entries.

**The modeled context has limits.** RESULTS.json explicitly says that glycans and the assay tag are absent from scored protein poses. Clearance from observed NAG atoms in the CSV does not certify compatibility with complete, flexible glycosylation. Backbone clash checks do not establish all-atom energetic compatibility or native receptor behavior.

**Eligibility remains unresolved.** The files flag novelty/type acceptance, Twin-Strep exclusion, and upload behavior as outstanding. I have not verified current organizer rules and do not assert that any entry is disqualified. These cannot be represented as completed checks.

There is no evidence here of fabricated experimental data: the current result files repeatedly identify the biological unknowns. The concern would arise if a submission narrative dropped those qualifications.

## 5. Rank order

**The rank order is not validated as an order of Challenge 4 performance.** All entries lack measurements on the primary criterion, mouse cross-reactivity, and affinity. The CSV's rank_score field is blank throughout; its listed clash and bad-bond metrics are uniformly zero, and hotspot obedience takes only two values. These columns cannot be interpreted as a calibrated probability of success.

Only the original top five received mutation scans, and RESULTS.json explicitly says the original order was preserved. The scan is therefore neither an independent validation of that ordering nor a comparable assessment of all 20. Favoring scanned entries simply because they have more annotations would confuse evidence availability with demonstrated function.

I cannot supply a candidate-specific alternative order for this application. Scientifically, the existing order should be labeled a computational selection order, with biological performance unknown.

## 6. Recommendation

**Do not present this shortlist as having established the challenge's defining property.** Preserve the sequences' status as exploratory computational hypotheses and report the geometric results within their actual scope. The files support neither a claim that all 20 are strong conditional-binder candidates nor a claim that all 20 are biologically impossible.

I am not recommending a submission subset or a new generation strategy. Within this evidence review, the concrete recommendation is to correct the composition summary, separate historical computations from current-shortlist evidence, and avoid converting mutation scores, cluster labels, or clash passes into claims of functional qualification. More candidates evaluated by the same proxies would not, by itself, resolve the missing biological evidence.

Files reviewed: `results/rfd4_his_scan_20261023/submission_for_ricardo_review.csv`, `METHOD.md`, `results/rfd4_his_scan_20261023/RESULTS.json`, and `results/rfd4_his_scan_20261023/MUTATION_ENERGIES.csv`. The CSV SHA-256 matches RESULTS.json. No designs, ranks, or source results were changed.
