# Mouse conservation of the current apical epitope

**Verdict: partially conserved only; as-is cross-reactivity is unestablished.**
Two of five identities survive. The other three are nonconserved, including two
losses of positive charge. Even the two identical positions have changed neighbours.
This is a sequence/structural conservation measurement, not a binding prediction.

## Alignment and numbering

Full-length human P02786 (760 aa) and mouse Q62351 (763 aa), read directly from
the local UniProt flatfiles, were globally aligned with Biopython PairwiseAligner,
BLOSUM62, affine gap opening -10 and extension -0.5. All six optimal alignments
(score 3100.5) are retained in `all_optimal_alignments.txt`; none was selected to
force an offset. `apical_alignment.csv` gives every human position in the project's
apical span, author 189–383, with every equally optimal mouse correspondence.
Mouse insertions and both complete sequences are visible in the full alignments.

Human author numbers equal UniProt positions. Every protein atom's author/label
pair in 11QC A was checked: author minus label = 121. Mouse AF PDB residue identities
were checked against Q62351, independently of the human alignment. The mouse offset
varies: for example human 189 maps to mouse 191, whereas human 205 maps to mouse 208.

| Human auth / 11QC label | Aligned mouse UniProt | Identity and chemistry |
|---|---|---|
| K205 / 84 | S208 | Substituted, nonconservative: long basic side chain → short polar hydroxyl; +1 → 0 |
| N206 / 85 | N209 | Identical neutral polar amide |
| R208 / 87 | N211 | Substituted, nonconservative: guanidinium → amide; +1 → 0 |
| L209 / 88 | L212 | Identical neutral aliphatic hydrophobe |
| V210 / 89 | D213 **or gap** | Nonconserved in all optima. If D: nonconservative hydrophobe → acidic carboxylate, 0 → -1. A deletion has no substitution chemistry. |

Human 210 is deleted in two of six optima and maps to D213 in four; the alternative
deletion falls at human 211 or 212. All other queried positions have unique mappings.
The sequences do not uniquely resolve that local indel. An AF domain fit on 192
unambiguous sequence-derived CA pairs gives RMSD 1.719 Å; it is a geometric check,
not a reason to relabel an ambiguous alignment as uniquely determined.

Across author 189–383, **146/195 (74.87%)** human positions are identical in every
optimum. Across human ectodomain 121–760, **499/640 (77.97%)** are identical using
human span length as denominator, counting gaps as nonidentical. This does not
reproduce the supplied approximate 86% over a modelled region. That region/identity
denominator is unspecified here; the measured canonical-sequence result is retained.

## Structural neighbourhoods

Reference: experimental 11QC chain A. A neighbour is another protein residue with
**any heavy atom within ≤8 Å of any heavy atom of the queried residue**. The central
residue is excluded; all target domains are included, so these are not merely
sequence-adjacent positions. Ranges reflect all six optimal sequence alignments.
Deletions are reported separately from substitutions.

| Human site | Other residues within 8 Å | Substituted neighbours | Deleted neighbours | Total nonidentical neighbours |
|---|---:|---:|---:|---:|
| K205 | 15 | 6–7 | 0–1 | 7 |
| N206 | 8 | 4 | 0 | 4 |
| R208 | 12 | 4–5 | 0–1 | 5 |
| L209 | 15 | 7 | 1 | 8 |
| V210 | 15 | 8–9 | 0–1 | 9 |

Substitution and deletion counts covary; do not add their independent maxima.
N206's four changed neighbours are D204→Q207, K205→S208, R208→N211 and
S370→N372. L209 additionally sits beside the ambiguous 210–212 deletion region.
Every changed neighbour and mouse correspondence is listed in `measurements.json`.
An independent KD-tree neighbour search reproduces all five neighbour totals and
finds **zero additional partner-chain neighbours** within this radius.

`targets/E1_mouse_his710.pdb` contains 110 residues in a discontinuous His710 crop;
it covers none of these five corresponding sites. The complete 763-residue
`data/AF-Q62351.pdb` was therefore used for mouse sequence/geometry checks.
Neighbour counts above describe the measured human reference environment; mouse
rotamers, loop conformations and complete assay glycans are not established.

## Best alternative by measured local conservation

The highest-ranked conserved **surface core** in this scan is centred on human E266:

| Human auth | Human 11QC label | Mouse UniProt | Identity |
|---|---:|---|---|
| K240 | 119 | K242 | Identical |
| I262 | 141 | I264 | Identical |
| E266 | 145 | E268 | Identical |
| A269 | 148 | A271 | Identical |
| L303 | 182 | L305 | Identical |

All five have relative solvent accessibility ≥0.20 in the human native dimer with
observed glycans, and CA distance ≤10 Å from E266. Their combined exposed area is
**364.7 Å²**. The union of their full 8 Å environments contains **52/59 identical
residues (88.14%)**, with seven changed residues. The same surface window includes
K261→E263 and E272→Q274: this is a conserved core in a partially substituted
surrounding surface, **not an entirely invariant recognition surface**.

The scan required ≥5 invariant exposed residues per window and ranked windows by
the invariant fraction of the union of their 8 Å environments, then conserved
exposed area. No qualifying window had a wholly invariant environment. This is the
best alternative under that explicit conservation criterion, not a globally best
binder epitope. Ranked alternatives and per-site neighbourhoods are saved in JSON.
No binding, mouse cross-reactivity, accessibility to a particular binder, or
functional compatibility was measured.

Reproduce locally, CPU only:

    .venv/bin/python scripts/mouse_epitope_conservation.py

All five input fingerprints are embedded in `measurements.json`. Human and mouse
coordinates are not changed by this script.
