# C4 mouse compatibility: fixed-pose audit

Written 2026-10-26T12:59:40Z. Evidence only. No ranking or submission was produced.

Twenty of twenty v4 sequences match unique actual binder chains by SHA256; every
target chain equals human P02786 residues 189-383. Mouse Q62351 is independently
sequence-verified. Retain all six full-length BLOSUM62 optimal alignments (gap open
-10, extension -0.5). Domain fits use 192 unambiguous CA pairs; local fits use the
unambiguous CA pairs within 12 A of any contacted human target CA, selected before
examining fit results. Only aligned apical mouse residues and human A189-383 enter
the transfer test. These are human-pose transfers onto a mouse monomer prediction,
not new mouse complex predictions or full native-context validation.

Human contacts use any heavy atom <=4.5 A. Conservation spans 60.87%-95.24% across
designs. All twenty have noninvariant contacts; fourteen contact K261 -> mouse E263
(a formal charge reversal). This does not establish a lost salt bridge, energetic
penalty, or pH response. F297 -> mouse L299 is frequent. site2_152 contacts the
ambiguous L212 -> gap/P214 region; every alignment is retained, never resolved to
improve the result. The old small hotspot summary does not describe every actual
v4 contact surface. See the full contact-level JSON and substitution lists in CSV.

At the prespecified <2.0 A threshold, 13/20 have mouse target-backbone overlaps in
both fits. Twelve of those also conflict with experimental human-reference apical
transfer in at least one fit: they cannot isolate a species effect. site1_185 has
the remaining differential backbone warning: mouse minima 1.915/1.784 A versus
human-reference 2.081/2.065 A, with only 2/1 mouse overlapping pairs. The human
reference still has side-chain overlaps, so this is not a human-pass/mouse-fail
binding test. Fits have different source conformations and are not independent.

Three designs are fit-sensitive. Four avoid target-backbone overlap under both fits:
site1_145, site2_172, site1_40 and site1_56. All four still have <2.0 A heavy-atom
overlaps in both fits, including local minima 0.826, 1.316, 0.661 and 0.452 A.
They remain backbone-compatible hypotheses with unresolved side-chain conflicts,
not clash-free designs. No design passes the full atom screen under BOTH fits.
site1_123 alone passes the local full-atom screen by 0.0015 A (minimum 2.0015 A),
while its domain fit has three backbone overlaps and minimum 1.565 A. This is a
fragile fit-dependent clearance, not robust compatibility. All other designs have
full-atom overlaps under both fits.

There is no evidence here that any sequence is biologically incapable of mouse
cross-reactivity. The current poses with persistent backbone overlaps cannot be
retained unchanged against this mouse model under this test. Resolving them would
require a pose or conformation change; neither rescue nor impossibility was tested.
Contact retention is proximity only: an overlapped atom can count as a contact.

Controls: twenty exact human self-transfers preserve contacts and minima; twenty
deliberately injected mouse backbone overlaps fire; twenty far-away controls have
no contacts; K205 -> S208 charge-loss mapping control fires. The matched experimental
human transfer is a separate scientific control, and exposes the transfer artifact.
These controls validate arithmetic and scope, not predictive accuracy.

Sensitivity: persistent mouse backbone conflicts are 9/20, 11/20 and 13/20 at
thresholds <1.5, <1.8 and <2.0 A. Of those, 8, 9 and 12 respectively also fail a human
transfer control. The mouse-only warning set changes with threshold; site1_185 is
not a cutoff-independent rejection. The companion table reports all 20 designs
at contact cutoffs 4.0/4.5/5.0 A and backbone thresholds 1.5/1.8/2.0 A (180 rows).
Do not choose a new cutoff based on which candidate it admits.

The table below preserves the input roster order; it is NOT a new ranking. M/H
are mouse/human-reference backbone overlap pair counts for domain/local fits.
Full molecule hashes and all changed contact residues are in per_design_assessment.csv.

| Design shorthand | Sequence SHA256 prefix | Invariant contacts | M domain/local | H domain/local | Interpretation |
|---|---|---:|---:|---:|---|
| site1_123 | `54e49124effc` | 17/21 | 3/0 | 3/2 | mouse fit sensitive conflict |
| site2_176 | `2527c2d491d0` | 25/27 | 1/0 | 0/0 | mouse fit sensitive conflict |
| site1_194 | `962e1a74ba7c` | 21/24 | 5/2 | 2/1 | mouse conflict with human transfer confounded |
| site2_120 | `6de4e987f029` | 28/35 | 5/5 | 5/5 | mouse conflict with human transfer confounded |
| site1_145 | `51b8f16f504e` | 14/21 | 0/0 | 0/0 | mouse backbone clear both fits |
| site2_108 | `c57083c527c4` | 36/40 | 2/4 | 2/1 | mouse conflict with human transfer confounded |
| site2_172 | `3c94a4ed1cb9` | 15/21 | 0/0 | 0/0 | mouse backbone clear both fits |
| site2_152 | `e9ab1a369927` | 37/47 | 5/9 | 9/10 | mouse conflict with human transfer confounded |
| site1_79 | `48dcf3612ea0` | 22/27 | 1/1 | 0/1 | mouse conflict with human transfer confounded |
| site2_12 | `7cfabe1fce28` | 32/40 | 8/8 | 9/9 | mouse conflict with human transfer confounded |
| site1_14 | `53ca8a6a3c7e` | 31/34 | 1/4 | 2/3 | mouse conflict with human transfer confounded |
| site1_73 | `a53ae434180e` | 14/23 | 4/0 | 2/1 | mouse fit sensitive conflict |
| site1_185 | `fb46aad4a7f2` | 24/30 | 2/1 | 0/0 | mouse conflict human reference clear both fits |
| site2_174 | `9a050f884378` | 22/27 | 3/2 | 3/2 | mouse conflict with human transfer confounded |
| site1_144 | `69b90475260b` | 20/21 | 27/46 | 24/33 | mouse conflict with human transfer confounded |
| site1_174 | `422533dc2cb8` | 21/24 | 9/16 | 8/10 | mouse conflict with human transfer confounded |
| site1_40 | `86959e1e6b45` | 24/27 | 0/0 | 0/0 | mouse backbone clear both fits |
| site2_1 | `477d52645685` | 25/31 | 2/2 | 1/1 | mouse conflict with human transfer confounded |
| site1_24 | `9749d13a5e70` | 26/30 | 11/18 | 15/12 | mouse conflict with human transfer confounded |
| site1_56 | `c1ce6ae493ac` | 15/22 | 0/0 | 0/0 | mouse backbone clear both fits |

Binding, KD, mouse cross-reactivity, protonation-dependent interface rearrangement,
pH7.4 binding/pH6.0 release, mouse pose relaxation, homodimer context, complete assay
glycans, and partner-chain accessibility remain unmeasured. Mouse/human reference
coordinates have different provenance; local fitting is a sensitivity check, not
independent validation. No new GPU job was launched.

Consultations: exact question and evidence digest are recorded under
results/decisions/20261026-c4-mouse-consult-*. Raw Astra and Kimi reviews are retained
separately. The evidence packet contains the full 20-row metrics and frozen v1 rule.

Recommendation: hold current C4 submissions. Keep v4 as the human-surrogate research
incumbent; do not promote the four backbone-clear rows into a new file. Geometry
and conservation do not resolve the required mouse/pH objectives. This is an
evidence limitation, not a claim of twenty nonbinders. Ricardo submits personally.

Reproduce in a fresh results directory, CPU only:

    OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 python3 scripts/c4_mouse_compatibility_20261026.py --out results/c4_mouse_reproduce_20261026
    OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 python3 scripts/c4_mouse_transfer_control_20261026.py --out results/c4_mouse_reproduce_20261026
    python3 scripts/c4_mouse_summary_20261026.py --out results/c4_mouse_reproduce_20261026

The scripts refuse to overwrite measurements or summary files. JSON input digests
bind all 20 original complex files, sequences, reference structures and the rule.
Biopython 1.85 and NumPy 2.0.2 were used with the installed Python 3.9; .venv does
not exist on this checkout. No historical results directory was rerun or overwritten.
