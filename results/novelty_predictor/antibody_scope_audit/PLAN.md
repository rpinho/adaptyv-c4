# CPU antigen-receptor annotation scope check

Registered before annotation outcomes. An actual full-query PLAbDab heavy-chain
match in the incoming 400-design batch makes the general-vs-antibody scale an
important unresolved issue. Keep the general-scale batch unchanged. This audit
annotates existing sequences; it does not generate, mutate or rank designs.

Run all 400 sequences, without score-based selection, through pinned ANARCII
2.0.8 in antibody/accuracy mode, CPU only, one numerical thread, fixed seed 0.
Retain raw annotations, per-query errors, model/source fingerprints and timing.
Check the publisher's classic ANARCI 12e8 FASTA as a positive control (H records
must identify heavy; L records must identify either light class), and the saved
195-residue native receptor control as a negative. Do not change published model
thresholds to pass a control. Require residue annotations to map to the actual
input sequence; never repair or synthesize sequence.

ANARCII is an explicitly different annotation model from the organizer's unnamed
structural typer. Its calls do not establish official molecule class, VHH-vs-VH
status or the applicable portal scale. IMGT CDR3 annotation alone does not measure
CDRH3 novelty. A full antibody-specific level additionally needs the actual
reference CDRH3 comparison; edit-distance normalization is unspecified in the
organizer's description. Leave that gap explicit rather than inventing a score.

The classic ANARCI source was pinned/read first, but it does not ship the HMM and
germline data: setup rebuilds them from live IMGT downloads. That dynamic build
was not executed. HMMER3.4 was installed from its Homebrew bottle while inspecting
that route; no classic ANARCI result is claimed. ANARCII is not silently labeled
as the same implementation.

Primary sources:
- https://www.adaptyvbio.com/blog/novelty
- https://github.com/oxpig/ANARCII
- https://pypi.org/project/anarcii/2.0.8/
- https://github.com/oxpig/ANARCI

New GPU/cloud compute: $0.00. Recorded Runpod $56.589988 + own Modal estimate
$10.718947 =~$67.31; separate generation lane and account-wide billing excluded.
