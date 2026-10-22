# Public patent-protein evidence added; no retired-pool local calls changed

NCBI pataa is now obtained and searched explicitly as a public patent source.
It is **not verified as the organizer's patent corpus**, and patent scope remains
incomplete in the manifest. The unchanged predictor returns **16 provisional
local level 2, 47 level 3, 4 level 4, and 6 coordinate-unavailable** retired
designs. Zero local calls changed. Official novelty levels established: **0**.

## What ran and measured numbers

The official NCBI pataa release dated 2026-10-21 contains **3,982,722 records and
787,058,342 residues**. Its 843,504,877-byte archive passed publisher MD5
verification. Official blastdbcmd 2.17.0+ for aarch64 macOS was obtained from the
versioned publisher archive, verified against its publisher MD5, and fingerprinted.
Its FASTA export contains exactly the metadata's record/residue counts and is
1,213,622,499 bytes. No records were deduplicated or filtered, and no sequences
were repaired. There are 3,694,389 fully canonical records; ambiguous letters
remain in the original export. Counts and alphabet are in ACQUISITION.json.

All five public sequence-source searches completed with production MMseqs2
18-8cc5c, four threads, 4G search memory, GPU disabled, sensitivity 7.5,
E-value <=0.001, 10,000 candidate cap, integer alignment identity, and no
identity/coverage prefilter. The patent search log confirms **3,982,722 target
records**, not a sampled database. Finite search and the logged default maximum
sequence length of 65,535 remain search-policy limitations; this experiment does
not establish exhaustive patent homology recall.

| Source | Records/sequences supplied to this search | Search seconds | Qualifying alignment rows |
|---|---:|---:|---:|
| SwissProt public snapshot | 487,476 deduplicated sequences | 10.384 | 32 |
| PDB public SEQRES protein snapshot | 188,017 deduplicated sequences | 3.629 | 34 |
| Thera-SAbDab public bulk snapshot | 2,211 deduplicated sequences | 1.037 | 0 |
| PLAbDab paired + unpaired public bulk | 372,631 deduplicated sequences | 3.769 | 0 |
| NCBI pataa public Patent-division proteins | 3,982,722 source records, duplicates retained | 42.822 | 231 |

Inputs: the unchanged 73 retired sequences and native TfR1 positive control,
plus one explicitly labeled patent-source control. Before searching, that control
was selected as the first canonical 60–250-aa record in source order. Its original
source identifier matches over **110/110 residues at 100% identity**, query
coverage **1.0**, E-value **2.794e-67**. It is a database control, never a design.
The 231 patent rows belong to controls: **0 of 73 retired designs** has a
qualifying patent hit under this policy. No-hit is not measured zero identity.

The full structural run reused the valid five-source sequence evidence and
repeated the unchanged structural pipeline on **all 68 frozen coordinate domains**
(67 retired + native control), including the same whole-chain assumptions,
historical PDB database, union discovery, full TM-align and exact score export.
It completed in **216.380 seconds**, with **42,718 raw rows / 24,653 distinct
query/target pairs**. Structural evidence for all 74 original inputs is identical
to the frozen exact-score report. Native control remains local level 1, qTM
**0.9928**. Source, input, database and raw-output fingerprints verify.

## Input-preparation error caught and corrected

The initial comparison manifest accidentally copied the small eight-domain
sensitivity subset rather than the complete 68-domain frozen input set. Its
full-pool structural comparison is invalid and is not a change in design quality.
The original inputs and misleading preliminary summary are preserved under
partial_structure_inputs/; the original run is sequence_and_partial_structure_run/.
PARTIAL_STRUCTURE_SCOPE_ERROR.json records the error and fingerprints.

The five sequence searches from that run remain valid. They were reused without
repeating the searches. Preparation now inherits every coordinate input directly
from the frozen exact-score report and requires exactly 68 domains before launch.
The final validator independently checks the complete coordinate-name set, full
sequence identities, input fingerprints, source control and all raw sequence
tables. Only run/predictions.json and VALIDATION.json are the final comparison.

The archive extraction guard also stopped on bundled taxonomy files before
their names were explicitly accounted for. The verified archive contains 15
regular files, including taxdb.btd, taxdb.bti and taxonomy4blast.sqlite3. Original
failure logs are retained; extraction resumed only after checking archive members
and every already-extracted file. No scientific threshold was changed.

## Remaining gaps and use

DATABASES.json identifies all five actual obtained sources and the partial patent
scope. It is the updated manifest for future sequence screens. Earlier manifests
and raw reports remain frozen as history. The organizer's patent corpus/provider,
all organizer snapshots/search knobs, exact structure predictor and domain
consensus versions, antibody typing/CDR scale and recall near boundaries remain
unestablished. Bulk antibody snapshots are older than this run. The historical
structural database and standalone TM-align residual disagreement remain caveats.
All 73 designs remain retired; this work produces no generator-accepted design.

New GPU/cloud CPU spend **$0.00**. Recorded Runpod **$56.589988** posted + our Modal
**$10.718947** estimated = **$67.308935 (~$67.31)**, excluding the separate RFD4
lane, storage/builds and consultations. No generation, cloud allocation, exporter,
ranker, account, portal/submission or main-branch action occurred.

Primary provenance: [NCBI metadata](https://ftp.ncbi.nlm.nih.gov/blast/db/pataa-prot-metadata.json),
[NCBI patent-data scope](https://www.ncbi.nlm.nih.gov/education/patent_and_ip_faqs/),
[official versioned exporter](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.17.0/),
[published novelty rules](https://www.adaptyvbio.com/blog/novelty).
