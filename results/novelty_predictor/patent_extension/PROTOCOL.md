# Explicit public patent-protein extension (CPU only)

Add NCBI's pataa collection as a specifically named public patent source. Do not
claim it is the organizer's unidentified patent corpus, or mark patent scope
complete. Do not silently replace the corpus with antibody databases. Existing
results remain frozen; a new manifest/report records this additional coverage.

Source metadata obtained before downloading: pataa, Protein, 3,982,722 sequences,
787,058,342 residues, release 2026-10-21. Archive expected 843,504,877 compressed
bytes and 2,308,834,559 unpacked bytes. Pin the observed metadata, HTTP headers,
publisher MD5 and local SHA256. Verify publisher MD5 before extraction. Export
every source record with official NCBI blastdbcmd 2.17.0+, aarch64 macOS binary
from the versioned publisher archive, itself verified against publisher MD5.
Record the exporter binary SHA/version, source-record count, exported count,
residue count and alphabet. Preserve duplicates and ambiguities; report any
exclusion or transformation explicitly. No proprietary corpus/account access.

Use unchanged production MMseqs2 settings and named boundaries. Search all 74
existing analysis sequences plus one patent-source positive control: the first
canonical 60–250-aa sequence in source order, selected before any search result.
The additional query is a database control only and never a design. A valid
control must have a complete 100%-identity source match; diagnose rather than
relax this requirement if it fails. Do not infer zero identity from no hits.

Recompute local novelty using the measured sequence maximum and frozen verified
exact structural evidence, without selecting a more favorable coordinate model.
For a complete production execution, rerun the unchanged structural half on the
same coordinates/settings and compare it to that frozen evidence. Any difference
must be reported; do not change structural search settings or thresholds.
Existing controls and all 73 retired designs stay labeled as such. Patent scope,
organizer snapshot/parameter, structure-prediction and domain-policy gaps remain
explicit; no official novelty level or qualified design is established.

Run locally with at most four CPU threads, 4G MMseqs search memory, finite stage
timeouts and a free-space check. No GPU, cloud CPU, generation, exporter/ranker
change, account, portal, submission or main-branch action. Recorded spend remains
Runpod $56.589988 posted + our Modal $10.718947 estimated = approximately $67.31;
new GPU/cloud CPU cost $0.00.

Primary sources:
- https://ftp.ncbi.nlm.nih.gov/blast/db/pataa-prot-metadata.json
- https://www.ncbi.nlm.nih.gov/books/NBK62345/
- https://www.ncbi.nlm.nih.gov/education/patent_and_ip_faqs/
- https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.17.0/
- https://www.adaptyvbio.com/blog/novelty
