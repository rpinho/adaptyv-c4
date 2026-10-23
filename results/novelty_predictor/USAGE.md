# CPU novelty workflow

Run a new CSV and optional saved coordinate package on this Mac:

```sh
.venv/bin/python scripts/run_cpu_novelty.py batch.csv \
  --structures coordinates.json --out results/novelty_batch_001
```

The output directory must be new. No account, upload, GPU, cloud service or new
protein structure is created. The workflow uses at most four numerical threads
and the existing pinned local reference databases/domain tools. Final JSON/CSV
paths are written to the output directory's STATUS.json. The CSV is an analysis
report, not a competition submission.

Input CSV headers: `name,sequence`. Names must be unique; sequences must use the
20 canonical uppercase amino acids. If molecule_class is present it must be
single_chain. The general novelty scale is supported; antibody typing and the
separate antibody/CDR scale remain unimplemented.

Coordinate manifest example:

```json
{
  "structures": [
    {"name": "design_001", "path": "models/design_001.cif", "chain": "B"},
    {"name": "design_002", "path": "models/design_002.pdb"},
    {"name": "design_003", "status": "unavailable", "reason": "No prediction yet"}
  ]
}
```

Paths are relative to that manifest. PDB and CIF are supported. Provide `sha256`
to enforce an expected fingerprint; otherwise the workflow acquires and records
the source hash. A supplied chain is sequence-verified. Without a chain, exactly
one chain must match the complete query sequence; ambiguity is rejected. An
explicit `author_residue_range: [first, last]` may select a known construct/control
region, and the complete resulting sequence must still match. No residue mapping
is inferred from a score, no sequence is repaired and no missing atoms are filled.

The pipeline searches all five public sequence families, runs a whole-chain
structural baseline, extracts verified complete backbones, runs CPU Chainsaw,
Merizo and source-built UniDoc, forms checked high/medium TED-style domains,
searches each domain, and adds every reference target found through the same
protein's whole-chain context to its actual full TM-align comparisons. Discontinuous
domains retain their actual residue sets and query-length normalization. All
named sequence/TM/coverage boundaries stay fixed. Model failures, missing backbone
atoms or coordinate/sequence mismatches stop the workflow with a failed status;
partial stage artifacts are retained for diagnosis.

Omit `--structures` for a sequence-only screen. Missing coordinates or accepted
domains leave the structural local level unestablished. If all three domain tools
assign no accepted domain, the workflow withholds numeric levels. It never copies
a stale whole-chain level into a missing-domain result. Structure inference from
sequence alone is not implemented on CPU in this workflow.

Default public-source manifest: patent_extension/DATABASES.json. The NCBI pataa
collection is obtained explicitly; it is not verified as the organizer's patent
corpus. Default structural reference: structural/DATABASE.json, a historical
1,562,678-chain PDB snapshot. Every source is fingerprinted; no missing database
is silently replaced. Alternative explicit manifests can be supplied with
`--databases` and `--structure-database`. Tool paths can be set with `--mmseqs`
and `--foldseek`.

The local domain runtime is under
`~/.cache/adaptyv-c4/novelty_domain_source_audit/`. Its native executables, weights,
source files and corrected consensus implementation are checked before inference.
Installation/calibration receipts are in domain_consensus/. The workflow does
not automatically install, update or rent a substitute if those resources differ.
Third-party wrappers receive safe, generated runtime paths and neutral query IDs.

All numeric levels are **provisional local predictions**. Organizer structure
prediction/domain versions, search settings, corpus snapshots, antibody policy,
finite retrieval recall and actual portal feedback remain gaps. Coverage-limited
local4 calls are not proof of a new fold. Current measured controls and limitations:
domain_consensus/RESULTS.md. A native-control end-to-end integration check of the
generic workflow is retained under workflow_control/.

Generation, entry ranking/export and submissions remain with their agreed owners.
This workflow does not start or stop their resources. Recorded tracked spend is
Runpod $56.589988 + our Modal estimate $10.718947 =~$67.31; new GPU/cloud CPU cost
for this local workflow is $0.00. Other lane costs/storage/builds are excluded.
