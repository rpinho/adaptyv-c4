# C4 rank-1 PDB extraction report

Final state verified: 2026-10-28T23:36:52.736022+00:00

**Outcome: incomplete. 0 of 342 requested PDBs recovered.** The destination `pdbs/` exists and is empty. All six pods are stopped (`EXITED`), and none have been terminated because extraction could not be verified.

| Pod | ID | Expected PDBs | Recovered | Result |
| --- | --- | ---: | ---: | --- |
| pod0 | ubrhji86dj3eqr | 61 | 0 | GPU capacity: all start attempts rejected |
| pod1 | v2d8dptjotvwoi | 61 | 0 | GPU capacity: all start attempts rejected |
| pod2 | h9jmb5eiu2hzwr | 60 | 0 | GPU capacity: all start attempts rejected |
| pod3 | qnfur45hmosdqu | 60 | 0 | GPU capacity: all start attempts rejected |
| pod4 | 63j0e3dmecmaxh | 50 | 0 | GPU capacity: all start attempts rejected |
| pod5 | injkt4l8tmdeex | 50 | 0 | Output directory absent; /workspace only contains .cache |

## Findings

All six initial authenticated start requests failed with HTTP 400: "There are not enough free GPUs on the host machine to start this pod." Six bounded retry rounds were performed. Pod5 became available in round 6; its start was accepted, and the retry loop exited for immediate inspection. A subsequent start attempt for each of pods 0–4 still failed. No retry process remains running.

Pod5 started at `2026-10-28T23:34:32.898Z`. SSH became available at `194.68.245.7:22078`, using the supplied key. Read-only inspection found that `/workspace/rescore_out` does not exist and `/workspace` contains only `.cache`. A bounded search of `/workspace`, `/root`, `/tmp`, `/mnt` and `/runpod-volume` found no `*rank_001*.pdb` files. Missing roots were tolerated; the search is not a whole-disk forensic recovery. The pod was promptly stopped and verified EXITED by `2026-10-28T23:35:25.775052+00:00`.

Every pod's API record reports `disk: 80` and `mounts: {}`. RunPod documents that container-disk data is cleared on stop and that persistent mounts are needed to preserve data. Pod5's absent output is consistent with the files having been on ephemeral container storage before the earlier stop. This is an inference; the other five workspaces remain uninspected. See [RunPod storage options](https://docs.runpod.io/pods/storage/types).

The original manifest and redistribution plan reconcile to **342 distinct designs** in groups of 61, 61, 60, 60, 50 and 50. Exact names are saved in `pdb_recovery/expected_designs.json`. A separate local audit found 341 per-design score files, with `tfr1_c4_acid_site2_96` absent from pod5's assigned designs. Score availability does not prove PDB availability.

## Access and next recovery step

The specified `~/.config/runpod/api_key` does not exist. The existing RunPod OAuth session was used through `scripts/runpod_control.py`; no credential was printed, copied into this repository, or changed.

For pods 0–4, RunPod documents a CPU Only mode for retrieving data when their GPUs are occupied. In the console, choose Start → CPU Only. With a RunPod API key, the GraphQL operation is:

```graphql
mutation {
  podResumeZeroGpu(input: { podId: "POD_ID" }) {
    id
    desiredStatus
  }
}
```

The current MCP action schema and REST start action cannot invoke this mode. An API key at the supplied path, or CPU Only console starts, is needed to inspect the remaining five pods without waiting for GPUs. A request for this missing access was sent during the work. See [RunPod zero-GPU recovery](https://docs.runpod.io/pods/troubleshooting/zero-gpus).

Pod5 requires another existing copy of the PDBs, or additional recovery evidence, because the requested files were absent. No scoring or regeneration was attempted, as instructed.

## Cost

Estimated additional compute for this recovery attempt: **at most $0.0087**, calculated from an upper bound of 52.877 seconds between pod5's reported start and the confirmed stopped-state observation, at the API-reported $0.59/h rate. All other start requests failed. The practical incremental estimate is **approximately $0.01**; exact disk charges and the final billed total are not yet available.

The final billing API snapshot for 2026-10-28 23:00–24:00 UTC reports **$0.00**, with no records for any of the six pods. This is not a final invoice and may reflect billing lag; it must not be read as proof of zero cost. Historical scoring spend is outside this recovery job.

## Preservation and evidence

No scoring, computation jobs, remote file writes, or termination requests were issued. All six pods were independently verified EXITED after the final start attempts. The conditional requirement to terminate after verified extraction could not be met. No local files were deleted.

Evidence in `pdb_recovery/`:

- `recovery_state.json`: consolidated final result and cost estimate.
- `final_pod_states.json`: all six final EXITED states and storage configuration.
- `pod5_remote_inspection.txt`: actual SSH filesystem listing.
- `pod5_stopped_state.json`: confirmed stopped state and timestamp.
- `pod0_recovery.json` through `pod4_recovery.json`: final rejected start attempts.
- `capacity_attempts.jsonl`: six timestamped retry rounds and pod5 start acceptance.
- `capacity_retry_state.json`: historical retry-loop snapshot.
- `final_billing_snapshot.json`: current-hour billing evidence.
- `expected_designs.json` and `local_score_audit.json`: expected assignments and score-file audit.
