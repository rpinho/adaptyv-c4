# C4 rank-1 PDB recovery status

Final check: 2026-10-28T23:36:52.736022+00:00

**Incomplete: 0 / 342 PDBs recovered. All six pods are EXITED; none terminated.**

| Pod | ID | Expected | Recovered | Final state | Failure |
| --- | --- | ---: | ---: | --- | --- |
| pod0 | ubrhji86dj3eqr | 61 | 0 | EXITED | GPU capacity: all start attempts rejected |
| pod1 | v2d8dptjotvwoi | 61 | 0 | EXITED | GPU capacity: all start attempts rejected |
| pod2 | h9jmb5eiu2hzwr | 60 | 0 | EXITED | GPU capacity: all start attempts rejected |
| pod3 | qnfur45hmosdqu | 60 | 0 | EXITED | GPU capacity: all start attempts rejected |
| pod4 | 63j0e3dmecmaxh | 50 | 0 | EXITED | GPU capacity: all start attempts rejected |
| pod5 | injkt4l8tmdeex | 50 | 0 | EXITED | Output directory absent; /workspace only contains .cache |

Pod5 started at 2026-10-28T23:34:32.898Z and was confirmed stopped at 2026-10-28T23:35:25.775052+00:00. Its new SSH endpoint was 194.68.245.7:22078. Read-only inspection found no `/workspace/rescore_out` directory and no rank-001 PDBs in a bounded search of /workspace, /root, /tmp, /mnt and /runpod-volume.

All six API records show `mounts: {}`. Pod5's empty workspace is consistent with loss of container-disk files on the earlier stop; the other five pods could not be inspected.

Normal starts were retried in six rounds, then attempted again for pods 0–4 after pod5 was stopped. Every start for pods 0–4 failed with insufficient host GPU capacity. No background retry process remains.

The specified API-key file is absent. Existing authenticated RunPod OAuth access was used through `scripts/runpod_control.py`. CPU Only recovery is documented by RunPod but unavailable through the current MCP action schema; it requires an API key for GraphQL or CPU Only starts in the console. A request for that missing access was sent to the user.

No scoring, computation jobs, remote file writes, or termination requests were issued. Termination remains gated on verified extraction.

Estimated incremental compute: at most $0.0087 (52.877s × $0.59/h), approximately $0.01. Current-hour billing API still reports $0.00 with no records; final billed total is not yet available.

See REPORT.md and `pdb_recovery/recovery_state.json` for final findings. `capacity_retry_state.json` is a historical snapshot of the retry loop, not the final pod state.
