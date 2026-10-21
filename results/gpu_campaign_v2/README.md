# Current GPU attempt and automatic continuation

The first full-dimer attempt could not fit in 80 GB even with narrower shards and
chunked attention. It was mirrored, the POD stopped and EXITED verified. Its
estimated compute was $0.3855; no accepted scaffold. The old supervisor files in
`../gpu_campaign/` describe that completed attempt.

This attempt keeps the **entire native human monomer, 639 residues**. No internal
crop. The other subunit is absent during generation: an explicit approximation.
Before CPU scoring, `postprocess_gpu_campaign.py` restores that subunit, checks
full native sequence and peptide connectivity, and rejects backbone/NAG overlap.
Glycan extensions, assay scar and tag are not modeled. No claim of verified mouse
cross-reactivity or experimental pH selectivity.

Active pod: `2p9x2qea1ietbw`, A100-SXM4-80GB. Fresh catalog and pod reads:
**$1.79/hour compute**. Overall cap **$3.58 across both attempts**, with the prior
$0.3855 deducted and three minutes reserved for shutdown/retries. Persistent
storage is additional; the pre-start account idle rate was $0.019/hour. Runtime
estimates are labeled; billing may lag by an hour.

- `supervisor.json`: live worker state, total estimated compute, deadline and pod stop verification.
- `worker_tail.log`: mirrored generation log.
- `mirrored/`: saved campaign outputs, without compile caches.
- `supervisor_config.json`: PID/connection/project configuration for this attempt.
- `continuation.json`: automatically starts CPU follow-up after verified pod shutdown.
- `../gpu_followup_v2/status.json`: appears when that CPU follow-up starts.

Do not launch duplicate workers/monitors. Do not terminate either original pod.
No account creation or challenge submission is delegated. All changes are off main.

Emergency stop and fresh verification:

```sh
python3 scripts/runpod_control.py --tool pod-action --arguments '{"id":"2p9x2qea1ietbw","body":{"action":"stop"}}'
python3 scripts/runpod_control.py --tool get-pod --arguments '{"id":"2p9x2qea1ietbw"}'
```

The default standalone control output redacts environment/authentication values.
The supervisor confirms actual EXITED state, including on worker exit or budget
expiry. `STOP_FAILED_REQUIRES_ATTENTION` is an honest shutdown failure requiring
action. It never treats an exited Python worker as proof that billing stopped.
