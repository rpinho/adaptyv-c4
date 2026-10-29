# Why no qualifying candidate yet — 2026-10-22 00:14 UTC

There is no candidate that has passed the current computational qualification
pipeline, and no experimentally validated binder has been established.

1. The original apparent candidate pool used a disconnected73-residue receptor
   crop. Restoring the full receptor produced backbone overlaps for all65
   well-aligned source poses. Their pKa scores cannot qualify them as binders
   of the native receptor. Evidence: ../full_context_audit/RESULTS.md.
2. Six completed native-target Modal trajectories produced no accepted designs.
   Five failed the initial fold-confidence screen; one passed that screen and
   failed refine. Consequently none reached sequence redesign, final complex
   qualification or CPU switch scoring. Evidence: ../modal_campaign_v3/RESULTS.md.
3. The subsequent coordinate-initialization pilot had an avoidable relative
   placement defect and was cancelled before a completed comparison. GPU
   release was verified. That defect does not explain the earlier v3 failures.
   Evidence: ../modal_campaign_v4/FRAME_AUDIT.md.

The observed immediate failure is low predicted binder fold confidence, not
an experimental observation of absent binding. The deeper cause has not been
established: limited stochastic sampling, generation settings and target/epitope
difficulty remain alternatives. Six unsuccessful trajectories do not prove
that the epitope is impossible. No measured affinity, mouse cross-reactivity or
pH-selective release has been demonstrated. More candidate generation is not
guaranteed improvement; these are proposals from pretrained models, not evidence
that model weights improve on every run.

At this checkpoint the corrected-frame pilot is active on one H10080GB, with
one running input, zero backlog and its first trajectory started. Its inputs
preserve the entire native target and pass the recorded geometry/identity
checks. It has not yet produced a new accepted design. The failure-dependent
policy, finite limits and unchanged acceptance criteria are in DECISION.md.

Estimated newly supervised Runpod runtime$1.7952. Modal runtime estimate through
the current submission approximately$3.60, including the prior completed and
cancelled attempts; setup, storage and earlier campaigns remain separate and
unreconciled. Live values supersede this snapshot: search.json and
batch_1/supervisor.json.
