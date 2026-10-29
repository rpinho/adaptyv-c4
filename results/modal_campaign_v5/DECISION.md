# Correct frame first; bounded failure-dependent search

The v4 pilot was cancelled after the CPU audit found an avoidable coordinate
placement error. Actual GPU release was verified. See ../modal_campaign_v4/FRAME_AUDIT.md.
No conclusion is drawn about v4 acceptance, or the cause of v3 failure, from
that aborted experiment.

The complete native639-residue target is rigidly translated so the unresolved
binder starts15A outside the intended patch, with7.955A nearest-atom clearance
against the full receptor dimer. Atom identities, sequence, numbering and all
CA pair distances are preserved. Exact audit: human_11qc_local_frame.audit.json.
No target crop, lysine replacement, mutation or peptide reconnection occurs.

`drive_modal_adaptive.py` replaces the previous pre-set repeats with a finite,
declared failure-dependent policy. It does not perform open-ended scientific
reasoning or guarantee improvement:

1. Corrected frame plus large-complex initialization, same65–75 length range.
2. If both trajectories stop at screen/refine with no accepted design, change
   to90–110 residues to test binder foldability. Keep all acceptance gates.
3. If that also fails at screen/refine, increase the pLDDT optimization weight
   from0.1 to0.5 with the same longer range; final confidence gates stay fixed.
4. A passing configuration gets an independent seed for scaffold diversity.
   A later-stage/interface failure stops this fold-focused policy for diagnosis.

At most6 calls or3 accepted computational scaffolds; each call has two
trajectories and a two-hour function cap, zero retries and one single-use GPU.
GPU release is independently verified before CPU analysis and any next call.
Each decision and its observed reason are written to search.json. No rejected
trajectory is promoted, and no pH model runs unless full-native receptor,
opposite-subunit and clash checks pass. Predicted acceptance is not affinity.

Spend envelope is **$27.90 cumulative Modal runtime upper estimate**, including
the already estimated$2.8760 v3 plus$0.2559 aborted-v4 runtime. A complete7800s
wall-time allowance at$4.649328/hour is reserved before each new call. Actual
submission-to-release wall time is charged conservatively to that estimate,
including queue/startup. Image build and storage are extra and not reconciled.
Previous newly supervised Runpod runtime is$1.7952, separate. Existing quoted
H100 GPU rate$3.9492/hour; A100 fallback$2.4984/hour before CPU/RAM. Ricardo's
existing higher-spend/continued-modeling authorization applies.

This policy has finite limits. It stops and releases GPUs on execution failure,
failed release, CPU follow-up failure, exhausted declared hypotheses or budget
reserve. It must never be described as continuous autonomous research beyond
these limits. An eventual stop requires a fresh scientific assessment, not
claiming an agent is still working or silently repeating the same failures.

Live pointer: search.json, with per-call state under batch_1/, batch_2/, etc.
Create CANCEL here to cancel the active call and its actual container.
No account creation, submission, repo main or KB main changes.
