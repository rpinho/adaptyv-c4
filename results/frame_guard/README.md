# Astra-confirmed initialization defect and minimal repair

Independent review: ../astra_review/REVIEW.md. The69-residue v5 binder collapsed
to one CA position. Zero N/CA/C coordinates produced a zero rotation matrix;
composing updates cannot escape that singular state. This explains v5, not v3.
The v5 search and successor length/weight variants remain cancelled, with GPU
release verified. Latest conservative cumulative Modal runtime estimate
$3.920816; earlier supervised Runpod$1.7952. Extra billing is unreconciled.

The version-locked patch substitutes canonical noncollinear initialization-only
N/CA/C positions for invalid/unobserved residue frames. Actual input positions,
template features and masks remain as supplied. Valid target frames are retained;
partial residues preserve a finite observed CA translation. Both first-input and
recycling paths are guarded, with sequence-derived masks for recycled predictions.

`CPU_VALIDATION.json` executes both actual pinned upstream frame conversions on
JAX0.11.2 CPU. Old zero-frame determinant/update norm are0/0. Repaired matrices
are finite and orthonormal with determinant+1; translation updates survive.
Valid, missing, partial, collinear, nonfinite, padded and recycled cases pass.
The same regression runs in the Linux image before GPU allocation. Source and
helper SHA256 hashes are recorded. These are numerical geometry tests, not a
full AlphaFold forward or evidence of binding. Astra's addendum records limits.

One subsequent integration diagnostic is authorized and launched, with state
in ../modal_frame_diagnostic_v8/. It uses the same initial seed/trajectory and
69-residue binder as failed v5, the local-frame complete639-residue target,
same objective parameters, one actual forward/sequence-gradient call and no
optimization campaign. It records finite gradients, noncollapse, physical
hotspot distance and patch/source hashes. It never enters candidate pH scoring.
No automatic design successor is queued. A software pass is not a binder pass.

Completion update: the real forward/gradient passed the narrow software check,
and its GPU was released. Full measurements and limits are in `GPU_RESULTS.md`.
The output is still off-patch and has abnormal adjacent-CA geometry. Follow-up
Opus/Fable consultation and stronger CPU regressions are recorded in
`../consultations/REVIEW.md`; no new helper patch was justified by that review.

H100 GPU$3.9492/hour; requested CPU/RAM runtime upper$4.649328/hour. One
25-minute subprocess timeout inside a30-minute function, with5-minute startup
allowance, single-use container, zero retries, exact-call cancellation and
independent zero-runner verification. Expected runtime cap about$2.32 plus
startup/teardown and CPU image build/storage; these are estimates, not a bill.

Before re-enabling any adaptive design search, enforce CPU native-context
admission rather than raw ranked rows, deduplicate trajectory/backbone hashes,
and route by actual geometry/patch acquisition and explicit failure reasons.
Astra found these controller defects; the cancelled controller is not yet safe
to call a qualifying-scaffold search. Do not restart it as though repaired.

CPU reproduction (fetch the pinned upstream source named in REVIEW.md first):

```sh
python3 scripts/patch_frame_guard.py --upstream /path/to/pristine/BindCraft2 --helper scripts/frame_guard.py
JAX_PLATFORMS=cpu .venv-frameguard/bin/python scripts/validate_frame_guard.py --upstream /path/to/patched/BindCraft2 --out results/frame_guard/CPU_VALIDATION.json
```

Existing local validation source cache: /tmp/c4-upstream. Patcher refuses an
unexpected/already-modified AF2 file. Normal local regressions:38 passed.
