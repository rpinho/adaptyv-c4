# Opus / Fable consultation — 2026-10-22 UTC

Ricardo explicitly authorized these consultations. Two independent, transparent
requests were sent through the existing authenticated Claude CLI, with tools and
MCP disabled, no session persistence, and a $2 per-request budget. The prompts
identify the target and legitimate biomedical challenge and request software /
numerical review. No credentials or designed sequences were included. Original
prompts and JSON responses are retained beside this file.

## Actual outcomes and usage

| Requested model | Outcome | Reported list-price usage |
|---|---|---:|
| `claude-opus-5-5` | Biology safeguard block (`[bio]`), `is_error=true`; no usable review | $0.5059474 |
| `claude-fable-5-1[1m]` | Review returned, `is_error=false` | $0.6863525 |

Total reported usage: **$1.1922999**. This is the CLI's `total_cost_usd`, not a
reconciled invoice or proof of an incremental subscription charge. Both responses
report usage by the requested family **and `claude-opus-5`**. Thus the successful
response is the Fable-request route's review, not an independently verified
pure-Fable completion. The blocked request was not retried, concealed or rephrased
to get around the block. CLI exit status alone is insufficient: Opus returned an
error inside a result envelope whose subtype is `success`.

## Review adjudication

The review's strongest useful point is experimental calibration: a finite,
noncollapsed initial prediction establishes a narrow software integration result;
it does not establish physical fold quality, patch acquisition, or a binder.
There is no matched guard-off GPU forward in this diagnostic. V5's collapsed
50-step trajectory is a different endpoint, and the comparisons must not be
presented as a paired causal experiment.

| Review concern | Check against actual code / evidence | Disposition |
|---|---|---|
| NaN gradients survive `where` on recycled coordinates | The pinned recycle wrapper explicitly applies `jax.lax.stop_gradient` before the guard. New CPU tests nevertheless differentiate through the guard and both actual upstream frame converters on zero, NaN-N, NaN-CA and valid inputs. All are finite. | Not reproduced; no speculative helper rewrite. This does not claim arbitrary downstream NaN computations are safe. |
| First call's `atom_mask` might be undefined or the wrong binding | AST checks confirm `atoms` and `atom_mask` are explicit arguments of `alphafold_input_features`. The completed real forward also exercised this path. | Not a demonstrated binding bug. |
| Invalid-frame fallback resets sidechain coordinates in the initialization copy | True of the copy. Both pinned `folding.py` and `folding_multimer.py` consume only N/CA/C from this field. The separate all-atom positions, template features and masks are untouched. | No sidechain consumer bug demonstrated in this pinned source. Document scope precisely. |
| CPU preservation / recycling tests were too weak | Added exact valid-frame rotation and translation preservation for both converters, direct CA translation assertions, and a degenerate recycled-input test under an all-present sequence mask. | Improved; `CPU_REVIEW_VALIDATION.json` passes on JAX 0.11.2 CPU. |
| Diagnostic geometry gate is not physical qualification | Agreed. Actual adjacent-CA minimum is 1.562 A and requested-patch distance is 38.541 A despite 69 distinct CA coordinates. | Keep status strictly software-only; no candidate admission. Suggested universal span bounds were not adopted as scientific gates. |
| Device check rejects mixed device listings | Existing diagnostic deliberately rejected any listed non-GPU backend; actual enumeration was CUDA only. | Conservative behavior, not an observed CPU fallback. Future checks can explicitly verify computation-device placement. |
| Runtime module provenance could be stronger | Manifest is build-time provenance; the successful run used that image, but it does not independently hash each imported module in the running process. | Useful follow-up hardening; no claim it is already implemented. |
| Patcher lacks atomic replacement / rollback | Accurate operational limitation. It is SHA-locked and used in a disposable image layer; already-patched input intentionally fails closed. | Record limitation; no urgent scientific result depends on idempotent repatching. |

The new regression also explicitly shows that the NaN-CA fallback has zero
coordinate gradient, as expected when its translation is replaced by zero. Zero
gradient there is not evidence of useful sequence optimization. The actual GPU
sequence-gradient norm is separately finite and nonzero (6.2586).

## Current project state

The one repaired forward/gradient completed, and its GPU was actually released;
see `../frame_guard/GPU_RESULTS.md`. No qualified candidate exists. Target
acquisition and the cancelled controller's native-context admission / unique
trajectory accounting remain unresolved. The old adaptive controller must stay
cancelled until those defects are repaired. No account, submission, main-branch
or KB-main action was taken.

Estimated newly supervised Runpod runtime: **$1.7952**. Cumulative Modal runtime
estimate: **$4.0531838**, including the completed $0.1323675 diagnostic.
Historical campaigns, unreconciled setup/storage, and this consultation usage
remain separate from those runtime estimates.
