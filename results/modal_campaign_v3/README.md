# Modal migration, 2026-10-21

Ricardo explicitly authorized other pods/providers, paying more, and using his
existing Modal account. Authentication to existing profile `pinho-ricardo` was
verified with read-only app/volume inventory. No account, payment method, or
submission was created. Unrelated existing Modal resources are untouched.

The original Runpod queue was cancelled. Its second read had failed with an
upstream Runpod authentication error (500), separate from host capacity. The three
known Runpod pods were last verified EXITED; fresh API reads are currently failing.
No new Runpod start was attempted during this provider migration.

Same BindCraft2 v1.0.4 revision `b6fb3ae15ab6a702f1156cffc11bcc3323c2fd06`,
same seven AF2 model checkpoints and negative ProteinMPNN checkpoint, checked by
SHA256 against the previous campaign metadata. CUDA12 extra follows the upstream
container recipe. Fresh dependency resolution is recorded in runtime-freeze.txt;
hardware and runtime differences prevent calling this a perfectly paired control.
The prepared native-target experiment and its finite quality gates are preserved.
No accepted new scaffold or experimental improvement is claimed.

Source/dependencies/AF2 weights build on CPU **before** a GPU is allocated.
One single-use function container, H100 preferred / A100-80GB fallback, maximum
one container, zero retries, no warm pool, two-hour function timeout, inner worker
deadline 115 minutes. Artifacts commit every 30 seconds. Both worker and function
terminate at completion; the local driver independently verifies zero runners
before starting CPU postprocessing. Cancellation uses exact call ID with
terminate_containers=True, not worker-only shutdown.

## Continued modeling without another user checkpoint

`drive_modal_series.py` runs three bounded batches, each two trajectories, with
independent campaign seeds 20261021, 20261022 and 20261023. The first preserves
the original native-target control; the next two estimate sampling variability.
No acceptance gates, target representation, epitope or losses are relaxed. This
is additional sampling, not evidence by itself of improved model quality.
Each batch releases its actual GPU container before local CPU postprocessing.
Execution failure, unverified release or failed CPU postprocessing halts further
allocation for diagnosis. Existing run/call IDs resume safely; never overwrite
or delete output directories to bypass duplicate-run protection.

The finite series is a six-hour GPU-runtime envelope, not an indefinite
allocation. This implements Ricardo's repeated request to keep modeling and
his explicit higher-spend authorization without asking him to restart each run.
After the series, retain results and stop all GPU containers. Further changes
must be justified by the saved model diagnostics, not lower acceptance gates.

## Cost and budget

Fresh [Modal prices](https://modal.com/pricing), read 2026-10-21:

- H100: $0.001097/sec = **$3.9492/hour** GPU.
- A100-80GB: $0.000694/sec = **$2.4984/hour** GPU.
- Four physical CPU cores: $0.18864/hour; 64 GiB RAM: $0.511488/hour.
- Requested runtime total: **$4.649328/hour H100**, **$3.198528/hour A100**.
- Two-hour function envelope: $9.30 H100 or $6.40 A100, plus image build,
  container startup/teardown and storage. This is a runtime estimate, not an
  enforced account-wide dollar cap. Startup timeout is five minutes.
- Results volume: $0.09/GiB-month. AF2 weights are in the image, not the volume.
- Three sequential batches: **$27.90 H100 / $19.20 A100 maximum requested
  function runtime**, plus build/startup/teardown/storage. Only one GPU at once.
- Previous newly supervised Runpod runtime estimate **$1.7952**, historical total
  and network-volume charges separate. The old $3.58 combined allowance is
  superseded for this pilot by Ricardo's explicit instruction to pay more.

GPU preference lists are supported by [Modal's official fallback documentation](https://modal.com/docs/examples/gpu_fallbacks).

## Live state

`deploy.log` is CPU image-build/deployment evidence, not proof of modeling.
`series.json` records the current batch, finite total and completed batch spend.
The first batch uses this directory; later batches use `replicate_2/` and
`replicate_3/`, each with the same state-file layout. `supervisor.json` records
the unique call ID, mirrored remote state, spend
estimate and independently checked runner count. `mirrored/worker.log` and
`mirrored/status.json` prove actual execution. Full heavy trajectories remain on
the dedicated `adaptyv-c4-results` volume; small scientific outputs mirror locally.
Deployed app definitions with zero runners do not keep GPUs reserved. All27
tests pass, including no next allocation after failed release, finite three-run
sequencing and restart without repeating completed batches.

Emergency stop: create `results/modal_campaign_v3/CANCEL`. The driver cancels
the exact function call and its containers, then checks zero runners. If the
driver is absent, load call_id from supervisor.json and invoke
`modal.FunctionCall.from_id(call_id).cancel(terminate_containers=True)` using the
existing authenticated `.venv-modal/bin/python`.
