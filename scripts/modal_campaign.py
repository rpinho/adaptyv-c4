"""Pinned C4 pilot on Modal; CPU image build precedes every GPU allocation.

Deploy: .venv-modal/bin/modal deploy scripts/modal_campaign.py
Launch/monitor series: .venv-modal/bin/python scripts/drive_modal_series.py
"""
from pathlib import Path

import modal

ROOT = Path(__file__).resolve().parents[1]
REVISION = "b6fb3ae15ab6a702f1156cffc11bcc3323c2fd06"
APP_NAME = "adaptyv-c4-native-control"
VOLUME_NAME = "adaptyv-c4-results"
RUN_NAME = "c4_native_apical_ablation_v3"
LARGE_INIT_RUN = "c4_native_largeinit_v4"
app = modal.App(APP_NAME)
volume = modal.Volume.from_name(VOLUME_NAME, create_if_missing=True)

# Mirrors upstream's container recipe, retaining the previous CUDA 12 family.
# Source, dependencies and AF2 weights are fetched on CPU, never a rented GPU.
image = (
    modal.Image.debian_slim(python_version="3.12")
    .apt_install("git", "curl", "build-essential", "libgomp1")
    .run_commands(
        "git init /opt/bindcraft && cd /opt/bindcraft && "
        "git remote add origin https://github.com/PacesaLab/BindCraft2.git && "
        f"git fetch --depth 1 origin {REVISION} && git checkout --detach FETCH_HEAD"
    )
    .run_commands("python -m pip install -e '/opt/bindcraft[cuda12]'")
    .env({
        "NVIDIA_VISIBLE_DEVICES": "all",
        "NVIDIA_DRIVER_CAPABILITIES": "compute,utility",
        "BINDCRAFT_WEIGHTS": "/opt/bindcraft-weights",
        "PYTHONUNBUFFERED": "1",
        "XLA_PYTHON_CLIENT_PREALLOCATE": "false",
        "OMP_NUM_THREADS": "8",
    })
    .run_commands(
        "python -c \"import nvidia,pathlib; print('\\n'.join(sorted(str(p) "
        "for r in nvidia.__path__ for p in pathlib.Path(r).glob('*/lib'))))\" "
        "> /etc/ld.so.conf.d/bindcraft-cuda.conf && ldconfig",
        "bindcraft fetch-weights",
        "python -m bindcraft.selfcheck cuda12 --shipped-only",
        "python -m pip freeze > /opt/bindcraft-runtime-freeze.txt",
    )
    .add_local_file(ROOT / "results/gpu_campaign_v3/campaign.json", "/opt/c4/campaign.json", copy=True)
    .add_local_file(ROOT / "results/gpu_campaign_v2/human_11qc_full_monomer.pdb", "/opt/c4/target.pdb", copy=True)
    .add_local_file(ROOT / "results/gpu_campaign_v2/mirrored/campaign_metadata.json", "/opt/c4/reference_metadata.json", copy=True)
    .add_local_file(ROOT / "results/modal_campaign_v5/human_11qc_local_frame.pdb", "/opt/c4/local-frame-target.pdb", copy=True)
    .add_local_file(ROOT / "scripts/modal_search_policy.py", "/opt/c4/search_policy.py", copy=True)
)


@app.function(
    image=image,
    gpu=["H100", "A100-80GB"],
    cpu=(4, 4),
    memory=(65536, 65536),
    volumes={"/results": volume},
    timeout=7200,
    startup_timeout=300,
    retries=0,
    max_containers=1,
    min_containers=0,
    scaledown_window=2,
    single_use_containers=True,
)
def run_campaign(run_name: str = RUN_NAME, overrides: dict | None = None):
    import hashlib
    import importlib.util
    import json
    import os
    import re
    import shutil
    import subprocess
    import time

    import jax
    from bindcraft.model_weights import alphafold_parameter_file, alphafold_parameters

    spec = importlib.util.spec_from_file_location('search_policy','/opt/c4/search_policy.py')
    policy = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(policy)
    policy.profile_for_name(run_name)
    overrides = overrides or {}
    if set(overrides) - {"campaign_seed"}:
        raise ValueError("Only a new independent seed is allowed in this series")
    started = time.time()
    run = Path("/results") / run_name
    if run.exists():
        raise RuntimeError("Run directory already exists; inspect it before any duplicate run")
    run.mkdir()

    def save(name, data):
        path = run / name
        tmp = path.with_suffix(path.suffix + ".tmp")
        tmp.write_text(json.dumps(data, indent=2) + "\n")
        tmp.replace(path)

    proc = None
    status = {"status": "preflight", "started": started, "revision": REVISION}
    try:
        devices = jax.devices()
        if not devices or any(d.platform != "gpu" for d in devices):
            raise RuntimeError(f"Refusing CPU fallback: {devices}")
        gpu = subprocess.check_output([
            "nvidia-smi", "--query-gpu=name,memory.total,driver_version", "--format=csv,noheader"
        ], text=True).strip()
        if int(gpu.split(",")[1].strip().split()[0]) < 79000:
            raise RuntimeError(f"GPU has insufficient memory: {gpu}")

        reference = json.loads(Path("/opt/c4/reference_metadata.json").read_text())
        checksums = {}
        params = alphafold_parameters(download=False)
        for old_path, expected in reference["checkpoints"].items():
            name = Path(old_path).name
            if name.startswith("params_"):
                path = alphafold_parameter_file(params, name.removeprefix("params_").removesuffix(".npz"))
            else:
                path = "/opt/bindcraft/bindcraft/weights/proteinmpnn/weights_negative/" + name
            if not path:
                raise RuntimeError(f"Missing checkpoint: {name}")
            actual = hashlib.sha256(Path(path).read_bytes()).hexdigest()
            if actual != expected:
                raise RuntimeError(f"Checkpoint identity mismatch: {name}")
            checksums[name] = actual

        config = json.loads(Path("/opt/c4/campaign.json").read_text())
        config = policy.make_config(config,run_name,overrides,str(run / 'campaign'))
        (run / "campaign.json").write_text(json.dumps(config, indent=2) + "\n")
        shutil.copy(config['targets'][0]['target_path'], run / "target.pdb")
        shutil.copy("/opt/bindcraft-runtime-freeze.txt", run / "runtime-freeze.txt")
        status.update(status="modeling", gpu=gpu, jax=jax.__version__, checkpoints=checksums,
                      container=os.environ.get("MODAL_TASK_ID"), max_runtime_seconds=7200)
        save("status.json", status)
        volume.commit()
        print(json.dumps(status), flush=True)
        with (run / "worker.log").open("w") as log:
            proc = subprocess.Popen(["bindcraft", "design", str(run / "campaign.json")],
                                    cwd="/opt/bindcraft", stdout=log, stderr=subprocess.STDOUT)
            deadline = started + 6900
            while proc.poll() is None:
                if time.time() >= deadline:
                    proc.terminate()
                    try:
                        proc.wait(timeout=30)
                    except subprocess.TimeoutExpired:
                        proc.kill()
                        proc.wait(timeout=10)
                    raise TimeoutError("Finite GPU pilot runtime reached")
                status.update(elapsed_seconds=time.time() - started, worker_pid=proc.pid)
                save("status.json", status)
                volume.commit()
                time.sleep(30)
        status.update(status="completed" if proc.returncode == 0 else "worker_failed", returncode=proc.returncode)
    except BaseException as exc:
        status.update(status="failed", error=repr(exc))
        if proc is not None and proc.poll() is None:
            proc.kill()
            proc.wait(timeout=10)
        raise
    finally:
        status.update(finished=time.time(), elapsed_seconds=time.time() - started,
                      shutdown="function returning; single-use container must terminate")
        save("status.json", status)
        volume.commit()
        print(json.dumps(status), flush=True)
    return status
