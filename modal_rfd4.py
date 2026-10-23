"""Run RFD4-Proteina on Modal for the TfR1 Challenge 4 entry.

Supports multiple specs: the original apical binder and the E1 (His707) patch.

Usage:
    modal run modal_rfd4.py::preflight_spec --spec e1_plain
    modal run modal_rfd4.py::generate --spec e1_plain --num-replicates 200
    modal run modal_rfd4.py::generate --spec e1_his015 --num-replicates 200
"""

import os
import pathlib
import subprocess

import modal

RFD4_REPO = pathlib.Path("/Users/ricardo/dev/RFD4-Proteina")
TARGET_CIF = pathlib.Path("/Users/ricardo/dev/adaptyv-c4/data/11qc.cif")
SPECS_DIR = pathlib.Path("/Users/ricardo/dev/adaptyv-c4/configs/rfd4")

INPUTS = "/workspace/inputs"
RUNS = "/runs"

GPU_CHEAP = "L40S"

SPEC_MAP = {
    "apical": "tfr1_apical_binder.yaml",
    "e1_plain": "tfr1_e1_plain.yaml",
    "e1_his015": "tfr1_e1_his015.yaml",
}

app = modal.App("rfd4-tfr1")

hf_cache = modal.Volume.from_name("rfd4-hf-cache", create_if_missing=True)
runs_vol = modal.Volume.from_name("rfd4-runs", create_if_missing=True)

image = (
    modal.Image.from_dockerfile(
        str(RFD4_REPO / "Dockerfile"),
        context_dir=str(RFD4_REPO),
        build_args={"RFD4_ACCELERATOR": "cu13"},
    )
    .env({"HF_HOME": "/cache/huggingface", "DATA_PATH": "/data"})
    .add_local_file(str(TARGET_CIF), f"{INPUTS}/11qc.cif")
    .add_local_file(str(SPECS_DIR / "tfr1_apical_binder.yaml"), f"{INPUTS}/tfr1_apical_binder.yaml")
    .add_local_file(str(SPECS_DIR / "tfr1_e1_plain.yaml"), f"{INPUTS}/tfr1_e1_plain.yaml")
    .add_local_file(str(SPECS_DIR / "tfr1_e1_his015.yaml"), f"{INPUTS}/tfr1_e1_his015.yaml")
)

VOLUMES = {"/cache/huggingface": hf_cache, RUNS: runs_vol}


def sh(cmd, cwd="/workspace", check=True):
    print("+ " + cmd, flush=True)
    p = subprocess.run(cmd, shell=True, cwd=cwd, text=True,
                       stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    print(p.stdout, flush=True)
    if check and p.returncode != 0:
        raise RuntimeError("command failed (%d): %s" % (p.returncode, cmd))
    return p.returncode, p.stdout


def resolve_spec(spec: str) -> str:
    fname = SPEC_MAP.get(spec, spec)
    path = f"{INPUTS}/{fname}"
    if not os.path.exists(path):
        raise FileNotFoundError(f"Spec not found: {path} (known specs: {list(SPEC_MAP)})")
    return path


@app.function(image=image, gpu=GPU_CHEAP, volumes=VOLUMES, timeout=3600)
def preflight_spec(spec: str = "apical"):
    """One design with prevalidation — validates the config on a real GPU for ~$2."""
    spec_path = resolve_spec(spec)
    out_name = f"preflight_{spec}"
    sh("nvidia-smi", cwd="/", check=False)
    sh(f"rfd4-generate inputs={spec_path} "
       f"recycle_strategy=medium/front_loaded "
       f"num_replicates=1 batch_size=1 seed=42 "
       f"save_chain_of_interest=true prevalidate_inputs=true "
       f"out_dir={RUNS}/{out_name} '~metrics'")
    sh(f"find {RUNS}/{out_name} -name '*.cif*' | head -5", cwd="/", check=False)
    runs_vol.commit()
    return out_name


@app.function(image=image, gpu=GPU_CHEAP, volumes=VOLUMES, timeout=3600)
def preflight():
    sh("nvidia-smi", cwd="/", check=False)
    sh("python -c \"import torch, rfproteina, atomworks; "
       "print('torch', torch.__version__); "
       "print('cuda available', torch.cuda.is_available()); "
       "print('device', torch.cuda.get_device_name(0) if torch.cuda.is_available() else 'NONE'); "
       "print('vram GiB', round(torch.cuda.get_device_properties(0).total_memory/2**30, 1) "
       "if torch.cuda.is_available() else 0)\"")
    sh("rfd4-generate --help | head -40", check=False)
    hf_cache.commit()


@app.function(image=image, gpu=GPU_CHEAP, volumes=VOLUMES, timeout=86400)
def generate(spec: str = "apical", num_replicates: int = 8, batch_size: int = 1,
             seed: int = 42, recycle_strategy: str = "medium/front_loaded",
             out_name: str = "", metrics: bool = False):
    """Design binders against the specified epitope."""
    spec_path = resolve_spec(spec)
    if not out_name:
        out_name = f"tfr1_{spec}"
    metrics_arg = ("metrics_tags='[ppi_refolding,hotspot_conditioning]'" if metrics
                   else "'~metrics'")
    sh(f"rfd4-generate inputs={spec_path} "
       f"recycle_strategy={recycle_strategy} "
       f"num_replicates={num_replicates} batch_size={batch_size} seed={seed} "
       f"save_chain_of_interest=true prevalidate_inputs=true "
       f"out_dir={RUNS}/{out_name} {metrics_arg}")
    sh("nvidia-smi --query-gpu=memory.used,memory.total --format=csv", cwd="/", check=False)
    sh(f"find {RUNS}/{out_name} -name '*.cif*' | wc -l", cwd="/", check=False)
    sh(f"find {RUNS}/{out_name} -name 'per_example_metrics.csv' -exec head -3 {{}} +",
       cwd="/", check=False)
    runs_vol.commit()


@app.function(image=image, volumes=VOLUMES, timeout=1800, cpu=2.0)
def fetch(path: str = ""):
    """List or dump what is on the runs volume."""
    sh(f"find {RUNS}/{path} -maxdepth 4 | head -100", cwd="/", check=False)
    sh(f"du -sh {RUNS}/*", cwd="/", check=False)


@app.function(image=image, volumes=VOLUMES, timeout=600, cpu=2.0)
def inspect_cif(path: str):
    """Print chain IDs and residue ranges from a generated CIF."""
    import biotite.structure.io.pdbx as pdbx
    import biotite.structure as struc
    full = f"{RUNS}/{path}"
    print(f"Reading {full}")
    f = pdbx.CIFFile.read(full)
    block = f.block
    arr = pdbx.get_structure(block, model=1)
    ca = arr[arr.atom_name == "CA"]
    chains = {}
    for a in ca:
        cid = a.chain_id
        rid = a.res_id
        if cid not in chains:
            chains[cid] = []
        chains[cid].append(rid)
    for cid in sorted(chains):
        rids = sorted(chains[cid])
        print(f"  chain {cid}: {len(rids)} residues, range {rids[0]}-{rids[-1]}")
    print(f"Total CA atoms: {len(ca)}")
