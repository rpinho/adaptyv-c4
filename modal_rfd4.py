"""Run RFD4-Proteina on Modal for the TfR1 Challenge 4 entry.

The image is built from the model repo's OWN Dockerfile rather than a hand-rolled
equivalent. That matters: their file pins the base image by digest, works around the
October 2026 eMudhra/emSign CA incident with a hash-pinned pip, installs only from
public PyPI with `--no-sources`, and ends with a build-time smoke test that asserts the
torch build, the TMol operators, Inductor and the CLI. Reimplementing that would be
strictly worse. CUDA libraries come from the PyTorch wheels, so no CUDA base image is
needed and Modal's driver is enough.

Cost control, in the order these should be run:

  1. dry_run   -- CPU only, no GPU. Materialises the annotated CIF the spec describes
                  so a wrong label_seq_id or an empty selection is caught for cents
                  rather than on a GPU. Upstream's own advice: "For anything
                  complicated, do this. It is the cheapest way to catch a mistake
                  before you spend GPU time on it."
  2. preflight -- cheapest supported GPU. Proves CUDA, the checkpoint pair and the
                  entry point, and populates the weight cache volume.
  3. smoke     -- one design from the shipped monomer benchmark.
  4. baseline  -- the shipped ppi_all_atom example UNCHANGED, to measure VRAM and
                  runtime per design on a known-good input. This is what separates
                  "the model works here" from "our spec is authored correctly".
  5. generate  -- our TfR1 spec.

Usage:
    modal run modal_rfd4.py::dry_run
    modal run modal_rfd4.py::preflight
    modal run modal_rfd4.py::smoke
    modal run modal_rfd4.py::baseline
    modal run modal_rfd4.py::generate --num-replicates 8
"""

import os
import pathlib
import subprocess

import modal

# --------------------------------------------------------------------------------------
# Paths
# --------------------------------------------------------------------------------------
RFD4_REPO = pathlib.Path("/Users/ricardo/dev/RFD4-Proteina")
TARGET_CIF = pathlib.Path("/Users/ricardo/dev/adaptyv-c4/data/11qc.cif")
SPEC = pathlib.Path(__file__).parent / "configs" / "rfd4" / "tfr1_apical_binder.yaml"

REPO_IN_IMAGE = "/opt/rfd4"      # WORKDIR in their Dockerfile is /workspace
INPUTS = "/workspace/inputs"     # spec + target CIF live together: `input:` is relative
RUNS = "/runs"                   # persisted across calls

# L40S is the cheapest GPU on RFD4's tested-hardware list (Ada Lovelace, 48 GB), at
# $0.000542/s ~ $1.95/hr. The docs state no VRAM number -- "Peak GPU memory depends on
# token count, batch size and recycle depth rather than on the weights" -- and the
# lowest card they list is an A100 40GB. Our cropped input is ~315 tokens, well under
# the 1024-token ceiling, so start cheap and move up only if it OOMs.
GPU_CHEAP = "L40S"
GPU_BIG = "A100-80GB"

app = modal.App("rfd4-tfr1")

# Weights (~2.8 GB: 1.74 GB model + 1.03 GB autoencoder) download once and persist.
# Volume storage is $0.09/GiB/month with 1 TiB free, so this is effectively free.
hf_cache = modal.Volume.from_name("rfd4-hf-cache", create_if_missing=True)
runs_vol = modal.Volume.from_name("rfd4-runs", create_if_missing=True)

image = (
    modal.Image.from_dockerfile(
        str(RFD4_REPO / "Dockerfile"),
        context_dir=str(RFD4_REPO),
        build_args={"RFD4_ACCELERATOR": "cu13"},
    )
    # Point HF_HOME at the volume so the checkpoint pair survives container teardown.
    .env({"HF_HOME": "/cache/huggingface", "DATA_PATH": "/data"})
    .add_local_file(str(TARGET_CIF), f"{INPUTS}/11qc.cif")
    .add_local_file(str(SPEC), f"{INPUTS}/tfr1_apical_binder.yaml")
)

VOLUMES = {"/cache/huggingface": hf_cache, RUNS: runs_vol}


def sh(cmd, cwd=REPO_IN_IMAGE, check=True):
    """Run a command, streaming output, and return (returncode, combined output)."""
    print("+ " + cmd, flush=True)
    p = subprocess.run(cmd, shell=True, cwd=cwd, text=True,
                       stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    print(p.stdout, flush=True)
    if check and p.returncode != 0:
        raise RuntimeError("command failed (%d): %s" % (p.returncode, cmd))
    return p.returncode, p.stdout


def sh_timed(cmd, cwd=REPO_IN_IMAGE, check=True):
    """Run a command while sampling GPU memory, and report peak use and wall time.

    Querying nvidia-smi after the command exits reports an idle GPU, which is how the
    first baseline run came back claiming 3 MiB. Upstream states no VRAM figure at all
    ("Peak GPU memory depends on token count, batch size and recycle depth rather than
    on the weights"), so this is the only way we get a real number to size a campaign.
    """
    import time
    sampler = subprocess.Popen(
        "nvidia-smi --query-gpu=memory.used --format=csv,noheader,nounits -l 2",
        shell=True, text=True, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
    t0 = time.time()
    try:
        rc, out = sh(cmd, cwd=cwd, check=check)
    finally:
        elapsed = time.time() - t0
        sampler.terminate()
        try:
            samples = [int(x) for x in (sampler.stdout.read() or "").split() if x.isdigit()]
        except Exception:
            samples = []
    peak = max(samples) if samples else None
    print("", flush=True)
    print("TIMING  wall=%.1fs" % elapsed, flush=True)
    print("VRAM    peak=%s MiB  (of 46068 on L40S)  samples=%d"
          % (peak if peak is not None else "UNMEASURED", len(samples)), flush=True)
    return rc, out, elapsed, peak


# --------------------------------------------------------------------------------------
# 1. CPU-only validation. Run this first, every time the spec changes.
# --------------------------------------------------------------------------------------
@app.function(image=image, volumes=VOLUMES, timeout=1800, cpu=2.0, memory=8192)
def dry_run():
    """Materialise the annotated CIF our spec describes. No GPU, no weights.

    `rfd4-generate` builds this in memory anyway; writing it out is how we confirm the
    conditioning landed on the atoms we meant. A wrong label_seq_id does not raise --
    it silently conditions the wrong residues, or reports "A condition matched 0 atoms".

    Note on the image: the repo's .dockerignore whitelists only pyproject.toml, src/,
    configs/ and data/benchmarks/**, so `scripts/` and `skills/` are NOT present. That
    is fine -- write_condition_spec_cifs.py is a thin wrapper over
    resolve_condition_spec_inputs, which the installed package exposes, and checking the
    emitted annotation columns directly is a sharper test than the generic inspector:
    it asserts the counts we specifically intended.
    """
    out = f"{RUNS}/cifs_dry_run"
    os.makedirs(out, exist_ok=True)
    sh(f"""python - <<'PY'
from pathlib import Path
from rfproteina.inference_utils.condition_spec import resolve_condition_spec_inputs
from rfproteina.utils.verbosity import configure_verbosity
configure_verbosity(True)
for p in resolve_condition_spec_inputs(
        Path("{INPUTS}/tfr1_apical_binder.yaml"), out_dir=Path("{out}"), validate=True):
    print("wrote", p)
PY""")
    sh(f"ls -la {out}", cwd="/", check=False)

    # Verify the conditioning landed where we meant. Expected, from the spec:
    #   C_HOT   -> exactly 10 atoms (2 tip-ward atoms on each of 5 epitope residues)
    #   C_SEQ   -> 195 residues of target sequence
    #   C_CRDGRP-> target backbone atoms only, no NaN-coordinate atoms
    sh(f"""python - <<'PY'
import glob, gzip, io, os
for path in sorted(glob.glob("{out}/*.cif*")):
    print("=" * 70); print(path)
    op = gzip.open if path.endswith(".gz") else open
    txt = op(path, "rt").read()
    lines = txt.splitlines()
    i = 0; cols = []; rows = []
    while i < len(lines):
        if lines[i].strip() == "loop_":
            j = i + 1; hdr = []
            while j < len(lines) and lines[j].lstrip().startswith("_"):
                hdr.append(lines[j].strip()); j += 1
            if hdr and hdr[0].startswith("_atom_site."):
                cols = [h.split(".", 1)[1] for h in hdr]
                while j < len(lines) and lines[j].strip() and not lines[j].lstrip().startswith(("_", "loop_", "#")):
                    rows.append(lines[j].split()); j += 1
                break
            i = j
        else:
            i += 1
    idx = {{c: k for k, c in enumerate(cols)}}
    print("atoms: %d   columns: %d" % (len(rows), len(cols)))
    cond = [c for c in cols if c.startswith("condition") or c.startswith("sa_") or c.startswith("SA_")]
    print("annotation columns present:"); [print("   ", c) for c in cond]
    # Values that mean "not set". 'False' and '-1' matter: a boolean column writes
    # False on every unmarked atom and CoordinateGroup writes -1, so counting those
    # as set reports the whole target as a hotspot. (It did, on the first run.)
    UNSET = {{".", "?", "0", "-1", "nan", "NaN", "False", "false", "UNSPECIFIED", "UNK"}}
    for c in cond:
        k = idx[c]
        vals = {{}}
        for r in rows:
            v = r[k] if k < len(r) else "."
            vals[v] = vals.get(v, 0) + 1
        set_atoms = sum(n for v, n in vals.items() if v not in UNSET)
        top = sorted(vals.items(), key=lambda x: -x[1])[:6]
        print("  %-36s set_on=%-6d values=%s" % (c, set_atoms, top))
        # Assert the two counts we can state exactly from the spec.
        if "hotspot" in c.lower():
            assert set_atoms == 10, "expected 10 hotspot atoms, got %d" % set_atoms
            print("       OK: 10 hotspot atoms (5 residues x 2 tip-ward atoms)")
        if "coordinategroup" in c.lower():
            assert set_atoms == 780, "expected 780 backbone atoms (195x4), got %d" % set_atoms
            print("       OK: 780 backbone atoms in one rigid group (195 residues x N,CA,C,O)")
    # which residues carry a hotspot
    hot = next((c for c in cond if "hotspot" in c.lower()), None)
    if hot:
        k = idx[hot]
        marked = {{}}
        for r in rows:
            v = r[k] if k < len(r) else "."
            if v not in UNSET:
                key = (r[idx["label_asym_id"]], r[idx["label_seq_id"]], r[idx["label_comp_id"]])
                marked.setdefault(key, []).append(r[idx["label_atom_id"]])
        print()
        print("HOTSPOT-MARKED residues (expect 5 residues / 10 atoms):")
        for kk, atoms in sorted(marked.items(), key=lambda x: int(x[0][1]) if x[0][1].isdigit() else 0):
            print("   chain %s  label_seq %-5s %-4s  atoms: %s" % (kk[0], kk[1], kk[2], ",".join(atoms)))
        print("   TOTAL residues=%d atoms=%d" % (len(marked), sum(len(a) for a in marked.values())))
PY""", check=False)
    runs_vol.commit()
    return out


# --------------------------------------------------------------------------------------
# 2. GPU readiness. Also the call that populates the weight cache.
# --------------------------------------------------------------------------------------
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
    # The release is a MATCHED PAIR; mixing a flow model with a different autoencoder
    # "produces nonsense rather than an error", so record exactly what resolved.
    sh("grep -nE 'ckpt|hf_ckpt' configs/inference.yaml | head -20", check=False)
    hf_cache.commit()


# --------------------------------------------------------------------------------------
# 3 & 4. Smoke test, then the shipped PPI example unchanged.
# --------------------------------------------------------------------------------------
@app.function(image=image, gpu=GPU_CHEAP, volumes=VOLUMES, timeout=3600)
def smoke():
    """One design, metrics off, no torch.compile warm-up."""
    sh(f"rfd4-generate inference_experiment=quickstart "
       f"inputs=data/benchmarks/monomer subset_to_examples=monomer-short "
       f"num_replicates=1 out_dir={RUNS}/smoke")
    sh(f"find {RUNS}/smoke -name '*.cif*' | head", cwd="/", check=False)
    runs_vol.commit()


@app.function(image=image, gpu=GPU_CHEAP, volumes=VOLUMES, timeout=7200)
def baseline(num_replicates: int = 4):
    """The shipped ppi_all_atom example, unmodified, to measure cost per design.

    Deliberately NOT our target. If this works and our spec does not, the fault is in
    our spec; if neither works, the fault is the environment. Keeping those separable
    is worth one short run.

    Metrics are left OFF: every *_refolding tag runs the AF3 oracle, which needs a
    separate container we have not set up. `'~metrics'` disables the metrics group.
    """
    inp = "data/benchmarks/ppi/ppi_all_atom/pdl1.cif.zst"
    sh(f"ls -la {inp}", check=False)
    sh(f"rfd4-generate inputs={inp} recycle_strategy=medium/front_loaded "
       f"num_replicates={num_replicates} batch_size=1 seed=42 "
       f"out_dir={RUNS}/baseline_pdl1 prevalidate_inputs=true '~metrics'")
    sh(f"nvidia-smi --query-gpu=memory.used,memory.total --format=csv", cwd="/", check=False)
    sh(f"find {RUNS}/baseline_pdl1 -name '*.cif*' | head -20", cwd="/", check=False)
    runs_vol.commit()


# --------------------------------------------------------------------------------------
# 5. Our target.
# --------------------------------------------------------------------------------------
@app.function(image=image, gpu=GPU_CHEAP, volumes=VOLUMES, timeout=86400)
def generate(num_replicates: int = 8, batch_size: int = 1, seed: int = 42,
             recycle_strategy: str = "medium/front_loaded", out_name: str = "tfr1",
             metrics: bool = False):
    """Design binders against the TfR1 apical epitope.

    batch_size defaults to 1 because upstream says to confirm the longest input fits at
    1 before raising it. save_chain_of_interest is on so designs can be foldseek-
    clustered later -- diversity is a selection criterion we intend to compete on.
    """
    metrics_arg = ("metrics_tags='[ppi_refolding,hotspot_conditioning]'" if metrics
                   else "'~metrics'")
    rc, out, elapsed, peak = sh_timed(
       f"rfd4-generate inputs={INPUTS}/tfr1_apical_binder.yaml "
       f"recycle_strategy={recycle_strategy} "
       f"num_replicates={num_replicates} batch_size={batch_size} seed={seed} "
       f"save_chain_of_interest=true prevalidate_inputs=true "
       f"out_dir={RUNS}/{out_name} {metrics_arg}")
    if num_replicates:
        print("PER DESIGN  %.1fs  -> $%.4f at L40S $1.95/hr"
              % (elapsed / num_replicates, elapsed / num_replicates / 3600 * 1.95), flush=True)
    sh(f"find {RUNS}/{out_name} -name '*.cif*' | head -20", cwd="/", check=False)
    sh(f"find {RUNS}/{out_name} -name 'per_example_metrics.csv' -exec head -3 {{}} +",
       cwd="/", check=False)
    runs_vol.commit()


@app.function(image=image, volumes=VOLUMES, timeout=1800, cpu=2.0)
def fetch(path: str = ""):
    """List or dump what is on the runs volume, so results can be pulled down."""
    sh(f"find {RUNS}/{path} -maxdepth 4 | head -100", cwd="/", check=False)
    sh(f"du -sh {RUNS}/*", cwd="/", check=False)
