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
    .add_local_file(str(SPEC.parent / "tfr1_apical_his.yaml"),
                    f"{INPUTS}/tfr1_apical_his.yaml")
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
def dry_run(spec: str = "tfr1_apical_binder"):
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
    out = f"{RUNS}/cifs_dry_run_{spec}"
    os.makedirs(out, exist_ok=True)
    sh(f"""python - <<'PY'
from pathlib import Path
from rfproteina.inference_utils.condition_spec import resolve_condition_spec_inputs
from rfproteina.utils.verbosity import configure_verbosity
configure_verbosity(True)
for p in resolve_condition_spec_inputs(
        Path("{INPUTS}/{spec}.yaml"), out_dir=Path("{out}"), validate=True):
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
             metrics: str = "geom", spec: str = "tfr1_apical_binder"):
    """Design binders against the TfR1 apical epitope.

    batch_size defaults to 1 because upstream says to confirm the longest input fits at
    1 before raising it. save_chain_of_interest is on so designs can be foldseek-
    clustered later -- diversity is a selection criterion we intend to compete on.
    """
    # metrics must run INLINE: hotspot_obedience reads the C_HOT conditioning, which is
    # in scope during generation but absent from the generated CIF, so offline
    # rfd4-evaluate cannot compute it. design+hotspot_conditioning are geometric and
    # need no oracle container; only *_refolding tags do.
    metrics_arg = {"none": "'~metrics'",
                   "geom": "metrics_tags='[design,hotspot_conditioning]'",
                   "full": "metrics_tags='[ppi_refolding,hotspot_conditioning]'"}[metrics]
    rc, out, elapsed, peak = sh_timed(
       f"rfd4-generate inputs={INPUTS}/{spec}.yaml "
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


@app.function(image=image, volumes=VOLUMES, timeout=1800, cpu=2.0)
def sequences(run: str = "tfr1_probe"):
    """Extract the designed chain's sequence from each generated CIF on the volume.

    The generated CIF holds both chains. The designed chain is the one whose residues
    were NOT sequence-conditioned, so it is identified by the absence of the target's
    conditioning rather than by chain letter.
    """
    sh(f"""python - <<'PY2'
import glob, gzip
THREE = {{'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E','GLY':'G',
          'HIS':'H','ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S',
          'THR':'T','TRP':'W','TYR':'Y','VAL':'V'}}
for path in sorted(glob.glob("{RUNS}/{run}/inference_samples/*/*_generated.cif*")):
    op = gzip.open if path.endswith('.gz') else open
    lines = op(path, 'rt').read().splitlines()
    i=0; cols=[]; rows=[]
    while i < len(lines):
        if lines[i].strip()=='loop_':
            j=i+1; hdr=[]
            while j<len(lines) and lines[j].lstrip().startswith('_'):
                hdr.append(lines[j].strip()); j+=1
            if hdr and hdr[0].startswith('_atom_site.'):
                cols=[h.split('.',1)[1] for h in hdr]
                while j<len(lines) and lines[j].strip() and not lines[j].lstrip().startswith(('_','loop_','#')):
                    rows.append(lines[j].split()); j+=1
                break
            i=j
        else: i+=1
    idx={{c:k for k,c in enumerate(cols)}}
    def g(r,c): return r[idx[c]] if c in idx else '?'
    chains={{}}
    for r in rows:
        ch=g(r,'label_asym_id'); rid=g(r,'label_seq_id'); comp=g(r,'label_comp_id')
        chains.setdefault(ch,{{}})[rid]=comp
    print('='*70); print(path.split('/')[-1])
    for ch, res in sorted(chains.items()):
        order=sorted(res, key=lambda x: int(x) if x.lstrip('-').isdigit() else 0)
        seq=''.join(THREE.get(res[k],'X') for k in order)
        print('  chain %s  len %d' % (ch, len(seq)))
        print('    %s' % seq)
PY2""", check=False)


@app.function(image=image, gpu=GPU_CHEAP, volumes=VOLUMES, timeout=14400)
def evaluate(run: str = "tfr1_c200",
             tags: str = "design,hotspot_conditioning,developability"):
    """Score an existing run offline, with NO refolding oracle.

    Every *_refolding tag runs the AF3/OpenFold3 oracle, which needs an apptainer
    container we have not built (and apptainer inside Modal's sandbox is doubtful). But
    the metrics that matter most for ranking are pure geometry or physics and need no
    container at all:

      hotspot_obedience   -- is the designed chain within 4.5 A of the conditioned
                             epitope. This is the exact failure the retired pipeline hid:
                             its binder sat 54-57 A away for 250 gradient steps.
                             rfproteina.metrics.hotspots imports only numpy/atomworks/
                             biotite -- no oracle.
      clash / bond_length -- geometric sanity on the generated structure
      developability      -- sequence and Rosetta developability, which maps onto the
                             competition's FIRST criterion, "a real chance to express
                             (cell-free) and bind"

    Note the known upstream bug: ppi.pass reads 0.0 for every C_CRDGRP-conditioned
    design because refolded_target_backbone_aligned_binder_rmsd is never produced
    (issue #59, PR #60 unmerged). We are not relying on ppi.pass.
    """
    out = f"{RUNS}/{run}_eval"
    staged = f"{RUNS}/{run}_eval_inputs"
    # Stage only the two-chain generated structures. A recursive scan of the run dir
    # also finds the *_chain_* single-chain exports written by save_chain_of_interest,
    # and every interface metric fails on those because the target is absent -- which is
    # what produced 200 empty "Unexpected error in context:" warnings.
    sh(f"rm -rf {staged} && mkdir -p {staged}", cwd="/")
    sh(f"for f in {RUNS}/{run}/inference_samples/*/*_generated.cif*; do "
       f"ln -s \"$f\" {staged}/$(basename $f); done", cwd="/", check=False)
    sh(f"ls {staged} | wc -l", cwd="/", check=False)
    rc, o, elapsed, peak = sh_timed(
        f"rfd4-evaluate inputs={staged} metrics_tags='[{tags}]' out_dir={out}",
        check=False)
    sh(f"find {out} -name '*.csv' | head", cwd="/", check=False)
    sh(f"""python - <<'PY2'
import csv, glob
for p in sorted(glob.glob("{out}/**/per_example_metrics.csv", recursive=True)):
    rows = list(csv.DictReader(open(p)))
    if not rows: continue
    print("=" * 70); print(p); print("rows:", len(rows))
    keys = [k for k in rows[0] if not k.startswith("generation.")]
    interesting = [k for k in keys if any(t in k.lower() for t in
                   ("hotspot", "clash", "ddg", "developab", "bond", "obedien", "rmsd"))]
    print("metric columns of interest:")
    for k in interesting: print("   ", k)
    hot = next((k for k in interesting if "obedien" in k.lower() or "hotspot" in k.lower()), None)
    if hot:
        def num(v):
            try: return float(v)
            except Exception: return None
        vals = [(num(r.get(hot)), r.get("name") or r.get("example") or "?") for r in rows]
        vals = [(v, n) for v, n in vals if v is not None]
        print()
        print("%s: %d of %d rows numeric" % (hot, len(vals), len(rows)))
        if vals:
            vals.sort(reverse=True)
            import statistics
            print("  max %.3f  median %.3f  min %.3f"
                  % (vals[0][0], statistics.median([v for v, _ in vals]), vals[-1][0]))
            print("  top 10:")
            for v, n in vals[:10]: print("    %.3f  %s" % (v, n))
PY2""", cwd="/", check=False)
    runs_vol.commit()


@app.function(image=image, gpu=GPU_CHEAP, volumes=VOLUMES, timeout=3600)
def debug_eval(run: str = "tfr1_c200"):
    """Score ONE design with full output captured to the volume.

    Modal truncates stdout ("Logs may not be continuous"), and rfd4-evaluate reports
    "Unexpected error in context:" with an empty message, so the real exception has to
    be captured to a file and read back.
    """
    import glob, os
    gen = sorted(glob.glob(f"{RUNS}/{run}/inference_samples/*/*_generated.cif*"))
    print("generated files: %d" % len(gen))
    chains = sorted(glob.glob(f"{RUNS}/{run}/inference_samples/*/*_chain_*.cif*"))
    print("chain-of-interest files: %d" % len(chains))
    if not gen:
        print("NO GENERATED FILES"); return
    one = os.path.dirname(gen[0])
    print("evaluating only:", one)
    log = f"{RUNS}/{run}_debug.log"
    sh(f"rfd4-evaluate inputs={one} metrics_tags='[design,hotspot_conditioning]' "
       f"out_dir={RUNS}/{run}_debug_eval --config-dir configs > {log} 2>&1 || true",
       check=False)
    runs_vol.commit()
    sh(f"tail -60 {log}", cwd="/", check=False)
    sh(f"grep -iE -A12 'traceback|error' {log} | head -40", cwd="/", check=False)


@app.function(image=image, volumes=VOLUMES, timeout=3600, cpu=2.0)
def export_sequences(run: str = "tfr1_c200", out: str = ""):
    """Write every designed chain's sequence to a CSV on the volume, for offline screening.

    The designed chain is the shorter one: the target is the preserved 195-residue apical
    domain, identical in every output, so the other chain is what was generated.
    Emits name,sequence so scripts/novelty_level_predictor.py and the diversity audit can
    consume it directly, and counts histidines so the switch arm's actual adherence can be
    checked against what the spec requested (upstream: these are model conditions, not
    hard output constraints).
    """
    dest = out or f"{RUNS}/{run}_sequences.csv"
    sh(f"""python - <<'PY2'
import csv, glob, gzip, os
THREE = {{'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E','GLY':'G',
          'HIS':'H','ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S',
          'THR':'T','TRP':'W','TYR':'Y','VAL':'V'}}
rows=[]
for path in sorted(glob.glob("{RUNS}/{run}/inference_samples/*/*_generated.cif*")):
    op = gzip.open if path.endswith('.gz') else open
    lines = op(path,'rt').read().splitlines()
    i=0; cols=[]; ats=[]
    while i < len(lines):
        if lines[i].strip()=='loop_':
            j=i+1; hdr=[]
            while j<len(lines) and lines[j].lstrip().startswith('_'):
                hdr.append(lines[j].strip()); j+=1
            if hdr and hdr[0].startswith('_atom_site.'):
                cols=[h.split('.',1)[1] for h in hdr]
                while j<len(lines) and lines[j].strip() and not lines[j].lstrip().startswith(('_','loop_','#')):
                    ats.append(lines[j].split()); j+=1
                break
            i=j
        else: i+=1
    idx={{c:k for k,c in enumerate(cols)}}
    def g(r,c): return r[idx[c]] if c in idx else '?'
    ch={{}}
    for r in ats:
        ch.setdefault(g(r,'label_asym_id'),{{}})[g(r,'label_seq_id')]=g(r,'label_comp_id')
    seqs={{}}
    for c,res in ch.items():
        order=sorted(res,key=lambda x:int(x) if x.lstrip('-').isdigit() else 0)
        seqs[c]=''.join(THREE.get(res[k],'X') for k in order)
    if not seqs: continue
    # designed chain = the one that is not the 195-residue preserved target
    cand = sorted(seqs.items(), key=lambda kv: len(kv[1]))
    designed = next((v for k,v in cand if len(v)!=195), cand[0][1])
    name = os.path.basename(path).split('_generated')[0]
    rows.append({{'name':name,'sequence':designed,'length':len(designed),
                  'n_his':designed.count('H'),
                  'pct_his':round(100.0*designed.count('H')/max(1,len(designed)),2),
                  'n_cys':designed.count('C')}})
with open("{dest}",'w',newline='') as fh:
    w=csv.DictWriter(fh,fieldnames=['name','sequence','length','n_his','pct_his','n_cys'],lineterminator=chr(10))
    w.writeheader(); w.writerows(rows)
print('wrote {dest}  rows=%d' % len(rows))
if rows:
    import statistics
    L=[r['length'] for r in rows]; H=[r['pct_his'] for r in rows]; C=[r['n_cys'] for r in rows]
    print('length   min %d  median %d  max %d' % (min(L),statistics.median(L),max(L)))
    print('%%His     min %.2f median %.2f max %.2f  (natural ~2.3%%)' % (min(H),statistics.median(H),max(H)))
    print('designs with >=3 His: %d/%d' % (sum(1 for r in rows if r['n_his']>=3), len(rows)))
    print('free Cys present in: %d/%d designs' % (sum(1 for c in C if c>0), len(rows)))
PY2""", check=False)
    runs_vol.commit()
    return dest


@app.function(image=image, gpu=GPU_CHEAP, volumes=VOLUMES, timeout=3600)
def trace_eval(run: str = "tfr1_c200"):
    """Get the ACTUAL exception. rfd4-evaluate swallows it into an empty
    'Unexpected error in context:' and Modal truncates stdout, so call the metric
    directly in-process on one structure and let Python raise."""
    import glob
    gen = sorted(glob.glob(f"{RUNS}/{run}/inference_samples/*/*_generated.cif*"))
    target = gen[0]
    print("target:", target)
    sh(f"""python - <<'PY2' 2>&1 | tail -60
import traceback
path = "{target}"
from atomworks.io.parser import parse
try:
    out = parse(path)
    aa = out["atom_array"] if isinstance(out, dict) else out
    print("parsed OK:", type(aa))
    import numpy as np
    for attr in ("chain_id", "res_id", "res_name", "atom_name"):
        if hasattr(aa, attr):
            v = getattr(aa, attr)
            print("  %-10s n=%d uniq=%s" % (attr, len(v), sorted(set(map(str, v)))[:6]))
    anns = [a for a in getattr(aa, "_annot", {{}}).keys()] if hasattr(aa, "_annot") else []
    print("annotations:", [a for a in anns if "ondition" in a or a.startswith("C_") or a.startswith("SA_")][:12])
except Exception:
    traceback.print_exc()

# now the metric itself
try:
    from rfproteina.metrics.hotspots import HotspotObedience, hotspot_obedience
    import inspect
    print()
    print("hotspot_obedience signature:", inspect.signature(hotspot_obedience))
    print("HotspotObedience.kwargs_to_compute_args:")
    print(inspect.getsource(HotspotObedience.kwargs_to_compute_args))
except Exception:
    traceback.print_exc()
PY2""", check=False)
