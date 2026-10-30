"""AF2-multimer complex screen for C3 (GLP1R ECD groove binders).

Authorized bounded screen (results/c3_improvement_20261024.md §6): predict the
binder+ECD complex for the current 20 (v6) AND a 24-design RFD4 shortlist against
the SAME 111-residue ECD target used by the ranked BindCraft CIFs, with ONE
predictor (ColabFold AlphaFold2-multimer v3) so iPTM is uniform and in the same
predictor family as BindCraft's own AF2 evidence. NOT Boltz (not cross-comparable).

Design choices for "smallest reliable, same-predictor":
  * ColabFold alphafold2_multimer_v3, num_recycle=3, 5 models, rank by multimer.
  * MSA via public mmseqs2 server: a de-novo binder has no homologs (-> effectively
    single-sequence), the ECD target gets a real MSA — this reproduces BindCraft's
    "target MSA + single-sequence binder" regime without custom a3m stitching.
  * Isolated persistent volume `c3-af2m-runs`, prefix `c3_ecd_complex_20261024/`.
    Does NOT touch rfd4-runs (live C4 job) or any submission CSV.
  * Per-design idempotent: a design with a valid scores.json is skipped on resume.
  * Remote driver (main) + remote consolidation so a Mini restart cannot interrupt.

Usage:
  # one-complex smoke test (debug image/weights cheaply, prints iPTM):
  modal run modal_af2m_c3_complex.py::smoke --name c3_bc_l57_fbcf76
  # full detached screen (survives client/Mini exit):
  modal run --detach modal_af2m_c3_complex.py::main
  # read-only consolidation to local JSON:
  modal run modal_af2m_c3_complex.py::fetch
"""
import csv
import json
import pathlib

import modal

CSV_PATH = pathlib.Path("results/c3_screen_input_20261024.csv")
PREFIX = "c3_ecd_complex_20261024"

app = modal.App("af2m-c3-ecd-complex")

runs_vol = modal.Volume.from_name("c3-af2m-runs", create_if_missing=True)
cf_cache = modal.Volume.from_name("colabfold-cache", create_if_missing=True)

image = (
    modal.Image.debian_slim(python_version="3.10")
    .apt_install("git", "wget", "aria2", "libgomp1", "tar", "gzip", "curl")
    .pip_install("colabfold[alphafold] @ git+https://github.com/sokrypton/ColabFold@v1.5.5")
    .run_commands(
        "pip install --upgrade 'jax[cuda12_pip]==0.4.23' "
        "-f https://storage.googleapis.com/jax-releases/jax_cuda_releases.html"
    )
    .env({
        "TF_FORCE_UNIFIED_MEMORY": "1",
        "XLA_PYTHON_CLIENT_MEM_FRACTION": "4.0",
        "XLA_PYTHON_CLIENT_PREALLOCATE": "false",
        "COLABFOLD_REBUILD": "v1",
    })
    .add_local_file(str(CSV_PATH), "/workspace/complexes.csv")
)

WEIGHTS = "/weights"


def _load_complexes():
    seqs = []
    with open("/workspace/complexes.csv") as f:
        for r in csv.DictReader(f):
            seqs.append(dict(name=r["name"], binder=r["binder"],
                             target=r["target"], source=r["source"]))
    return seqs


def _valid_scores(path):
    try:
        with open(path) as f:
            d = json.load(f)
        return isinstance(d.get("iptm"), (int, float)) and "ptm" in d
    except Exception:
        return False


@app.function(
    image=image,
    gpu="L40S",
    volumes={"/runs": runs_vol, WEIGHTS: cf_cache},
    timeout=10800,
)
def predict(entry: dict) -> dict:
    """Predict one binder+target complex with ColabFold AF2-multimer. Idempotent.

    Reuse this same remote entry point for one-off retries with the longer wall.
    """
    import glob
    import os
    import shutil
    import subprocess

    name = entry["name"]
    out_base = pathlib.Path(f"/runs/{PREFIX}/{name}")
    scores_path = out_base / "scores.json"

    if _valid_scores(scores_path):
        with open(scores_path) as f:
            return json.load(f)

    out_base.mkdir(parents=True, exist_ok=True)
    # Multimer FASTA: binder chain : target chain
    work = pathlib.Path(f"/tmp/{name}")
    work.mkdir(parents=True, exist_ok=True)
    fasta = work / f"{name}.fasta"
    fasta.write_text(f">{name}\n{entry['binder']}:{entry['target']}\n")
    cf_out = work / "out"
    cf_out.mkdir(exist_ok=True)

    cmd = [
        "colabfold_batch",
        str(fasta), str(cf_out),
        "--model-type", "alphafold2_multimer_v3",
        "--num-recycle", "3",
        "--num-models", "5",
        "--rank", "multimer",
        "--msa-mode", "mmseqs2_uniref_env",
        "--data", WEIGHTS,
    ]
    print(f"[{name}] {' '.join(cmd)}", flush=True)
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=9900)
    print(f"[{name}] rc={r.returncode}\n{r.stdout[-1500:]}", flush=True)
    if r.returncode != 0:
        print(f"[{name}] STDERR {r.stderr[-1500:]}", flush=True)
        cf_cache.commit()
        return {"name": name, "error": r.stderr[-400:] or "colabfold failed"}

    # commit downloaded weights so later designs reuse them
    cf_cache.commit()

    score_files = sorted(glob.glob(str(cf_out / "*_scores_rank_001_*.json")))
    pdb_files = sorted(glob.glob(str(cf_out / "*_rank_001_*.pdb")))
    if not score_files:
        print(f"[{name}] no rank_001 scores; dir:", flush=True)
        for p in sorted(cf_out.rglob("*")):
            print("   ", p.name, flush=True)
        return {"name": name, "error": "no rank_001 scores"}

    with open(score_files[0]) as f:
        sc = json.load(f)
    plddt = sc.get("plddt", [])
    mean_plddt = sum(plddt) / len(plddt) if plddt else None
    result = {
        "name": name,
        "source": entry["source"],
        "binder_len": len(entry["binder"]),
        "target_len": len(entry["target"]),
        "iptm": sc.get("iptm"),
        "ptm": sc.get("ptm"),
        "mean_plddt": mean_plddt,
        "model": os.path.basename(score_files[0]),
        "predictor": "colabfold_alphafold2_multimer_v3",
        "num_recycle": 3,
        "num_models": 5,
    }
    # persist raw outputs + distilled scores
    shutil.copy2(score_files[0], out_base / "scores_rank_001.json")
    if pdb_files:
        shutil.copy2(pdb_files[0], out_base / "complex_rank_001.pdb")
    with open(scores_path, "w") as f:
        json.dump(result, f, indent=2)
    runs_vol.commit()
    print(f"[{name}] iptm={result['iptm']} ptm={result['ptm']} plddt={mean_plddt}", flush=True)
    return result


@app.local_entrypoint()
def smoke(name: str = "c3_bc_l57_fbcf76"):
    """Run a single complex to validate the image/weights/parse path cheaply.
    local_entrypoint → one GPU container (predict, 10800s), no nested-timeout."""
    import csv as _csv
    entries = {}
    with open(CSV_PATH) as f:
        for r in _csv.DictReader(f):
            entries[r["name"]] = dict(name=r["name"], binder=r["binder"],
                                      target=r["target"], source=r["source"])
    if name not in entries:
        print(f"name {name} not in input; options: {list(entries)[:5]}...")
        return
    res = predict.remote(entries[name])
    print("SMOKE RESULT:", json.dumps(res, indent=2))
    return res


# Allow multiple waves of long predictions plus consolidation within the driver wall.
@app.function(image=image, volumes={"/runs": runs_vol}, timeout=86400)
def main():
    """Remote driver: idempotent fan-out over all 44 complexes, then consolidate.
    Launch with `modal run --detach ...::main` so it survives client/Mini exit."""
    entries = _load_complexes()
    runs_vol.reload()
    pending = []
    done = 0
    for e in entries:
        sp = pathlib.Path(f"/runs/{PREFIX}/{e['name']}/scores.json")
        if _valid_scores(sp):
            done += 1
        else:
            pending.append(e)
    print(f"total={len(entries)} done={done} pending={len(pending)}", flush=True)

    results = []
    for res in predict.map(pending, order_outputs=False, return_exceptions=True):
        if isinstance(res, Exception):
            print("task exception:", repr(res), flush=True)
        else:
            results.append(res)
            print("done:", res.get("name"), "iptm=", res.get("iptm"), res.get("error", ""), flush=True)

    _consolidate()
    return {"pending_run": len(pending), "previously_done": done}


def _consolidate():
    import glob
    runs_vol.reload()
    rows = []
    for sp in sorted(glob.glob(f"/runs/{PREFIX}/*/scores.json")):
        try:
            with open(sp) as f:
                rows.append(json.load(f))
        except Exception as ex:
            print("skip", sp, ex, flush=True)
    out = pathlib.Path(f"/runs/{PREFIX}/consolidated_scores.json")
    with open(out, "w") as f:
        json.dump(rows, f, indent=2)
    runs_vol.commit()
    scored = [r for r in rows if isinstance(r.get("iptm"), (int, float))]
    print(f"consolidated {len(rows)} ({len(scored)} scored) -> {out}", flush=True)
    return rows


@app.function(image=image, volumes={"/runs": runs_vol}, timeout=600)
def consolidate():
    return _consolidate()


@app.local_entrypoint()
def fetch():
    """Read-only: pull consolidated scores to local results/."""
    rows = consolidate.remote()
    local = pathlib.Path("results/c3_screen_af2m_scores_20261024.json")
    with open(local, "w") as f:
        json.dump(rows, f, indent=2)
    scored = [r for r in rows if isinstance(r.get("iptm"), (int, float))]
    print(f"Wrote {len(rows)} rows ({len(scored)} scored) -> {local}")
    for r in sorted(scored, key=lambda d: -(d.get("iptm") or 0)):
        print(f"  {r['name']:<24s} {r.get('source',''):<14s} iptm={r['iptm']:.3f} "
              f"ptm={r.get('ptm',0):.3f} plddt={r.get('mean_plddt') or 0:.1f}")
