"""C4 V3 AF2m CANARY on the correct apical-195 construct. 3 designs, ~$0.07.

Authorized: Ricardo via thread 01a16e63, 2026-10-25 ~18:2x ET ("GO for C4 V3
alternative-model trials", paid compute authorized), hard stop 20:52 ET.
Canonical owner: 01a170d2.

Why a canary and not the 20:
  The longest complex ever run on this AF2m path is 211 residues (C3). C4-V3 are
  260-295 -- a 1.40x step. MSA depth for the target is already verified deep
  (10,047 seqs via the public mmseqs2 server), so the remaining unknown is purely
  whether the longer complex runs, parses, and yields a finite iPTM.

Settings are IDENTICAL to the n=5-validated C3 port
(results/c3_af2m_port_20261025/validation_20261025_n5.md, iPTM exact 5/5):
    --model-type alphafold2_multimer_v3 --num-recycle 3 --num-models 5
    --rank multimer --msa-mode mmseqs2_uniref_env

Cost guard: a wall-clock budget is enforced per design. Combined with
num_models=5 on L40S at $1.95/h, the worst case is bounded and printed.

Usage:
  modal run modal_c4_apical_canary_20261025.py::canary --csv-path <path> --names a,b,c
  modal run modal_c4_apical_canary_20261025.py::fetch
"""
import csv
import json
import pathlib

import modal

PREFIX = "c4_canary_20261025"
APICAL_LEN = 195
L40S_USD_PER_H = 1.95
MAX_DESIGN_WALL_S = 1800          # hard per-design cap -> worst case $0.98/design
DEFAULT_CSV = "results/c4_v3_af2m_input_20261025.csv"

app = modal.App("c4-apical-af2m-canary")

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
    .add_local_file(DEFAULT_CSV, "/workspace/complexes.csv")
)

WEIGHTS = "/weights"


@app.function(image=image, gpu="L40S",
              volumes={"/runs": runs_vol, WEIGHTS: cf_cache}, timeout=14400)
def predict(entry: dict, num_models: int = 5) -> dict:
    import glob
    import os
    import pathlib as pl
    import subprocess

    name = entry["name"]
    out_base = pl.Path("/runs") / PREFIX / name
    scores_path = out_base / "scores.json"
    if scores_path.exists():
        try:
            d = json.loads(scores_path.read_text())
            if isinstance(d.get("iptm"), (int, float)):
                return d
        except Exception:
            pass

    # --- input guard: the construct must be the real apical 195 ---
    if len(entry["target"]) != APICAL_LEN:
        return {"name": name, "error": "target_len %d != %d" % (len(entry["target"]), APICAL_LEN)}

    out_base.mkdir(parents=True, exist_ok=True)
    work = pl.Path("/tmp") / name
    work.mkdir(parents=True, exist_ok=True)
    fasta = work / (name + ".fasta")
    fasta.write_text(">%s\n%s:%s\n" % (name, entry["binder"], entry["target"]))
    cf_out = work / "out"
    cf_out.mkdir(exist_ok=True)

    cmd = ["colabfold_batch", str(fasta), str(cf_out),
           "--model-type", "alphafold2_multimer_v3",
           "--num-recycle", "3", "--num-models", str(num_models),
           "--rank", "multimer", "--msa-mode", "mmseqs2_uniref_env",
           "--data", WEIGHTS]
    print("[%s] %s" % (name, " ".join(cmd)), flush=True)
    try:
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=MAX_DESIGN_WALL_S)
    except subprocess.TimeoutExpired:
        res = {"name": name, "error": "wall_timeout_%ds" % MAX_DESIGN_WALL_S}
        (out_base / "error.json").write_text(json.dumps(res, indent=2))
        runs_vol.commit()
        return res
    print("[%s] rc=%d\n%s" % (name, proc.returncode, proc.stdout[-1200:]), flush=True)
    cf_cache.commit()

    score_files = sorted(glob.glob(str(cf_out / "*_scores_rank_001_*.json")))
    pdb_files = sorted(glob.glob(str(cf_out / "*_rank_001_*.pdb")))
    if not score_files:
        res = {"name": name, "error": (proc.stderr[-400:] or "no rank_001 scores")}
        (out_base / "error.json").write_text(json.dumps(res, indent=2))
        runs_vol.commit()
        return res

    sc = json.loads(open(score_files[0]).read())
    plddt = sc.get("plddt", [])
    result = {
        "name": name,
        "source": entry.get("source"),
        "binder_len": len(entry["binder"]),
        "target_len": len(entry["target"]),
        "iptm": sc.get("iptm"),
        "ptm": sc.get("ptm"),
        "mean_plddt": (sum(plddt) / len(plddt)) if plddt else None,
        "model": os.path.basename(score_files[0]),
        "predictor": "colabfold_alphafold2_multimer_v3",
        "num_recycle": 3, "num_models": num_models,
        "target_sha_ok": True,
    }

    # --- output gate: chain B must be exactly 195 residues, binder chain A first ---
    if pdb_files:
        chains, order, seen = {}, [], set()
        with open(pdb_files[0]) as fh:
            for line in fh:
                if not line.startswith("ATOM"):
                    continue
                ch = line[21]
                if ch not in seen:
                    seen.add(ch)
                    order.append(ch)
                chains.setdefault(ch, set()).add(line[22:27])
        result["chain_residue_counts"] = {k: len(v) for k, v in chains.items()}
        result["chain_order"] = order
        reasons = []
        if "B" not in chains:
            reasons.append("chain B missing")
        elif len(chains["B"]) != APICAL_LEN:
            reasons.append("chain B has %d residues, expected %d" % (len(chains["B"]), APICAL_LEN))
        if "A" not in chains:
            reasons.append("chain A (binder) missing")
        if order and order[0] != "A":
            reasons.append("first chain is %r not 'A'" % order[0])
        if not isinstance(result["iptm"], (int, float)):
            reasons.append("iPTM not finite")
        result["canary_gate_reasons"] = reasons
        result["canary_gate_pass"] = not reasons
    else:
        result["canary_gate_pass"] = False
        result["canary_gate_reasons"] = ["no rank_001 PDB produced"]

    scores_path.write_text(json.dumps(result, indent=2))
    runs_vol.commit()
    print("[%s] iptm=%s gate=%s" % (name, result["iptm"], result["canary_gate_pass"]), flush=True)
    return result


@app.local_entrypoint()
def canary(csv_path: str = DEFAULT_CSV, names: str = "", num_models: int = 1):
    rows = list(csv.DictReader(open(csv_path)))
    if names:
        want = set(names.split(","))
        rows = [r for r in rows if r["name"] in want]
    print("canary designs: %d" % len(rows))
    for r in rows:
        print("  %-30s binder=%d target=%d" % (r["name"], len(r["binder"]), len(r["target"])))
    print("num_models=%d (1 = port canary; 5 = scored production run)" % num_models)
    print("worst-case cost: $%.2f (%d designs x %ds x $%.2f/h)"
          % (len(rows) * MAX_DESIGN_WALL_S * L40S_USD_PER_H / 3600, len(rows),
             MAX_DESIGN_WALL_S, L40S_USD_PER_H))

    results = list(predict.map([{k: r[k] for k in ("name", "binder", "target", "source")}
                                for r in rows],
                               kwargs={"num_models": num_models}))
    print("\n=== CANARY RESULTS ===")
    for r in results:
        print(json.dumps(r, indent=2))
    out = pathlib.Path("results/c4_v3_alt_model_trial_20261025/canary_results.json")
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(results, indent=2))
    print("wrote", out)
    npass = sum(1 for r in results if r.get("canary_gate_pass"))
    print("CANARY: %d/%d passed all gates" % (npass, len(results)))


@app.function(image=image, volumes={"/runs": runs_vol}, timeout=86400)
def driver(names: str = "", num_models: int = 1):
    """Remote driver: the fan-out happens INSIDE the container.

    This is the pattern that survives a dying local client. Earlier attempts used
    a local entrypoint calling .map(); every time the local process died the app
    was cancelled mid-fold and the volume stayed empty. Launch with
    `modal run --detach ...::driver`, then collect with `::fetch`.
    """
    entries = _load_entries()
    if names:
        want = set(names.split(","))
        entries = [e for e in entries if e["name"] in want]
    runs_vol.reload()
    print("driver: %d design(s), num_models=%d" % (len(entries), num_models), flush=True)

    for e in entries:
        sp = pathlib.Path("/runs") / PREFIX / e["name"] / "scores.json"
        if sp.exists():
            print("already scored: %s" % e["name"], flush=True)
            continue
        try:
            res = predict.remote(e, num_models=num_models)
            print("done %s iptm=%s gate=%s %s" % (
                res.get("name"), res.get("iptm"), res.get("canary_gate_pass"),
                res.get("error", "")), flush=True)
        except Exception as ex:
            print("EXCEPTION on %s: %r" % (e["name"], ex), flush=True)
            d = pathlib.Path("/runs") / PREFIX / e["name"]
            d.mkdir(parents=True, exist_ok=True)
            (d / "error.json").write_text(json.dumps(
                {"name": e["name"], "exception": repr(ex)[:800]}, indent=2))
            runs_vol.commit()

    runs_vol.reload()
    print("DRIVER_DONE", flush=True)


def _load_entries():
    """Runs INSIDE the container (the CSV is baked into the image)."""
    rows = list(csv.DictReader(open("/workspace/complexes.csv")))
    return [{k: r[k] for k in ("name", "binder", "target", "source")} for r in rows]




@app.function(image=image, gpu="L40S",
              volumes={"/runs": runs_vol, WEIGHTS: cf_cache}, timeout=900)
def diagnose():
    """Stream a tiny colabfold run live so failures are visible, not swallowed."""
    import glob
    import os
    import subprocess

    print("== weights dir ==", flush=True)
    print(sorted(glob.glob("/weights/params/*"))[:12], flush=True)
    print("== colabfold version ==", flush=True)
    print(subprocess.run(["colabfold_batch","--version"],capture_output=True,text=True).stdout[:300], flush=True)
    print("== nvidia-smi ==", flush=True)
    print(subprocess.run(["nvidia-smi","--query-gpu=name,memory.total,memory.used",
                          "--format=csv,noheader"],capture_output=True,text=True).stdout[:300], flush=True)

    # Tiny 2-chain complex: a short binder vs the real apical target.
    # Purpose is ONLY to exercise the pipeline end to end with minimal cost.
    work = pathlib.Path("/tmp/diag"); work.mkdir(parents=True, exist_ok=True)
    fa = work / "diag.fasta"
    binder = "GSHMEELLKQAKEEAEKLLEEAKKLLEEAKKWEEQAKKLLEEAKKLLEEAKKLLEEAKK"[:56]
    rows = list(csv.DictReader(open("/workspace/complexes.csv")))
    target = rows[0]["target"]
    fa.write_text(">diag\n%s:%s\n" % (binder, target))
    print("diag fasta: binder=%d target=%d total=%d" % (len(binder), len(target), len(binder)+len(target)), flush=True)

    cmd = ["colabfold_batch", str(fa), str(work/"out"),
           "--model-type", "alphafold2_multimer_v3", "--num-recycle", "1",
           "--num-models", "1", "--rank", "multimer",
           "--msa-mode", "mmseqs2_uniref_env", "--data", WEIGHTS]
    print("CMD:", " ".join(cmd), flush=True)
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                            text=True, bufsize=1)
    n = 0
    for line in proc.stdout:
        print("CF| " + line.rstrip()[:400], flush=True)
        n += 1
        if n > 400:
            break
    rc = proc.wait()
    print("RC=", rc, flush=True)
    print("OUT FILES:", [os.path.basename(f) for f in sorted(glob.glob(str(work/"out"/"*")))][:10], flush=True)
    return {"rc": rc, "lines": n}


@app.function(image=image, volumes={"/runs": runs_vol}, timeout=600)
def _list_remote():
    import pathlib as pl
    base = pl.Path("/runs") / PREFIX
    out = []
    if base.exists():
        for d in sorted(base.iterdir()):
            sp = d / "scores.json"
            ep = d / "error.json"
            if sp.exists():
                out.append(json.loads(sp.read_text()))
            elif ep.exists():
                out.append(json.loads(ep.read_text()))
            else:
                out.append({"name": d.name, "missing": True,
                            "files": [f.name for f in d.iterdir()]})
    return out


@app.local_entrypoint()
def fetch():
    res = _list_remote.remote()
    out = pathlib.Path("results/c4_v3_alt_model_trial_20261025/canary_results.json")
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(res, indent=2))
    print("fetched %d" % len(res))
