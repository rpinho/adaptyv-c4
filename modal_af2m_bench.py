"""GPU cost-per-design benchmark for the C3 AF2M screen (Ricardo approval #2).

Reuses the EXACT image and AF2M settings from modal_af2m_c3_complex.py
(alphafold2_multimer_v3, num_recycle=3, num_models=5, rank=multimer, same 111-res
ECD target). Only the GPU type changes. Writes to a separate `bench_` prefix on the
same volume so it can never satisfy a real design's idempotency check.

Measurement: full function wall (container start -> scores written) and the
colabfold subprocess wall. cost_per_design = wall_h * modal_hourly_rate.

  modal run --detach modal_af2m_bench.py::main
"""
import json
import pathlib
import time

import modal

import modal_af2m_c3_complex as _base

# The remote container must also contain the frozen AF2M module so the image and
# the scoring/validation helpers are byte-identical to the real screen.
image = _base.image.add_local_file(
    str(pathlib.Path(__file__).parent / "modal_af2m_c3_complex.py"),
    "/root/modal_af2m_c3_complex.py",
)
WEIGHTS = _base.WEIGHTS
cf_cache = _base.cf_cache
runs_vol = _base.runs_vol
_valid_scores = _base._valid_scores

app = modal.App("af2m-c3-gpu-bench")

PREFIX = "bench_20261024"
RATES = {"T4": 0.59, "L4": 0.80, "A10G": 1.10, "L40S": 1.95}


def _bench_core(entry: dict, gpu_label: str) -> dict:
    import glob
    import os
    import shutil
    import subprocess

    name = entry["name"]
    out_base = pathlib.Path(f"/runs/{PREFIX}/{gpu_label}/{name}")
    scores_path = out_base / "scores.json"
    if _valid_scores(scores_path):
        with open(scores_path) as f:
            d = json.load(f)
        d["cached"] = True
        return d

    t0 = time.time()
    out_base.mkdir(parents=True, exist_ok=True)
    work = pathlib.Path(f"/tmp/{name}")
    work.mkdir(parents=True, exist_ok=True)
    fasta = work / f"{name}.fasta"
    fasta.write_text(f">{name}\n{entry['binder']}:{entry['target']}\n")
    cf_out = work / "out"
    cf_out.mkdir(exist_ok=True)

    cmd = [
        "colabfold_batch", str(fasta), str(cf_out),
        "--model-type", "alphafold2_multimer_v3",
        "--num-recycle", "3",
        "--num-models", "5",
        "--rank", "multimer",
        "--msa-mode", "mmseqs2_uniref_env",
        "--data", WEIGHTS,
    ]
    t_start = time.time()
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=5000)
    sub_wall = time.time() - t_start
    if r.returncode != 0:
        print(f"[{gpu_label}] STDERR {r.stderr[-1200:]}", flush=True)
        return {"name": name, "gpu": gpu_label, "error": r.stderr[-300:] or "failed",
                "subprocess_wall_s": sub_wall}

    score_files = sorted(glob.glob(str(cf_out / "*_scores_rank_001_*.json")))
    pdb_files = sorted(glob.glob(str(cf_out / "*_rank_001_*.pdb")))
    if not score_files:
        return {"name": name, "gpu": gpu_label, "error": "no rank_001 scores",
                "subprocess_wall_s": sub_wall}

    with open(score_files[0]) as f:
        sc = json.load(f)
    plddt = sc.get("plddt", [])
    wall = time.time() - t0
    result = {
        "name": name,
        "gpu": gpu_label,
        "binder_len": len(entry["binder"]),
        "target_len": len(entry["target"]),
        "iptm": sc.get("iptm"),
        "ptm": sc.get("ptm"),
        "mean_plddt": sum(plddt) / len(plddt) if plddt else None,
        "subprocess_wall_s": round(sub_wall, 1),
        "function_wall_s": round(wall, 1),
        "cached": False,
    }
    shutil.copy2(score_files[0], out_base / "scores_rank_001.json")
    if pdb_files:
        shutil.copy2(pdb_files[0], out_base / "complex_rank_001.pdb")
    with open(scores_path, "w") as f:
        json.dump(result, f, indent=2)
    runs_vol.commit()
    print(f"[{gpu_label}] {name} iptm={result['iptm']} sub={sub_wall:.0f}s wall={wall:.0f}s",
          flush=True)
    return result


@app.function(
    image=image,
    gpu="T4",
    volumes={"/runs": runs_vol, WEIGHTS: cf_cache},
    timeout=5400,
    retries=0,
)
def bench_t4(entry: dict) -> dict:
    return _bench_core(entry, "T4")


@app.function(
    image=image,
    gpu="L4",
    volumes={"/runs": runs_vol, WEIGHTS: cf_cache},
    timeout=5400,
    retries=0,
)
def bench_l4(entry: dict) -> dict:
    return _bench_core(entry, "L4")


@app.function(
    image=image,
    gpu="A10G",
    volumes={"/runs": runs_vol, WEIGHTS: cf_cache},
    timeout=5400,
    retries=0,
)
def bench_a10g(entry: dict) -> dict:
    return _bench_core(entry, "A10G")


@app.function(
    image=image,
    gpu="L40S",
    volumes={"/runs": runs_vol, WEIGHTS: cf_cache},
    timeout=5400,
    retries=0,
)
def bench_l40s(entry: dict) -> dict:
    return _bench_core(entry, "L40S")


BENCH = {"T4": bench_t4, "L4": bench_l4, "A10G": bench_a10g, "L40S": bench_l40s}


def _entry(name):
    import csv
    with open("results/c3_screen_input_20261024.csv") as f:
        for r in csv.DictReader(f):
            if r["name"] == name:
                return dict(name=r["name"], binder=r["binder"], target=r["target"],
                            source=r["source"])
    raise SystemExit(f"{name} not in input csv")


@app.local_entrypoint()
def main(name: str = "c3_glp1r_ecd_40", gpus: str = "T4,L4,A10G"):
    e = _entry(name)
    print(f"bench design {name} len={len(e['binder'])}+{len(e['target'])}", flush=True)
    rows = []
    for g in [x.strip() for x in gpus.split(",") if x.strip()]:
        t0 = time.time()
        try:
            res = BENCH[g].remote(e)
        except Exception as ex:
            print(f"[{g}] FAILED {ex!r}", flush=True)
            rows.append({"gpu": g, "error": repr(ex)[:200]})
            continue
        rem = time.time() - t0
        if res.get("error"):
            print(f"[{g}] error {res['error'][:200]}", flush=True)
            rows.append({"gpu": g, "error": res["error"][:200]})
            continue
        h = (res["subprocess_wall_s"] or res["function_wall_s"]) / 3600.0
        res["usd_per_design_subprocess"] = round(h * RATES[g], 4)
        hf = res["function_wall_s"] / 3600.0
        res["usd_per_design_function"] = round(hf * RATES[g], 4)
        res["call_wall_s"] = round(rem, 1)
        rows.append(res)
        print(f"[{g}] iptm={res['iptm']:.3f} sub={res['subprocess_wall_s']:.0f}s "
              f"${res['usd_per_design_subprocess']:.3f}/design (func ${res['usd_per_design_function']:.3f})",
              flush=True)
    out = pathlib.Path("results/c3_gpu_bench_20261024.json")
    with open(out, "w") as f:
        json.dump(rows, f, indent=2)
    print(f"wrote {out}", flush=True)
    for r in rows:
        if "usd_per_design_subprocess" in r:
            print(f"{r['gpu']:>5s}  ${r['usd_per_design_subprocess']:.3f}/design  "
                  f"(${r['usd_per_design_function']:.3f} incl. startup)  "
                  f"iptm={r['iptm']:.3f} plddt={r['mean_plddt']:.1f}", flush=True)
