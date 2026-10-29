"""Boltz-2 complex prediction — PARALLEL version (one container per design).

Same logic as modal_boltz2_acid.py but fans out to parallel containers.
Faster at the cost of more concurrent GPU allocation.

Usage:
    modal run modal_boltz2_acid_parallel.py::main --refold-json results/c4_acid_combined_refold.json --max-designs 30
"""
import json
import pathlib

import modal

TFR1_APICAL = (
    "GKITFAEKVHAHLGTGDPYTPGFPSPPSRSSGLPVEKLTLDNAAFPFLAYSGIPSQLLSFVRD"
    "LNQYRAKEMGLSLQWLYSVMRVEYHFLSPYVSPKESPFRHVFWGSGSHTLPALLENLKLRKQN"
    "NGAFNETLFRNQLALATWTIQGAANALS"
)

app = modal.App("boltz2-c4-acid-par")

boltz_cache = modal.Volume.from_name("boltz2-cache", create_if_missing=True)
runs_vol = modal.Volume.from_name("rfd4-runs", create_if_missing=True)

image = (
    modal.Image.debian_slim(python_version="3.10")
    .pip_install(
        "boltz[cuda]",
        "pyyaml",
        "propka",
        "biopython==1.84",
    )
    .env({"BOLTZ_CACHE": "/cache/boltz", "REBUILD": "v2"})
)


@app.function(
    image=image,
    gpu="L40S",
    volumes={"/cache/boltz": boltz_cache, "/runs": runs_vol},
    timeout=3600,
)
def predict_one(name: str, binder_seq: str):
    """Predict one complex with Boltz-2 and run PROPKA."""
    import subprocess
    import shutil
    import yaml

    print(f"Predicting {name} ({len(binder_seq)} res)", flush=True)

    yaml_dir = pathlib.Path(f"/tmp/boltz2/{name}")
    yaml_dir.mkdir(parents=True, exist_ok=True)
    yaml_path = yaml_dir / f"{name}.yaml"

    config = {
        "version": 1,
        "sequences": [
            {"protein": {"id": "A", "sequence": binder_seq}},
            {"protein": {"id": "B", "sequence": TFR1_APICAL}},
        ],
    }
    with open(yaml_path, "w") as f:
        yaml.dump(config, f)

    out_dir = yaml_dir / "output"
    out_dir.mkdir(exist_ok=True)

    cmd = [
        "boltz", "predict",
        str(yaml_path),
        "--use_msa_server",
        "--out_dir", str(out_dir),
        "--output_format", "pdb",
        "--diffusion_samples", "1",
        "--recycling_steps", "3",
    ]

    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=1800)
        if result.returncode != 0:
            print(f"  ERROR: {result.stderr[-300:]}", flush=True)
            return {"name": name, "error": result.stderr[-200:]}
    except subprocess.TimeoutExpired:
        print(f"  TIMEOUT", flush=True)
        return {"name": name, "error": "timeout"}

    pred_dir = out_dir / f"boltz_results_{name}" / "predictions" / name
    if not pred_dir.exists():
        return {"name": name, "error": "no prediction dir"}

    conf_files = list(pred_dir.rglob("confidence*.json"))
    pdb_files = list(pred_dir.rglob("*.pdb"))

    if not conf_files:
        return {"name": name, "error": "no confidence output"}

    with open(conf_files[0]) as f:
        conf = json.load(f)

    iptm = conf.get("iptm", 0)
    ptm = conf.get("ptm", 0)
    confidence = conf.get("confidence_score", 0)
    print(f"  ipTM={iptm:.3f}, pTM={ptm:.3f}, conf={confidence:.3f}", flush=True)

    r = {"name": name, "sequence": binder_seq, "iptm": iptm, "ptm": ptm,
         "confidence": confidence, "complex_plddt": conf.get("complex_plddt", 0)}

    if pdb_files:
        pdb_path = str(pdb_files[0])
        out_base = pathlib.Path("/runs/boltz2_c4_acid")
        out_base.mkdir(parents=True, exist_ok=True)
        dest = out_base / f"{name}_complex.pdb"
        shutil.copy2(pdb_path, dest)

        try:
            from propka.run import single as propka_single
            mol = propka_single(pdb_path)
            best_conf = max(mol.conformations.values(), key=lambda c: len(c.groups))

            his_bound = {}
            for group in best_conf.groups:
                if type(group).__name__ == "HISGroup":
                    chain = group.atom.chain_id
                    resnum = group.atom.res_num
                    his_bound[f"{chain}{resnum}"] = round(group.pka_value, 2)

            r["his_complex_pkas"] = his_bound
            if his_bound:
                binder_his = {k: v for k, v in his_bound.items() if k.startswith("A")}
                r["avg_complex_his_pka"] = round(
                    sum(binder_his.values()) / len(binder_his), 2) if binder_his else 0
                r["n_elevated_complex"] = sum(1 for v in binder_his.values() if v > 6.5)
                print(f"  Complex His pKa (binder): {binder_his}", flush=True)
        except Exception as e:
            print(f"  PROPKA error: {e}", flush=True)
            r["propka_error"] = str(e)

    runs_vol.commit()
    return r


@app.local_entrypoint()
def main(refold_json: str = "", max_designs: int = 30):
    if not refold_json:
        print("Usage: modal run modal_boltz2_acid_parallel.py::main --refold-json <path>")
        return

    with open(refold_json) as f:
        refold = json.load(f)

    if isinstance(refold, dict):
        refold = refold.get("results", [])

    passed = [r for r in refold if r.get("status") == "PASS"]
    passed.sort(key=lambda x: (-x.get("avg_his_pka", 0), -x.get("mean_plddt", 0)))

    selected = passed[:max_designs]
    print(f"Refold: {len(refold)} total, {len(passed)} PASS, selected top {len(selected)}")

    results = list(predict_one.starmap(
        [(r["name"], r["sequence"]) for r in selected]
    ))

    tag = refold_json.replace(".json", "").split("_")[-1] if "_batch" in refold_json else "batch1"
    out = f"/Users/ricardo/dev/adaptyv-c4/results/c4_acid_boltz2_{tag}.json"
    with open(out, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\nSaved {len(results)} results to {out}")

    good = [r for r in results if "error" not in r and r.get("iptm", 0) > 0.5]
    print(f"\n{len(good)}/{len(results)} with ipTM>0.5")
    for r in sorted(good, key=lambda x: -x.get("iptm", 0)):
        print(f"  {r['name']}: ipTM={r.get('iptm',0):.3f}, "
              f"complex_His_elevated={r.get('n_elevated_complex', '?')}")
