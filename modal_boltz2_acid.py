"""Boltz-2 complex prediction for C4 acid-site designs that pass ESMFold refold.

Reads refold results JSON, selects PASS designs, predicts binder+TfR1 complexes.
After complex prediction, runs PROPKA for interface His pKa analysis.

Usage:
    modal run modal_boltz2_acid.py::main --refold-json results/c4_acid_combined_refold.json --max-designs 30
"""
import json
import pathlib

import modal

TFR1_APICAL = (
    "GKITFAEKVHAHLGTGDPYTPGFPSPPSRSSGLPVEKLTLDNAAFPFLAYSGIPSQLLSFVRD"
    "LNQYRAKEMGLSLQWLYSVMRVEYHFLSPYVSPKESPFRHVFWGSGSHTLPALLENLKLRKQN"
    "NGAFNETLFRNQLALATWTIQGAANALS"
)

app = modal.App("boltz2-c4-acid")

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
    timeout=14400,
)
def predict_and_score(seqs_json: str):
    """Predict complexes with Boltz-2 and run PROPKA on results."""
    import subprocess
    import shutil
    import yaml

    seqs = json.loads(seqs_json)
    print(f"Predicting {len(seqs)} complexes with Boltz-2", flush=True)

    results = []
    out_base = pathlib.Path("/runs/boltz2_c4_acid")
    out_base.mkdir(parents=True, exist_ok=True)

    for i, entry in enumerate(seqs):
        name = entry["name"]
        binder_seq = entry["sequence"]
        print(f"\n[{i+1}/{len(seqs)}] {name} ({len(binder_seq)} res)", flush=True)

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
        print(f"  Running Boltz-2...", flush=True)

        try:
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=900)
            if result.returncode != 0:
                print(f"  ERROR: {result.stderr[-300:]}", flush=True)
                results.append({"name": name, "error": result.stderr[-200:]})
                continue
        except subprocess.TimeoutExpired:
            print(f"  TIMEOUT after 900s", flush=True)
            results.append({"name": name, "error": "timeout"})
            continue

        pred_dir = out_dir / f"boltz_results_{name}" / "predictions" / name
        if not pred_dir.exists():
            for p in sorted(out_dir.rglob("*")):
                if p.is_file():
                    print(f"    {p.relative_to(out_dir)}", flush=True)
            results.append({"name": name, "error": "no prediction dir"})
            continue

        conf_files = list(pred_dir.rglob("confidence*.json"))
        pdb_files = list(pred_dir.rglob("*.pdb"))

        if not conf_files:
            results.append({"name": name, "error": "no confidence output"})
            continue

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
            dest = out_base / f"{name}_complex.pdb"
            shutil.copy2(pdb_path, dest)

            # PROPKA on complex
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

        results.append(r)

    runs_vol.commit()

    # Summary
    good = [r for r in results if "error" not in r and r.get("iptm", 0) > 0.5]
    print(f"\n{'='*60}")
    print(f"SUMMARY: {len(good)}/{len(results)} with ipTM>0.5")
    for r in sorted(good, key=lambda x: -x.get("iptm", 0)):
        print(f"  {r['name']}: ipTM={r.get('iptm', 0):.3f}, "
              f"complex_His_elevated={r.get('n_elevated_complex', '?')}")

    return results


@app.local_entrypoint()
def main(refold_json: str = "", max_designs: int = 30):
    if not refold_json:
        print("Usage: modal run modal_boltz2_acid.py::main --refold-json <path> [--max-designs 30]")
        return

    with open(refold_json) as f:
        refold = json.load(f)

    if isinstance(refold, dict):
        refold = refold.get("results", [])

    passed = [r for r in refold if r.get("status") == "PASS"]
    passed.sort(key=lambda x: (-x.get("avg_his_pka", 0), -x.get("mean_plddt", 0)))

    selected = passed[:max_designs]
    print(f"Refold results: {len(refold)} total, {len(passed)} PASS, "
          f"selected top {len(selected)} for Boltz-2")

    seqs = [{"name": r["name"], "sequence": r["sequence"]} for r in selected]
    results = predict_and_score.remote(seqs_json=json.dumps(seqs))

    out = "/Users/ricardo/dev/adaptyv-c4/results/c4_acid_boltz2.json"
    with open(out, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\nSaved {len(results)} results to {out}")
