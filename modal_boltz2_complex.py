"""Boltz-2 complex prediction for top C4 v8 designs.

Predicts binder-TfR1 complex structures, extracts ipTM scores,
and saves complex PDBs for PROPKA interface analysis.

Usage:
    modal run modal_boltz2_complex.py::predict_complexes
"""
import csv
import json
import pathlib

import modal

CSV_PATH = pathlib.Path("/Users/ricardo/dev/adaptyv-c4/results/submission_v8_20261023.csv")

# TfR1 apical domain sequence (from E1_his707.pdb chain A, res 260-750)
TFR1_APICAL = (
    "GKITFAEKVHAHLGTGDPYTPGFPSPPSRSSGLPVEKLTLDNAAFPFLAYSGIPSQLLSFVRD"
    "LNQYRAKEMGLSLQWLYSVMRVEYHFLSPYVSPKESPFRHVFWGSGSHTLPALLENLKLRKQN"
    "NGAFNETLFRNQLALATWTIQGAANALS"
)

# Top 5 by composite score (from PROPKA + spatial analysis)
TOP5 = [
    "c4_v_58_4x",
    "tfr1_e1_his015_131",
    "tfr1_e1_his015_145",
    "c4_v_33_4x_v2",
    "c4_v_62_max",
]

app = modal.App("boltz2-c4-complex")

boltz_cache = modal.Volume.from_name("boltz2-cache", create_if_missing=True)
runs_vol = modal.Volume.from_name("rfd4-runs", create_if_missing=True)

image = (
    modal.Image.debian_slim(python_version="3.10")
    .pip_install(
        "boltz[cuda]",
        "pyyaml",
    )
    .env({"BOLTZ_CACHE": "/cache/boltz"})
    .add_local_file(str(CSV_PATH), "/workspace/sequences.csv")
)


@app.function(
    image=image,
    gpu="L40S",
    volumes={"/cache/boltz": boltz_cache, "/runs": runs_vol},
    timeout=3600,  # 1 hour — MSA + folding for 5 complexes
)
def predict_complexes():
    import subprocess
    import os
    import yaml

    # Load sequences
    seqs = {}
    with open("/workspace/sequences.csv") as f:
        for r in csv.DictReader(f):
            seqs[r["name"]] = r["sequence"]

    # Filter to top 5
    targets = {name: seqs[name] for name in TOP5 if name in seqs}
    print(f"Predicting {len(targets)} complexes with Boltz-2", flush=True)
    print(f"TfR1 target: {len(TFR1_APICAL)} residues", flush=True)

    results = []
    out_base = pathlib.Path("/runs/boltz2_complexes")
    out_base.mkdir(parents=True, exist_ok=True)

    for name, binder_seq in targets.items():
        print(f"\n{'='*60}", flush=True)
        print(f"Processing {name} ({len(binder_seq)} res binder + {len(TFR1_APICAL)} res target)", flush=True)

        # Write YAML input
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

        # Run Boltz-2
        cmd = [
            "boltz", "predict",
            str(yaml_path),
            "--use_msa_server",
            "--out_dir", str(out_dir),
            "--output_format", "pdb",
            "--diffusion_samples", "1",
            "--recycling_steps", "3",
        ]
        print(f"Running: {' '.join(cmd)}", flush=True)

        try:
            result = subprocess.run(
                cmd, capture_output=True, text=True, timeout=600,
            )
            print(f"stdout: {result.stdout[-500:]}", flush=True)
            if result.returncode != 0:
                print(f"stderr: {result.stderr[-500:]}", flush=True)
                results.append({"name": name, "error": result.stderr[-200:]})
                continue
        except subprocess.TimeoutExpired:
            print(f"TIMEOUT after 600s", flush=True)
            results.append({"name": name, "error": "timeout"})
            continue

        # Parse confidence output
        pred_dir = out_dir / "predictions" / name
        conf_files = list(pred_dir.glob("confidence_*.json")) if pred_dir.exists() else []
        pdb_files = list(pred_dir.glob("*.pdb")) if pred_dir.exists() else []

        if conf_files:
            with open(conf_files[0]) as f:
                conf = json.load(f)
            print(f"  ipTM={conf.get('iptm', 'N/A')}", flush=True)
            print(f"  pTM={conf.get('ptm', 'N/A')}", flush=True)
            print(f"  confidence={conf.get('confidence_score', 'N/A')}", flush=True)
            print(f"  complex_plddt={conf.get('complex_plddt', 'N/A')}", flush=True)

            entry = {"name": name, **conf}
            results.append(entry)

            # Copy PDB and confidence to volume
            if pdb_files:
                import shutil
                dest = out_base / f"{name}_complex.pdb"
                shutil.copy2(pdb_files[0], dest)
                print(f"  Saved complex PDB to {dest}", flush=True)

            dest_conf = out_base / f"{name}_confidence.json"
            shutil.copy2(conf_files[0], dest_conf)
        else:
            print(f"  No confidence output found", flush=True)
            results.append({"name": name, "error": "no confidence output"})

    # Summary
    print(f"\n\n{'='*80}", flush=True)
    print("BOLTZ-2 COMPLEX PREDICTION SUMMARY", flush=True)
    print(f"{'='*80}", flush=True)
    print(f"{'Design':<35s} {'ipTM':>6s} {'pTM':>6s} {'Conf':>6s} {'cpLDDT':>7s}", flush=True)
    print("-" * 80, flush=True)

    for r in results:
        if "error" in r:
            print(f"{r['name']:<35s} ERROR: {r['error']}", flush=True)
        else:
            print(f"{r['name']:<35s} {r.get('iptm', 0):>6.3f} {r.get('ptm', 0):>6.3f} "
                  f"{r.get('confidence_score', 0):>6.3f} {r.get('complex_plddt', 0):>7.1f}", flush=True)

    # Save results
    out_json = out_base / "boltz2_results.json"
    with open(out_json, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\nSaved results to {out_json}", flush=True)

    runs_vol.commit()
    return results
