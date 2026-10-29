"""Boltz-2 complex prediction for C4 replacement candidates.

Evaluating 4 pool designs to replace weak ΔΔG_switch designs in v8.
"""
import json
import pathlib

import modal

TFR1_APICAL = (
    "GKITFAEKVHAHLGTGDPYTPGFPSPPSRSSGLPVEKLTLDNAAFPFLAYSGIPSQLLSFVRD"
    "LNQYRAKEMGLSLQWLYSVMRVEYHFLSPYVSPKESPFRHVFWGSGSHTLPALLENLKLRKQN"
    "NGAFNETLFRNQLALATWTIQGAANALS"
)

REPLACEMENTS = {
    "tfr1_e1_his015_105": "SLSEDITEHVSKFAKRLFLSDTQTDELISTVLSLFNENSDISFDDFEHDIKSKLSHLGLSHHDIEELYHNVQHIIESH",
    "tfr1_e1_his015_118": "STFEQALTEAHTLLEEHQPSTLSELLGIINKVMQKHHIPSEHQLRFSNMVIHWYMEHHPEFKFPHF",
    "tfr1_e1_his015_176": "SRILILTVTSTTTDIALEVANHHTDVDVMFAAGQLALSSVIDKVLSHHYDGVLVVSLVGPNYHHDESTAHFFHHALHEKLPSHVFVDLLLLNATDSRDKQLHELKEAVHNLKSKIASKKSS",
    "tfr1_e1_his015_60": "SVSGDPLDSFAEELFEELTSHLLGKTRAEAKTLLHEILDEHHLSHQLLEETDDKDVYRIYDGLRHFDLHLHYENDVVTKIELKDNTNGNTKTAHL",
}

app = modal.App("boltz2-c4-replacements")
boltz_cache = modal.Volume.from_name("boltz2-cache", create_if_missing=True)
runs_vol = modal.Volume.from_name("rfd4-runs", create_if_missing=True)

image = (
    modal.Image.debian_slim(python_version="3.10")
    .pip_install("boltz[cuda]", "pyyaml")
    .env({"BOLTZ_CACHE": "/cache/boltz"})
)


@app.function(
    image=image,
    gpu="L40S",
    volumes={"/cache/boltz": boltz_cache, "/runs": runs_vol},
    timeout=3600,
)
def predict_replacements():
    import subprocess
    import shutil
    import yaml

    out_base = pathlib.Path("/runs/boltz2_replacements")
    out_base.mkdir(parents=True, exist_ok=True)

    print(f"Predicting complexes for {len(REPLACEMENTS)} replacement candidates", flush=True)

    results = []
    for name, binder_seq in REPLACEMENTS.items():
        print(f"\n{'='*60}", flush=True)
        print(f"Processing {name} ({len(binder_seq)} res)", flush=True)

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
        print(f"Running: {' '.join(cmd)}", flush=True)

        try:
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
            print(f"stdout tail: {result.stdout[-200:]}", flush=True)
            if result.returncode != 0:
                print(f"stderr tail: {result.stderr[-200:]}", flush=True)
                results.append({"name": name, "error": result.stderr[-200:]})
                continue
        except subprocess.TimeoutExpired:
            print("TIMEOUT", flush=True)
            results.append({"name": name, "error": "timeout"})
            continue

        pred_dir = out_dir / f"boltz_results_{name}" / "predictions" / name
        if not pred_dir.exists():
            results.append({"name": name, "error": "pred_dir not found"})
            continue

        conf_file = sorted(pred_dir.glob("confidence_*.json"))[0]
        pdb_file = sorted(pred_dir.glob("*.pdb"))[0]

        with open(conf_file) as f:
            conf = json.load(f)

        entry = {
            "name": name,
            "iptm": conf.get("iptm", 0),
            "ptm": conf.get("ptm", 0),
            "confidence_score": conf.get("confidence_score", 0),
            "complex_plddt": conf.get("complex_plddt", 0),
        }
        results.append(entry)
        print(f"  ipTM: {entry['iptm']:.3f}, pTM: {entry['ptm']:.3f}", flush=True)

        dest_pdb = out_base / f"{name}_complex.pdb"
        dest_conf = out_base / f"{name}_confidence.json"
        shutil.copy2(pdb_file, dest_pdb)
        shutil.copy2(conf_file, dest_conf)

    print(f"\n\n{'='*80}", flush=True)
    print("REPLACEMENT CANDIDATES - BOLTZ-2 RESULTS", flush=True)
    print(f"{'='*80}", flush=True)
    print(f"{'Design':<35s} {'ipTM':>6s} {'pTM':>6s} {'Conf':>6s}", flush=True)
    for r in results:
        if "error" in r:
            print(f"{r['name']:<35s} ERROR: {r['error']}", flush=True)
        else:
            print(f"{r['name']:<35s} {r['iptm']:>6.3f} {r['ptm']:>6.3f} {r['confidence_score']:>6.3f}", flush=True)

    out_json = out_base / "replacement_results.json"
    with open(out_json, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\nSaved to {out_json}", flush=True)

    runs_vol.commit()
    return results
