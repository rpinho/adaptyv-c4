"""Boltz-2 multi-seed validation for top C4 designs.

Runs 5 diffusion samples per design to check ipTM variance.
Fable recommendation: if ipTM variance >0.15, score is noise.

Usage:
    modal run modal_boltz2_multiseed.py::validate_top5
"""
import csv
import json
import pathlib

import modal

CSV_PATH = pathlib.Path("/Users/ricardo/dev/adaptyv-c4/results/submission_v8_20261023.csv")

TFR1_APICAL = (
    "GKITFAEKVHAHLGTGDPYTPGFPSPPSRSSGLPVEKLTLDNAAFPFLAYSGIPSQLLSFVRD"
    "LNQYRAKEMGLSLQWLYSVMRVEYHFLSPYVSPKESPFRHVFWGSGSHTLPALLENLKLRKQN"
    "NGAFNETLFRNQLALATWTIQGAANALS"
)

# Top 5 by ipTM — validate these
TOP5 = [
    "c4_v_185_4x",       # ipTM 0.929
    "tfr1_e1_his015_148", # ipTM 0.924
    "tfr1_e1_plain_98",   # ipTM 0.923
    "c4_v_58_4x",         # ipTM 0.916 + ΔΔG 3.958 (best combined)
    "tfr1_e1_his015_115", # ipTM 0.895
]

app = modal.App("boltz2-c4-multiseed")

boltz_cache = modal.Volume.from_name("boltz2-cache", create_if_missing=True)
runs_vol = modal.Volume.from_name("rfd4-runs", create_if_missing=True)

image = (
    modal.Image.debian_slim(python_version="3.10")
    .pip_install("boltz[cuda]", "pyyaml")
    .env({"BOLTZ_CACHE": "/cache/boltz", "REBUILD": "v3"})
    .add_local_file(str(CSV_PATH), "/workspace/sequences.csv")
)


@app.function(
    image=image,
    gpu="L40S",
    volumes={"/cache/boltz": boltz_cache, "/runs": runs_vol},
    timeout=7200,
)
def validate_top5():
    import subprocess
    import shutil
    import yaml

    seqs = {}
    with open("/workspace/sequences.csv") as f:
        for r in csv.DictReader(f):
            seqs[r["name"]] = r["sequence"]

    targets = {name: seqs[name] for name in TOP5 if name in seqs}
    print(f"Multi-seed validation: {len(targets)} designs × 5 seeds", flush=True)

    results = []
    out_base = pathlib.Path("/runs/boltz2_multiseed")
    out_base.mkdir(parents=True, exist_ok=True)

    for name, binder_seq in targets.items():
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
            "--diffusion_samples", "5",
            "--recycling_steps", "3",
        ]
        print(f"Running: {' '.join(cmd)}", flush=True)

        try:
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=900)
            print(f"stdout tail: {result.stdout[-300:]}", flush=True)
            if result.returncode != 0:
                print(f"stderr tail: {result.stderr[-300:]}", flush=True)
                results.append({"name": name, "error": result.stderr[-200:]})
                continue
        except subprocess.TimeoutExpired:
            print("TIMEOUT after 900s", flush=True)
            results.append({"name": name, "error": "timeout"})
            continue

        pred_dir = out_dir / f"boltz_results_{name}" / "predictions" / name
        if not pred_dir.exists():
            for p in sorted(out_dir.rglob("*"))[:20]:
                print(f"  {p.relative_to(out_dir)}", flush=True)
            results.append({"name": name, "error": "pred_dir not found"})
            continue

        conf_files = sorted(pred_dir.glob("confidence_*.json"))
        pdb_files = sorted(pred_dir.glob("*.pdb"))

        seed_scores = []
        for cf in conf_files:
            with open(cf) as f:
                conf = json.load(f)
            seed_scores.append({
                "file": cf.name,
                "iptm": conf.get("iptm", 0),
                "ptm": conf.get("ptm", 0),
                "confidence_score": conf.get("confidence_score", 0),
            })

        if seed_scores:
            iptms = [s["iptm"] for s in seed_scores]
            mean_iptm = sum(iptms) / len(iptms)
            var_iptm = max(iptms) - min(iptms)
            print(f"  Seeds: {len(seed_scores)}", flush=True)
            print(f"  ipTM: mean={mean_iptm:.3f}, min={min(iptms):.3f}, max={max(iptms):.3f}, range={var_iptm:.3f}", flush=True)
            print(f"  RELIABLE: {'YES' if var_iptm < 0.15 else 'NO (variance > 0.15)'}", flush=True)

            entry = {
                "name": name,
                "seeds": seed_scores,
                "mean_iptm": mean_iptm,
                "min_iptm": min(iptms),
                "max_iptm": max(iptms),
                "range_iptm": var_iptm,
                "reliable": var_iptm < 0.15,
            }
            results.append(entry)

            # Save best PDB
            best_idx = iptms.index(max(iptms))
            if best_idx < len(pdb_files):
                dest = out_base / f"{name}_best.pdb"
                shutil.copy2(pdb_files[best_idx], dest)
        else:
            results.append({"name": name, "error": "no confidence files"})

    # Summary
    print(f"\n\n{'='*80}", flush=True)
    print("MULTI-SEED VALIDATION SUMMARY", flush=True)
    print(f"{'='*80}", flush=True)
    print(f"{'Design':<35s} {'Mean':>6s} {'Min':>6s} {'Max':>6s} {'Range':>6s} {'OK?':>5s}", flush=True)
    for r in results:
        if "error" in r:
            print(f"{r['name']:<35s} ERROR: {r['error']}", flush=True)
        else:
            ok = "YES" if r["reliable"] else "NO"
            print(f"{r['name']:<35s} {r['mean_iptm']:>6.3f} {r['min_iptm']:>6.3f} "
                  f"{r['max_iptm']:>6.3f} {r['range_iptm']:>6.3f} {ok:>5s}", flush=True)

    out_json = out_base / "multiseed_results.json"
    with open(out_json, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\nSaved to {out_json}", flush=True)

    runs_vol.commit()
    return results
