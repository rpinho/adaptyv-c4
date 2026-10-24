"""ESMFold + PROPKA analysis for C4 v8 designs.

Folds all 20 designs, saves PDBs, runs PROPKA on each monomer to get
His pKa values. Reports surface vs buried His and pKa environment effects.

Usage:
    modal run modal_esmfold_c4_propka.py::fold_and_propka
"""
import csv
import json
import pathlib

import modal

CSV_PATH = pathlib.Path("/Users/ricardo/dev/adaptyv-c4/results/submission_v8_20261023.csv")

app = modal.App("esmfold-c4-propka")

hf_cache = modal.Volume.from_name("rfd4-hf-cache", create_if_missing=True)
runs_vol = modal.Volume.from_name("rfd4-runs", create_if_missing=True)

image = (
    modal.Image.debian_slim(python_version="3.10")
    .pip_install(
        "numpy<2",
        "torch==2.1.2",
        "transformers==4.36.2",
        "accelerate==0.25.0",
        "biopython==1.84",
        "propka==3.5.1",
    )
    .env({"HF_HOME": "/cache/huggingface"})
    .add_local_file(str(CSV_PATH), "/workspace/sequences.csv")
)


@app.function(
    image=image,
    gpu="L40S",
    volumes={"/cache/huggingface": hf_cache, "/runs": runs_vol},
    timeout=900,
)
def fold_and_propka():
    import io
    import tempfile
    import torch
    from transformers import AutoTokenizer, EsmForProteinFolding

    print("Loading ESMFold model...", flush=True)
    tokenizer = AutoTokenizer.from_pretrained("facebook/esmfold_v1")
    model = EsmForProteinFolding.from_pretrained(
        "facebook/esmfold_v1", low_cpu_mem_usage=True
    )
    model = model.cuda()
    model.eval()
    model.esm = model.esm.half()

    seqs = {}
    with open("/workspace/sequences.csv") as f:
        for r in csv.DictReader(f):
            seqs[r["name"]] = r["sequence"]
    print(f"Folding {len(seqs)} sequences...", flush=True)

    results = []
    pdb_texts = {}

    for i, (name, seq) in enumerate(seqs.items()):
        print(f"\n[{i+1}/{len(seqs)}] {name} ({len(seq)} res)...", flush=True)

        tokenized = tokenizer([seq], return_tensors="pt", add_special_tokens=False)
        tokenized = {k: v.cuda() for k, v in tokenized.items()}

        with torch.no_grad():
            output = model(**tokenized)

        pdb_text = model.output_to_pdb(output)[0]
        pdb_texts[name] = pdb_text

        plddt_vals = []
        for line in pdb_text.splitlines():
            if line.startswith("ATOM") and line[12:16].strip() == "CA":
                bfactor = float(line[60:66].strip())
                plddt_vals.append(bfactor * 100 if bfactor <= 1.0 else bfactor)

        mean_plddt = sum(plddt_vals) / len(plddt_vals) if plddt_vals else 0
        print(f"  pLDDT={mean_plddt:.1f}", flush=True)

        # Run PROPKA on this monomer
        propka_result = run_propka(name, pdb_text, seq)

        results.append({
            "name": name,
            "length": len(seq),
            "mean_plddt": round(mean_plddt, 1),
            "his_count": seq.count("H"),
            "his_fraction": round(seq.count("H") / len(seq) * 100, 1),
            **propka_result,
        })

        del output, tokenized
        torch.cuda.empty_cache()

    print("\n" + "=" * 80)
    print("PROPKA SUMMARY")
    print("=" * 80)
    print(f"{'Design':<30s} {'pLDDT':>6s} {'His':>4s} {'pKa_mean':>8s} {'pKa_range':>12s} {'surf':>5s} {'bur':>5s} {'ΔΔG_est':>8s}")
    print("-" * 80)
    for r in results:
        pka_range = f"{r.get('pka_min', 0):.1f}-{r.get('pka_max', 0):.1f}" if r.get("his_propka_count", 0) > 0 else "N/A"
        print(f"{r['name']:<30s} {r['mean_plddt']:>6.1f} {r['his_count']:>4d} "
              f"{r.get('pka_mean', 0):>8.2f} {pka_range:>12s} "
              f"{r.get('surface_his', 0):>5d} {r.get('buried_his', 0):>5d} "
              f"{r.get('ddg_monomer', 0):>8.2f}")

    # Save PDBs to volume
    pdb_dir = pathlib.Path("/runs/c4_v8_pdbs")
    pdb_dir.mkdir(parents=True, exist_ok=True)
    for name, pdb in pdb_texts.items():
        (pdb_dir / f"{name}.pdb").write_text(pdb)
    runs_vol.commit()
    print(f"\nSaved {len(pdb_texts)} PDB files to /runs/c4_v8_pdbs/")

    return {"results": results, "pdb_texts": pdb_texts}


def run_propka(name, pdb_text, sequence):
    """Run PROPKA on a monomer PDB and extract His pKa values."""
    import tempfile
    import os

    try:
        from propka.run import single

        with tempfile.NamedTemporaryFile(mode="w", suffix=".pdb", delete=False) as f:
            f.write(pdb_text)
            pdb_path = f.name

        try:
            mc = single(pdb_path)

            his_pkas = []
            for conformation in mc.conformations:
                for group in mc.conformations[conformation].groups:
                    if type(group).__name__ == "HISGroup":
                        his_pkas.append({
                            "res_num": group.res_num,
                            "pka": group.pka_value,
                            "model_pka": group.model_pka,
                            "shift": group.pka_value - group.model_pka,
                        })
                break  # first conformation only

            if not his_pkas:
                print(f"  PROPKA: no His residues found", flush=True)
                return {
                    "his_propka_count": 0,
                    "pka_mean": 0,
                    "pka_min": 0,
                    "pka_max": 0,
                    "surface_his": 0,
                    "buried_his": 0,
                    "ddg_monomer": 0,
                    "his_pkas": [],
                }

            pkas = [h["pka"] for h in his_pkas]
            # Surface His: pKa close to model (6.50 ± 0.5)
            surface = sum(1 for p in pkas if 5.5 <= p <= 7.5)
            buried = len(pkas) - surface

            # Wyman-linkage monomer estimate:
            # ΔΔG_monomer = -RT * Σ ln[(1 + 10^(pKa_i - 5.5)) / (1 + 10^(pKa_i - 7.4))]
            import math
            RT = 0.592  # kcal/mol at 25°C
            ddg = 0
            for h in his_pkas:
                pka = h["pka"]
                ratio = (1 + 10**(pka - 5.5)) / (1 + 10**(pka - 7.4))
                ddg += -RT * math.log(ratio)

            print(f"  PROPKA: {len(his_pkas)} His, pKa mean={sum(pkas)/len(pkas):.2f} "
                  f"(range {min(pkas):.1f}-{max(pkas):.1f}), "
                  f"surface={surface}, buried={buried}, "
                  f"ΔΔG_monomer={ddg:.2f} kcal/mol", flush=True)

            for h in his_pkas:
                env = "surface" if 5.5 <= h["pka"] <= 7.5 else ("buried" if h["pka"] > 7.5 else "shifted-low")
                print(f"    H{h['res_num']:3d}: pKa={h['pka']:.2f} (model={h['model_pka']:.2f}, shift={h['shift']:+.2f}) [{env}]", flush=True)

            return {
                "his_propka_count": len(his_pkas),
                "pka_mean": round(sum(pkas) / len(pkas), 2),
                "pka_min": round(min(pkas), 2),
                "pka_max": round(max(pkas), 2),
                "surface_his": surface,
                "buried_his": buried,
                "ddg_monomer": round(ddg, 2),
                "his_pkas": his_pkas,
            }
        finally:
            os.unlink(pdb_path)

    except Exception as e:
        print(f"  PROPKA ERROR: {e}", flush=True)
        return {
            "his_propka_count": 0,
            "pka_mean": 0,
            "pka_min": 0,
            "pka_max": 0,
            "surface_his": 0,
            "buried_his": 0,
            "ddg_monomer": 0,
            "his_pkas": [],
            "propka_error": str(e),
        }
