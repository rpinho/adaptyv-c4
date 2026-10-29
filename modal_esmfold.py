"""ESMFold single-sequence refold on Modal GPU.

Folds all submission candidates and computes pLDDT as a quality gate.

Usage:
    modal run modal_esmfold.py::refold_all
"""
import csv
import io
import json
import os
import pathlib

import modal

RESULTS_DIR = pathlib.Path("/Users/ricardo/dev/adaptyv-c4/results")
CSV_PATH = RESULTS_DIR / "submission_reviewed_20261023.csv"

app = modal.App("esmfold-refold")

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
    )
    .env({"HF_HOME": "/cache/huggingface"})
    .add_local_file(str(CSV_PATH), "/workspace/sequences.csv")
)


@app.function(
    image=image,
    gpu="L40S",
    volumes={"/cache/huggingface": hf_cache, "/runs": runs_vol},
    timeout=1800,
)
def refold_all():
    """Fold all submission sequences with ESMFold, return pLDDT stats."""
    import torch
    from transformers import AutoTokenizer, EsmForProteinFolding
    from transformers.models.esm.openfold_utils.protein import Protein, to_pdb

    print(f"Loading ESMFold model...", flush=True)
    tokenizer = AutoTokenizer.from_pretrained("facebook/esmfold_v1")
    model = EsmForProteinFolding.from_pretrained(
        "facebook/esmfold_v1", low_cpu_mem_usage=True
    )
    model = model.cuda()
    model.eval()
    model.esm = model.esm.half()

    results = []
    seqs = {}
    with open("/workspace/sequences.csv") as f:
        for r in csv.DictReader(f):
            seqs[r["name"]] = r["sequence"]
    print(f"Folding {len(seqs)} sequences...", flush=True)

    for i, (name, seq) in enumerate(seqs.items()):
        print(f"[{i+1}/{len(seqs)}] {name} ({len(seq)} res)...", end=" ", flush=True)

        tokenized = tokenizer([seq], return_tensors="pt", add_special_tokens=False)
        tokenized = {k: v.cuda() for k, v in tokenized.items()}

        with torch.no_grad():
            output = model(**tokenized)

        # Get PDB with pLDDT in B-factors via model's built-in method
        pdb_text = model.output_to_pdb(output)[0]
        pdb_dir = "/runs/esmfold_refolds"
        os.makedirs(pdb_dir, exist_ok=True)
        with open(f"{pdb_dir}/{name}.pdb", "w") as f:
            f.write(pdb_text)

        # Parse pLDDT from B-factor column of CA atoms
        # HF ESMFold writes pLDDT on 0-1 scale; multiply by 100
        plddt_vals = []
        for line in pdb_text.splitlines():
            if line.startswith("ATOM") and line[12:16].strip() == "CA":
                bfactor = float(line[60:66].strip())
                plddt_vals.append(bfactor * 100 if bfactor <= 1.0 else bfactor)

        if not plddt_vals:
            print("NO CA ATOMS", flush=True)
            results.append({"name": name, "length": len(seq), "mean_plddt": -1,
                            "min_plddt": -1, "plddt_gt70_frac": 0, "status": "PARSE_FAILED"})
            del output, tokenized
            torch.cuda.empty_cache()
            continue

        mean_plddt = sum(plddt_vals) / len(plddt_vals)
        min_plddt = min(plddt_vals)
        gt70_frac = sum(1 for p in plddt_vals if p > 70) / len(plddt_vals)

        status = "PASS" if mean_plddt >= 70 else ("MARGINAL" if mean_plddt >= 50 else "FAIL")
        print(f"pLDDT={mean_plddt:.1f} (min={min_plddt:.1f}, >70: {gt70_frac*100:.0f}%) {status}", flush=True)

        results.append({
            "name": name,
            "length": len(seq),
            "mean_plddt": round(mean_plddt, 1),
            "min_plddt": round(min_plddt, 1),
            "plddt_gt70_frac": round(gt70_frac, 3),
            "status": status,
        })

        # Clear GPU cache between folds
        del output, tokenized
        torch.cuda.empty_cache()

    # Save results CSV to volume
    csv_out = "/runs/esmfold_refold_results.csv"
    buf = io.StringIO()
    w = csv.DictWriter(buf, fieldnames=["name", "length", "mean_plddt", "min_plddt",
                                         "plddt_gt70_frac", "status"])
    w.writeheader()
    w.writerows(results)
    with open(csv_out, "w") as f:
        f.write(buf.getvalue())

    # Summary
    passed = [r for r in results if r["status"] == "PASS"]
    marginal = [r for r in results if r["status"] == "MARGINAL"]
    failed = [r for r in results if r["status"] in ("FAIL",)]
    print(f"\nSummary: {len(passed)} PASS, {len(marginal)} MARGINAL, {len(failed)} FAIL", flush=True)
    print(f"Saved to {csv_out}", flush=True)

    return {"results": results, "csv": buf.getvalue()}


def convert_to_pdb(output, sequence):
    """Convert ESMFold output to minimal PDB text with pLDDT in B-factor."""
    import torch
    positions = output["positions"][-1].squeeze().cpu().numpy()
    plddt = output["plddt"].squeeze().cpu().numpy() * 100

    atom_names = ["N", "CA", "C", "O", "CB"]
    lines = []
    atom_idx = 1

    for res_i, aa in enumerate(sequence):
        n_atoms = min(len(atom_names), positions.shape[1]) if positions.ndim == 3 else 1
        for a_i in range(n_atoms):
            if a_i == 4 and aa == "G":
                continue
            aname = atom_names[a_i] if positions.ndim == 3 else "CA"
            if positions.ndim == 3:
                x, y, z = positions[res_i, a_i, :3]
            else:
                x, y, z = positions[res_i, :3]
            b = plddt[res_i] if plddt.ndim > 0 else float(plddt)
            lines.append(
                f"ATOM  {atom_idx:5d}  {aname:<3s} {'ALA':>3s} A{res_i+1:4d}    "
                f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00{b:6.2f}"
            )
            atom_idx += 1
    lines.append("END")
    return "\n".join(lines)


@app.function(image=image, volumes={"/runs": runs_vol}, timeout=120)
def fetch_results():
    """Fetch refold results CSV from the volume."""
    csv_path = "/runs/esmfold_refold_results.csv"
    if os.path.exists(csv_path):
        with open(csv_path) as f:
            return f.read()
    return "NOT_FOUND"
