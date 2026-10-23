"""ESMFold single-sequence refold for C3 (GLP1R) designs on Modal GPU.

Usage:
    modal run modal_esmfold_c3.py::refold_all
"""
import csv
import io
import os
import pathlib

import modal

RESULTS_DIR = pathlib.Path("/Users/ricardo/dev/adaptyv-c4/results")
CSV_PATH = RESULTS_DIR / "c3_sequences.csv"

app = modal.App("esmfold-refold-c3")

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
    """Fold all C3 submission sequences with ESMFold, return pLDDT stats."""
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

        pdb_text = model.output_to_pdb(output)[0]
        pdb_dir = "/runs/esmfold_c3_refolds"
        os.makedirs(pdb_dir, exist_ok=True)
        with open(f"{pdb_dir}/{name}.pdb", "w") as f:
            f.write(pdb_text)

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

        del output, tokenized
        torch.cuda.empty_cache()

    csv_out = "/runs/esmfold_c3_refold_results.csv"
    buf = io.StringIO()
    w = csv.DictWriter(buf, fieldnames=["name", "length", "mean_plddt", "min_plddt",
                                         "plddt_gt70_frac", "status"])
    w.writeheader()
    w.writerows(results)
    with open(csv_out, "w") as f:
        f.write(buf.getvalue())

    passed = [r for r in results if r["status"] == "PASS"]
    marginal = [r for r in results if r["status"] == "MARGINAL"]
    failed = [r for r in results if r["status"] in ("FAIL",)]
    print(f"\nSummary: {len(passed)} PASS, {len(marginal)} MARGINAL, {len(failed)} FAIL", flush=True)
    print(f"Saved to {csv_out}", flush=True)

    return {"results": results, "csv": buf.getvalue()}
