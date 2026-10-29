"""ESMFold refold for C4 v8 replacement candidate (his015_165).

Usage:
    modal run modal_esmfold_c4_v8.py::refold_all
"""
import csv
import pathlib

import modal

CSV_PATH = pathlib.Path("/Users/ricardo/dev/adaptyv-c4/results/c4_v8_candidates.csv")

app = modal.App("esmfold-c4-v8")

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
    timeout=600,
)
def refold_all():
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

    for r in results:
        print(f"  {r['name']}: pLDDT={r['mean_plddt']}, {r['status']}")

    return {"results": results}
