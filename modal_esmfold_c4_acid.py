"""ESMFold refold + monomer PROPKA for C4 acid-site RFD4 designs.

Reads binder sequences from a CSV (output of extract_sequences),
refolds with ESMFold, checks pLDDT, and runs PROPKA for free-state His pKa.

Usage:
    # First extract sequences from the RFD4 volume:
    modal run modal_rfd4.py::extract_sequences --run-name tfr1_c4_acid1
    # Then filter locally and run refold:
    python3 scripts/filter_his_sequences.py results/c4_acid1_sequences.csv > results/c4_acid1_filtered.csv
    modal run modal_esmfold_c4_acid.py::main --csv results/c4_acid1_filtered.csv --tag acid1
"""
import csv
import json
import pathlib

import modal

app = modal.App("esmfold-c4-acid")

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
        "propka",
    )
    .env({"HF_HOME": "/cache/huggingface"})
)


@app.function(
    image=image,
    gpu="L40S",
    volumes={"/cache/huggingface": hf_cache, "/runs": runs_vol},
    timeout=3600,
)
def refold_batch(seqs_json: str):
    """Refold sequences and run monomer PROPKA. Input: JSON list of {name, sequence}."""
    import torch
    from transformers import AutoTokenizer, EsmForProteinFolding

    seqs = json.loads(seqs_json)
    print(f"Refolding {len(seqs)} sequences...", flush=True)

    print("Loading ESMFold model...", flush=True)
    tokenizer = AutoTokenizer.from_pretrained("facebook/esmfold_v1")
    model = EsmForProteinFolding.from_pretrained(
        "facebook/esmfold_v1", low_cpu_mem_usage=True
    )
    model = model.cuda().eval()
    model.esm = model.esm.half()
    print("Model loaded.", flush=True)

    results = []
    for i, entry in enumerate(seqs):
        name = entry["name"]
        seq = entry["sequence"]
        n_his = seq.count("H")
        print(f"\n[{i+1}/{len(seqs)}] {name} ({len(seq)} aa, {n_his} His)", flush=True)

        tokenized = tokenizer([seq], return_tensors="pt", add_special_tokens=False)
        tokenized = {k: v.cuda() for k, v in tokenized.items()}

        with torch.no_grad():
            output = model(**tokenized)

        pdb_text = model.output_to_pdb(output)[0]

        plddt_vals = []
        for line in pdb_text.splitlines():
            if line.startswith("ATOM") and line[12:16].strip() == "CA":
                bf = float(line[60:66].strip())
                plddt_vals.append(bf * 100 if bf <= 1.0 else bf)

        if not plddt_vals:
            print("  NO CA ATOMS", flush=True)
            results.append({"name": name, "sequence": seq, "mean_plddt": -1, "status": "PARSE_FAILED"})
            del output, tokenized
            torch.cuda.empty_cache()
            continue

        mean_plddt = sum(plddt_vals) / len(plddt_vals)
        gt70_frac = sum(1 for p in plddt_vals if p > 70) / len(plddt_vals)
        status = "PASS" if mean_plddt >= 70 else ("MARGINAL" if mean_plddt >= 50 else "FAIL")

        # PROPKA on monomer
        pdb_path = f"/tmp/{name}.pdb"
        with open(pdb_path, "w") as f:
            f.write(pdb_text)

        his_pkas = {}
        try:
            from propka.run import single as propka_single
            mol = propka_single(pdb_path)
            best_conf = max(mol.conformations.values(), key=lambda c: len(c.groups))
            for group in best_conf.groups:
                if type(group).__name__ == "HISGroup":
                    resnum = group.atom.res_num
                    his_pkas[resnum] = round(group.pka_value, 2)
        except Exception as e:
            print(f"  PROPKA error: {e}", flush=True)

        avg_pka = sum(his_pkas.values()) / len(his_pkas) if his_pkas else 0

        r = {
            "name": name,
            "sequence": seq,
            "length": len(seq),
            "n_his": n_his,
            "mean_plddt": round(mean_plddt, 1),
            "plddt_gt70_frac": round(gt70_frac, 3),
            "status": status,
            "his_pkas": his_pkas,
            "avg_his_pka": round(avg_pka, 2),
        }
        results.append(r)
        print(f"  pLDDT={mean_plddt:.1f} ({status}), His pKa: {his_pkas}", flush=True)

        del output, tokenized
        torch.cuda.empty_cache()

    passed = [r for r in results if r["status"] == "PASS"]
    print(f"\n{'='*60}")
    print(f"SUMMARY: {len(passed)}/{len(results)} PASS (pLDDT≥70)")
    return results


@app.local_entrypoint()
def main(csv_path: str = "", tag: str = "acid"):
    if not csv_path:
        print("Usage: modal run modal_esmfold_c4_acid.py::main --csv <path> [--tag <tag>]")
        return

    seqs = []
    with open(csv_path) as f:
        for r in csv.DictReader(f):
            seqs.append({"name": r["name"], "sequence": r["sequence"]})
    print(f"Read {len(seqs)} sequences from {csv_path}")

    if not seqs:
        print("No sequences to refold.")
        return

    results = refold_batch.remote(seqs_json=json.dumps(seqs))

    out = f"/Users/ricardo/dev/adaptyv-c4/results/c4_{tag}_refold.json"
    with open(out, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\nSaved {len(results)} results to {out}")

    passed = [r for r in results if r["status"] == "PASS"]
    print(f"\nPASS: {len(passed)} designs with pLDDT≥70")
    for r in sorted(passed, key=lambda x: -x.get("avg_his_pka", 0)):
        print(f"  {r['name']}: pLDDT={r['mean_plddt']}, His={r['n_his']}, "
              f"avg_pKa={r.get('avg_his_pka', '?')}, pKas={r.get('his_pkas', {})}")
