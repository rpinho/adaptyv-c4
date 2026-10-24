"""ESMFold refold + monomer PROPKA for C4 acid-site RFD4 designs.

Reads binder (chain C) sequences from the RFD4 runs volume,
refolds with ESMFold, checks pLDDT, and runs PROPKA for free-state His pKa.

Usage:
    modal run modal_esmfold_c4_acid.py::refold_all --run-name tfr1_c4_acid1
    modal run modal_esmfold_c4_acid.py::refold_all --run-name tfr1_c4_acid2
"""
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

AA3TO1 = {
    'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F',
    'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L',
    'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R',
    'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y',
}


def extract_binder_sequences(run_name: str, min_his: int = 2):
    """Extract chain C binder sequences from RFD4 CIFs on the volume."""
    import glob
    import biotite.structure.io.pdbx as pdbx

    base = f"/runs/{run_name}/inference_samples"
    cifs = sorted(glob.glob(f"{base}/*/{run_name}_*_generated.cif"))
    print(f"Found {len(cifs)} CIF files in {base}")

    seqs = {}
    for cif_path in cifs:
        design_id = pathlib.Path(cif_path).parent.name
        try:
            f = pdbx.CIFFile.read(cif_path)
            arr = pdbx.get_structure(f.block, model=1)
            ca = arr[arr.atom_name == "CA"]
            chain_c = [a for a in ca if a.chain_id == "C"]
            if not chain_c:
                continue
            seq = "".join(AA3TO1.get(a.res_name, "X") for a in chain_c)
            n_his = seq.count("H")
            if n_his >= min_his:
                seqs[design_id] = seq
                print(f"  {design_id}: {len(seq)} aa, {n_his} His ✓")
            else:
                print(f"  {design_id}: {len(seq)} aa, {n_his} His (skip, <{min_his})")
        except Exception as e:
            print(f"  ERROR {design_id}: {e}")

    print(f"\nKept {len(seqs)}/{len(cifs)} designs with ≥{min_his} His")
    return seqs


@app.function(
    image=image,
    gpu="L40S",
    volumes={"/cache/huggingface": hf_cache, "/runs": runs_vol},
    timeout=3600,
)
def refold_all(run_name: str = "tfr1_c4_acid1", min_his: int = 2):
    import torch
    from transformers import AutoTokenizer, EsmForProteinFolding

    # Step 1: extract binder sequences
    seqs = extract_binder_sequences(run_name, min_his)
    if not seqs:
        print("No sequences with enough His. Exiting.")
        return {"results": [], "run_name": run_name}

    # Step 2: load ESMFold
    print("\nLoading ESMFold model...", flush=True)
    tokenizer = AutoTokenizer.from_pretrained("facebook/esmfold_v1")
    model = EsmForProteinFolding.from_pretrained(
        "facebook/esmfold_v1", low_cpu_mem_usage=True
    )
    model = model.cuda().eval()
    model.esm = model.esm.half()
    print("Model loaded.", flush=True)

    results = []
    for i, (name, seq) in enumerate(seqs.items()):
        print(f"\n[{i+1}/{len(seqs)}] {name} ({len(seq)} aa, {seq.count('H')} His)", flush=True)

        tokenized = tokenizer([seq], return_tensors="pt", add_special_tokens=False)
        tokenized = {k: v.cuda() for k, v in tokenized.items()}

        with torch.no_grad():
            output = model(**tokenized)

        pdb_text = model.output_to_pdb(output)[0]

        # Parse pLDDT from b-factors
        plddt_vals = []
        for line in pdb_text.splitlines():
            if line.startswith("ATOM") and line[12:16].strip() == "CA":
                bfactor = float(line[60:66].strip())
                plddt_vals.append(bfactor * 100 if bfactor <= 1.0 else bfactor)

        if not plddt_vals:
            print("  NO CA ATOMS", flush=True)
            results.append({"name": name, "sequence": seq, "mean_plddt": -1, "status": "PARSE_FAILED"})
            del output, tokenized
            torch.cuda.empty_cache()
            continue

        mean_plddt = sum(plddt_vals) / len(plddt_vals)
        gt70_frac = sum(1 for p in plddt_vals if p > 70) / len(plddt_vals)
        status = "PASS" if mean_plddt >= 70 else ("MARGINAL" if mean_plddt >= 50 else "FAIL")

        # Step 3: PROPKA on monomer (free-state His pKa)
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
        elevated = sum(1 for p in his_pkas.values() if p > 6.5)

        entry = {
            "name": name,
            "sequence": seq,
            "length": len(seq),
            "n_his": seq.count("H"),
            "mean_plddt": round(mean_plddt, 1),
            "plddt_gt70_frac": round(gt70_frac, 3),
            "status": status,
            "his_pkas": his_pkas,
            "avg_his_pka": round(avg_pka, 2),
            "elevated_his": elevated,
        }
        results.append(entry)
        print(f"  pLDDT={mean_plddt:.1f} ({status}), His pKa: {his_pkas}, "
              f"avg={avg_pka:.2f}, {elevated} elevated", flush=True)

        del output, tokenized
        torch.cuda.empty_cache()

    # Summary
    passed = [r for r in results if r["status"] == "PASS"]
    print(f"\n{'='*60}")
    print(f"SUMMARY: {len(passed)}/{len(results)} PASS (pLDDT≥70)")
    if passed:
        for r in sorted(passed, key=lambda x: -x["avg_his_pka"]):
            print(f"  {r['name']}: pLDDT={r['mean_plddt']}, "
                  f"His={r['n_his']}, avg_pKa={r['avg_his_pka']}, "
                  f"elevated={r['elevated_his']}")

    # Save results to volume
    out_path = f"/runs/{run_name}_refold_results.json"
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2)
    runs_vol.commit()
    print(f"\nSaved to {out_path}")

    return {"results": results, "run_name": run_name}


@app.local_entrypoint()
def main(run_name: str = "tfr1_c4_acid1", min_his: int = 2):
    result = refold_all.remote(run_name=run_name, min_his=min_his)
    out = f"/Users/ricardo/dev/adaptyv-c4/results/{run_name}_refold.json"
    with open(out, "w") as f:
        json.dump(result, f, indent=2)
    print(f"\nSaved {len(result.get('results', []))} results to {out}")

    passed = [r for r in result.get("results", []) if r["status"] == "PASS"]
    print(f"\n{'='*60}")
    print(f"PASS: {len(passed)} designs with pLDDT≥70")
    for r in sorted(passed, key=lambda x: -x.get("avg_his_pka", 0)):
        print(f"  {r['name']}: pLDDT={r['mean_plddt']}, His={r['n_his']}, "
              f"avg_pKa={r.get('avg_his_pka', '?')}")
