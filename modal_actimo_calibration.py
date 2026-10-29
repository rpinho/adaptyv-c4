"""Calibration: ESMFold + PROPKA on actimo-labs C2 (TNF-α) sequences.

Folds each monomer with ESMFold (HuggingFace), runs PROPKA on the structure,
reports His pKa. C2 designs place His near binder carboxylates → free-state pKa
should be ELEVATED (>6.5). This validates our PROPKA methodology against a
known-working competition entry.
"""
import json
import modal

app = modal.App("actimo-calibration")
hf_cache = modal.Volume.from_name("rfd4-hf-cache", create_if_missing=True)

ACTIMO_SEQS = {
    "tnf3_01": "MKKEIEEKKKEAEKLHEEAKKLLIKADELSHRYHELKEQAETEELKKEAEKAKEELTKVLKEHHEAVEKYEKLYAEAHELEHKL",
    "tnf3_02": "MKKEIEEKKKEAEKLHEEAKKLLIKAHELSERYHELKEQAETEELKKEAEKAKEELTKVLKEHHEAVEKYEKLYAEAHELEHKL",
    "tnf3_03": "KEEEERERHRLEDLRALAELYRETGEGHLAVARKAAEEGLRNPEKGEENLHIAEGEGAHAVAAYEAAAHAYREAGDEERAREAEEERERVREEVERLIREVREAL",
    "tnf3_04": "SDEAESLLHQAKAGADHLAGHTIARSHEAAETKGPEAARAAAAEGEAEMDALAAAHPEIAAYVAAQKAAARAAVDAHIAEHE",
    "tnf3_06": "LEKELEEVKEELKEVTEEVARFAERAAREVPEHAELLLHERDHVVAAAKHAHTASEIIQAAHHAAGLAAIAAHVASEELKEEAEKLKEKAEEIAKRAKELEEKIAK",
    "tnf3_08": "SEEIKRAEEEAKLNLEAAKARGEAAKHTFSEEIGEAAKKYNEALGKAVAAVHLAEAYKDDKEKTEELLEEAKKAGEEANKYHDEFMKLLKEEH",
    "tnf3_12": "SEEELLEILRKTLEHHVAAAEGLAEKLREAAERAAEIPELRPLAEEGLRVAERAREVARRAREVAEHGTAEEVREAAREVQAAIDELHAKHEEIHRGIEEAE",
    "tnf3_15": "KEEVERLLREAEELARRFEEKAHELAERSTPEGAHALDAAAHAARAAVEAARLGDRERAEHHIESARALLEKAREIEPELAEGIEEAEELFERLAELVRKALEL",
    "tnf3_17": "EVLELLERAIELAEEALRRAREEGLHEETIHALEEALRLARKALELARRALETGSEEDRIHAERYLGGALNHLATAAKRAEEEGSEETAELAKRAAEEVSKALEL",
    "tnf3_18": "SEEELLEILRKTLEHHVAAAEGLAEKLREAAERAAEIPELRPLAEEGLRVAERAREVARRAREVAEHGTAEEVREAAREVQAAIDELHAKHEEILHGIEEAE",
}

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
    volumes={"/cache/huggingface": hf_cache},
    timeout=1800,
)
def fold_and_score():
    import torch
    from transformers import AutoTokenizer, EsmForProteinFolding

    print("Loading ESMFold model...", flush=True)
    tokenizer = AutoTokenizer.from_pretrained("facebook/esmfold_v1")
    model = EsmForProteinFolding.from_pretrained(
        "facebook/esmfold_v1", low_cpu_mem_usage=True
    )
    model = model.cuda().eval()
    model.esm = model.esm.half()
    print("Model loaded.", flush=True)

    results = []
    for name, seq in ACTIMO_SEQS.items():
        print(f"\n{'='*60}", flush=True)
        print(f"Folding {name} ({len(seq)} aa, {seq.count('H')} His)", flush=True)

        tokenized = tokenizer([seq], return_tensors="pt", add_special_tokens=False)
        tokenized = {k: v.cuda() for k, v in tokenized.items()}

        with torch.no_grad():
            output = model(**tokenized)

        pdb_str = model.output_to_pdb(output)[0]

        # Parse pLDDT from CA b-factors
        plddt_vals = []
        for line in pdb_str.splitlines():
            if line.startswith("ATOM") and line[12:16].strip() == "CA":
                bf = float(line[60:66].strip())
                plddt_vals.append(bf * 100 if bf <= 1.0 else bf)
        plddt = sum(plddt_vals) / len(plddt_vals) if plddt_vals else 0

        pdb_path = f"/tmp/{name}.pdb"
        with open(pdb_path, "w") as f:
            f.write(pdb_str)

        print(f"  pLDDT = {plddt:.1f}", flush=True)

        # Run PROPKA
        from propka.run import single as propka_single
        mol = propka_single(pdb_path)

        best_conf = max(mol.conformations.values(), key=lambda c: len(c.groups))
        his_pkas = {}
        for group in best_conf.groups:
            if type(group).__name__ == "HISGroup":
                resnum = group.atom.res_num
                his_pkas[resnum] = group.pka_value
                print(f"  H{resnum}: pKa = {group.pka_value:.2f}", flush=True)

        avg_pka = sum(his_pkas.values()) / len(his_pkas) if his_pkas else 0
        elevated = sum(1 for p in his_pkas.values() if p > 6.5)

        entry = {
            "name": name,
            "length": len(seq),
            "n_his": seq.count('H'),
            "plddt": round(float(plddt), 1),
            "his_pkas": {str(k): round(v, 2) for k, v in sorted(his_pkas.items())},
            "avg_his_pka": round(avg_pka, 2),
            "elevated_count": elevated,
        }
        results.append(entry)
        print(f"  avg pKa = {avg_pka:.2f}, {elevated}/{len(his_pkas)} elevated", flush=True)

    return results


@app.local_entrypoint()
def main():
    results = fold_and_score.remote()
    out = "/Users/ricardo/dev/adaptyv-c4/results/actimo_calibration.json"
    with open(out, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\nSaved {len(results)} results to {out}")

    print("\n" + "=" * 80)
    print("ACTIMO-LABS FREE-STATE PROPKA CALIBRATION")
    print("=" * 80)
    print(f"{'Name':10s} {'Len':>4s} {'#H':>3s} {'pLDDT':>6s} {'Avg pKa':>8s} {'Elevated':>9s}")
    print("-" * 80)
    for r in results:
        print(f"{r['name']:10s} {r['length']:4d} {r['n_his']:3d} "
              f"{r['plddt']:6.1f} {r['avg_his_pka']:8.2f} {r['elevated_count']:>9d}")

    avg_all = sum(r["avg_his_pka"] for r in results) / len(results) if results else 0
    elev_all = sum(r["elevated_count"] for r in results)
    total_his = sum(r["n_his"] for r in results)
    print(f"\nOverall: avg His pKa = {avg_all:.2f}, {elev_all}/{total_his} elevated (>6.5)")
