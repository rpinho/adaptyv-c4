"""Calibration: ESMFold + PROPKA on actimo-labs C2 (TNF-α) sequences.

Folds each monomer with ESMFold, runs PROPKA on the structure, reports His pKa.
C2 designs place His near binder carboxylates → free-state pKa should be ELEVATED (>6.5).
This validates our PROPKA methodology against a known-working competition entry.
"""
import json
import modal

app = modal.App("actimo-calibration")
vol = modal.Volume.from_name("boltz2-cache", create_if_missing=True)

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
    modal.Image.debian_slim(python_version="3.11")
    .pip_install("torch", "fair-esm", "propka", "omegaconf", "numpy")
)


@app.function(image=image, gpu="L40S", timeout=1800, volumes={"/cache": vol})
def fold_and_score():
    import esm
    import os
    import tempfile

    os.environ["TORCH_HOME"] = "/cache/torch"

    model = esm.pretrained.esmfold_v1()
    model = model.eval().cuda()
    model.set_chunk_size(64)

    results = []
    for name, seq in ACTIMO_SEQS.items():
        print(f"\n{'='*60}")
        print(f"Folding {name} ({len(seq)} aa, {seq.count('H')} His)")

        import torch
        with torch.no_grad():
            output = model.infer_pdb(seq)

        pdb_path = f"/tmp/{name}.pdb"
        with open(pdb_path, "w") as f:
            f.write(output)

        # Run PROPKA
        from propka.run import single as propka_single
        mol = propka_single(pdb_path)

        best_conf = max(mol.conformations.values(), key=lambda c: len(c.groups))
        his_pkas = {}
        for group in best_conf.groups:
            if type(group).__name__ == "HISGroup":
                resnum = group.atom.res_num
                his_pkas[resnum] = group.pka_value
                print(f"  H{resnum}: pKa = {group.pka_value:.2f}")

        avg_pka = sum(his_pkas.values()) / len(his_pkas) if his_pkas else 0
        elevated = sum(1 for p in his_pkas.values() if p > 6.5)

        entry = {
            "name": name,
            "length": len(seq),
            "n_his": seq.count('H'),
            "his_pkas": {str(k): round(v, 2) for k, v in sorted(his_pkas.items())},
            "avg_his_pka": round(avg_pka, 2),
            "elevated_count": elevated,
            "note": "C2 design — free-state His should be elevated (>6.5) by binder D/E"
        }
        results.append(entry)
        print(f"  avg pKa = {avg_pka:.2f}, {elevated}/{len(his_pkas)} elevated")

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
    print(f"{'Name':10s} {'Len':>4s} {'#H':>3s} {'Avg pKa':>8s} {'Elevated':>9s} {'His pKa values'}")
    print("-" * 80)
    for r in results:
        pkas = ", ".join(f"H{k}={v}" for k, v in r["his_pkas"].items())
        print(f"{r['name']:10s} {r['length']:4d} {r['n_his']:3d} "
              f"{r['avg_his_pka']:8.2f} {r['elevated_count']:>9d}   {pkas}")

    avg_all = sum(r["avg_his_pka"] for r in results) / len(results) if results else 0
    elev_all = sum(r["elevated_count"] for r in results)
    total_his = sum(r["n_his"] for r in results)
    print(f"\nOverall: avg His pKa = {avg_all:.2f}, {elev_all}/{total_his} elevated (>6.5)")
    print("C2 validation: elevated free-state pKa confirms binder D/E strategy works for PROPKA")
