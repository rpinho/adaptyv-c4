"""Run PROPKA on ESMFold monomer PDBs for C4 v8 designs.

Reads PDB files from results/c4_v8_pdbs/, runs PROPKA on each,
reports His pKa values and Wyman-linkage ΔΔG estimates.

Usage:
    python3 -I scripts/propka_analysis.py
"""
import csv
import json
import math
import os
import sys
from pathlib import Path

sys.path.insert(0, "/Users/ricardo/Library/Python/3.9/lib/python/site-packages")

from propka.run import single

RESULTS = Path(__file__).resolve().parent.parent / "results"
CSV = RESULTS / "submission_v8_20261023.csv"
PDB_DIR = RESULTS / "c4_v8_pdbs"

RT = 0.592  # kcal/mol at 25°C


def analyze_pdb(name, pdb_path):
    try:
        mc = single(str(pdb_path))
    except Exception as e:
        return {"name": name, "error": str(e)}

    his_pkas = []
    for conf_key in mc.conformations:
        for group in mc.conformations[conf_key].groups:
            if type(group).__name__ == "HISGroup":
                his_pkas.append({
                    "res_num": group.atom.res_num,
                    "pka": round(group.pka_value, 2),
                    "model_pka": round(group.model_pka, 2),
                    "shift": round(group.pka_value - group.model_pka, 2),
                })
        break  # first conformation only

    if not his_pkas:
        return {"name": name, "his_count": 0, "ddg_monomer": 0}

    pkas = [h["pka"] for h in his_pkas]

    # Wyman-linkage monomer estimate
    ddg = 0
    for h in his_pkas:
        pka = h["pka"]
        ratio = (1 + 10**(pka - 5.5)) / (1 + 10**(pka - 7.4))
        ddg += -RT * math.log(ratio)

    # Surface vs buried classification
    surface = sum(1 for p in pkas if 5.0 <= p <= 7.5)
    buried_low = sum(1 for p in pkas if p < 5.0)
    buried_high = sum(1 for p in pkas if p > 7.5)

    # His with pKa in the switching window (5.5-7.4)
    switching = sum(1 for p in pkas if 5.5 <= p <= 7.4)

    return {
        "name": name,
        "his_count": len(his_pkas),
        "pka_mean": round(sum(pkas) / len(pkas), 2),
        "pka_min": round(min(pkas), 2),
        "pka_max": round(max(pkas), 2),
        "surface_his": surface,
        "buried_low": buried_low,
        "buried_high": buried_high,
        "switching_his": switching,
        "ddg_monomer": round(ddg, 2),
        "his_pkas": his_pkas,
    }


def main():
    if not PDB_DIR.exists():
        print(f"PDB directory not found: {PDB_DIR}")
        print("Run the Modal ESMFold script first, then download PDBs.")
        sys.exit(1)

    # Load design info
    designs = {}
    with open(CSV) as f:
        for r in csv.DictReader(f):
            designs[r["name"]] = r["sequence"]

    pdb_files = sorted(PDB_DIR.glob("*.pdb"))
    if not pdb_files:
        print(f"No PDB files in {PDB_DIR}")
        sys.exit(1)

    print(f"PROPKA analysis for {len(pdb_files)} C4 v8 PDB files")
    print("=" * 130)

    results = []
    for pdb_path in pdb_files:
        name = pdb_path.stem
        seq = designs.get(name, "")
        print(f"\nProcessing {name}...", end=" ", flush=True)
        r = analyze_pdb(name, pdb_path)
        results.append(r)

        if "error" in r:
            print(f"ERROR: {r['error']}")
            continue

        if r["his_count"] == 0:
            print("no His residues")
            continue

        print(f"{r['his_count']} His, pKa {r['pka_min']:.1f}-{r['pka_max']:.1f} "
              f"(mean {r['pka_mean']:.2f}), "
              f"switching={r['switching_his']}, ΔΔG_mono={r['ddg_monomer']:.2f} kcal/mol")

        for h in r["his_pkas"]:
            env = "SWITCHING" if 5.5 <= h["pka"] <= 7.4 else (
                "ALWAYS-ON" if h["pka"] > 7.4 else "LOW-PKA")
            print(f"  H{h['res_num']:3d}: pKa={h['pka']:.2f} "
                  f"(model={h['model_pka']:.2f}, shift={h['shift']:+.2f}) [{env}]")

    # Summary table
    print("\n\n" + "=" * 130)
    print("SUMMARY TABLE")
    print("=" * 130)
    print(f"{'Design':<40s} {'His':>4s} {'pKa_mean':>8s} {'pKa_range':>12s} "
          f"{'Switch':>6s} {'AlwOn':>5s} {'LowPk':>5s} {'ΔΔG_mono':>9s}")
    print("-" * 130)

    for r in sorted(results, key=lambda x: -x.get("ddg_monomer", 0)):
        if r.get("his_count", 0) == 0 or "error" in r:
            continue
        pka_range = f"{r['pka_min']:.1f}-{r['pka_max']:.1f}"
        always_on = r.get("buried_high", 0)
        low_pka = r.get("buried_low", 0)
        print(f"{r['name']:<40s} {r['his_count']:>4d} {r['pka_mean']:>8.2f} "
              f"{pka_range:>12s} {r['switching_his']:>6d} {always_on:>5d} "
              f"{low_pka:>5d} {r['ddg_monomer']:>9.2f}")

    # Save JSON
    out_path = RESULTS / "c4_propka_analysis.json"
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\nSaved to {out_path}")

    # Interpretation
    print("\n\nINTERPRETATION:")
    print("-" * 80)
    print("ΔΔG_monomer: Wyman-linkage estimate from MONOMER pKa values only.")
    print("  More negative = stronger pH-dependent charge change = better switch potential.")
    print("  This is a LOWER BOUND — actual ΔΔG_bind includes interface pKa shifts.")
    print()
    print("His categories:")
    print("  SWITCHING (pKa 5.5-7.4): ideal for pH switch — protonated at 5.5, neutral at 7.4")
    print("  ALWAYS-ON (pKa >7.4): protonated at BOTH pH values — no switching contribution")
    print("  LOW-PKA (<5.5): neutral at BOTH pH values — no switching contribution")
    print()
    print("For a full analysis, need COMPLEX structures (Boltz-2) to compare bound vs free pKa.")


if __name__ == "__main__":
    main()
