"""Score C4 acid-site designs: combine ESMFold refold + Boltz-2 complex + PROPKA.

Computes:
  - ΔpKa = complex_pKa - free_pKa (per His)
  - Wyman linkage L = ΔΔG / (RT * ln10) at pH 5.5/7.4
  - Gate: L ≥ +2

Usage:
    python3 -I scripts/score_c4_acid.py \
        --refold results/c4_acid_combined_refold.json \
        --boltz results/c4_acid_boltz2.json
"""
import argparse
import json
import math

RT = 0.592  # kcal/mol at 25°C
PH_LOW = 5.5
PH_HIGH = 7.4


def wyman_L(his_pkas_free: dict, his_pkas_complex: dict):
    """Compute Wyman-linkage ΔΔG from free and complex His pKa values.

    For C4: positive L = more binding at pH 5.5 = CORRECT direction.
    """
    total_ddg = 0.0
    per_his = {}

    for key in his_pkas_complex:
        pka_bound = his_pkas_complex[key]
        pka_free = his_pkas_free.get(key, 6.0)  # default model pKa

        # ΔΔG = -RT * ln[(1 + 10^(pKa_bound - pH_low)) / (1 + 10^(pKa_bound - pH_high))]
        #      + RT * ln[(1 + 10^(pKa_free - pH_low)) / (1 + 10^(pKa_free - pH_high))]
        ddg_bound = -RT * math.log(
            (1 + 10**(pka_bound - PH_LOW)) / (1 + 10**(pka_bound - PH_HIGH))
        )
        ddg_free = -RT * math.log(
            (1 + 10**(pka_free - PH_LOW)) / (1 + 10**(pka_free - PH_HIGH))
        )
        ddg = ddg_bound - ddg_free
        total_ddg += ddg
        per_his[key] = {
            "pka_free": round(pka_free, 2),
            "pka_bound": round(pka_bound, 2),
            "delta_pka": round(pka_bound - pka_free, 2),
            "ddg": round(ddg, 3),
        }

    L = total_ddg / (RT * math.log(10)) if total_ddg != 0 else 0
    return round(L, 2), round(total_ddg, 3), per_his


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--refold", required=True)
    p.add_argument("--boltz", required=True)
    p.add_argument("--out", default="")
    p.add_argument("--gate", type=float, default=2.0)
    args = p.parse_args()

    with open(args.refold) as f:
        refold_data = json.load(f)
    if isinstance(refold_data, dict):
        refold_data = refold_data.get("results", [])

    with open(args.boltz) as f:
        boltz_data = json.load(f)

    refold_map = {r["name"]: r for r in refold_data}

    scored = []
    for b in boltz_data:
        name = b["name"]
        if "error" in b:
            continue

        r = refold_map.get(name, {})
        free_pkas = r.get("his_pkas", {})
        complex_pkas = b.get("his_complex_pkas", {})

        # Map free pKa keys (resnum) to complex keys (chain+resnum)
        # Free state is monomer → chain A in complex
        free_mapped = {}
        for k, v in free_pkas.items():
            free_mapped[f"A{k}"] = v

        # Only score binder His (chain A)
        binder_complex = {k: v for k, v in complex_pkas.items() if k.startswith("A")}

        L, total_ddg, per_his = wyman_L(free_mapped, binder_complex)

        entry = {
            "name": name,
            "sequence": b.get("sequence", r.get("sequence", "")),
            "length": len(b.get("sequence", r.get("sequence", ""))),
            "n_his": b.get("sequence", "").count("H"),
            "mean_plddt": r.get("mean_plddt", -1),
            "iptm": b.get("iptm", 0),
            "ptm": b.get("ptm", 0),
            "confidence": b.get("confidence", 0),
            "avg_free_pka": r.get("avg_his_pka", 0),
            "avg_complex_pka": b.get("avg_complex_his_pka", 0),
            "L": L,
            "total_ddg": total_ddg,
            "switch_ok": L >= args.gate,
            "per_his": per_his,
        }
        scored.append(entry)

    scored.sort(key=lambda x: -x["L"])

    print(f"{'Name':<40s} {'Len':>4s} {'#H':>3s} {'pLDDT':>6s} {'ipTM':>6s} "
          f"{'L':>6s} {'ΔΔG':>7s} {'Gate':>5s}")
    print("-" * 90)

    passing = []
    for s in scored:
        gate = "✓" if s["switch_ok"] else "✗"
        print(f"{s['name']:<40s} {s['length']:4d} {s['n_his']:3d} "
              f"{s['mean_plddt']:6.1f} {s['iptm']:6.3f} "
              f"{s['L']:6.2f} {s['total_ddg']:7.3f} {gate:>5s}")
        if s["switch_ok"]:
            passing.append(s)

    print(f"\n{len(passing)}/{len(scored)} pass L ≥ {args.gate} gate")

    if args.out:
        with open(args.out, "w") as f:
            json.dump(scored, f, indent=2)
        print(f"Saved to {args.out}")

    return scored


if __name__ == "__main__":
    main()
