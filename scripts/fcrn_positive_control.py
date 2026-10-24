"""FcRn/IgG positive control for C4 PROPKA Wyman-linkage pipeline.

FcRn–IgG Fc is a well-characterized pH-dependent interaction:
  - Binds at pH ~6.0 (endosomal), releases at pH 7.4 (physiological)
  - Same direction as C4: acid-preferential binding
  - Key His residues on IgG Fc: His310, His435, His436 (EU numbering)
  - These His form salt bridges with FcRn Glu at low pH

PDB structures:
  - 1I1A: Human FcRn + IgG1 Fc heterodimer (2.8Å)
  - 4N0U: High-resolution FcRn–Fc complex (2.6Å)
  - 1OQO: Free IgG1 Fc fragment (2.4Å, no FcRn)

We run PROPKA on:
  1. Free Fc (1OQO) → baseline His pKa
  2. FcRn–Fc complex (1I1A) → bound His pKa
Then compute Wyman linkage L from the pKa shift.

Expected: L > 0 (proton uptake on binding = acid-preferential).
"""
import math
import os
import sys
import tempfile
import urllib.request

sys.path.insert(0, "/Users/ricardo/Library/Python/3.9/lib/python/site-packages")

RT = 0.592       # kcal/mol at 25°C
PH_LOW = 5.5
PH_HIGH = 7.4

# Known pH-switch His residues on IgG Fc (EU numbering)
# In PDB structures the numbering may differ by chain/offset
KNOWN_SWITCH_HIS = {310, 435, 436}


def download_pdb(pdb_id: str, out_dir: str) -> str:
    url = f"https://files.rcsb.org/download/{pdb_id.upper()}.pdb"
    out_path = os.path.join(out_dir, f"{pdb_id.lower()}.pdb")
    if os.path.exists(out_path):
        print(f"  Already downloaded: {out_path}")
        return out_path
    print(f"  Downloading {url} ...")
    urllib.request.urlretrieve(url, out_path)
    print(f"  Saved to {out_path}")
    return out_path


def run_propka(pdb_path: str) -> dict:
    """Run PROPKA and return {chain_resnum: pka} for all HIS groups."""
    from propka.run import single as propka_single

    mol = propka_single(pdb_path)
    best_conf = max(mol.conformations.values(), key=lambda c: len(c.groups))

    his_pkas = {}
    for group in best_conf.groups:
        if type(group).__name__ == "HISGroup":
            chain = group.atom.chain_id
            resnum = group.atom.res_num
            key = f"{chain}{resnum}"
            his_pkas[key] = round(group.pka_value, 2)
    return his_pkas


def wyman_L(free_pkas: dict, bound_pkas: dict):
    """Compute Wyman linkage from free and bound His pKa values.

    Returns (L, total_ddg, per_his_details).
    Positive L = proton uptake on binding = acid-preferential.
    """
    total_ddg = 0.0
    per_his = {}

    for key in bound_pkas:
        pka_bound = bound_pkas[key]
        pka_free = free_pkas.get(key, 6.04)  # model His pKa if no free match

        ddg_bound = -RT * math.log(
            (1 + 10**(pka_bound - PH_LOW)) / (1 + 10**(pka_bound - PH_HIGH))
        )
        ddg_free = -RT * math.log(
            (1 + 10**(pka_free - PH_LOW)) / (1 + 10**(pka_free - PH_HIGH))
        )
        ddg = ddg_bound - ddg_free
        total_ddg += ddg
        per_his[key] = {
            "pka_free": pka_free,
            "pka_bound": pka_bound,
            "delta_pka": round(pka_bound - pka_free, 2),
            "ddg": round(ddg, 3),
        }

    L = total_ddg / (RT * math.log(10)) if total_ddg != 0 else 0
    return round(L, 2), round(total_ddg, 3), per_his


def identify_fc_chains(pdb_path: str) -> list:
    """Identify Fc chains by looking for known Fc His residues."""
    chains_with_his = {}
    with open(pdb_path) as f:
        for line in f:
            if line.startswith("ATOM") and line[17:20].strip() == "HIS":
                chain = line[21]
                resnum = int(line[22:26].strip())
                chains_with_his.setdefault(chain, set()).add(resnum)

    # Report what we find
    print(f"\n  Chains with HIS residues:")
    for chain in sorted(chains_with_his):
        resnums = sorted(chains_with_his[chain])
        print(f"    Chain {chain}: HIS at positions {resnums}")

    return chains_with_his


def main():
    out_dir = "/private/tmp/claude-501/-Users-ricardo--worktrees-personal-kb-claude-adaptyv-c4-v2/e08d870a-a286-4846-bb11-f5d7bbc95128/scratchpad/fcrn_control"
    os.makedirs(out_dir, exist_ok=True)

    # ----- 1. Download structures -----
    print("=" * 60)
    print("FcRn/IgG Positive Control — PROPKA Wyman Linkage Calibration")
    print("=" * 60)

    print("\n1. Downloading PDB structures...")
    complex_pdb = download_pdb("1I1A", out_dir)  # FcRn–Fc complex
    free_fc_pdb = download_pdb("1OQO", out_dir)   # Free IgG1 Fc

    # ----- 2. Identify chains -----
    print("\n2. Identifying chains in complex (1I1A)...")
    complex_chains = identify_fc_chains(complex_pdb)

    print("\n   Identifying chains in free Fc (1OQO)...")
    free_chains = identify_fc_chains(free_fc_pdb)

    # ----- 3. Run PROPKA -----
    print("\n3. Running PROPKA on complex (1I1A: FcRn + Fc)...")
    complex_pkas = run_propka(complex_pdb)
    print(f"   All His pKas in complex: {complex_pkas}")

    print("\n   Running PROPKA on free Fc (1OQO)...")
    free_pkas = run_propka(free_fc_pdb)
    print(f"   All His pKas in free Fc: {free_pkas}")

    # ----- 4. Match Fc His between free and complex -----
    # 1I1A: chains A,B = FcRn α-chain and β2m; chains C,D = Fc (two heavy chains)
    # 1OQO: chains A,B = Fc (two heavy chains)
    # We need to map free Fc His to complex Fc His

    # Find which chains are Fc by His count and known positions
    print("\n4. Mapping Fc His residues between free and complex...")

    # In complex 1I1A, Fc chains are typically C and D (or E,F depending on structure)
    # Let's use all His from all chains and report per-chain L
    print("\n   Computing Wyman linkage for ALL His across all chains:")

    L_all, ddg_all, per_his_all = wyman_L(free_pkas, complex_pkas)
    print(f"\n   ALL chains combined:")
    print(f"     L = {L_all:.2f}  (ΔΔG = {ddg_all:.3f} kcal/mol)")
    for k in sorted(per_his_all, key=lambda x: -abs(per_his_all[x]["delta_pka"])):
        h = per_his_all[k]
        print(f"     {k}: free={h['pka_free']:.2f} → bound={h['pka_bound']:.2f} "
              f"(Δ={h['delta_pka']:+.2f}, ΔΔG={h['ddg']:+.3f})")

    # ----- 5. Per-chain analysis -----
    print("\n5. Per-chain Wyman linkage:")
    chain_ids = sorted(set(k[0] for k in complex_pkas.keys()))
    for cid in chain_ids:
        chain_complex = {k: v for k, v in complex_pkas.items() if k[0] == cid}
        if not chain_complex:
            continue
        # Try to match free-state His by residue number
        # For Fc chains (C,D in 1I1A), match to free Fc chains (A,B in 1OQO)
        # Try same chain letter first, then try matching by residue number
        chain_free = {k: v for k, v in free_pkas.items() if k[0] == cid}
        if not chain_free:
            # Try matching by resnumber only (strip chain ID)
            resnums_free = {k[1:]: v for k, v in free_pkas.items()}
            chain_free_mapped = {}
            for k, v in chain_complex.items():
                resnum = k[1:]
                if resnum in resnums_free:
                    chain_free_mapped[k] = resnums_free[resnum]
            chain_free = chain_free_mapped

        L_chain, ddg_chain, per_his_chain = wyman_L(chain_free, chain_complex)
        his_str = ", ".join(f"{k}={v:.2f}" for k, v in sorted(chain_complex.items()))
        print(f"   Chain {cid}: L = {L_chain:+.2f} (ΔΔG = {ddg_chain:+.3f})")
        print(f"     Complex pKas: {his_str}")
        if chain_free:
            free_str = ", ".join(f"{k}={v:.2f}" for k, v in sorted(chain_free.items()))
            print(f"     Free pKas:    {free_str}")

    # ----- 6. Correct chain mapping: complex C,D = Fc ↔ free A,B = Fc -----
    print("\n6. Corrected chain mapping (complex C,D = Fc ↔ free A,B = Fc):")
    # 1I1A: A = FcRn α-chain, B = β2m, C = Fc heavy chain 1, D = Fc heavy chain 2
    # 1OQO: A = Fc heavy chain 1, B = Fc heavy chain 2, C,D = small peptides
    fc_complex_chains = ["C", "D"]
    fc_free_chains = ["A", "B"]

    # Map by residue number: complex chain C → free chain A, complex chain D → free chain B
    chain_map = {"C": "A", "D": "B"}
    fc_free_mapped = {}
    fc_complex_only = {}
    for cc, fc in chain_map.items():
        for k, v in complex_pkas.items():
            if k[0] != cc:
                continue
            resnum = k[1:]
            free_key = f"{fc}{resnum}"
            bound_key = f"{cc}{resnum}"
            fc_complex_only[bound_key] = v
            if free_key in free_pkas:
                fc_free_mapped[bound_key] = free_pkas[free_key]

    print(f"  Complex Fc His: {fc_complex_only}")
    print(f"  Matched free:   {fc_free_mapped}")

    L_fc, ddg_fc, per_his_fc = wyman_L(fc_free_mapped, fc_complex_only)
    print(f"\n  Fc-only Wyman L = {L_fc:+.2f} (ΔΔG = {ddg_fc:+.3f} kcal/mol)")
    for k in sorted(per_his_fc, key=lambda x: -abs(per_his_fc[x]["delta_pka"])):
        h = per_his_fc[k]
        is_switch = int(k[1:]) in KNOWN_SWITCH_HIS
        marker = " *** KNOWN SWITCH" if is_switch else ""
        print(f"    {k} (His{k[1:]}): free={h['pka_free']:.2f} → bound={h['pka_bound']:.2f} "
              f"(Δ={h['delta_pka']:+.2f}, ΔΔG={h['ddg']:+.3f}){marker}")

    # ----- 7. Focus on known switch His -----
    print("\n7. Known pH-switch His (EU: 310, 435, 436) on Fc chains C,D:")
    for cid in fc_complex_chains:
        chain_pkas = {k: v for k, v in complex_pkas.items() if k[0] == cid}
        for k, v in sorted(chain_pkas.items()):
            resnum = int(k[1:])
            if resnum not in KNOWN_SWITCH_HIS:
                continue
            free_chain = chain_map[cid]
            free_key = f"{free_chain}{resnum}"
            free_match = free_pkas.get(free_key, None)
            delta = (v - free_match) if free_match is not None else None
            delta_str = f"Δ={delta:+.2f}" if delta is not None else "no free match"
            free_str = f"{free_match:.2f}" if free_match is not None else "?"
            print(f"   {k} (His{resnum}): complex={v:.2f}, free={free_str} ({delta_str})")

    # ----- 8. Summary -----
    print("\n" + "=" * 60)
    print("SUMMARY — FcRn/IgG Positive Control")
    print("=" * 60)
    print(f"Complex: 1I1A (FcRn–Fc), Free: 1OQO (Fc alone)")
    print(f"Fc-only Wyman L = {L_fc:+.2f} (ΔΔG = {ddg_fc:+.3f} kcal/mol)")
    print(f"(All-chains L = {L_all:+.2f}, but includes FcRn His which aren't relevant)")
    n_positive = sum(1 for h in per_his_fc.values() if h["delta_pka"] > 0)
    n_negative = sum(1 for h in per_his_fc.values() if h["delta_pka"] < 0)
    print(f"Fc His with positive ΔpKa (proton uptake): {n_positive}")
    print(f"Fc His with negative ΔpKa (proton release): {n_negative}")
    print(f"\nInterpretation:")
    if L_fc > 0:
        print(f"  L > 0 confirms acid-preferential binding (correct C4 direction)")
        print(f"  This calibrates our pipeline: known switch gives L = {L_fc:+.2f}")
    elif L_fc < 0:
        print(f"  L < 0 — opposite of expected for FcRn/IgG acid-preferential binding.")
        print(f"  Possible reasons:")
        print(f"    - PROPKA errors on crystal structures with missing atoms")
        print(f"    - His pKa shifts dominated by non-interface His")
        print(f"    - Known switch His may have unusual PROPKA behavior")
        print(f"  Focus on known switch His (310, 435, 436) for targeted calibration.")
    print(f"\nFor C4 gate calibration:")
    print(f"  Raw L from full Fc is noisy with non-interface His")
    print(f"  Better metric: L computed on interface His only (within 5Å of FcRn)")


if __name__ == "__main__":
    main()
