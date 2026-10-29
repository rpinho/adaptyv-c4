#!/usr/bin/env python3 -I
"""Complex PROPKA analysis: compute true ΔΔG_switch from Boltz-2 complexes.

For each design:
1. Run PROPKA on the complex PDB → pKa_bound for each His
2. Extract binder chain → run PROPKA on free binder → pKa_free for each His
3. ΔpKa = pKa_bound - pKa_free per His
4. ΔΔG_switch = -RT * Σ_interface_His ln[(1+10^(pKa_b-5.5))/(1+10^(pKa_b-7.4))]
                + RT * Σ_interface_His ln[(1+10^(pKa_f-5.5))/(1+10^(pKa_f-7.4))]

Also identifies interface His (within 5Å of target) and His707 contacts.

Usage:
    python3 -I scripts/complex_propka_switch.py
"""
import json
import math
import pathlib
import sys

RESULTS_DIR = pathlib.Path(__file__).resolve().parent.parent / "results"
COMPLEX_DIR = RESULTS_DIR / "boltz2_complexes"

RT = 0.592  # kcal/mol at 25°C
PH_LOW = 5.5
PH_HIGH = 7.4


def run_propka_on_pdb(pdb_path):
    """Run PROPKA on a PDB file, return dict of {(chain, resnum): pka} for HIS."""
    try:
        from propka.run import single as propka_single
    except ImportError:
        print("PROPKA not installed: pip install propka")
        sys.exit(1)

    mol = propka_single(str(pdb_path))
    his_pkas = {}
    # Use the conformation with the most groups (multi-chain PDBs put HIS in "1A", not "AVR")
    best_conf = max(mol.conformations.values(), key=lambda c: len(c.groups))
    for group in best_conf.groups:
        if type(group).__name__ == "HISGroup":
            chain = group.atom.chain_id
            resnum = group.atom.res_num
            his_pkas[(chain, resnum)] = group.pka_value
    return his_pkas


def extract_chain(pdb_path, chain_id, out_path):
    """Extract a single chain from a PDB file."""
    with open(pdb_path) as f:
        lines = f.readlines()
    with open(out_path, "w") as f:
        for line in lines:
            if line.startswith(("ATOM", "HETATM")):
                if len(line) > 21 and line[21] == chain_id:
                    f.write(line)
            elif line.startswith("END"):
                f.write(line)


def find_interface_his(pdb_path, binder_chain="A", target_chain="B", cutoff=5.0):
    """Find binder His residues within cutoff of target atoms. Returns list of resnum."""
    try:
        from Bio.PDB import PDBParser
    except ImportError:
        print("BioPython not installed — using distance-free fallback")
        return None

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("complex", str(pdb_path))
    model = structure[0]

    if binder_chain not in model or target_chain not in model:
        return []

    target_atoms = list(model[target_chain].get_atoms())
    interface_his = []

    for res in model[binder_chain]:
        if res.get_resname() != "HIS":
            continue
        resnum = res.get_id()[1]
        min_dist = min(
            (ra - ta for ra in res.get_atoms() for ta in target_atoms),
            default=999
        )
        if min_dist < cutoff:
            interface_his.append(resnum)

    return interface_his


def find_his707_contacts(pdb_path, binder_chain="A", target_chain="B", cutoff=5.0):
    """Find binder residues within cutoff of target His residues (His707 proxy)."""
    try:
        from Bio.PDB import PDBParser
    except ImportError:
        return []

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("complex", str(pdb_path))
    model = structure[0]

    if target_chain not in model:
        return []

    contacts = []
    for target_res in model[target_chain]:
        if target_res.get_resname() != "HIS":
            continue
        target_his_num = target_res.get_id()[1]
        target_atoms = list(target_res.get_atoms())

        for binder_res in model[binder_chain]:
            binder_atoms = list(binder_res.get_atoms())
            min_dist = min(
                (ba - ta for ba in binder_atoms for ta in target_atoms),
                default=999
            )
            if min_dist < cutoff:
                contacts.append({
                    "target_his": target_his_num,
                    "binder_res": binder_res.get_id()[1],
                    "binder_resname": binder_res.get_resname(),
                    "distance": round(min_dist, 1),
                })

    return contacts


def wyman_linkage_ddg(pka, ph_low=PH_LOW, ph_high=PH_HIGH):
    """Wyman-linkage ΔΔG contribution for a single His at given pKa."""
    return -RT * math.log(
        (1 + 10 ** (pka - ph_low)) / (1 + 10 ** (pka - ph_high))
    )


def analyze_design(name, complex_pdb):
    """Full analysis for one design."""
    import tempfile

    result = {"name": name}

    # 1. PROPKA on complex
    complex_pkas = run_propka_on_pdb(complex_pdb)
    binder_his_bound = {r: pka for (c, r), pka in complex_pkas.items() if c == "A"}
    target_his_bound = {r: pka for (c, r), pka in complex_pkas.items() if c == "B"}

    # 2. Extract free binder, PROPKA on it
    with tempfile.NamedTemporaryFile(suffix=".pdb", delete=False) as tmp:
        free_binder_path = tmp.name
    extract_chain(complex_pdb, "A", free_binder_path)
    free_pkas = run_propka_on_pdb(free_binder_path)
    binder_his_free = {r: pka for (c, r), pka in free_pkas.items() if c == "A"}

    # 3. Interface His (None = BioPython missing; [] = no interface His found)
    interface_his_result = find_interface_his(complex_pdb)
    if interface_his_result is None:
        interface_his = list(binder_his_bound.keys())
    else:
        interface_his = interface_his_result

    # 4. His707 contacts
    his707_contacts = find_his707_contacts(complex_pdb)

    # 5. Compute ΔΔG_switch for interface His only
    his_analysis = []
    total_ddg = 0.0
    common_his = sorted(set(binder_his_bound.keys()) & set(binder_his_free.keys()))
    for resnum in common_his:
        pka_b = binder_his_bound[resnum]
        pka_f = binder_his_free[resnum]
        dpka = pka_b - pka_f
        is_interface = resnum in interface_his
        ddg_bound = wyman_linkage_ddg(pka_b)
        ddg_free = wyman_linkage_ddg(pka_f)
        ddg_switch = ddg_bound - ddg_free  # contribution to binding switch

        entry = {
            "resnum": resnum,
            "pka_bound": round(pka_b, 2),
            "pka_free": round(pka_f, 2),
            "dpka": round(dpka, 2),
            "ddg_switch": round(ddg_switch, 3),
            "interface": is_interface,
        }
        his_analysis.append(entry)

        if is_interface:
            total_ddg += ddg_switch

    result["his_analysis"] = his_analysis
    result["interface_his_count"] = sum(1 for h in his_analysis if h["interface"])
    result["total_his_count"] = len(his_analysis)
    # C4 is acid-OFF (bind at pH 7.4, release at pH 6.0). Positive = pKa DEPRESSED on binding
    # = the correct direction. Verified numerically: a His going 7.0 -> 5.0 on binding gives
    # +1.70; a His going 6.5 -> 8.0 gives -1.11.
    result["ddg_switch_interface"] = round(total_ddg, 3)
    result["his707_contacts"] = his707_contacts
    result["target_his_bound"] = {str(k): round(v, 2) for k, v in target_his_bound.items()}

    return result


def main():
    # Load Boltz-2 results for ipTM scores
    boltz_results_file = COMPLEX_DIR / "boltz2_results.json"
    if boltz_results_file.exists():
        with open(boltz_results_file) as f:
            boltz_results = {r["name"]: r for r in json.load(f) if "error" not in r}
    else:
        boltz_results = {}

    # Find all complex PDBs
    complex_pdbs = sorted(COMPLEX_DIR.glob("*_complex.pdb"))
    if not complex_pdbs:
        print(f"No complex PDBs found in {COMPLEX_DIR}")
        sys.exit(1)

    print(f"Analyzing {len(complex_pdbs)} complex structures")
    print("=" * 100)

    all_results = []
    for pdb_path in complex_pdbs:
        name = pdb_path.stem.replace("_complex", "")
        print(f"\n{'='*80}")
        print(f"Design: {name}")
        boltz = boltz_results.get(name, {})
        iptm = boltz.get("iptm", "N/A")
        print(f"  ipTM: {iptm}")

        try:
            result = analyze_design(name, pdb_path)
            result["iptm"] = iptm
            all_results.append(result)

            print(f"  Interface His: {result['interface_his_count']}/{result['total_his_count']}")
            print(f"  ΔΔG_switch (interface): {result['ddg_switch_interface']:.3f} kcal/mol")
            if result["his707_contacts"]:
                print(f"  His707 contacts: {len(result['his707_contacts'])}")
                for c in result["his707_contacts"]:
                    print(f"    Target His{c['target_his']} ← {c['distance']}Å → "
                          f"Binder {c['binder_resname']}{c['binder_res']}")
            else:
                print(f"  His707 contacts: NONE")

            print(f"\n  Per-His analysis:")
            print(f"  {'Res':>5s} {'pKa_b':>7s} {'pKa_f':>7s} {'ΔpKa':>7s} {'ΔΔG':>8s} {'Iface':>6s}")
            for h in result["his_analysis"]:
                iface = "YES" if h["interface"] else "no"
                print(f"  H{h['resnum']:>4d} {h['pka_bound']:>7.2f} {h['pka_free']:>7.2f} "
                      f"{h['dpka']:>7.2f} {h['ddg_switch']:>8.3f} {iface:>6s}")

        except Exception as e:
            print(f"  ERROR: {e}")
            all_results.append({"name": name, "error": str(e)})

    # Summary ranking
    print(f"\n\n{'='*100}")
    print("COMPLEX PROPKA RANKING — C4 ΔΔG_switch (positive = pKa depressed on binding = acid-OFF = correct)")
    print(f"{'='*100}")
    print(f"{'Design':<40s} {'ipTM':>6s} {'iHis':>5s} {'ΔΔG_sw':>8s} {'His707':>7s}")
    print("-" * 100)

    ranked = sorted(
        [r for r in all_results if "error" not in r],
        key=lambda x: x.get("ddg_switch_interface", 0),
        reverse=True,
    )
    for r in ranked:
        iptm_str = f"{r['iptm']:.3f}" if isinstance(r["iptm"], float) else str(r["iptm"])
        his707 = len(r.get("his707_contacts", []))
        print(f"{r['name']:<40s} {iptm_str:>6s} {r['interface_his_count']:>5d} "
              f"{r['ddg_switch_interface']:>8.3f} {his707:>7d}")

    # Save results
    out_file = RESULTS_DIR / "c4_complex_propka_results.json"
    with open(out_file, "w") as f:
        json.dump(all_results, f, indent=2, default=str)
    print(f"\nSaved to {out_file}")


if __name__ == "__main__":
    main()
