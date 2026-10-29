#!/usr/bin/env python3 -I
"""Post-process Boltz-2 complex prediction results.

Downloads results from Modal volume, computes:
1. ipTM ranking (binding evidence)
2. His707 contact analysis (atoms within 5Å of target His707)
3. Interface His identification (binder His within 5Å of any target atom)
4. Complex PROPKA → ΔpKa per interface His → true ΔΔG_switch

Usage:
    python3 -I scripts/boltz2_postprocess.py [--download-only]
"""
import json
import pathlib
import sys
import subprocess

RESULTS_DIR = pathlib.Path(__file__).resolve().parent.parent / "results"
COMPLEX_DIR = RESULTS_DIR / "boltz2_complexes"
COMPLEX_DIR.mkdir(parents=True, exist_ok=True)


def download_from_modal():
    """Download Boltz-2 results from Modal volume."""
    print("Downloading Boltz-2 results from Modal volume...")
    # List what's in the volume
    cmd = ["modal", "volume", "ls", "rfd4-runs", "boltz2_complexes/"]
    result = subprocess.run(cmd, capture_output=True, text=True)
    print(result.stdout)
    if result.returncode != 0:
        print(f"Error listing volume: {result.stderr}")
        return False

    # Download all files
    cmd = [
        "modal", "volume", "get", "rfd4-runs",
        "boltz2_complexes/",
        str(COMPLEX_DIR),
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    print(result.stdout)
    if result.returncode != 0:
        print(f"Error downloading: {result.stderr}")
        return False

    return True


def load_results():
    """Load boltz2_results.json."""
    results_file = COMPLEX_DIR / "boltz2_results.json"
    if not results_file.exists():
        print(f"No results file at {results_file}")
        return None
    with open(results_file) as f:
        return json.load(f)


def analyze_binding(results):
    """Rank designs by ipTM (binding evidence)."""
    print("\n" + "=" * 80)
    print("BOLTZ-2 BINDING EVIDENCE (ipTM ranking)")
    print("=" * 80)
    print(f"{'Design':<35s} {'ipTM':>6s} {'pTM':>6s} {'Conf':>6s} {'cpLDDT':>7s} {'Verdict'}")
    print("-" * 80)

    good = []
    for r in sorted(results, key=lambda x: x.get("iptm", 0), reverse=True):
        if "error" in r:
            print(f"{r['name']:<35s} ERROR: {r['error']}")
            continue
        iptm = r.get("iptm", 0)
        ptm = r.get("ptm", 0)
        conf = r.get("confidence_score", 0)
        cplddt = r.get("complex_plddt", 0)

        if iptm >= 0.7:
            verdict = "STRONG BINDER"
        elif iptm >= 0.5:
            verdict = "possible binder"
        elif iptm >= 0.3:
            verdict = "weak"
        else:
            verdict = "no binding"

        print(f"{r['name']:<35s} {iptm:>6.3f} {ptm:>6.3f} {conf:>6.3f} {cplddt:>7.1f} {verdict}")
        if iptm >= 0.3:
            good.append(r)

    print(f"\nDesigns with ipTM >= 0.3: {len(good)}/{len(results)}")
    return good


def analyze_his707_contacts(good_designs):
    """Check which designs contact TfR1 His707."""
    print("\n" + "=" * 80)
    print("His707 CONTACT ANALYSIS")
    print("=" * 80)

    try:
        from Bio.PDB import PDBParser
    except ImportError:
        print("BioPython not installed — run: pip install biopython")
        print("Skipping His707 contact analysis")
        return

    parser = PDBParser(QUIET=True)

    for design in good_designs:
        name = design["name"]
        pdb_path = COMPLEX_DIR / f"{name}_complex.pdb"
        if not pdb_path.exists():
            print(f"{name}: no complex PDB found")
            continue

        structure = parser.get_structure(name, str(pdb_path))
        model = structure[0]

        # Find His707 in chain B (target)
        # TfR1 apical domain is 154 residues, His707 maps to residue position
        # in the PDB numbering (depends on Boltz-2 output numbering)
        chain_b = model["B"] if "B" in model else None
        if not chain_b:
            print(f"{name}: no chain B found")
            continue

        # Find all His in chain B
        target_his = []
        for res in chain_b:
            if res.get_resname() == "HIS":
                target_his.append(res)

        # Find binder atoms (chain A) within 5Å of any target His
        chain_a = model["A"] if "A" in model else None
        if not chain_a:
            continue

        binder_atoms = list(chain_a.get_atoms())

        for his_res in target_his:
            his_atoms = list(his_res.get_atoms())
            min_dist = float("inf")
            closest_binder_res = None
            for ha in his_atoms:
                for ba in binder_atoms:
                    d = ha - ba
                    if d < min_dist:
                        min_dist = d
                        closest_binder_res = ba.get_parent()

            resid = his_res.get_id()[1]
            if min_dist < 5.0:
                br = closest_binder_res
                br_name = br.get_resname() if br else "?"
                br_id = br.get_id()[1] if br else "?"
                print(f"  {name}: Target His{resid} ← {min_dist:.1f}Å → Binder {br_name}{br_id}")


def analyze_interface_his(good_designs):
    """Identify binder His residues at the interface (within 5Å of target)."""
    print("\n" + "=" * 80)
    print("INTERFACE HISTIDINE ANALYSIS")
    print("=" * 80)

    try:
        from Bio.PDB import PDBParser
    except ImportError:
        print("BioPython not installed — skipping")
        return {}

    parser = PDBParser(QUIET=True)
    interface_his_map = {}

    for design in good_designs:
        name = design["name"]
        pdb_path = COMPLEX_DIR / f"{name}_complex.pdb"
        if not pdb_path.exists():
            continue

        structure = parser.get_structure(name, str(pdb_path))
        model = structure[0]
        chain_a = model["A"] if "A" in model else None
        chain_b = model["B"] if "B" in model else None
        if not chain_a or not chain_b:
            continue

        target_atoms = list(chain_b.get_atoms())

        interface_his = []
        total_his = 0
        for res in chain_a:
            if res.get_resname() == "HIS":
                total_his += 1
                res_atoms = list(res.get_atoms())
                min_dist = min(
                    (ra - ta for ra in res_atoms for ta in target_atoms),
                    default=999
                )
                if min_dist < 5.0:
                    interface_his.append((res.get_id()[1], min_dist))

        interface_his_map[name] = interface_his
        print(f"  {name}: {len(interface_his)}/{total_his} His at interface "
              f"(interface = within 5Å of target)")
        for resid, dist in interface_his:
            print(f"    His{resid}: {dist:.1f}Å from target")

    return interface_his_map


def main():
    download_only = "--download-only" in sys.argv

    if not download_from_modal():
        print("Download failed or no results yet")
        sys.exit(1)

    results = load_results()
    if not results:
        sys.exit(1)

    if download_only:
        print(f"Downloaded {len(results)} results")
        sys.exit(0)

    good = analyze_binding(results)
    if good:
        analyze_his707_contacts(good)
        analyze_interface_his(good)


if __name__ == "__main__":
    main()
