"""Score C4 designs by geometric His–carboxylate contacts in Boltz-2 complex PDBs.

More reliable than PROPKA-based Wyman linkage for pH-switch ranking.
Counts binder His within salt-bridge distance (4.0Å) of target Asp/Glu.

Usage:
    python3 scripts/score_geometric_contacts.py \
        --pdb-dir results/boltz2_complexes/ \
        --scored results/c4_acid_scored_combined.json \
        --out results/c4_acid_geom_scored.json
"""
import argparse
import json
import os
import math

# Target TfR1 acidic residues (Boltz-2 chain B numbering)
# D538→B63, E561→B86, E573→B98, E593→B118, E607→B132
TARGET_ACIDIC_RESIDS = {63, 86, 98, 118, 132}
CUTOFF = 4.0
CUTOFF2 = CUTOFF * CUTOFF

HIS_ATOMS = ("ND1", "NE2")
ACID_ATOMS = ("OD1", "OD2", "OE1", "OE2")
ACID_RESNAMES = ("ASP", "GLU")


def count_contacts(pdb_path, binder_chain="A", target_chain="B",
                   target_resids=None, cutoff=4.0):
    from Bio.PDB import PDBParser
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("cx", pdb_path)
    model = next(structure.get_models())

    if target_chain not in model:
        return 0, []

    tgt_oxygens = []
    for res in model[target_chain]:
        if res.id[0] != " ":
            continue
        if target_resids and res.id[1] not in target_resids:
            continue
        if res.get_resname().upper() in ACID_RESNAMES:
            for aname in ACID_ATOMS:
                if aname in res:
                    tgt_oxygens.append((res.id[1], res.get_resname(), aname,
                                       res[aname].get_vector().get_array()))

    if not tgt_oxygens:
        return 0, []

    c2 = cutoff * cutoff
    contacts = []

    for res in model[binder_chain]:
        if res.id[0] != " " or res.get_resname().upper() != "HIS":
            continue
        resnum = res.id[1]
        for haname in HIS_ATOMS:
            if haname not in res:
                continue
            hcoord = res[haname].get_vector().get_array()
            for tresnum, tresname, taname, tcoord in tgt_oxygens:
                d2 = sum((a - b) ** 2 for a, b in zip(hcoord, tcoord))
                if d2 <= c2:
                    contacts.append({
                        "his_resnum": resnum,
                        "his_atom": haname,
                        "target_resnum": tresnum,
                        "target_resname": tresname,
                        "target_atom": taname,
                        "distance": round(math.sqrt(d2), 2),
                    })

    unique_his = len(set(c["his_resnum"] for c in contacts))
    return unique_his, contacts


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--pdb-dir", required=True)
    p.add_argument("--scored", required=True, help="Wyman-scored JSON for metadata")
    p.add_argument("--out", default="results/c4_acid_geom_scored.json")
    p.add_argument("--iptm-gate", type=float, default=0.6)
    args = p.parse_args()

    with open(args.scored) as f:
        scored = json.load(f)

    scored_map = {s["name"]: s for s in scored}

    results = []
    for s in scored:
        name = s["name"]
        pdb_file = os.path.join(args.pdb_dir, f"{name}_complex.pdb")
        if not os.path.exists(pdb_file):
            continue

        n_contacts, contact_details = count_contacts(
            pdb_file, target_resids=TARGET_ACIDIC_RESIDS, cutoff=CUTOFF)

        entry = {
            "name": name,
            "sequence": s.get("sequence", ""),
            "length": s.get("length", 0),
            "n_his": s.get("n_his", 0),
            "mean_plddt": s.get("mean_plddt", 0),
            "iptm": s.get("iptm", 0),
            "L": s.get("L", 0),
            "n_his_contacts": n_contacts,
            "contacts": contact_details,
        }
        results.append(entry)

    # Sort: ipTM gate, then by n_contacts desc, then ipTM desc
    gated = [r for r in results if r["iptm"] >= args.iptm_gate]
    gated.sort(key=lambda x: (-x["n_his_contacts"], -x["iptm"]))

    print(f"{'Name':<40s} {'ipTM':>6s} {'#H':>3s} {'Cont':>4s} {'L':>6s} {'pLDDT':>6s} {'Contacts'}")
    print("-" * 110)
    for r in gated:
        contact_str = "; ".join(
            f"H{c['his_resnum']}.{c['his_atom']}-{c['target_resname']}{c['target_resnum']}.{c['target_atom']}({c['distance']})"
            for c in r["contacts"]
        )
        print(f"{r['name']:<40s} {r['iptm']:6.3f} {r['n_his']:3d} {r['n_his_contacts']:4d} "
              f"{r['L']:6.2f} {r['mean_plddt']:6.1f}  {contact_str}")

    n_with_contacts = sum(1 for r in gated if r["n_his_contacts"] > 0)
    n_multi = sum(1 for r in gated if r["n_his_contacts"] >= 2)
    print(f"\n{len(gated)} designs pass ipTM ≥ {args.iptm_gate}")
    print(f"{n_with_contacts} have ≥1 His–carboxylate contact")
    print(f"{n_multi} have ≥2 His–carboxylate contacts")

    with open(args.out, "w") as f:
        json.dump(gated, f, indent=2)
    print(f"Saved to {args.out}")