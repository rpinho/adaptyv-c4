"""Fetch all Boltz-2 complex PDBs from Modal volume and score geometrically.

Downloads any new PDB files from Modal volume boltz2_c4_acid/ that aren't
already local, then runs geometric His-D/E contact scoring on all of them.
"""
import json
import os
import subprocess
import math
from Bio.PDB import PDBParser

TARGET_ACIDIC_RESIDS = {63, 86, 98, 118, 132}

def count_contacts(pdb_path, cutoff=4.0):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("cx", pdb_path)
    model = next(structure.get_models())

    if 'B' not in model:
        return 0, []

    tgt_oxygens = []
    for res in model['B']:
        if res.id[0] != ' ':
            continue
        if res.id[1] not in TARGET_ACIDIC_RESIDS:
            continue
        if res.get_resname() in ('ASP', 'GLU'):
            for an in ('OD1', 'OD2', 'OE1', 'OE2'):
                if an in res:
                    tgt_oxygens.append((res.id[1], res.get_resname(), an,
                                       res[an].get_vector().get_array()))

    c2 = cutoff * cutoff
    contacts = []
    for res in model['A']:
        if res.id[0] != ' ' or res.get_resname() != 'HIS':
            continue
        resnum = res.id[1]
        for ha in ('ND1', 'NE2'):
            if ha not in res:
                continue
            hc = res[ha].get_vector().get_array()
            for tresnum, tresname, taname, tc in tgt_oxygens:
                d2 = sum((a - b) ** 2 for a, b in zip(hc, tc))
                if d2 <= c2:
                    contacts.append({
                        "his_resnum": resnum,
                        "his_atom": ha,
                        "target_resnum": tresnum,
                        "target_resname": tresname,
                        "target_atom": taname,
                        "distance": round(math.sqrt(d2), 2),
                    })

    unique_his = len(set(c["his_resnum"] for c in contacts))
    return unique_his, contacts


def main():
    base = "/Users/ricardo/dev/adaptyv-c4"
    pdb_dir = os.path.join(base, "results/boltz2_complexes/boltz2_c4_acid")
    os.makedirs(pdb_dir, exist_ok=True)

    # 1. Fetch new PDBs from Modal volume
    print("Fetching PDB list from Modal volume...")
    result = subprocess.run(
        ["/Users/ricardo/Library/Python/3.9/bin/modal", "volume", "ls",
         "rfd4-runs", "boltz2_c4_acid/"],
        capture_output=True, text=True
    )

    remote_files = []
    for line in result.stdout.strip().split('\n'):
        line = line.strip()
        if line.endswith('_complex.pdb'):
            remote_files.append(line)

    local_files = set(os.listdir(pdb_dir)) if os.path.isdir(pdb_dir) else set()
    new_files = [f for f in remote_files if f not in local_files]

    if new_files:
        print("Downloading {} new PDB files...".format(len(new_files)))
        subprocess.run(
            ["/Users/ricardo/Library/Python/3.9/bin/modal", "volume", "get",
             "rfd4-runs", "boltz2_c4_acid/", pdb_dir + "/"],
            capture_output=True, text=True
        )
    else:
        print("No new PDB files to download ({} already local)".format(len(local_files)))

    # 2. Load all refold data for metadata
    refold_path = os.path.join(base, "results/c4_acid_combined_refold.json")
    with open(refold_path) as f:
        all_refold = json.load(f)
    refold_map = {r["name"]: r for r in all_refold}

    # 3. Load existing Boltz-2 scores for ipTM
    boltz_data = {}
    for bf in ["c4_acid_boltz2_batch1.json", "c4_acid_boltz2_batch2.json",
               "c4_acid_boltz2_batch3.json"]:
        bpath = os.path.join(base, "results", bf)
        if os.path.exists(bpath):
            with open(bpath) as f:
                for r in json.load(f):
                    boltz_data[r["name"]] = r

    # 4. Score all PDB files
    pdb_files = sorted(f for f in os.listdir(pdb_dir) if f.endswith("_complex.pdb"))
    print("\nScoring {} complex PDBs...".format(len(pdb_files)))

    results = []
    for pf in pdb_files:
        name = pf.replace("_complex.pdb", "")
        pdb_path = os.path.join(pdb_dir, pf)

        n4, contacts4 = count_contacts(pdb_path, 4.0)
        n6, _ = count_contacts(pdb_path, 6.0)

        refold = refold_map.get(name, {})
        boltz = boltz_data.get(name, {})

        results.append({
            "name": name,
            "sequence": refold.get("sequence", boltz.get("sequence", "")),
            "length": refold.get("length", 0),
            "n_his": refold.get("n_his", 0),
            "mean_plddt": refold.get("mean_plddt", 0),
            "iptm": boltz.get("iptm", 0),
            "ptm": boltz.get("ptm", 0),
            "confidence": boltz.get("confidence", 0),
            "n_contacts_4A": n4,
            "n_contacts_6A": n6,
            "contacts": contacts4,
        })

    # 5. Sort and display
    results.sort(key=lambda x: (-x["n_contacts_4A"], -x["iptm"]))

    fmt = "{:<40s} {:>6s} {:>4s} {:>4s} {:>6s} {:>3s}  {}"
    print("\n" + fmt.format('Name', 'ipTM', '<4A', '<6A', 'pLDDT', '#H', 'Contacts'))
    print("-" * 110)

    n_with_contacts = 0
    for r in results:
        if r["n_contacts_4A"] > 0:
            n_with_contacts += 1
        contact_str = "; ".join(
            "H{}.{}-{}{}.{}({})".format(
                c['his_resnum'], c['his_atom'],
                c['target_resname'], c['target_resnum'],
                c['target_atom'], c['distance'])
            for c in r["contacts"]
        )
        print(fmt.format(
            r['name'],
            "{:.3f}".format(r['iptm']),
            str(r['n_contacts_4A']),
            str(r['n_contacts_6A']),
            "{:.1f}".format(r['mean_plddt']),
            str(r['n_his']),
            contact_str))

    print("\n{} total scored, {} with His-D/E contacts <4A".format(
        len(results), n_with_contacts))

    # 6. Save
    out_path = os.path.join(base, "results/c4_acid_all_geom_scored.json")
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2)
    print("Saved to {}".format(out_path))

    return results


if __name__ == "__main__":
    main()
