"""Build C4 submission v8 using geometric contact + ipTM ranking.

PROPKA Wyman linkage L is unreliable (FcRn control: L=+0.13 vs real ~+2-3).
Instead rank by: (1) His-to-target-D/E contacts within 4A, (2) ipTM.

Tier 1: designs with His contacts AND ipTM >= 0.75
Tier 2: designs with His contacts AND ipTM >= 0.5
Tier 3: high ipTM unconditional binders (no contacts, controls)
"""
import argparse
import csv
import json
import os
import math
from Bio.PDB import PDBParser

TARGET_ACIDIC_RESIDS = {63, 86, 98, 118, 132}
CUTOFF_STRICT = 4.0
CUTOFF_LOOSE = 6.0


def count_contacts(pdb_path, cutoff=4.0):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("cx", pdb_path)
    model = next(structure.get_models())

    tgt_oxygens = []
    for res in model['B']:
        if res.id[0] != ' ':
            continue
        if res.id[1] not in TARGET_ACIDIC_RESIDS:
            continue
        if res.get_resname() in ('ASP', 'GLU'):
            for an in ('OD1', 'OD2', 'OE1', 'OE2'):
                if an in res:
                    tgt_oxygens.append(res[an].get_vector().get_array())

    c2 = cutoff * cutoff
    his_hit = set()
    for res in model['A']:
        if res.id[0] != ' ' or res.get_resname() != 'HIS':
            continue
        for ha in ('ND1', 'NE2'):
            if ha not in res:
                continue
            hc = res[ha].get_vector().get_array()
            for tc in tgt_oxygens:
                d2 = sum((a - b) ** 2 for a, b in zip(hc, tc))
                if d2 <= c2:
                    his_hit.add(res.id[1])
                    break
            if res.id[1] in his_hit:
                break

    return len(his_hit)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--scored", required=True)
    p.add_argument("--pdb-dir", required=True)
    p.add_argument("--max", type=int, default=20)
    p.add_argument("--out", default="results/submission_v8.csv")
    args = p.parse_args()

    with open(args.scored) as f:
        scored = json.load(f)

    designs = []
    for s in scored:
        pdb = os.path.join(args.pdb_dir, "{}_complex.pdb".format(s['name']))
        n4 = count_contacts(pdb, CUTOFF_STRICT) if os.path.exists(pdb) else 0
        n6 = count_contacts(pdb, CUTOFF_LOOSE) if os.path.exists(pdb) else 0
        s['n_contacts_4A'] = n4
        s['n_contacts_6A'] = n6
        designs.append(s)

    # Tier 1: contacts within 4A + ipTM >= 0.75
    tier1 = [d for d in designs if d['n_contacts_4A'] > 0 and d.get('iptm', 0) >= 0.75]
    tier1.sort(key=lambda x: (-x['n_contacts_4A'], -x.get('iptm', 0)))

    # Tier 2: contacts within 4A + ipTM >= 0.5 (not in tier 1)
    tier2 = [d for d in designs if d['n_contacts_4A'] > 0 and d.get('iptm', 0) >= 0.5 and d not in tier1]
    tier2.sort(key=lambda x: (-x['n_contacts_4A'], -x.get('iptm', 0)))

    # Tier 2b: contacts within 6A + ipTM >= 0.75 (not in tier 1/2)
    tier2b = [d for d in designs if d['n_contacts_6A'] > 0 and d['n_contacts_4A'] == 0
              and d.get('iptm', 0) >= 0.75]
    tier2b.sort(key=lambda x: (-x['n_contacts_6A'], -x.get('iptm', 0)))

    # Tier 3: high ipTM unconditional binders (controls)
    used = set(d['name'] for d in tier1 + tier2 + tier2b)
    tier3 = [d for d in designs if d['name'] not in used and d.get('iptm', 0) >= 0.75]
    tier3.sort(key=lambda x: -x.get('iptm', 0))

    selected = []
    for tier_name, tier_list in [("T1:contact+ipTM", tier1), ("T2:contact", tier2),
                                  ("T2b:near+ipTM", tier2b), ("T3:unconditional", tier3)]:
        for d in tier_list:
            if len(selected) >= args.max:
                break
            d['tier'] = tier_name
            selected.append(d)

    fmt = "{:<40s} {:>6s} {:>4s} {:>4s} {:>6s} {:>6s} {:>3s}  {}"
    print(fmt.format('Name', 'ipTM', '<4A', '<6A', 'L', 'pLDDT', '#H', 'Tier'))
    print("-" * 100)
    for s in selected:
        print(fmt.format(
            s['name'], "{:.3f}".format(s.get('iptm', 0)),
            str(s['n_contacts_4A']), str(s['n_contacts_6A']),
            "{:.2f}".format(s.get('L', 0)),
            "{:.1f}".format(s.get('mean_plddt', 0)),
            str(s.get('n_his', 0)),
            s.get('tier', '')))

    print("\nTier counts: T1={}, T2={}, T2b={}, T3={}".format(
        sum(1 for s in selected if 'T1' in s.get('tier', '')),
        sum(1 for s in selected if s.get('tier', '') == 'T2:contact'),
        sum(1 for s in selected if 'T2b' in s.get('tier', '')),
        sum(1 for s in selected if 'T3' in s.get('tier', ''))))

    rows = []
    for s in selected:
        rows.append({
            "name": s["name"],
            "sequence": s["sequence"],
            "molecule_class": "single_chain",
        })

    with open(args.out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["name", "sequence", "molecule_class"])
        w.writeheader()
        w.writerows(rows)
    print("\nSaved {} designs to {}".format(len(rows), args.out))


if __name__ == "__main__":
    main()
