"""Score C4 designs for ACID-OFF pH switching (bind at pH 7.4, release at pH 6.0).

Correct direction for Challenge 4: TfR1 binder that binds at neutral pH
and releases in the acidifying endosome.

Scoring criteria for acid-OFF His mechanisms:
1. His near target Arg/Lys (cationic) — protonation creates charge repulsion
2. His in buried interface positions — protonation penalty in hydrophobic environment
3. His making interface H-bonds as neutral acceptor — lost upon protonation
4. PENALTY: His near target Asp/Glu — stabilizes low-pH binding (wrong direction)

Usage:
    python3 scripts/score_acid_off.py \
        --pdb-dir results/boltz2_complexes/boltz2_c4_acid/ \
        --refold results/c4_acid_combined_refold.json \
        --out results/c4_acid_off_scored.json
"""
import argparse
import json
import math
import os
from Bio.PDB import PDBParser, NeighborSearch, DSSP

CATIONIC_ATOMS = {
    'ARG': ('NH1', 'NH2', 'NE', 'CZ'),
    'LYS': ('NZ',),
}
ANIONIC_ATOMS = {
    'ASP': ('OD1', 'OD2'),
    'GLU': ('OE1', 'OE2'),
}
HIS_N_ATOMS = ('ND1', 'NE2')

CUTOFF_CONTACT = 4.5
CUTOFF_BURIAL = 5.0


def get_interface_his(structure, binder_chain='A', target_chain='B'):
    """Find binder His residues at the interface (any atom within 5A of target)."""
    model = next(structure.get_models())
    if binder_chain not in model or target_chain not in model:
        return []

    target_atoms = list(model[target_chain].get_atoms())
    ns = NeighborSearch(target_atoms)

    interface_his = []
    for res in model[binder_chain]:
        if res.id[0] != ' ' or res.get_resname() != 'HIS':
            continue
        for atom in res:
            nearby = ns.search(atom.get_vector().get_array(), CUTOFF_BURIAL)
            if nearby:
                interface_his.append(res)
                break

    return interface_his


def score_design(pdb_path):
    """Score a single design for acid-OFF potential."""
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('cx', pdb_path)
    model = next(structure.get_models())

    if 'A' not in model or 'B' not in model:
        return None

    # Collect target cationic atom coords
    tgt_cationic = []
    for res in model['B']:
        if res.id[0] != ' ':
            continue
        rn = res.get_resname()
        if rn in CATIONIC_ATOMS:
            for aname in CATIONIC_ATOMS[rn]:
                if aname in res:
                    tgt_cationic.append({
                        'resnum': res.id[1], 'resname': rn, 'atom': aname,
                        'coord': res[aname].get_vector().get_array()
                    })

    # Collect target anionic atom coords (penalty)
    tgt_anionic = []
    for res in model['B']:
        if res.id[0] != ' ':
            continue
        rn = res.get_resname()
        if rn in ANIONIC_ATOMS:
            for aname in ANIONIC_ATOMS[rn]:
                if aname in res:
                    tgt_anionic.append({
                        'resnum': res.id[1], 'resname': rn, 'atom': aname,
                        'coord': res[aname].get_vector().get_array()
                    })

    # Collect ALL target atoms for burial calculation
    tgt_all_coords = []
    for atom in model['B'].get_atoms():
        tgt_all_coords.append(atom.get_vector().get_array())

    # Collect binder non-His atoms for burial
    binder_non_his_coords = []
    for res in model['A']:
        if res.id[0] != ' ':
            continue
        if res.get_resname() != 'HIS':
            for atom in res:
                binder_non_his_coords.append(atom.get_vector().get_array())

    # Score each binder His
    his_scores = []
    for res in model['A']:
        if res.id[0] != ' ' or res.get_resname() != 'HIS':
            continue

        resnum = res.id[1]
        h = {'resnum': resnum, 'cationic_contacts': [], 'anionic_contacts': [],
             'burial_count': 0, 'min_cationic_dist': 999.0, 'min_anionic_dist': 999.0}

        for ha in HIS_N_ATOMS:
            if ha not in res:
                continue
            hc = res[ha].get_vector().get_array()

            # Cationic contacts (GOOD for acid-OFF)
            for tc in tgt_cationic:
                d = math.sqrt(sum((a - b) ** 2 for a, b in zip(hc, tc['coord'])))
                if d < h['min_cationic_dist']:
                    h['min_cationic_dist'] = d
                if d <= CUTOFF_CONTACT:
                    h['cationic_contacts'].append({
                        'his_atom': ha, 'target': '{}{}'.format(tc['resname'], tc['resnum']),
                        'target_atom': tc['atom'], 'distance': round(d, 2)
                    })

            # Anionic contacts (BAD for acid-OFF — stabilizes low-pH binding)
            for tc in tgt_anionic:
                d = math.sqrt(sum((a - b) ** 2 for a, b in zip(hc, tc['coord'])))
                if d < h['min_anionic_dist']:
                    h['min_anionic_dist'] = d
                if d <= CUTOFF_CONTACT:
                    h['anionic_contacts'].append({
                        'his_atom': ha, 'target': '{}{}'.format(tc['resname'], tc['resnum']),
                        'target_atom': tc['atom'], 'distance': round(d, 2)
                    })

            # Burial: count heavy atoms within 5A (target + binder non-His)
            burial = 0
            for tc in tgt_all_coords:
                d2 = sum((a - b) ** 2 for a, b in zip(hc, tc))
                if d2 <= CUTOFF_BURIAL * CUTOFF_BURIAL:
                    burial += 1
            for bc in binder_non_his_coords:
                d2 = sum((a - b) ** 2 for a, b in zip(hc, bc))
                if d2 <= CUTOFF_BURIAL * CUTOFF_BURIAL:
                    burial += 1
            h['burial_count'] = max(h['burial_count'], burial)

        h['min_cationic_dist'] = round(h['min_cationic_dist'], 2)
        h['min_anionic_dist'] = round(h['min_anionic_dist'], 2)

        # Composite score for this His:
        # + for cationic contacts (protonation → repulsion)
        # + for high burial (protonation penalty in buried site)
        # - for anionic contacts (protonation → attraction, wrong direction)
        n_cat = len(h['cationic_contacts'])
        n_ani = len(h['anionic_contacts'])
        burial_score = 1.0 if h['burial_count'] >= 30 else 0.5 if h['burial_count'] >= 15 else 0.0

        h['switch_score'] = n_cat * 2.0 + burial_score - n_ani * 2.0
        his_scores.append(h)

    # Design-level aggregation
    n_cationic_his = sum(1 for h in his_scores if len(h['cationic_contacts']) > 0)
    n_anionic_his = sum(1 for h in his_scores if len(h['anionic_contacts']) > 0)
    n_buried_his = sum(1 for h in his_scores if h['burial_count'] >= 15)
    total_switch = sum(h['switch_score'] for h in his_scores)

    return {
        'n_interface_his': len(his_scores),
        'n_cationic_contacts': n_cationic_his,
        'n_anionic_contacts': n_anionic_his,
        'n_buried_his': n_buried_his,
        'total_switch_score': round(total_switch, 2),
        'his_details': his_scores,
    }


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--pdb-dir', required=True)
    p.add_argument('--refold', required=True)
    p.add_argument('--boltz-jsons', nargs='*', default=[])
    p.add_argument('--out', default='results/c4_acid_off_scored.json')
    p.add_argument('--top', type=int, default=30)
    args = p.parse_args()

    with open(args.refold) as f:
        all_refold = json.load(f)
    refold_map = {r['name']: r for r in all_refold}

    boltz_map = {}
    for bf in args.boltz_jsons:
        if os.path.exists(bf):
            with open(bf) as f:
                for r in json.load(f):
                    boltz_map[r['name']] = r

    # Also load default boltz files
    for bf in ['c4_acid_boltz2_batch1.json', 'c4_acid_boltz2_batch2.json',
               'c4_acid_boltz2_batch3.json']:
        bpath = os.path.join(os.path.dirname(args.out), bf)
        if os.path.exists(bpath) and bf not in [os.path.basename(x) for x in args.boltz_jsons]:
            with open(bpath) as f:
                for r in json.load(f):
                    if r['name'] not in boltz_map:
                        boltz_map[r['name']] = r

    pdb_files = sorted(f for f in os.listdir(args.pdb_dir) if f.endswith('_complex.pdb'))
    print('Scoring {} complex PDBs for acid-OFF...'.format(len(pdb_files)))

    results = []
    for pf in pdb_files:
        name = pf.replace('_complex.pdb', '')
        pdb_path = os.path.join(args.pdb_dir, pf)

        scores = score_design(pdb_path)
        if scores is None:
            continue

        refold = refold_map.get(name, {})
        boltz = boltz_map.get(name, {})

        entry = {
            'name': name,
            'sequence': refold.get('sequence', boltz.get('sequence', '')),
            'length': refold.get('length', 0),
            'n_his': refold.get('n_his', 0),
            'mean_plddt': refold.get('mean_plddt', 0),
            'iptm': boltz.get('iptm', 0),
        }
        entry.update(scores)
        results.append(entry)

    # Rank: high ipTM + positive switch score + no anionic contacts
    results.sort(key=lambda x: (-x['total_switch_score'], -x['iptm']))

    fmt = '{:<40s} {:>6s} {:>4s} {:>4s} {:>4s} {:>5s} {:>6s}'
    print('\n' + fmt.format('Name', 'ipTM', 'Cat', 'Ani', 'Bur', 'Sw', 'pLDDT'))
    print('-' * 85)
    for r in results[:args.top]:
        print(fmt.format(
            r['name'],
            '{:.3f}'.format(r['iptm']),
            str(r['n_cationic_contacts']),
            str(r['n_anionic_contacts']),
            str(r['n_buried_his']),
            '{:.1f}'.format(r['total_switch_score']),
            '{:.1f}'.format(r['mean_plddt'])))

    n_cat = sum(1 for r in results if r['n_cationic_contacts'] > 0)
    n_ani = sum(1 for r in results if r['n_anionic_contacts'] > 0)
    n_bur = sum(1 for r in results if r['n_buried_his'] > 0)
    print('\nSummary: {} designs scored'.format(len(results)))
    print('  With His-cationic contacts: {}'.format(n_cat))
    print('  With His-anionic contacts (penalty): {}'.format(n_ani))
    print('  With buried interface His: {}'.format(n_bur))

    with open(args.out, 'w') as f:
        json.dump(results, f, indent=2)
    print('Saved to {}'.format(args.out))


if __name__ == '__main__':
    main()
