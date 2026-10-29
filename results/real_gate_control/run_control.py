"""CPU-only deposited-complex control; execute the existing audit unchanged."""
import csv
import hashlib
import json
from pathlib import Path
import sys

import numpy as np
from Bio.Data.PDBData import protein_letters_3to1 as AA
from Bio.PDB import MMCIFParser, Model, Chain
from Bio.PDB.MMCIF2Dict import MMCIF2Dict

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'scripts'))
from run_context_ensemble import fit, residues, write_protein, sha, min_dist
import audit_full_context as audit

OUT = Path(__file__).resolve().parent


def save(name, value):
    (OUT / name).write_text(json.dumps(value, indent=2) + '\n')


def table(cif, prefix):
    keys = [k for k in cif if k.startswith(prefix + '.')]
    return [{k.split('.', 1)[1]: v for k, v in zip(keys, row)}
            for row in zip(*(cif[k] for k in keys))]


def main():
    prereg = OUT / 'PREREGISTRATION.md'
    assert sha(prereg) == (OUT / 'PREREGISTRATION.sha256').read_text().strip()
    protected = [ROOT / 'scripts' / name for name in
                 ['audit_full_context.py', 'run_context_ensemble.py',
                  'postprocess_gpu_campaign.py']]
    before = {str(p.relative_to(ROOT)): sha(p) for p in protected}
    x = MMCIF2Dict(str(ROOT / 'data/6w3h.cif'))
    refs = {r['id']: r for r in table(x, '_struct_ref')}
    segments = [r for r in table(x, '_struct_ref_seq')
                if r['pdbx_db_accession'] == 'P02786']
    entities = {refs[r['ref_id']]['entity_id'] for r in segments}
    assert len(entities) == 1
    entity = entities.pop()
    scheme = [r for r in table(x, '_pdbx_poly_seq_scheme')
              if r['entity_id'] == entity]
    declared = [r for r in table(x, '_entity_poly_seq') if r['entity_id'] == entity]
    assert len(declared) == 166
    differences = table(x, '_struct_ref_seq_dif')
    native_map = {}
    for s in segments:
        lo, hi = int(s['seq_align_beg']), int(s['seq_align_end'])
        dblo, dbhi = int(s['db_align_beg']), int(s['db_align_end'])
        alo, ahi = int(s['pdbx_auth_seq_align_beg']), int(s['pdbx_auth_seq_align_end'])
        assert hi-lo == dbhi-dblo == ahi-alo
        for label, auth, native in zip(range(lo, hi+1), range(alo, ahi+1), range(dblo, dbhi+1)):
            native_map[s['pdbx_strand_id'], label] = (auth, native)

    # Parse label identifiers explicitly; author identifiers come from the scheme.
    source = MMCIFParser(QUIET=True, auth_chains=False, auth_residues=False).get_structure(
        '6w3h', ROOT / 'data/6w3h.cif')[0]
    template = MMCIFParser(QUIET=True).get_structure('11qc', ROOT / 'data/11qc.cif')[0]
    chain_info = {}
    mapping_rows = []
    for s in scheme:
        label, author = s['asym_id'], s['pdb_strand_id']
        pos = int(s['seq_id'])
        entry = native_map.get((author, pos))
        diff = [d for d in differences if d['pdbx_pdb_strand_id'] == author
                and int(d['seq_num']) == pos]
        observed = (' ', pos, ' ') in source[label] and 'CA' in source[label][pos]
        if entry:
            assert entry[0] == int(s['pdb_seq_num'])
        else:
            assert len(diff) == 1, 'Every unmapped construct residue needs annotation'
        mapping_rows.append(dict(label_chain=label, auth_chain=author,
            label_seq_id=pos, declared_author_number=int(s['pdb_seq_num']),
            observed_author_number=s['auth_seq_num'], monomer=s['mon_id'],
            uniprot_position=entry[1] if entry else None, observed_CA=observed,
            unmapped_reason=diff[0]['details'] if not entry else None))
        chain_info[label] = dict(auth_chain=author, entity_id=entity)

    # Extend label/auth identity to other polymer entities without chain-letter assumptions.
    for r in table(x, '_pdbx_poly_seq_scheme'):
        chain_info[r['asym_id']] = dict(auth_chain=r['pdb_strand_id'], entity_id=r['entity_id'])
    assemblies = table(x, '_pdbx_struct_assembly_gen')
    pairs = []
    for a in assemblies:
        assert a['oper_expression'] == '1', 'Need explicit symmetry expansion otherwise'
        labels = a['asym_id_list'].split(',')
        receptor = [c for c in labels if c in chain_info and chain_info[c]['entity_id'] == entity]
        partner = [c for c in labels if c in chain_info and chain_info[c]['entity_id'] != entity]
        assert len(receptor) == len(partner) == 1
        pairs.append(dict(assembly=a['assembly_id'], receptor_label=receptor[0],
                          partner_label=partner[0], receptor_author=chain_info[receptor[0]]['auth_chain'],
                          partner_author=chain_info[partner[0]]['auth_chain']))
    save('MAPPING.json', dict(receptor_entity_id=entity, declared_residues=len(declared),
                             segments=segments, all_positions=mapping_rows, pairs=pairs))
    with (OUT / 'MAPPING.csv').open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(mapping_rows[0])); w.writeheader(); w.writerows(mapping_rows)

    sandbox = OUT / 'audit_fixture'
    (sandbox / 'data').mkdir(parents=True, exist_ok=True)
    inputs = sandbox / 'results/v3_all_pdb'; inputs.mkdir(parents=True, exist_ok=True)
    link = sandbox / 'data/11qc.cif'
    if not link.exists(): link.symlink_to(ROOT / 'data/11qc.cif')
    fits = []
    full_clearance = []
    for pair in pairs:
        rows = [r for r in mapping_rows if r['label_chain'] == pair['receptor_label']
                and r['uniprot_position'] is not None and r['observed_CA']]
        for native_chain in 'AB':
            moving = [source[r['label_chain']][r['label_seq_id']]['CA'].coord for r in rows]
            fixed = [template[native_chain][r['uniprot_position']]['CA'].coord for r in rows]
            rotation, translation, rmsd = fit(moving, fixed)
            reverse_r, reverse_t, partner_fit = fit(fixed, moving)
            assert abs(rmsd-partner_fit) < 1e-5
            assert np.max(np.abs(rotation @ reverse_r - np.eye(3))) < 1e-5
            name = f"assembly{pair['assembly']}_Fc{pair['partner_author']}_TfR{pair['receptor_author']}_native{native_chain}"
            model = Model.Model(0); target = Chain.Chain(native_chain); model.add(target)
            mismatches = []
            for row in rows:
                r = source[row['label_chain']][row['label_seq_id']].copy()
                r.id = (' ', row['uniprot_position'], ' ')
                if AA[r.resname] != AA[template[native_chain][r.id].resname]:
                    mismatches.append(dict(label=row['label_seq_id'], native=row['uniprot_position'],
                                           source=r.resname, template=template[native_chain][r.id].resname))
                for atom in r.get_atoms(): atom.coord = atom.coord @ rotation + translation
                target.add(r)
            fc = source[pair['partner_label']].copy(); fc.id = 'C'
            for atom in fc.get_atoms(): atom.coord = atom.coord @ rotation + translation
            model.add(fc); write_protein(model, inputs / (name + '.pdb'))
            restored = Model.Model(0)
            for cid in 'AB': restored.add(template[cid].copy())
            restored.add(fc.copy()); write_protein(restored, OUT / (name + '_restored.pdb'))
            bb = lambda ids: np.array([a.coord for cid in ids for r in residues(restored[cid])
                                       for a in r if a.name in ('N', 'CA', 'C', 'O')])
            full_clearance.append(dict(design=name, full_native_backbone_minimum_A=min_dist(bb('C'), bb('AB')),
                method='Existing min_dist() helper on full restored artifact; supplementary, not audit reimplementation'))
            fits.append(dict(design=name, **pair, native_chain=native_chain,
                observed_mapped_CA=len(rows), source_to_native_CA_rmsd=rmsd,
                partner_fit=partner_fit, sequence_mismatches=mismatches,
                rotation=rotation.tolist(), translation=translation.tolist(),
                partner_observed_residues=len(residues(fc))))
    save('FITS.json', fits)
    save('FULL_NATIVE_CLEARANCE.json', full_clearance)

    # This calls the original function, not a copied/reimplemented clash calculation.
    audit.ROOT = sandbox
    audit.main()
    raw = json.loads((sandbox / 'results/full_context_audit/candidates.json').read_text())
    by_name = {r['design']: r for r in fits}
    classified = []
    for r in raw:
        f = by_name[r['design']]
        fit_error = max(f['partner_fit'], r['target_ca_rmsd'])
        n, minimum = r['restored_backbone_pair_clashes_under2A'], r['restored_backbone_min_distance']
        if fit_error > 2: band = 'INCONCLUSIVE'
        elif n <= 9 and minimum >= 1.5: band = 'GATE VINDICATED'
        elif n >= 10 and minimum < 1: band = 'GATE FALSIFIED'
        else: band = 'INCONCLUSIVE'
        classified.append(dict(**r, partner_fit=f['partner_fit'], band=band))
    bands = {r['band'] for r in classified}
    overall = ('GATE FALSIFIED' if 'GATE FALSIFIED' in bands else
               'GATE VINDICATED' if 'GATE VINDICATED' in bands else 'INCONCLUSIVE')
    after = {str(p.relative_to(ROOT)): sha(p) for p in protected}
    assert before == after
    save('RESULTS.json', dict(overall_band=overall, controls=classified,
        protected_script_hashes_before=before, protected_script_hashes_after=after,
        source_sha256=sha(ROOT / 'data/6w3h.cif'), template_sha256=sha(link),
        preregistration_sha256=sha(prereg), new_gpu_spend_usd=0,
        audit_execution='Unchanged audit_full_context.main(), ROOT redirected to isolated fixture'))
    print(json.dumps(dict(overall_band=overall, controls=[dict(design=r['design'],
        band=r['band'], rmsd=r['partner_fit'], clashes=r['restored_backbone_pair_clashes_under2A'],
        minimum=r['restored_backbone_min_distance']) for r in classified])))


if __name__ == '__main__': main()
