"""Read deposited 6W3H mappings; CPU geometry calibration only, never design input."""
import hashlib
import json
from pathlib import Path
import numpy as np
from Bio.PDB import MMCIFParser, MMCIF2Dict

ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).parent
BB = {'N', 'CA', 'C', 'O'}
EXPECTED = {'6w3h': '8b8886ba4abb68f5206d721a5c29c4613170ab8d79a80b13952687107fa3a155',
            '11qc': '4eaa5db67ed83a896713c4b42d1dc214c778bdac240fdf3a9a82e86a66f79495'}

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def protein(chain):
    return [r for r in chain if r.id[0] == ' ']

def atom_records(residues, backbone=False):
    return [(a.coord.astype(float), [r.parent.id, r.id[1], r.resname, a.name])
            for r in residues for a in r
            if a.element != 'H' and (not backbone or a.name in BB)]

def nearest(a, b, cutoff=2.0):
    if not a or not b:
        raise ValueError('Missing geometry')
    ax = np.array([x for x, _ in a]); bx = np.array([x for x, _ in b])
    if not np.isfinite(ax).all() or not np.isfinite(bx).all():
        raise ValueError('Nonfinite geometry')
    best = (float('inf'), None, None); count = 0
    for start in range(0, len(ax), 100):
        d = np.linalg.norm(ax[start:start+100, None] - bx[None], axis=-1)
        i, j = np.unravel_index(d.argmin(), d.shape)
        if d[i, j] < best[0]:
            best = (float(d[i, j]), a[start+i][1], b[j][1])
        count += int((d < cutoff).sum())
    return dict(minimum_A=best[0], Fc_atom=best[1], context_atom=best[2],
                atom_pairs_below_2A=count)

def fit(x, y):
    xc, yc = x.mean(0), y.mean(0)
    u, _, vt = np.linalg.svd((x-xc).T @ (y-yc))
    correction = np.eye(3); correction[2, 2] = np.linalg.det(u @ vt)
    rotation = u @ correction @ vt
    return rotation, yc-xc @ rotation

def residuals(x, y, rotation, translation):
    distances = np.linalg.norm(x @ rotation + translation-y, axis=-1)
    return dict(RMSD_A=float(np.sqrt(np.mean(distances**2))),
                median_A=float(np.median(distances)), p95_A=float(np.quantile(distances, .95)),
                maximum_A=float(distances.max()))

def main():
    models = {}; dictionaries = {}
    for name, expected in EXPECTED.items():
        path = ROOT/'data'/(name+'.cif')
        if sha(path) != expected:
            raise ValueError('Input identity differs')
        structure = MMCIFParser(QUIET=True).get_structure(name, path)
        if len(structure) != 1:
            raise ValueError('Unexpected model multiplicity')
        models[name] = structure[0]
        dictionaries[name] = MMCIF2Dict.MMCIF2Dict(path)
    source, native = models['6w3h'], models['11qc']
    d = dictionaries['6w3h']
    if d['_pdbx_struct_oper_list.type'] != ['identity operation']:
        raise ValueError('Assembly operator needs additional handling')
    if d['_pdbx_struct_assembly_gen.asym_id_list'] != ['A,C,E', 'B,D']:
        raise ValueError('Deposited assembly pairing differs')
    label_author = {}
    for label, author in zip(d['_atom_site.label_asym_id'], d['_atom_site.auth_asym_id']):
        if label in label_author and label_author[label] != author:
            raise ValueError('Ambiguous assembly label/author mapping')
        label_author[label] = author
    assembly_pairs = []
    for labels in d['_pdbx_struct_assembly_gen.asym_id_list']:
        authors = {label_author[label] for label in labels.split(',')}
        receptor = authors & set('CD'); fc = authors & set('AB')
        if len(receptor) != 1 or len(fc) != 1:
            raise ValueError('Cannot resolve deposited therapeutic complex pairing')
        assembly_pairs.append((next(iter(receptor)), next(iter(fc))))
    maps = {c: {} for c in 'CD'}; segments = {c: [] for c in 'CD'}
    pre = '_struct_ref_seq.'
    for index, chain in enumerate(d[pre+'pdbx_strand_id']):
        if chain not in maps or d[pre+'pdbx_db_accession'][index] != 'P02786':
            continue
        lo = int(d[pre+'pdbx_auth_seq_align_beg'][index]); hi = int(d[pre+'pdbx_auth_seq_align_end'][index])
        nlo = int(d[pre+'db_align_beg'][index]); nhi = int(d[pre+'db_align_end'][index])
        if hi-lo != nhi-nlo:
            raise ValueError('Nonlinear deposited alignment')
        segments[chain].append(dict(author_begin=lo, author_end=hi, native_begin=nlo, native_end=nhi))
        for old, new in zip(range(lo, hi+1), range(nlo, nhi+1)):
            if old in maps[chain]:
                raise ValueError('Overlapping deposited mappings')
            maps[chain][old] = new
    report = dict(scope='Therapeutic experimental geometry calibration only; never design seed or admission',
                  input_sha256=EXPECTED, script_sha256=sha(__file__),
                  assembly_label_to_author=label_author,
                  local_fit_definition='Mapped receptor residues with any heavy atom within5A of paired Fc heavy atom in deposited complex',
                  all_mapped_CA_used_for_global_fit=True, artificial_join_treated_as_native_bond=False,
                  cases=[], assemblies=[])
    for receptor, fc in assembly_pairs:
        rr = protein(source[receptor]); fr = protein(source[fc])
        if any(r.id[2] != ' ' or r.id[1] not in maps[receptor] or 'CA' not in r for r in rr):
            raise ValueError('Observed receptor cannot be mapped unambiguously')
        mapping = [(r, maps[receptor][r.id[1]]) for r in rr]
        for r, n in mapping:
            if any(r.resname != native[c][n].resname for c in 'AB'):
                raise ValueError('Mapped receptor identity differs from native')
        fx = atom_records(fr)
        interface = [i for i, (r, _) in enumerate(mapping)
                     if nearest(fx, atom_records([r]))['minimum_A'] <= 5]
        if len(interface) < 3:
            raise ValueError('Insufficient observed interface for local diagnostic')
        report['assemblies'].append(dict(receptor=receptor, Fc=fc, mapped_observed_receptor_residues=len(rr),
            observed_Fc_residues=len(fr), declared_Fc_residues=227, declared_receptor_residues=166,
            deposited_segments=segments[receptor], observed_mapping=[[r.id[1], n] for r, n in mapping],
            interface_native_positions=[mapping[i][1] for i in interface],
            deposited_backbone_clearance=nearest(atom_records(fr, True), atom_records(rr, True))))
        x = np.array([r['CA'].coord for r, _ in mapping], dtype=float)
        for native_chain in 'AB':
            y = np.array([native[native_chain][n]['CA'].coord for _, n in mapping], dtype=float)
            mapped_positions = {n for _, n in mapping}
            for mode, indices in [('global_all_mapped', list(range(len(mapping)))), ('local_interface_diagnostic', interface)]:
                rotation, translation = fit(x[indices], y[indices])
                fb = [(xyz @ rotation+translation, ident) for xyz, ident in atom_records(fr, True)]
                fh = [(xyz @ rotation+translation, ident) for xyz, ident in fx]
                opposite = 'B' if native_chain == 'A' else 'A'
                groups = {
                    'mapped_same_subunit_backbone': atom_records([native[native_chain][n] for n in sorted(mapped_positions)], True),
                    'omitted_same_subunit_backbone': atom_records([r for r in protein(native[native_chain]) if r.id[1] not in mapped_positions], True),
                    'opposite_subunit_backbone': atom_records(protein(native[opposite]), True),
                    'full_native_AB_backbone': atom_records(protein(native['A'])+protein(native['B']), True)}
                clearances = {name: nearest(fb, atoms) for name, atoms in groups.items()}
                nag = [r for c in native for r in c if r.resname == 'NAG']
                clearances['native_NAG_to_Fc_heavy'] = nearest(fh, atom_records(nag))
                report['cases'].append(dict(receptor=receptor, Fc=fc, fitted_native_chain=native_chain,
                    fit_mode=mode, fit_CA_count=len(indices), global_residual=residuals(x,y,rotation,translation),
                    interface_residual=residuals(x[interface],y[interface],rotation,translation),
                    receptor_CA_residuals=[dict(source_author_id=r.id[1], native_author_id=n,
                        residual_A=float(distance), source_CA_B_factor=float(r['CA'].bfactor))
                        for (r,n),distance in zip(mapping,np.linalg.norm(x@rotation+translation-y,axis=-1))],
                    rotation=rotation.tolist(), translation=translation.tolist(), clearances=clearances))
    (OUT/'CPU_CALIBRATION.json').write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
    for c in report['cases']:
        print(c['receptor'],c['Fc'],c['fitted_native_chain'],c['fit_mode'],
              'globalRMS',round(c['global_residual']['RMSD_A'],3),
              'localRMS',round(c['interface_residual']['RMSD_A'],3),
              'fullBBmin',round(c['clearances']['full_native_AB_backbone']['minimum_A'],3),
              'NAGmin',round(c['clearances']['native_NAG_to_Fc_heavy']['minimum_A'],3))

if __name__ == '__main__':
    main()
