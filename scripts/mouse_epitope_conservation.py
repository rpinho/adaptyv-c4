#!/usr/bin/env python3
"""Sequence-derived mouse mapping and measured human structural neighbourhoods."""
import csv
import hashlib
import json
from pathlib import Path

import numpy as np
from Bio import Align, SwissProt
from Bio.Align import substitution_matrices
from Bio.Data.PDBData import protein_letters_3to1 as AA
from Bio.PDB import MMCIFParser, PDBParser, ShrakeRupley
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from Bio.PDB.DSSP import residue_max_acc

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'results/mouse_epitope_conservation'
PATCH = (205, 206, 208, 209, 210)
DOMAIN = range(189, 384)  # Current project apical-domain boundary, human auth.


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def measure():
    OUT.mkdir(parents=True, exist_ok=True)
    paths = [ROOT/p for p in ('data/P02786.txt', 'data/Q62351.txt',
             'data/11qc.cif', 'data/AF-Q62351.pdb', 'targets/E1_mouse_his710.pdb')]
    seq = []
    for p in paths[:2]:
        with p.open() as f:
            seq.append(SwissProt.read(f).sequence)
    human, mouse = seq
    aligner = Align.PairwiseAligner(mode='global')
    aligner.substitution_matrix = substitution_matrices.load('BLOSUM62')
    aligner.open_gap_score, aligner.extend_gap_score = -10, -.5
    alignments = aligner.align(human, mouse)
    maps = [{int(h)+1: int(m)+1 for h, m in zip(*a.indices) if h >= 0 and m >= 0}
            for a in alignments]
    options = {h: sorted({m.get(h) for m in maps}, key=lambda x: -1 if x is None else x)
               for h in range(1, len(human)+1)}

    def row(h):
        opts = options[h]
        aas = [mouse[m-1] if m else '-' for m in opts]
        return dict(human_auth=h, human_label=h-121, human_aa=human[h-1],
                    mouse_uniprot_options=opts, mouse_aa_options=aas,
                    identical=all(a == human[h-1] for a in aas),
                    mapping_ambiguous=len(opts) > 1)

    rows = [row(h) for h in DOMAIN]
    with (OUT/'apical_alignment.csv').open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    (OUT/'all_optimal_alignments.txt').write_text('\n\n'.join(
        f'Optimal alignment {i+1}; score {a.score}\n{a}' for i, a in enumerate(alignments)))
    native = MMCIFParser(QUIET=True).get_structure('human', paths[2])[0]
    cif = MMCIF2Dict(str(paths[2]))
    label_pairs = {(int(a), int(l)) for a, l, c, g in zip(
        cif['_atom_site.auth_seq_id'], cif['_atom_site.label_seq_id'],
        cif['_atom_site.auth_asym_id'], cif['_atom_site.group_PDB'])
        if c == 'A' and g == 'ATOM'}
    assert all(a-l == 121 for a, l in label_pairs)
    af = PDBParser(QUIET=True).get_structure('mouse', paths[3])[0]['A']
    trim = PDBParser(QUIET=True).get_structure('trim', paths[4])[0]['A']
    for chain, sequence in ((native['A'], human), (af, mouse), (trim, mouse)):
        for r in chain:
            if r.id[0] == ' ' and r.resname in AA:
                assert AA[r.resname] == sequence[r.id[1]-1], (r, 'sequence mismatch')
    rr = {r.id[1]: r for r in native['A'] if r.id[0] == ' ' and r.resname in AA}
    coords = {n: np.array([a.coord for a in r if a.element not in ('H', 'D')])
              for n, r in rr.items()}

    def neighbours(n, radius=8):
        return [p for p in rr if p != n and
                np.linalg.norm(coords[n][:, None]-coords[p][None, :], axis=-1).min() <= radius]

    def counts(ns):
        # A deletion is counted separately, never mislabeled a substitution.
        substituted = [sum(m.get(h) is not None and mouse[m[h]-1] != human[h-1]
                           for h in ns) for m in maps]
        deleted = [sum(m.get(h) is None for h in ns) for m in maps]
        return dict(neighbour_count=len(ns), substitution_count_range=[min(substituted), max(substituted)],
                    deletion_count_range=[min(deleted), max(deleted)],
                    nonidentical_count=sum(not row(h)['identical'] for h in ns),
                    changed_residues=[row(h) for h in ns if not row(h)['identical']])

    chemistry = {
        205: 'Nonconservative K->S: long basic cationic side chain -> short polar hydroxyl; formal charge +1->0',
        206: 'Identical N: polar amide, neutral',
        208: 'Nonconservative R->N: basic guanidinium -> neutral polar amide; formal charge +1->0',
        209: 'Identical L: aliphatic hydrophobe, neutral',
        210: 'Ambiguous gap or nonconservative V->D: aliphatic hydrophobe -> acidic carboxylate; formal charge 0->-1 if D'}
    hotspots = [row(h) | counts(neighbours(h)) | dict(chemistry=chemistry[h]) for h in PATCH]
    # Fit only unambiguous sequence-derived domain pairs as a geometric check.
    paired = [(h, options[h][0]) for h in DOMAIN if len(options[h]) == 1 and options[h][0]]
    x = np.array([af[m]['CA'].coord for h, m in paired]); y = np.array([rr[h]['CA'].coord for h, m in paired])
    xc, yc = x.mean(0), y.mean(0)
    u, _, vt = np.linalg.svd((x-xc).T @ (y-yc))
    fix = np.eye(3); fix[2, 2] = np.linalg.det(u@vt); rotation = u@fix@vt
    deviations = np.linalg.norm((x-xc)@rotation+yc-y, axis=-1)
    # Solvent exposure in native dimer incl. observed glycans, 1.4 A probe.
    ShrakeRupley(n_points=200).compute(native, level='R')
    maxima = residue_max_acc['Wilke']
    rsa = {n: rr[n].sasa/maxima[rr[n].resname] for n in rr}
    candidates = []
    for n in DOMAIN:
        if rsa[n] < .2 or not row(n)['identical']:
            continue
        surface = [h for h in DOMAIN if rsa[h] >= .2 and
                   np.linalg.norm(rr[h]['CA'].coord-rr[n]['CA'].coord) <= 10]
        invariant = [h for h in surface if row(h)['identical']]
        if len(invariant) < 5:
            continue
        env = sorted(set(h for p in invariant for h in [p]+neighbours(p)))
        changed = sum(not row(h)['identical'] for h in env)
        candidates.append(dict(center=row(n), conserved_surface=[row(h) | dict(rsa=rsa[h]) for h in invariant],
            substituted_surface=[row(h) for h in surface if h not in invariant],
            environment_changed=changed, environment_size=len(env),
            environment_conservation=1-changed/len(env),
            conserved_surface_sasa_A2=sum(rr[h].sasa for h in invariant),
            per_residue_neighbourhood=[dict(human_auth=h, **counts(neighbours(h))) for h in invariant]))
    candidates.sort(key=lambda p: (-p['environment_conservation'], -p['conserved_surface_sasa_A2'], p['center']['human_auth']))
    report = dict(inputs=[dict(path=str(p.relative_to(ROOT)), sha256=sha(p)) for p in paths],
        method=dict(alignment='Full-length global affine-gap BLOSUM62 from local UniProt flatfile sequences',
            gap_open=-10, gap_extend=-.5, score=alignments.score, optimal_alignments=len(maps),
            numbering='Human auth=UniProt; all 11qc A protein atom author/label pairs verified auth-label=121; mouse=UniProt/AF PDB',
            domain_human_auth=[189,383], neighbourhood='Other same-chain protein residues with any heavy atom within <=8 A of any hotspot heavy atom, human 11qc A; all domains included; hotspot excluded',
            patch_scan='Centres: identical apical residues RSA>=0.2; patch: identical RSA>=0.2 residues within CA10 A; >=5 residues. Rank by invariant fraction of union of their 8 A environments, then exposed conserved area. RSA from native dimer+observed glycans, 1.4 A probe, 200 points, Wilke maxima. Conservation-only ranking; no binding claim.'),
        apical_identity=dict(identical=sum(r['identical'] for r in rows), human_length=len(rows),
                            fraction=sum(r['identical'] for r in rows)/len(rows)),
        ectodomain_identity=dict(human_auth=[121,760],mouse_assay=[123,763],
            human_positions=640,
            identical=sum(row(h)['identical'] for h in range(121,761)),
            fraction_human_denominator=sum(row(h)['identical'] for h in range(121,761))/640,
            note='Sequence-derived identity on this full ectodomain span, not the challenge approximately86% claim for its modelled region'),
        hotspots=hotspots, fit=dict(unambiguous_pairs=len(paired), CA_rmsd_A=float(np.sqrt(np.mean(deviations**2)))),
        mouse_trim=dict(residues=len(trim), hotspot_coverage={str(h): [m for m in options[h] if m and m in trim] for h in PATCH},
            warning='His710-centred discontinuous crop does not cover this epitope; full AF mouse model used'),
        top_conserved_surface_patches=candidates[:10],
        verdict='Partially conserved only: 2/5 identities; two charge-losing substitutions and one ambiguous nonconserved indel/substitution. As-is cross-reactivity remains unestablished.',
        limitations=['Sequence/geometry and exposure only, no binding or cross-reactivity measurement',
                     'Human local neighbourhood, not a mouse side-chain recognition model',
                     'Observed glycan atoms only, no complete assay glycans or tags',
                     'All equally optimal alignments retained; human210 cannot be assigned uniquely from these sequences/settings'])
    (OUT/'measurements.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(dict(hotspots=[{k:r[k] for k in ('human_auth','mouse_uniprot_options','mouse_aa_options','neighbour_count','substitution_count_range','deletion_count_range')} for r in hotspots], best_patch=candidates[0], apical_identity=report['apical_identity']), indent=2))


if __name__ == '__main__':
    measure()
