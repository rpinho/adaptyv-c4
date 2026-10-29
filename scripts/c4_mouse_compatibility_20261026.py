#!/usr/bin/env python3
"""CPU-only fixed-pose mouse compatibility, joined on SHA256(sequence).

Does not predict binding, optimize a pose, run a model, or rank designs.
All optimal canonical-sequence alignments and both prespecified fits are retained.
"""
import argparse
import csv
import hashlib
import json
import subprocess
from pathlib import Path

import numpy as np
from Bio import Align, SwissProt
from Bio.Align import substitution_matrices
from Bio.Data.PDBData import protein_letters_3to1 as AA
from Bio.PDB import PDBParser

ROOT = Path(__file__).resolve().parents[1]
CSV = ROOT / 'results/c4_ranked_v4_apical_L6_20261025/c4_acid_ranked_v4_apical_L6_20261025.csv'
PDBS = ROOT / 'results/c4_apical_pdb_audit_20261025/all343'
DOMAIN = range(189, 384)
BACKBONE = {'N', 'CA', 'C', 'O'}
CHARGE = {'K': 1, 'R': 1, 'D': -1, 'E': -1}


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def seqsha(s):
    return hashlib.sha256(s.encode('ascii')).hexdigest()


def residues(chain):
    return [r for r in chain if r.id[0] == ' ' and r.resname in AA]


def sequence(chain):
    return ''.join(AA[r.resname] for r in residues(chain))


def atoms(rs):
    return [(r.id[1], a.name, a.coord.astype(float)) for r in rs for a in r
            if a.element not in ('H', 'D')]


def fit(x, y):
    """Proper least-squares row-vector rigid transformation, x -> y."""
    xc, yc = x.mean(0), y.mean(0)
    u, _, vt = np.linalg.svd((x-xc).T @ (y-yc))
    fix = np.eye(3)
    fix[2, 2] = np.linalg.det(u @ vt)
    rot = u @ fix @ vt
    shift = yc - xc @ rot
    dev = np.linalg.norm(x @ rot + shift - y, axis=1)
    return rot, shift, float(np.sqrt(np.mean(dev**2)))


def distances(binder, target):
    """Chunk atom pairs to bound memory; return exact pair-level minima/counts."""
    bxyz = np.array([a[2] for a in binder])
    txyz = np.array([a[2] for a in target])
    tmin = np.full(len(target), np.inf)
    bmin = np.full(len(binder), np.inf)
    clash_pairs, backbone_pairs, clashes = 0, 0, []
    mind = np.inf
    bbdmin = np.inf
    bbmask = np.array([a[1] in BACKBONE for a in target])
    for start in range(0, len(binder), 128):
        d2 = ((bxyz[start:start+128, None, :] - txyz[None, :, :])**2).sum(-1)
        d = np.sqrt(d2)
        mind = min(mind, float(d.min()))
        bbdmin = min(bbdmin, float(d[:, bbmask].min()))
        tmin = np.minimum(tmin, d.min(0))
        bmin[start:start+len(d)] = d.min(1)
        ii, jj = np.where(d < 2.)
        clash_pairs += len(ii)
        backbone_pairs += int(bbmask[jj].sum())
        for i, j in zip(ii, jj):
            clashes.append(dict(binder_res=binder[start+i][0], binder_atom=binder[start+i][1],
                                target_res=target[j][0], target_atom=target[j][1],
                                distance_A=float(d[i,j]), target_backbone=bool(bbmask[j])))
    perres = {}
    for (n, _, _), d in zip(target, tmin):
        perres[n] = min(perres.get(n, np.inf), float(d))
    return dict(min_heavy_A=float(mind), min_target_backbone_A=float(bbdmin),
                pairs_lt2A=clash_pairs, target_backbone_pairs_lt2A=backbone_pairs,
                contact_residues={str(c): sorted(n for n,d in perres.items() if d<=c)
                                  for c in (4.,4.5,5.)},
                per_target_residue_min_A=perres, clash_pairs=clashes,
                binder_contact_residues=sorted(set(binder[i][0] for i,d in enumerate(bmin) if d<=4.5)))


def run(out):
    out.mkdir(parents=True, exist_ok=True)
    assert not (out/'measurements.json').exists(), 'Refuse to overwrite a measurement'
    now = subprocess.check_output(['date','-u','+%Y-%m-%dT%H:%M:%SZ'],text=True).strip()
    paths = [CSV, ROOT/'data/P02786.txt', ROOT/'data/Q62351.txt', ROOT/'data/AF-Q62351.pdb',
             ROOT/'results/decisions/20261026-c4-selection-rule-v1.md', Path(__file__)]
    fingerprints = {str(p.relative_to(ROOT)): sha(p) for p in paths}
    with paths[1].open() as f:
        human = SwissProt.read(f).sequence
    with paths[2].open() as f:
        mouse = SwissProt.read(f).sequence
    aligner = Align.PairwiseAligner(mode='global')
    matrix = substitution_matrices.load('BLOSUM62')
    aligner.substitution_matrix = matrix
    aligner.open_gap_score, aligner.extend_gap_score = -10, -.5
    aligns = aligner.align(human, mouse)
    maps = [{int(h)+1:int(m)+1 for h,m in zip(*a.indices) if h>=0 and m>=0} for a in aligns]
    opts = {h:sorted({m.get(h) for m in maps}, key=lambda x:-1 if x is None else x) for h in DOMAIN}
    parser = PDBParser(QUIET=True)
    mc = parser.get_structure('mouse', paths[3])[0]['A']
    assert sequence(mc) == mouse
    unambig = [h for h in DOMAIN if len(opts[h])==1 and opts[h][0] is not None]
    # Restrict mouse coordinates to the alignment-derived apical span, no full-domain extension.
    mouse_ids = sorted({m[h] for m in maps for h in DOMAIN if h in m})
    mt = atoms([mc[n] for n in mouse_ids])
    rows = list(csv.DictReader(CSV.open()))
    requested = {seqsha(r['sequence']):r for r in rows}
    assert len(rows)==len(requested)==20
    # Discover molecules by extracting actual chain A; names never join evidence.
    candidates = {s:[] for s in requested}
    for p in sorted(PDBS.glob('*.pdb')):
        model = parser.get_structure('complex',p)[0]
        if 'A' not in model:
            continue
        s = seqsha(sequence(model['A']))
        if s in candidates:
            candidates[s].append(p)
    assert all(len(v)==1 for v in candidates.values()), {k:len(v) for k,v in candidates.items()}
    raw, table = [], []
    for row in rows:
        sid = seqsha(row['sequence'])
        p = candidates[sid][0]
        fingerprints[str(p.relative_to(ROOT))] = sha(p)
        model = parser.get_structure('complex',p)[0]
        bc, hc = model['A'], model['B']
        assert sequence(bc) == row['sequence']
        assert sequence(hc) == human[188:383]
        hr = {r.id[1]+188:r for r in residues(hc)}
        assert set(hr)==set(DOMAIN)
        bt = atoms(residues(bc))
        ht = atoms(residues(hc))
        baseline = distances(bt, ht)
        contact = [n+188 for n in baseline['contact_residues']['4.5']]
        assert contact, sid
        contactrows=[]
        for h in contact:
            aa = human[h-1]
            maa = [mouse[n-1] if n is not None else '-' for n in opts[h]]
            contactrows.append(dict(human_uniprot=h,human_aa=aa,mouse_options=opts[h],
                mouse_aa_options=maa, invariant=all(a==aa for a in maa),
                mapping_ambiguous=len(opts[h])>1,
                blosum62_options=[float(matrix[aa,a]) if a!='-' else None for a in maa],
                charge_reversal=any(CHARGE.get(aa,0)*CHARGE.get(a,0)<0 for a in maa),
                charge_loss=any(CHARGE.get(aa,0)!=0 and CHARGE.get(a,0)==0 for a in maa),
                min_human_distance_A=baseline['per_target_residue_min_A'][h-188]))
        # Local fitting uses only sequence-unambiguous pairs, selected by human geometry.
        # 12 A to any contacted human target CA, never selected by resulting mouse clash count.
        cc = np.array([hr[h]['CA'].coord for h in contact])
        local = [h for h in unambig if np.linalg.norm(hr[h]['CA'].coord-cc,axis=1).min()<=12.]
        assert len(local)>=10
        fits = {}
        for label, hs in (('domain',unambig),('local',local)):
            x = np.array([mc[opts[h][0]]['CA'].coord for h in hs],dtype=float)
            y = np.array([hr[h]['CA'].coord for h in hs],dtype=float)
            rot,shift,rms = fit(x,y)
            fitted = [(n,a,xyz@rot+shift) for n,a,xyz in mt]
            metrics = distances(bt,fitted)
            retention=[]
            for mapping in maps:
                paired=[h for h in contact if h in mapping]
                retained=sum(metrics['per_target_residue_min_A'][mapping[h]]<=4.5 for h in paired)
                retention.append(dict(mapped=len(paired),retained=retained,total_human_contacts=len(contact),
                                      fraction=retained/len(contact)))
            fits[label] = dict(n_fit_pairs=len(hs),fit_CA_rmsd_A=rms,
                               fit_human_positions=hs, human_contacts_retained_by_alignment=retention,**metrics)
        # Self-transfer control: same receptor to itself; all atom distances preserved.
        ca = np.array([hr[h]['CA'].coord for h in unambig],dtype=float)
        r,t,self_rms = fit(ca,ca)
        selfmetrics = distances(bt, [(n,a,xyz@r+t) for n,a,xyz in ht])
        assert self_rms<1e-6
        assert abs(selfmetrics['min_heavy_A']-baseline['min_heavy_A'])<1e-6
        assert selfmetrics['contact_residues']==baseline['contact_residues']
        # Positive control: place one binder heavy atom EXACTLY onto a mouse backbone atom.
        mn = next(a for a in mt if a[1]=='CA')
        synthetic = distances([(1,'CA',mn[2].copy())],mt)
        assert synthetic['target_backbone_pairs_lt2A']>=1 and synthetic['min_heavy_A']==0
        far = distances([(1,'CA',mn[2]+np.array([10000.,10000.,10000.]))],mt)
        assert not far['contact_residues']['5.0'] and not far['pairs_lt2A']
        invariant=sum(c['invariant'] for c in contactrows)
        reversed_n=sum(c['charge_reversal'] for c in contactrows)
        both = all(f['target_backbone_pairs_lt2A']>0 for f in fits.values())
        # Mechanical labels only. Panel-reviewed interpretation is written separately.
        label = ('fixed_pose_backbone_conflict_both_fits' if both else
                 'fit_sensitive_backbone_conflict' if any(f['target_backbone_pairs_lt2A'] for f in fits.values()) else
                 'no_backbone_conflict_detected')
        rec=dict(sequence_sha256=sid,name=row['name'],source_pdb=str(p.relative_to(ROOT)),
                 human_contact_positions=contact,contact_residue_evidence=contactrows,
                 human_baseline=baseline, fits=fits, mechanical_label=label,
                 controls=dict(human_self_fit_rmsd_A=self_rms,
                    human_self_min_distance_delta_A=abs(selfmetrics['min_heavy_A']-baseline['min_heavy_A']),
                    human_self_contacts_exact=True,positive_injected_overlap_pairs=synthetic['pairs_lt2A'],
                    positive_backbone_overlap_pairs=synthetic['target_backbone_pairs_lt2A'],
                    far_control_contacts=len(far['contact_residues']['5.0'])))
        raw.append(rec)
        tr=dict(sequence_sha256=sid,name=row['name'],human_contact_count=len(contact),
                invariant_contact_count=invariant,invariant_fraction=round(invariant/len(contact),4),
                changed_contact_count=len(contact)-invariant,
                ambiguous_contact_count=sum(c['mapping_ambiguous'] for c in contactrows),
                charge_reversal_count=reversed_n,charge_loss_count=sum(c['charge_loss'] for c in contactrows),
                human_min_heavy_A=round(baseline['min_heavy_A'],3),
                human_backbone_pairs_lt2A=baseline['target_backbone_pairs_lt2A'],
                mechanical_label=label,method='sequence_all_optima_plus_fixed_pose_domain_and_local_CA_fits',
                binding_claim='unmeasured',source_pdb=str(p.relative_to(ROOT)))
        for label,f in fits.items():
            for key in ('fit_CA_rmsd_A','min_heavy_A','min_target_backbone_A'):
                tr[label+'_'+key]=round(f[key],3)
            for key in ('pairs_lt2A','target_backbone_pairs_lt2A'):
                tr[label+'_'+key]=f[key]
            ps=[v['fraction'] for v in f['human_contacts_retained_by_alignment']]
            tr[label+'_contact_retention_min']=round(min(ps),4)
            tr[label+'_contact_retention_max']=round(max(ps),4)
        table.append(tr)
    # Explicit mapping/chemistry control at the previously flagged species-sensitive patch.
    assert opts[205]==[208] and human[204]=='K' and mouse[207]=='S'
    assert float(matrix['K','S'])==0 and CHARGE['K']!=CHARGE.get('S',0)
    assert len(maps)==6 and len(unambig)==192
    for pstr,expected in fingerprints.items():
        assert sha(ROOT/pstr)==expected, 'Input drift: '+pstr
    report=dict(created_at_utc=now,inputs=fingerprints,
        method=dict(scope='195-aa human apical surrogate only; mouse aligned apical domain only',
          join='sha256(actual extracted binder chain A sequence), no name join',
          alignment='global BLOSUM62, open -10, extension -0.5; retain all six optimal mappings',
          optimal_alignment_count=len(maps),score=aligns.score,unambiguous_CA_pairs=len(unambig),
          domain_fit='all 192 unambiguous apical CA pairs, proper Kabsch',
          local_fit='unambiguous CA pairs within 12 A of any human contact target CA; >=10 pairs',
          contacts='any heavy atom <=4.5 A; also raw 4.0/5.0 A sets; clashes strictly <2.0 A',
          backbone_clash='any binder heavy atom vs mouse target N,CA,C,O; no hydrogen atoms',
          chemistry='invariance across every optimal alignment, BLOSUM62, formal D/E=-1 K/R=+1; His unassigned',
          caveat='No docking, relaxation, protonation energy, KD, or binding prediction; static mouse AF rotamers'),
        mapping_options={str(k):v for k,v in opts.items()},
        controls=dict(all_passed=True,exact_human_self_transfers=20,positive_overlap_controls_fired=20,
                      far_controls_silent=20,K205_to_S208_charge_loss_detected=True,
                      hashes_unique=20,structures_unique_by_actual_sequence=20), designs=raw)
    (out/'measurements.json').write_text(json.dumps(report,indent=2)+'\n')
    with (out/'mouse_compatibility.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(table[0]));w.writeheader();w.writerows(table)
    (out/'alignment_optima.txt').write_text('\n\n'.join(str(a) for a in aligns)+'\n')
    print(json.dumps(dict(rows=len(table),controls=report['controls'],
        labels={s:sum(r['mechanical_label']==s for r in table) for s in sorted({r['mechanical_label'] for r in table})},
        invariant_fraction_range=[min(r['invariant_fraction'] for r in table),max(r['invariant_fraction'] for r in table)]),indent=2))


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--out',type=Path,required=True)
    run(ap.parse_args().out)
