#!/usr/bin/env python
"""Bounded CPU robustness screen in resolved 11QC dimer context.

Fixed backbone, local target/all-binder side-chain sampling; separately repacked
unbound binder and receptor. These replicas are sensitivity checks, not an
equilibrium ensemble or a calibrated affinity predictor. No GPU allocation.
"""
import argparse
import csv
import hashlib
import io
import json
import math
import os
from pathlib import Path
import sys
import time

import numpy as np
from Bio import Align
from Bio.Data.PDBData import protein_letters_3to1 as AA
from Bio.PDB import MMCIFParser, PDBParser, PDBIO, Select

import install_switch as install
from ph_linkage import run_propka, linkage

ROOT = Path(__file__).resolve().parents[1]


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def dump(path, obj):
    install.save_json(Path(path), obj)


def residues(chain):
    return [r for r in chain if r.id[0] == ' ' and r.resname in AA]


def atoms(rr):
    return np.array([a.coord for r in rr for a in r if a.element != 'H'])


def fit(mobile, fixed):
    x, y = np.asarray(mobile), np.asarray(fixed)
    xc, yc = x.mean(0), y.mean(0)
    u, _, vt = np.linalg.svd((x-xc).T @ (y-yc))
    correction = np.eye(3)
    correction[2, 2] = np.linalg.det(u @ vt)
    rot = u @ correction @ vt
    trans = yc - xc @ rot
    rmsd = float(np.sqrt(np.mean(np.sum((x @ rot + trans-y)**2, axis=1))))
    return rot, trans, rmsd


def min_dist(x, y):
    best = math.inf
    for chunk in np.array_split(x, max(1, len(x)//200)):
        best = min(best, float(np.linalg.norm(chunk[:, None]-y[None], axis=-1).min()))
    return best


def compare_binders(a, b):
    aligner = Align.PairwiseAligner(mode='global', match_score=2, mismatch_score=-1,
                                    open_gap_score=-3, extend_gap_score=-0.5)
    sa, sb = [''.join(AA[r.resname] for r in rr) for rr in (a,b)]
    al = aligner.align(sa, sb)[0]
    pairs = [(i,j) for (lo,hi),(lo2,hi2) in zip(*al.aligned)
             for i,j in zip(range(lo,hi),range(lo2,hi2))]
    identity = sum(sa[i] == sb[j] for i,j in pairs)/max(len(sa),len(sb))
    coverage = len(pairs)/max(len(sa),len(sb))
    _,_,rmsd = fit([a[i]['CA'].coord for i,j in pairs], [b[j]['CA'].coord for i,j in pairs])
    return dict(identity=identity, aligned_coverage=coverage, ca_rmsd=rmsd,
                similar=bool((identity >= .7 and coverage >= .7) or
                             (rmsd < 2 and coverage >= .8)))


class Protein(Select):
    def accept_residue(self, r):
        return r.id[0] == ' ' and r.resname in AA
    def accept_atom(self, a):
        return a.element != 'H' and a.altloc in (' ', 'A')


def write_protein(model, path):
    writer=PDBIO(); writer.set_structure(model); writer.save(str(path), Protein())


def prepare_inputs(out):
    priorities = json.loads((ROOT/'results/next_campaign/parent_priorities.json').read_text())
    old = json.loads((ROOT/'results/switch_install/campaign/campaign.json').read_text())
    parser=PDBParser(QUIET=True)
    models={p['backbone']: parser.get_structure('x',p['source_path'])[0] for p in priorities}
    comparisons=[]; chosen=[]
    for p in priorities:
        duplicate=False
        for other in chosen:
            comparison=compare_binders(residues(models[p['backbone']]['C']),
                                       residues(models[other['backbone']]['C']))
            comparisons.append(dict(a=p['backbone'],b=other['backbone'],**comparison))
            duplicate |= comparison['similar']
        if not duplicate and len(chosen)<5:
            chosen.append(dict(p))
    # A separate native-H24 calibration, using its best complete matched pair.
    calibration=next(d for d in old['designs'] if 'l58_9f94b53776050aad' in d['backbone'])
    vv=[v for v in calibration['variants'] if v['kind']=='switch' and v.get('linkage')
        and any(q['pair']==v['pair'] and q['kind']=='control' and q.get('linkage')
                for q in calibration['variants'])]
    if vv:
        v=min(vv,key=lambda x:x['linkage']['L'])
        chosen.append(dict(backbone=calibration['backbone'],pair=v['pair'],
                           source_path=calibration['source_pdb'],calibration=True))
    template=MMCIFParser(QUIET=True).get_structure('11qc',ROOT/'data/11qc.cif')[0]
    human_text=(ROOT/'data/P02786.txt').read_text()
    human=''.join(human_text.split('\nSQ   ',1)[1].split('\n',1)[1].split('//')[0].split())
    target_sequence_mismatches=[]; topology=[]
    for c in template:
        rr=residues(c)
        for r in rr:
            if AA[r.resname] != human[r.id[1]-1]:
                target_sequence_mismatches.append([c.id,r.id[1],r.resname,human[r.id[1]-1]])
        gaps=[]
        for r,s in zip(rr,rr[1:]):
            d=float(np.linalg.norm(r['C'].coord-s['N'].coord))
            if s.id[1]!=r.id[1]+1 or d>1.9:
                gaps.append([r.id[1],s.id[1],d])
        topology.append(dict(chain=c.id,residues=len(rr),first=rr[0].id[1],last=rr[-1].id[1],gaps=gaps))
    if target_sequence_mismatches or any(t['gaps'] for t in topology):
        raise ValueError('11QC sequence/connectivity needs explicit repair before modelling')
    # Map the transferrin contact footprint by receptor sequence, not chain names.
    tf=MMCIFParser(QUIET=True).get_structure('tf',ROOT/'data/3s9l.cif')[0]
    receptors=[]; partners=[]
    for c in tf:
        rr=residues(c); seq=''.join(AA[r.resname] for r in rr)
        if len(seq)<8: continue
        probe=seq[len(seq)//2:len(seq)//2+25]
        (receptors if probe in human else partners).append(rr)
    footprint=set()
    if receptors and partners:
        pts=atoms([r for rr in partners for r in rr])
        for rr in receptors:
            for r in rr:
                if min_dist(atoms([r]),pts)<5:
                    footprint.add(r.id[1])
    target_path=out/'human_dimer.pdb';write_protein(template,target_path)
    target_lines=[l for l in target_path.read_text().splitlines(True) if l.startswith(('ATOM','TER'))]
    rows=[]
    for index,p in enumerate(chosen):
        d=next(d for d in old['designs'] if d['backbone']==p['backbone'])
        manifest=next(Path(d['directory']).glob('*variants.json'))
        original=json.loads(manifest.read_text())
        sites=original['pairs'][p['pair']-1]['sites']
        m=parser.get_structure('s',p['source_path'])[0]
        moving=[]; fixed=[]
        for cid in 'AB':
            for r in residues(m[cid]):
                f=template[cid][r.id]
                if r.resname!=f.resname: raise ValueError('target mutation in source')
                moving.append(r['CA'].coord);fixed.append(f['CA'].coord)
        rot,trans,rmsd=fit(moving,fixed)
        if rmsd>2: raise ValueError('target alignment exceeds 2 A; context transfer requires refinement')
        binder=m['C'].copy()
        for a in binder.get_atoms(): a.coord=a.coord @ rot+trans
        bp=out/f'{index+1:02d}_binder.pdb';write_protein(binder,bp)
        context=out/f'{index+1:02d}_context.pdb'
        context.write_text(''.join(target_lines)+''.join(l for l in bp.read_text().splitlines(True)
                                                       if l.startswith(('ATOM','TER')))+'END\n')
        bpts=atoms(residues(binder)); contacts=[]; added_clashes=[]
        original_keys={(c.id,r.id[1]) for c in m if c.id in 'AB' for r in residues(c)}
        for c in template:
            for r in residues(c):
                dist=min_dist(atoms([r]),bpts)
                if dist<5: contacts.append(dict(chain=c.id,resnum=r.id[1],resname=r.resname,distance=dist))
                if dist<2 and (c.id,r.id[1]) not in original_keys:
                    added_clashes.append(dict(chain=c.id,resnum=r.id[1],distance=dist))
        glycan=[]
        for c in template:
            for n in [251,317,727]:
                glycan.append(dict(chain=c.id,site=n,distance=min_dist(atoms([c[n]]),bpts)))
        observed_glycans=[a.coord for c in template for r in c if r.resname=='NAG' for a in r]
        cend=residues(binder)[-1]
        cterminal=min_dist(atoms([cend]),atoms([r for c in template for r in residues(c)]))
        audit=dict(alignment_ca_rmsd=rmsd,aligned_residues=len(moving),contacts=contacts,
                   alignment_policy='all source target CAs; <=2 A transfer; full template replaces refolded target',
                   newly_restored_target_clashes_under_2A=added_clashes,glycan_attachment_distances=glycan,
                   observed_NAG_min_distance=min_dist(np.array(observed_glycans),bpts),
                   transferrin_footprint_overlap=sorted({r['resnum'] for r in contacts}&footprint),
                   binder_Cterminal_residue_target_min_distance=cterminal,
                   glycan_model='attachment/observed-NAG proximity audit only; glycans omitted from Rosetta/PROPKA',
                   terminal_model='resolved receptor starts122; assay121+GP scar absent; binder tag absent; native charged termini',
                   source_sha256=sha(p['source_path']),context_sha256=sha(context))
        rows.append(dict(p,parent_index=index+1,sites=sites,context_pdb=str(context),audit=audit))
    result=dict(parents=rows,comparisons=comparisons,target_topology=topology,
                target_sequence_mismatches=target_sequence_mismatches,transferrin_footprint=sorted(footprint),
                template_sha256=sha(ROOT/'data/11qc.cif'),script_sha256=sha(__file__),
                selection='greedy priority with aligned sequence/CA redundancy checks; not universal fold clustering')
    dump(out/'inputs.json',result)
    return result


def local_prepare(pose,sf,indices):
    """Repack only declared sidechains; backbone/jumps and other residues fixed."""
    pr=install.pyrosetta
    work=pose.clone();task=pr.standard_packer_task(work);task.restrict_to_repacking()
    task.or_include_current(True)
    selected=set(indices)
    mm=pr.rosetta.core.kinematics.MoveMap();mm.set_bb(False);mm.set_chi(False);mm.set_jump(False)
    for j in range(1,work.total_residue()+1):
        if j not in selected: task.nonconst_residue_task(j).prevent_repacking()
        else: mm.set_chi(j,True)
    pr.rosetta.protocols.minimization_packing.PackRotamersMover(sf,task).apply(work)
    pr.rosetta.protocols.minimization_packing.MinMover(mm,sf,'lbfgs_armijo_nonmonotone',.001,True).apply(work)
    return work


def extract(pose,keep):
    work=pose.clone()
    for j in range(work.total_residue(),0,-1):
        if work.pdb_info().chain(j) not in keep: work.delete_residue_slow(j)
    return work


def mutate_fixed_backbone(pose,index,aa):
    """MutateResidue idealizes carbonyl O; restore all input backbone heavy atoms."""
    pr=install.pyrosetta
    names=[n for n in ['N','CA','C','O','OXT'] if pose.residue(index).has(n)]
    coords={n:pose.residue(index).xyz(n) for n in names}
    name3=next(k for k,v in install.AA1.items() if v==aa)
    pr.rosetta.protocols.simple_moves.MutateResidue(index,name3).apply(pose)
    for name,xyz in coords.items():
        if not pose.residue(index).has(name):raise ValueError('mutation removed backbone atom '+name)
        pose.set_xyz(pr.rosetta.core.id.AtomID(pose.residue(index).atom_index(name),index),xyz)


def active_indices(pose,cutoff=10):
    binder=[i for i in range(1,pose.total_residue()+1) if pose.pdb_info().chain(i)=='C']
    centres=np.array([[pose.residue(i).nbr_atom_xyz()[k] for k in range(3)] for i in binder])
    active=[]
    for i in range(1,pose.total_residue()+1):
        xyz=np.array([pose.residue(i).nbr_atom_xyz()[k] for k in range(3)])
        if pose.pdb_info().chain(i)=='C' or np.linalg.norm(centres-xyz,axis=1).min()<cutoff:
            active.append(i)
    return active


def model_triplet(parent,seed,out,deadline):
    pr=install.pyrosetta;sf=pr.get_fa_scorefxn()
    directory=out/f"parent{parent['parent_index']:02d}_seed{seed}"
    directory.mkdir(exist_ok=True)
    result=dict(parent_index=parent['parent_index'],backbone=parent['backbone'],seed=seed,
                status='running',variants=[],started=time.time(),sites=parent['sites'])
    path=directory/'triplet.json';dump(path,result)
    base=pr.pose_from_pdb(parent['context_pdb'])
    if base.num_chains()!=3 or base.total_residue()!=1278+sum(
            base.pdb_info().chain(i)=='C' for i in range(1,base.total_residue()+1)):
        raise ValueError('full-context import changed chain/residue counts')
    positions=install.map_sites(base,parent['sites'])
    active=active_indices(base)
    # The same separately prepared free receptor is shared across WT/HH/QQ.
    pr.rosetta.numeric.random.rg().set_seed(seed)
    free_target=extract(base,set('AB'))
    target_active=[i for i in active if base.pdb_info().chain(i) in 'AB']
    free_target=local_prepare(free_target,sf,target_active)
    target_path=directory/'free_target.pdb';free_target.dump_pdb(str(target_path))
    et=float(sf(free_target));ft=run_propka(str(target_path),str(directory))
    if not ft: raise ValueError('free receptor PROPKA failed')
    result.update(target_active_residues=len(target_active),scorefunction=sf.get_name(),
                  pyrosetta_version=pr.version(),target_groups=len(ft))
    dump(path,result)
    for label,aas in [('WT',[s['wt'] for s in parent['sites']]),('HH',['H']*len(positions)),
                       ('QQ',['Q']*len(positions))]:
        if time.monotonic()>=deadline:
            result['status']='time_cap_partial';dump(path,result);return result
        pr.rosetta.numeric.random.rg().set_seed(seed)
        work=base.clone()
        for i,aa in zip(positions,aas):
            mutate_fixed_backbone(work,i,aa)
        work=local_prepare(work,sf,active)
        cxpath=directory/f'{label}_complex.pdb';work.dump_pdb(str(cxpath))
        fbpose=extract(work,{'C'})
        pr.rosetta.numeric.random.rg().set_seed(seed)
        fbpose=local_prepare(fbpose,sf,range(1,fbpose.total_residue()+1))
        fbpath=directory/f'{label}_free_binder.pdb';fbpose.dump_pdb(str(fbpath))
        dg=float(sf(work))-et-float(sf(fbpose))
        cb=run_propka(str(cxpath),str(directory));fb=run_propka(str(fbpath),str(directory))
        if not cb or not fb: raise ValueError(label+' PROPKA failed')
        free=dict(ft);free.update(fb)
        missing_bound=sorted(set(free)-set(cb));missing_free=sorted(set(cb)-set(free))
        total,terms=linkage(cb,free,6,7.4,0)
        missing_his=[s['resnum'] for s in parent['sites'] if label=='HH'
                     and ('C',s['resnum'],'HIS') not in (set(cb)&set(free))]
        expected_sidechains={(work.pdb_info().chain(i),work.pdb_info().number(i),
                              work.residue(i).name3())
                             for i in range(1,work.total_residue()+1)
                             if work.residue(i).name3() in {'HIS','ASP','GLU','LYS','ARG','TYR','CYS'}}
        missing_sidechains=sorted(expected_sidechains-(set(cb)&set(free)))
        # All titratable-group coverage must match, not only installed His.
        complete=not (missing_bound or missing_free or missing_his or missing_sidechains)
        row=dict(label=label,complex=str(cxpath),free_binder=str(fbpath),
                 complex_sha256=sha(cxpath),free_binder_sha256=sha(fbpath),
                 dG_local_repack_REU=dg,coverage_complete=complete,
                 missing_bound=missing_bound,missing_free=missing_free,missing_selected_his=missing_his,
                 missing_expected_sidechains=missing_sidechains,
                 matched_group_L=total,linkage_L=total if complete else None,rows=terms,
                 positive_budget=sum(t['contribution'] for t in terms if t['contribution']>0),
                 negative_budget=sum(t['contribution'] for t in terms if t['contribution']<0))
        result['variants'].append(row);dump(path,result)
        print(json.dumps(dict(event='model_scored',parent=parent['parent_index'],seed=seed,label=label,
                              dG=dg,L=row['linkage_L'],coverage_complete=complete)),flush=True)
    wt,hh,qq=result['variants']
    result.update(status='completed',finished=time.time(),
                  hh_minus_qq=(hh['linkage_L']-qq['linkage_L']) if hh['coverage_complete'] and qq['coverage_complete'] else None,
                  neutral_retention=bool(wt['dG_local_repack_REU']<0 and hh['dG_local_repack_REU']<0
                                        and hh['dG_local_repack_REU']-wt['dG_local_repack_REU']<=3))
    dump(path,result);return result


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--out',default='results/context_ensemble')
    ap.add_argument('--hours',type=float,default=2);ap.add_argument('--seeds',type=int,default=5)
    ap.add_argument('--parents',type=int,default=6);ap.add_argument('--prepare-only',action='store_true')
    args=ap.parse_args();out=Path(args.out).resolve();out.mkdir(exist_ok=True,parents=True)
    if (out/'status.json').exists(): raise SystemExit('fresh output directory required')
    started=time.time();deadline=time.monotonic()+args.hours*3600
    status=dict(status='preparing',pid=os.getpid(),started=started,wall_cap_seconds=args.hours*3600,
                max_parents=args.parents,seeds=args.seeds,models_scored=0,triplets=[],errors=[],
                new_gpu_spend_usd=0,script_sha256=sha(__file__),
                limitations=['fixed backbone','no modeled glycans or tag','empirical independent-site pKas',
                             'local repack binding-energy proxy, not IAM or experimental affinity',
                             'replica robustness, not equilibrium ensemble'])
    dump(out/'status.json',status)
    inputs=prepare_inputs(out)
    if args.prepare_only:
        status['status']='inputs_prepared';dump(out/'status.json',status);return
    install.init();status['status']='running';dump(out/'status.json',status)
    for s in range(args.seeds):
        for p in inputs['parents'][:args.parents]:
            if time.monotonic()>=deadline: break
            if p['audit']['newly_restored_target_clashes_under_2A']:
                status['errors'].append(dict(parent=p['parent_index'],reason='restored target steric clash; deferred'))
                dump(out/'status.json',status);continue
            try:
                r=model_triplet(p,202610210+s,out,deadline)
                status['triplets'].append(dict(parent_index=p['parent_index'],seed=r['seed'],
                                               status=r['status'],hh_minus_qq=r.get('hh_minus_qq'),
                                               neutral_retention=r.get('neutral_retention')))
                status['models_scored']+=len(r['variants'])
            except Exception as e:
                status['errors'].append(dict(parent=p['parent_index'],seed=s,error=repr(e)))
                print('ERROR',repr(e),flush=True)
            status['updated']=time.time();dump(out/'status.json',status)
        if time.monotonic()>=deadline: break
    planned=min(args.parents,len(inputs['parents']))*args.seeds
    status.update(status='completed' if len(status['triplets'])==planned and
                  all(t['status']=='completed' for t in status['triplets']) else 'partial',
                  finished=time.time(),elapsed_seconds=time.time()-started)
    dump(out/'status.json',status)


if __name__=='__main__': main()
