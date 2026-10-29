#!/usr/bin/env python
"""Continue accepted GPU scaffolds into bounded CPU switch validation.

This queue runs after the pod is verified stopped. Only native-sequence, full
human dimer complexes are admitted; no lysinated trajectory can enter pKa scoring.
"""
import argparse
import csv
import json
from pathlib import Path
import time
import numpy as np
from Bio.PDB import PDBParser,MMCIFParser
from Bio.Data.PDBData import protein_letters_3to1 as AA
from cif2pdb import convert
import install_switch as install
from campaign_admission import unique_parents
from run_context_ensemble import ROOT, residues, atoms, fit, min_dist, write_protein, model_triplet, dump, sha


def expected_target():
    text=(ROOT/'data/P02786.txt').read_text()
    full=''.join(text.split('\nSQ   ',1)[1].split('\n',1)[1].split('//')[0].split())
    return full[121:760]


def canonicalize(source,destination):
    model=PDBParser(QUIET=True).get_structure('x',str(source))[0]
    restored=False;partner_fit=None
    if set(c.id for c in model)==set('AB'):
        expected=expected_target()
        aa=residues(model['A']);bb=residues(model['B'])
        if ''.join(AA[r.resname] for r in aa)!=expected or not 30<=len(bb)<=150:
            raise ValueError('monomer fallback needs complete native target A and binder B')
        template=MMCIFParser(QUIET=True).get_structure('11qc',ROOT/'data/11qc.cif')[0]
        rotation,translation,partner_fit=fit([r['CA'].coord for r in residues(template['A'])],
                                             [r['CA'].coord for r in aa])
        if partner_fit>2:raise ValueError('opposite-subunit restoration exceeds 2A target fit')
        binder=model['B'];model.detach_child('B');binder.id='C'
        partner=template['B'].copy()
        for atom in partner.get_atoms():atom.coord=atom.coord@rotation+translation
        model.add(partner);model.add(binder);restored=True
    if set(c.id for c in model)!=set('ABC'):
        raise ValueError('expected full target A/B and one binder C')
    expected=expected_target()
    for cid in 'AB':
        rr=residues(model[cid]);sequence=''.join(AA[r.resname] for r in rr)
        if sequence!=expected:
            raise ValueError('target sequence differs from full native 122-760; lysinated/cropped/mutated input forbidden')
        # Rebuild IDs without temporary collisions in Biopython's child dictionary.
        chain=model[cid]
        for r in rr:chain.detach_child(r.id)
        for n,r in enumerate(rr,122):r.id=(' ',n,' ');chain.add(r)
        for a,b in zip(rr,rr[1:]):
            if np.linalg.norm(a['C'].coord-b['N'].coord)>1.9:
                raise ValueError('native target has a broken peptide bond')
    binder=residues(model['C'])
    if not 30<=len(binder)<=150:raise ValueError('unexpected binder architecture')
    bb=lambda rr: np.array([a.coord for r in rr for a in r if a.name in ('N','CA','C','O')])
    target=[r for cid in 'AB' for r in residues(model[cid])]
    dimer_min=min_dist(bb(residues(model['A'])),bb(residues(model['B'])))
    if dimer_min<1.8:raise ValueError('reconstructed target subunits overlap at backbone level')
    bb_min=min_dist(bb(binder),bb(target))
    if bb_min<2:raise ValueError(f'binder/full-target backbone clash {bb_min:.3f} A')
    template=MMCIFParser(QUIET=True).get_structure('glycans',ROOT/'data/11qc.cif')[0]
    rot,trans,glycan_fit=fit([r['CA'].coord for r in residues(template['A'])],
                            [r['CA'].coord for r in residues(model['A'])])
    glycan_pts=np.array([atom.coord for chain in template for residue in chain
                         if residue.resname=='NAG' for atom in residue])@rot+trans
    glycan_distance=min_dist(atoms(binder),glycan_pts)
    if glycan_distance<2:raise ValueError('binder overlaps resolved NAG atoms after context alignment')
    write_protein(model,destination)
    return dict(source=str(source),source_sha256=sha(source),canonical_sha256=sha(destination),
                binder_residues=len(binder),minimum_binder_target_backbone_distance=bb_min,
                opposite_subunit_restored=restored,opposite_subunit_fit_rmsd=partner_fit,
                target_subunit_backbone_min_distance=dimer_min,
                observed_NAG_min_distance=glycan_distance,glycan_template_fit_rmsd=glycan_fit,
                target_sequence_verified=True,
                limitations=['glycans and assay tag absent','human only; mouse cross-reactivity unverified'])


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--out',default='results/gpu_followup')
    ap.add_argument('--source',default='results/gpu_campaign/mirrored');ap.add_argument('--seeds',type=int,default=5)
    args=ap.parse_args();out=Path(args.out).resolve();out.mkdir(exist_ok=True,parents=True)
    if (out/'status.json').exists():raise SystemExit('fresh follow-up output required')
    status=dict(status='preflight',started=time.time(),models_scored=0,designs=[],errors=[],new_gpu_spend_usd=0,
                script_sha256=sha(__file__),max_parents=2,max_pairs_per_parent=1,seeds=args.seeds,
                interpretation='CPU structural/PROPKA proxies; not experimentally validated binding')
    dump(out/'status.json',status)
    ranked=Path(args.source).resolve()/'3_Ranked'
    table=ranked/'!_Ranked.csv'
    if not table.exists():
        status['status']='no_accepted_scaffolds';dump(out/'status.json',status);return
    with table.open() as fh:records=unique_parents(list(csv.DictReader(fh)),limit=2)
    if not records:
        status['status']='no_accepted_scaffolds';dump(out/'status.json',status);return
    deadline=time.monotonic()+7200
    install.init()
    for index,row in enumerate(records,1):
        design=row['design'];record=dict(design=design,triplets=[])
        try:
            source=ranked/(design+'.cif');raw=out/(design+'_raw.pdb');convert(str(source),str(raw))
            canonical=out/(design+'_canonical.pdb');record['audit']=canonicalize(raw,canonical)
            plan=install.make_plan(canonical,'C','AB',max_pairs=1)
            dump(out/(design+'_plan.json'),plan)
            if not plan['pairs']:
                record['status']='no_eligible_pair';status['designs'].append(record)
                dump(out/'status.json',status);continue
            parent=dict(parent_index=index,backbone=design,sites=plan['pairs'][0]['sites'],context_pdb=str(canonical))
            for seed in range(args.seeds):
                if time.monotonic()>=deadline:break
                triplet=model_triplet(parent,202610220+seed,out,deadline)
                record['triplets'].append(triplet);status['models_scored']+=len(triplet['variants'])
                dump(out/'status.json',status|{'active_design':record})
            valid=[t for t in record['triplets'] if t.get('hh_minus_qq') is not None]
            supporting=[t for t in valid if t['hh_minus_qq']<0 and t.get('neutral_retention')]
            record.update(status='modeled',valid_replicas=len(valid),supporting_replicas=len(supporting),
                          engineering_advance=bool(len(valid)==5 and len(supporting)>=4))
        except Exception as exc:
            record.update(status='failed',error=repr(exc));status['errors'].append(dict(design=design,error=repr(exc)))
        status['designs'].append(record);dump(out/'status.json',status)
    status.update(status='completed' if time.monotonic()<deadline else 'time_cap_partial',finished=time.time())
    dump(out/'status.json',status)


if __name__=='__main__':main()
