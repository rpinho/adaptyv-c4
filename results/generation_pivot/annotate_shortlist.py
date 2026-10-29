"""Attach measured context/composition and descriptive His-contact evidence.

Preserves the canonical ranker's order. No affinity, pH or eligibility score.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np
from Bio.PDB import PDBParser

CONTACT_DISTANCE_A=4.5  # Descriptive geometry convention, not a pH pass cutoff.
ROOT=Path(__file__).resolve().parents[2]
PRIORS={'tfr1_c200m_20261023':None,'tfr1_his200m_20261023':.1,
        'tfr1_his020_200m_20261023':.2,'tfr1_his030_200m_20261023':.3}


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def contacts(path):
    model=PDBParser(QUIET=True).get_structure('native_frame',str(path))[0]
    target=[a for cid in 'AB' for a in model[cid].get_atoms() if a.element!='H']
    groups=[target,[a for a in target if (a.parent.resname=='LYS' and a.name=='NZ') or
                     (a.parent.resname=='ARG' and a.name in ('NE','NH1','NH2'))],
                    [a for a in target if (a.parent.resname=='ASP' and a.name in ('OD1','OD2')) or
                     (a.parent.resname=='GLU' and a.name in ('OE1','OE2'))]]
    results=[]
    for r in model['C']:
        if r.resname!='HIS':continue
        names=('CG','ND1','CD2','CE1','NE2');missing=[name for name in names if name not in r]
        item=dict(position=r.id[1],missing_ring_atoms=missing)
        if not missing:
            for label,points,atoms in [('ring_to_target',names,groups[0]),
                                      ('nitrogen_to_basic',('ND1','NE2'),groups[1]),
                                      ('nitrogen_to_acidic',('ND1','NE2'),groups[2])]:
                distances=np.linalg.norm(np.array([r[name].coord for name in points])[:,None]-np.array([a.coord for a in atoms])[None],axis=-1)
                _,j=np.unravel_index(distances.argmin(),distances.shape);a=atoms[j]
                item[label+'_minimum_A']=float(distances.min())
                item[label+'_nearest_residue']=dict(chain=a.parent.parent.id,number=a.parent.id[1],name=a.parent.resname,atom=a.name)
        results.append(item)
    return results


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--shortlist',type=Path,required=True);parser.add_argument('--candidates',type=Path,nargs='+',required=True)
    parser.add_argument('--context',type=Path,required=True);parser.add_argument('--out',type=Path,required=True);args=parser.parse_args()
    if args.out.exists():raise FileExistsError('Preserve the existing review table')
    pool={}
    for path in args.candidates:
        with path.open(newline='') as stream:
            for r in csv.DictReader(stream):
                if r['name'] in pool:raise ValueError('Duplicate source identity')
                pool[r['name']]=r
    context=json.loads((args.context/'RESULTS.json').read_text());audits={r['source_id']:r for r in context['designs']}
    prior=json.loads((ROOT/'results/novelty_predictor/adaptive_completion/predictions.json').read_text())
    prior_by_hash={r['sequence_sha256']:r for r in prior['designs']}
    with args.shortlist.open(newline='') as stream:rows=list(csv.DictReader(stream))
    evidence=[]
    for row in rows:
        source=pool[row['source_id']];audit=audits[row['source_id']];assert row['sequence']==source['sequence']
        previous=prior_by_hash.get(hashlib.sha256(row['sequence'].encode()).hexdigest())
        row.update(source_arm=source['source_arm'],his_prior=PRIORS[source['source_arm']],n_his=row['sequence'].count('H'),
                   n_cys=row['sequence'].count('C'),hotspot_obedience=source['hotspot_obedience.hotspot_obedience'],
                   clash_fraction_1p5A=source['clash_metrics.frac_atoms_clashing_1p5A'],
                   bad_bond_fraction=source['bond_length_metrics.bond_length.frac_unreasonable'],
                   native_context_staging_pass=audit.get('full_native_staging_pass'),
                   target_fit_CA_RMSD_A=audit.get('unchanged_gate',{}).get('target_ca_rmsd'),
                   omitted_backbone_pairs_under2A=audit.get('unchanged_gate',{}).get('restored_backbone_pair_clashes_under2A'),
                   native_backbone_minimum_A=audit.get('staging',{}).get('minimum_binder_full_target_backbone_distance_A'),
                   observed_NAG_minimum_A=audit.get('staging',{}).get('observed_NAG_min_distance_A'),
                   prior_local_general_novelty_same_sequence=(previous or {}).get('local_predicted_level'),
                   novelty_evidence_scope='Earlier supplied-conformer local general-scale call for an identical sequence; current conformer and organizer route not rescored' if previous else 'New sequence; organizer novelty unmeasured',
                   pH_binding_evidence='unmeasured',organizer_acceptance='unmeasured',
                   his_rings_near_target_4p5A=None,his_nitrogens_near_basic_4p5A=None,his_nitrogens_near_acidic_4p5A=None)
        details=[]
        if audit.get('full_native_staging_pass'):
            pdb=args.context/'restored_staging'/(audit['fixture_name']+'.pdb');assert sha(pdb)==audit['staging']['output_sha256']
            details=contacts(pdb);assert len(details)==row['n_his']
            if all(not r['missing_ring_atoms'] for r in details):
                for field,key in [('his_rings_near_target_4p5A','ring_to_target_minimum_A'),
                                  ('his_nitrogens_near_basic_4p5A','nitrogen_to_basic_minimum_A'),
                                  ('his_nitrogens_near_acidic_4p5A','nitrogen_to_acidic_minimum_A')]:
                    row[field]=sum(r[key]<=CONTACT_DISTANCE_A for r in details)
        evidence.append(dict(source_id=row['source_id'],histidines=details,prior_matching_design=(previous or {}).get('name')))
    args.out.parent.mkdir(parents=True,exist_ok=True)
    with args.out.open('w',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(rows[0]),lineterminator='\n');writer.writeheader();writer.writerows(rows)
    receipt=dict(state='complete',rows=len(rows),order_preserved=True,descriptive_contact_cutoff_A=CONTACT_DISTANCE_A,
                 source_shortlist_sha256=sha(args.shortlist),context_report_sha256=sha(args.context/'RESULTS.json'),
                 output_sha256=sha(args.out),script_sha256=sha(__file__),designs=evidence,
                 limitations=['Contact counts and His priors do not establish protonation, affinity or pH selectivity',
                              'Uses observed native sidechain rotamers in rigid staging; no refolding, relaxation or pKa calculation',
                              'Prior novelty is conditional evidence from an earlier conformer, not an organizer score'],new_gpu_spend_usd=0)
    args.out.with_suffix('.evidence.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps(dict(rows=len(rows),order_preserved=True,affinity_or_pH_score_computed=False)))


if __name__=='__main__':main()
