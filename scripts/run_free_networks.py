#!/usr/bin/env python
"""CPU scout of reusable binder pH modules; no claim of receptor binding.

The prior docked poses fail the full-context audit. Here only isolated binder
sequences/conformations are modeled. Favorable free His pKa does not establish a
binding switch. This is a small mechanistic substrate for later full-context design.
"""
import argparse
import json
import os
from pathlib import Path
import time
import install_switch as install
from ph_linkage import run_propka
from run_context_ensemble import local_prepare, extract, mutate_fixed_backbone, sha, dump, ROOT


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--out',default='results/free_networks')
    ap.add_argument('--limit',type=int,default=24);ap.add_argument('--seeds',type=int,default=3)
    a=ap.parse_args();out=Path(a.out).resolve();out.mkdir(parents=True,exist_ok=True)
    if (out/'status.json').exists():raise SystemExit('fresh output required')
    queue=json.loads((ROOT/'results/next_campaign/parent_priorities.json').read_text())
    parent=queue[0];source=Path(parent['variant_path'])
    old=json.loads(next(source.parent.glob('*variants.json')).read_text())
    sites=old['pairs'][int(parent['pair'])-1]['sites']
    started=time.time();deadline=time.monotonic()+3600
    state=dict(status='running',pid=os.getpid(),started=started,models=[],errors=[],new_gpu_spend_usd=0,
               source=str(source),source_sha256=sha(source),script_sha256=sha(__file__),
               context_status='old docking rejected; isolated binder scout only',
               interpretation='free-state pKa and fixed-backbone energy only; no receptor affinity or pH selectivity claim')
    dump(out/'status.json',state)
    install.init();pr=install.pyrosetta;sf=pr.get_fa_scorefxn()
    base=extract(pr.pose_from_pdb(str(source)),{'C'})
    centres=[base.pdb_info().pdb2pose('C',s['resnum']) for s in sites]
    positions=[]
    for i in range(1,base.total_residue()+1):
        if i in centres or base.residue(i).name1() in 'GP':continue
        d=min(base.residue(i).nbr_atom_xyz().distance(base.residue(c).nbr_atom_xyz()) for c in centres)
        if d<=8:positions.append((d,i))
    mutations=[('HH',[])]
    for c in centres:
        mutations.append((f'H{base.pdb_info().number(c)}Q',[(c,'Q')]))
    mutations.append(('QQ',[(i,'Q') for i in centres]))
    for _,i in sorted(positions):
        for aa in 'DE':
            if base.residue(i).name1()!=aa:
                mutations.append((f'{base.residue(i).name1()}{base.pdb_info().number(i)}{aa}',[(i,aa)]))
    mutations=mutations[:a.limit]
    dump(out/'queue.json',dict(parent=parent,sites=sites,variants=mutations,seeds=a.seeds,
                               selection='nearest binder neighbor atoms <=8A from selected His; D/E support scout',
                               min_neutral_retention_REU='ddE <=3 vs same-seed HH; engineering filter only'))
    for seed_index in range(a.seeds):
        baseline=None
        for label,changes in mutations:
            if time.monotonic()>=deadline:break
            seed=202610210+seed_index
            try:
                pr.rosetta.numeric.random.rg().set_seed(seed)
                work=base.clone()
                for i,aa in changes:
                    mutate_fixed_backbone(work,i,aa)
                work=local_prepare(work,sf,range(1,work.total_residue()+1))
                path=out/f'seed{seed_index}_{label}.pdb';work.dump_pdb(str(path))
                energy=float(sf(work))
                if label=='HH':baseline=energy
                groups=run_propka(str(path),str(out))
                if not groups:raise ValueError('PROPKA failed')
                his=[dict(resnum=s['resnum'],pka=groups.get(('C',s['resnum'],'HIS')))
                     for s in sites if work.residue(work.pdb_info().pdb2pose('C',s['resnum'])).name1()=='H']
                row=dict(seed=seed,label=label,pdb=str(path),sha256=sha(path),energy_REU=energy,
                         ddE_vs_HH_REU=energy-baseline,neutral_energy_filter=energy-baseline<=3,
                         selected_his=his,missing_selected_his=any(h['pka'] is None for h in his),
                         rows=[dict(chain=c,resnum=n,resname=r,pka=p) for (c,n,r),p in groups.items()])
                state['models'].append(row)
                print(json.dumps({k:v for k,v in row.items() if k not in ['rows','pdb','sha256']}),flush=True)
            except Exception as exc:
                state['errors'].append(dict(seed=seed,label=label,error=repr(exc)))
                print('ERROR',repr(exc),flush=True)
            state['updated']=time.time();dump(out/'status.json',state)
        if time.monotonic()>=deadline:break
    state.update(status='completed' if len(state['models'])==len(mutations)*a.seeds else 'partial',finished=time.time())
    dump(out/'status.json',state)


if __name__=='__main__':main()
