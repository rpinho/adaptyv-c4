"""Verify Fable 5's Q4 on the FULL 342 apical pool (it measured a 66-design proxy).
Claim: L6 is carried almost entirely by a handful of interface histidines, so the
L6<0 gate selects for 'has His near the interface', not for a pH switch."""
import importlib.util, os, json
import numpy as np
spec=importlib.util.spec_from_file_location('la','scripts/c4_linkage_audit_20261025.py')
LA=importlib.util.module_from_spec(spec); spec.loader.exec_module(LA)
P='results/c4_apical_pdb_audit_20261025/all343'
tot={}; cnt={}; pka={}; perdesign=[]
for f in sorted(x for x in os.listdir(P) if x.endswith('.pdb')):
    r=LA.score(os.path.join(P,f),'B','A',6.0,7.4)
    if r is None: continue
    his=0.0; nhis=0
    for x in r['rows']:
        rn=x['resname']; c=x['contribution']
        tot[rn]=tot.get(rn,0.0)+c; cnt[rn]=cnt.get(rn,0)+1
        pka.setdefault(rn,[]).append(x['pka_free'])
        if rn=='HIS': his+=c; nhis+=1
    perdesign.append({'name':f.replace('_complex.pdb',''),'L6':r['L'],
                      'his_contrib':his,'n_his':nhis,'n_groups':r['n_groups'],
                      'frac_abs':abs(his)/abs(r['L']) if abs(r['L'])>1e-9 else None})
L=np.array([d['L6'] for d in perdesign]); H=np.array([d['his_contrib'] for d in perdesign])
fr=np.array([d['frac_abs'] for d in perdesign if d['frac_abs'] is not None])
nh=np.array([d['n_his'] for d in perdesign])
out={'n_designs':len(perdesign),
 'by_residue_type':{k:{'n':cnt[k],'total_contribution':round(tot[k],3),
    'mean':round(tot[k]/cnt[k],5),'median_pka_free':round(float(np.median(pka[k])),2)}
    for k in sorted(tot,key=lambda k:tot[k])},
 'pool_total_L6':round(float(L.sum()),2),
 'pool_total_HIS_contribution':round(float(H.sum()),2),
 'median_frac_absL6_carried_by_His':round(float(np.median(fr)),3),
 'median_n_contributing_His_per_design':float(np.median(nh)),
 'median_n_contributing_groups':float(np.median([d['n_groups'] for d in perdesign])),
 'designs_where_His_alone_flips_the_sign':int(sum(1 for d in perdesign
    if d['frac_abs'] is not None and np.sign(d['L6'])==np.sign(d['his_contrib']) and abs(d['his_contrib'])>abs(d['L6']-d['his_contrib']))),
 'corr_L6_vs_his_contrib':round(float(np.corrcoef(L,H)[0,1]),4),
 'corr_L6_vs_n_his':round(float(np.corrcoef(L,nh)[0,1]),4)}
json.dump({'summary':out,'per_design':perdesign},open('/tmp/hiscensus.json','w'),indent=1)
print(json.dumps(out,indent=1))
