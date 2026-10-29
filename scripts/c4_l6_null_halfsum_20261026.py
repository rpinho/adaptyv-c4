"""Astra's requested statistic: per-design half-sum of L(pf+d,pf) and L(pf-d,pf) =
the expected L6 under equiprobable shift SIGNS, holding the observed |shift| magnitudes
and the observed pKa composition fixed. Plus the zero-gate split under that null."""
import importlib.util, os, json
import numpy as np
spec=importlib.util.spec_from_file_location('la','scripts/c4_linkage_audit_20261025.py')
LA=importlib.util.module_from_spec(spec); spec.loader.exec_module(LA)
P='results/c4_apical_pdb_audit_20261025/all343'
rng=np.random.default_rng(20261026)
names,orig,half,rand=[],[],[],[]
res_contrib={}
for f in sorted(x for x in os.listdir(P) if x.endswith('.pdb')):
    r=LA.score(os.path.join(P,f),'B','A',6.0,7.4)
    if r is None: continue
    names.append(f.replace('_complex.pdb','')); orig.append(r['L'])
    hs=0.0; rs=0.0
    for x in r['rows']:
        pf=x['pka_free']; d=abs(x['pka_bound']-pf)
        up=LA.term(pf+d,pf,6.0)-LA.term(pf+d,pf,7.4)
        dn=LA.term(pf-d,pf,6.0)-LA.term(pf-d,pf,7.4)
        hs+=0.5*(up+dn)
        rs+=up if rng.random()<0.5 else dn
        res_contrib.setdefault(x['resname'],[]).append(0.5*(up+dn))
    half.append(hs); rand.append(rs)
o=np.array(orig); h=np.array(half); rd=np.array(rand)
out={'n':len(o),
 'L6_observed':{'mean':float(o.mean()),'median':float(np.median(o)),'n_below_zero':int((o<0).sum())},
 'NULL_half_sum_expected_L6':{'mean':float(h.mean()),'median':float(np.median(h)),
   'p10':float(np.percentile(h,10)),'p90':float(np.percentile(h,90)),
   'n_below_zero':int((h<0).sum())},
 'NULL_single_random_sign_draw':{'mean':float(rd.mean()),'n_below_zero':int((rd<0).sum())},
 'observed_minus_null_per_design':{'mean':float((o-h).mean()),'median':float(np.median(o-h)),
   'n_below_zero':int(((o-h)<0).sum())},
 'null_contribution_by_residue_type_mean':{k:float(np.mean(v)) for k,v in
   sorted(res_contrib.items(),key=lambda kv:-abs(np.mean(kv[1])))},
 'zero_gate_split':{'observed_L6_lt_0':int((o<0).sum()),
   'null_half_sum_lt_0':int((h<0).sum()),
   'null_random_sign_lt_0':int((rd<0).sum())}}
json.dump({'summary':out,'per_design':[{'name':n,'L6':float(a),'null_half_sum':float(b),
  'excess_over_null':float(a-b)} for n,a,b in zip(names,o,h)]},open('/tmp/l6null.json','w'),indent=1)
print(json.dumps(out,indent=1))
