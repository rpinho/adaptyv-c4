import importlib.util, json, os, tempfile
import numpy as np
import types
# reuse helpers without running main()
src=open('scripts/c4_l6_testretest_20261026.py').read().split('def main(')[0]
src=src.replace('os.path.abspath(__file__)', repr('/Users/ricardo/dev/adaptyv-c4/scripts/x.py'))
m=types.ModuleType('trh'); exec(compile(src,'trh','exec'), m.__dict__)
eight=['tfr1_c4_acid_site1_13','tfr1_c4_acid_site1_147','tfr1_c4_acid_site1_113',
 'tfr1_c4_acid_site2_24','tfr1_c4_acid_site1_143','tfr1_c4_acid_site2_128',
 'tfr1_c4_acid_site1_125','tfr1_c4_acid_site2_55']
seats=[s['name'] for s in json.load(open('results/c4_ranking_v4_acid_off_20261026/v4_ranking_acid_off.json'))['seats']]
targets=sorted(set(eight+seats))
rng=np.random.default_rng(20261027); out=[]
with tempfile.TemporaryDirectory() as td:
    for nm in targets:
        p=os.path.join(m.POOL, nm+'_complex.pdb')
        if not os.path.exists(p): continue
        lines,idx,xyz=m.read_atoms(p)
        base=m.score_coords(lines,idx,xyz,td,'b')
        if base is None: continue
        vals=[m.score_coords(lines,idx,xyz+rng.normal(0,0.10,size=xyz.shape),td,f'j{k}') for k in range(10)]
        vals=[v for v in vals if v is not None]; a=np.array(vals)
        out.append(dict(design=nm, L6=round(base,4), sd=round(float(a.std(ddof=1)),4),
                        lo=round(float(a.min()),4), hi=round(float(a.max()),4),
                        flips=int((np.sign(a)!=np.sign(base)).sum()), n=len(a),
                        in_eight=nm in eight, in_seats=nm in seats))
json.dump(out,open('/tmp/tr_focus.json','w'),indent=1)
for grp,lab in ((True,'THE EIGHT'),(False,'v4 SEATS (not in the eight)')):
    sub=[r for r in out if (r['in_eight'] if grp else (r['in_seats'] and not r['in_eight']))]
    if not sub: continue
    print(f"\n=== {lab} ===")
    print(f"{'design':32} {'L6':>8} {'sd':>7} {'range':>17} {'flips/10':>9}")
    for r in sorted(sub,key=lambda x:x['L6']):
        print(f"{r['design']:32} {r['L6']:8.3f} {r['sd']:7.3f} {r['lo']:8.3f}..{r['hi']:7.3f} {r['flips']:>9}")
    s=[r for r in sub if r['flips']>0]
    print(f"  sign undetermined (>=1 flip in 10): {len(s)} / {len(sub)}")
    u=[r for r in sub if r['sd']>abs(r['L6'])]
    print(f"  sd exceeds |L6|: {len(u)} / {len(sub)}")
