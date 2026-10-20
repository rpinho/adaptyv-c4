#!/usr/bin/env python
"""Per-residue check: is a candidate epitope residue actually reachable in the dimer?

For each residue of interest, report relative SASA in the A+B dimer, the gain in SASA
when chain B is removed (how much of its surface the partner subunit covers) and its
closest approach to chain B. A residue whose surface only appears on the isolated
monomer is not a usable epitope residue.
"""
import sys, os, json
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from epitope_scan2 import (uniprot_seq, conservation_map, model_of, aa_res, MAXASA, T3, DATA)
import freesasa
from collections import defaultdict


def sasa_of(res_lists, chain_of_interest):
    cls = freesasa.Classifier()
    radii, coords, keys = [], [], []
    for cid, rl in res_lists:
        for r in rl:
            for at in r:
                if at.element == 'H':
                    continue
                try:
                    rad = cls.radius(r.get_resname(), at.get_name())
                except Exception:
                    rad = 1.8
                radii.append(rad if rad > 0 else 1.8)
                coords.extend(float(x) for x in at.coord)
                keys.append((cid, r.id[1]))
    res = freesasa.calcCoord(coords, radii)
    per = defaultdict(float)
    for i, k in enumerate(keys):
        per[k] += res.atomArea(i)
    return {n: v for (c, n), v in per.items() if c == chain_of_interest}


m = model_of("11qc")
A, B = aa_res(m['A']), aa_res(m['B'])
dimer = sasa_of([('A', A), ('B', B)], 'A')
mono = sasa_of([('A', A)], 'A')
bpts = np.array([at.coord for r in B for at in r if at.element != 'H'])
aa = {r.id[1]: T3.get(r.get_resname(), 'X') for r in A}
human = uniprot_seq(os.path.join(DATA, "P02786.txt"))
mouse = uniprot_seq(os.path.join(DATA, "Q62351.txt"))
cons = conservation_map(human, mouse)

want = [int(x) for x in sys.argv[1].split(",")]
print("%-6s %-3s %-9s %-9s %-9s %-8s %s" % ("res", "aa", "relSASA", "monoSASA", "occluded", "d_to_B", "mouse-identical"))
rows = []
for n in want:
    if n not in aa:
        print("%-6d missing from the template" % n); continue
    pts = np.array([at.coord for at in m['A'][n] if at.element != 'H'])
    d = float(np.linalg.norm(pts[:, None, :] - bpts[None, :, :], axis=-1).min())
    rel_d = dimer.get(n, 0) / MAXASA[aa[n]]
    rel_m = mono.get(n, 0) / MAXASA[aa[n]]
    occ = rel_m - rel_d
    rows.append(dict(res=n, aa=aa[n], rel_dimer=round(rel_d, 3), rel_mono=round(rel_m, 3),
                     occluded=round(occ, 3), dist_to_B=round(d, 1), conserved=bool(cons.get(n, 0))))
    print("%-6d %-3s %-9.3f %-9.3f %-9.3f %-8.1f %s" % (n, aa[n], rel_d, rel_m, occ, d,
          "yes" if cons.get(n, 0) else "NO"))
json.dump(rows, open(os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "results", "dimer_occlusion.json"), "w"), indent=1)
clean = [r['res'] for r in rows if r['dist_to_B'] >= 8 and r['rel_dimer'] >= 0.15 and r['conserved']]
print("\nclean, conserved, reachable (>=8 A from chain B, relSASA>=0.15): %s" % (",".join(map(str, clean)) or "none"))
