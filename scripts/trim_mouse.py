#!/usr/bin/env python
"""Trim the AlphaFold mouse TfR1 model to an epitope neighbourhood, chain A.

The mouse ortholog has no experimental structure, so the multi-target arm uses the
AlphaFold model. Only well-predicted residues are kept: pLDDT (B-factor column) >= 70.
"""
import sys
import numpy as np
from Bio.PDB import PDBParser, PDBIO, Select
from Bio.PDB.Polypeptide import is_aa

hot = [int(x) for x in sys.argv[1].split(",")]
out = sys.argv[2]
radius = float(sys.argv[3]) if len(sys.argv) > 3 else 12.0
m = PDBParser(QUIET=True).get_structure("mouse", "data/AF-Q62351.pdb")[0]
ch = list(m)[0]
res = [r for r in ch if is_aa(r, standard=True)]
plddt = {r.id[1]: float(np.mean([a.get_bfactor() for a in r])) for r in res}
hot_atoms = np.array([a.coord for r in res if r.id[1] in hot for a in r])
keep = []
for r in res:
    pts = np.array([a.coord for a in r])
    d = float(np.linalg.norm(pts[:, None, :] - hot_atoms[None, :, :], axis=-1).min())
    if d <= radius and plddt[r.id[1]] >= 70:
        keep.append(r.id[1])
nums = sorted(keep)
ivs, start, prev = [], nums[0], nums[0]
for n in nums[1:]:
    if n - prev <= 2: prev = n
    else: ivs.append((start, prev)); start = prev = n
ivs.append((start, prev))
final = {n for lo, hi in ivs if hi - lo + 1 >= 4 for n in range(lo, hi + 1) if n in plddt}

class K(Select):
    def accept_residue(self, r): return is_aa(r, standard=True) and r.id[1] in final
    def accept_atom(self, a): return a.element != "H"

ch.id = "A"
io = PDBIO(); io.set_structure(m); io.save(out, select=K())
print("mouse hotspots %s" % ",".join(map(str, hot)))
print("kept %d residues, fragments: %s" % (len(final), ", ".join("%d-%d" % i for i in ivs)))
print("mean pLDDT over hotspots: %.1f" % np.mean([plddt[h] for h in hot if h in plddt]))
print("wrote %s" % out)
