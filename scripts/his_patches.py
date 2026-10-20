#!/usr/bin/env python
"""Native-His epitope options (mechanism family 6: the target supplies half the linkage).

For each native histidine in the assay construct, report whether it is usable as the
centre of an epitope: exposure, human/mouse identity of the surrounding patch, whether
it falls in an excluded site, and the patch's residue list.
"""
import json, os, sys
import numpy as np
HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(HERE, "scripts"))
from epitope_scan2 import uniprot_seq, conservation_map, dimer_sasa, DATA

scan = json.load(open(os.path.join(HERE, "results", "epitope_scan2.json")))
ks = scan["meta"]["known_sites"]
S = {lab: set(d["residues"]) for lab, d in ks.items()}
tag = set(scan["meta"]["tag_excluded"]); glyc = set(scan["meta"]["glycan_sites"])
block = S["transferrin"] | S["HFE"] | tag | glyc

human = uniprot_seq(os.path.join(DATA, "P02786.txt"))
mouse = uniprot_seq(os.path.join(DATA, "Q62351.txt"))
cons = conservation_map(human, mouse)
sasa, cb = dimer_sasa("11qc")
HIS = [n for n in sorted(sasa) if sasa[n][0] == 'H']
RAD = 10.0
rows = []
print("%-6s %-9s %-9s %-8s %-6s %-5s %s" % ("His", "relSASA", "mouseHis?", "excluded", "ident", "hyd", "patch"))
for h in HIS:
    if h not in cb:
        continue
    members = [n for n in sasa if n in cb and n not in block
               and sasa[n][2] >= 0.15 and float(np.linalg.norm(cb[n] - cb[h])) <= RAD]
    if h in block:
        status = "YES(" + ",".join(lab for lab, s in S.items() if h in s and lab in ("transferrin", "HFE")) + ")"
        if h in tag: status += "tag"
        if h in glyc: status += "glycan"
    else:
        status = "no"
    ident = float(np.mean([cons.get(n, 0) for n in members])) if members else 0.0
    seq = "".join(sasa[n][0] for n in sorted(members))
    hyd = sum(1 for c in seq if c in "AVLIMFWYC")
    rows.append(dict(his=h, rel_sasa=round(sasa[h][2], 3), conserved_his=bool(cons.get(h, 0)),
                     excluded=status, patch_identity=round(ident, 3), hydrophobic=hyd,
                     patch=sorted(members), seq=seq,
                     bindcraft_hotspots=",".join("A%d" % n for n in sorted(members))))
    print("%-6d %-9.3f %-9s %-8s %-6.2f %-5d %s" % (
        h, sasa[h][2], "yes" if cons.get(h, 0) else "NO", status, ident, hyd,
        ",".join(str(n) for n in sorted(members))))
json.dump(rows, open(os.path.join(HERE, "results", "his_patches.json"), "w"), indent=1)
print("\nwrote results/his_patches.json")
