#!/usr/bin/env python
"""Derive the final hotspot sets for the generators.

Three epitope classes, each emitted as a BindCraft/RFdiffusion-ready residue list:
  A. de-risked apical: the consensus of the three designed TfR1 binders in the PDB
     (6WRV/6WRW/6WRX) MINUS the transferrin and HFE footprints, glycans and the tag
     vicinity. Proven designable, and outside the transferrin site.
  B. ligandable secondary: the bicyclic-peptide site (9GH7) under the same subtraction.
  C. novel: the top-scoring scanned patch with no overlap to any known site.
Each set is reported with human/mouse identity, exposure and native-His content.
"""
import json, os, sys
import numpy as np

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
scan = json.load(open(os.path.join(HERE, "results", "epitope_scan2.json")))
ks = scan["meta"]["known_sites"]
S = {lab: set(d["residues"]) for lab, d in ks.items()}
tag = set(scan["meta"]["tag_excluded"])
glyc = set(scan["meta"]["glycan_sites"])
block = S["transferrin"] | S["HFE"] | tag | glyc

designed = (S["designed_3DS18"] & S["designed_2DS25.5"]) | (S["designed_2DS25.5"] & S["designed_2DS25.1"]) | (S["designed_3DS18"] & S["designed_2DS25.1"])
sets = {
    "A_apical_derisked": sorted(designed - block),
    "A2_apical_designed_union": sorted((S["designed_3DS18"] | S["designed_2DS25.5"] | S["designed_2DS25.1"]) - block),
    "A3_machupo_only": sorted(S["MachupoGP1"] - block - designed),
    "B_bicyclic_site": sorted(S["bicyclic_peptide"] - block),
}
best_novel = next((p for p in scan["patches"] if not p["known_site_overlap"]), None)
if best_novel:
    sets["C_novel_" + str(best_novel["centre"])] = best_novel["residues"]

# per-residue annotation from the scan meta is not stored, so recompute identity/exposure
sys.path.insert(0, os.path.join(HERE, "scripts"))
from epitope_scan2 import uniprot_seq, conservation_map, dimer_sasa, DATA
human = uniprot_seq(os.path.join(DATA, "P02786.txt"))
mouse = uniprot_seq(os.path.join(DATA, "Q62351.txt"))
cons = conservation_map(human, mouse)
sasa, cb = dimer_sasa("11qc")

out = {}
print("%-28s %-4s %-7s %-8s %-5s %-10s %s" % ("set", "n", "ident", "relSASA", "hyd", "nativeHis", "residues"))
for name, res in sets.items():
    res = [r for r in res if r in sasa]
    if not res:
        continue
    ident = float(np.mean([cons.get(r, 0) for r in res]))
    expo = float(np.mean([sasa[r][2] for r in res]))
    seq = "".join(sasa[r][0] for r in sorted(res))
    hyd = sum(1 for c in seq if c in "AVLIMFWYC")
    his = [r for r in res if sasa[r][0] == 'H']
    # spatial compactness: max CB-CB distance
    pts = np.array([cb[r] for r in res if r in cb])
    span = float(np.max(np.linalg.norm(pts[:, None, :] - pts[None, :, :], axis=-1))) if len(pts) > 1 else 0.0
    out[name] = dict(residues=res, n=len(res), identity=round(ident, 3),
                     mean_rel_sasa=round(expo, 3), hydrophobic=hyd, native_his=his,
                     span_a=round(span, 1), seq=seq,
                     bindcraft_hotspots=",".join("A%d" % r for r in res))
    print("%-28s %-4d %-7.2f %-8.3f %-5d %-10s %s" % (
        name, len(res), ident, expo, hyd, ",".join(str(h) for h in his) or "-",
        ",".join(str(r) for r in res)))
    print("%-28s span %.1f A   seq %s" % ("", span, seq))

json.dump(out, open(os.path.join(HERE, "results", "hotspot_sets.json"), "w"), indent=1)
print("\nwrote results/hotspot_sets.json")
