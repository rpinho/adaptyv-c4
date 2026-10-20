#!/usr/bin/env python
"""Epitope scan v2 for the Challenge 4 target.

Upgrades over v1:
  * template is 11QC (2.4 A cryo-EM, TfR1 ectodomain homodimer, residues 122-760 =
    exactly the assay construct) instead of 1CX8 (3.2 A)
  * the transferrin footprint comes from 3S9L (3.22 A X-ray TfR1-transferrin complex)
    instead of 1SUV (7.5 A)
  * every known binding site on human TfR1 that is in the PDB is mapped and each
    candidate patch is annotated against it: transferrin (3S9L), HFE (1DE4),
    three computationally designed binders (6WRV/6WRW/6WRX), Machupo virus GP1
    (3KAS) and a bicyclic peptide (9GH7)

Masks applied: inside the assay construct, solvent-exposed in the physiological dimer,
not in the transferrin or HFE footprint, not an N-glycosylation site, not within 15 A of
the C-terminal Twin-Strep anchor. Known designed/viral/peptide epitopes are ANNOTATED,
never silently excluded: overlap there means a proven-designable site, which is a
de-risking choice, while a patch far from all of them is a novelty choice.

Usage:  python epitope_scan2.py [--window LO-HI] [--min-ident F] [--min-hyd N] [--top N]
"""
import argparse, json, os
from collections import defaultdict
import numpy as np
import freesasa
from Bio.PDB import MMCIFParser
from Bio.PDB.Polypeptide import is_aa

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(HERE, "data")
RESULTS = os.path.join(HERE, "results")
CONSTRUCT = (121, 760)
TAG_EXCLUSION_A = 15.0
CONTACT_A = 5.0
HYDROPHOBIC = set("AVLIMFWYC")
CHARGED = set("DEKR")
GLYCAN_SITES = [251, 317, 727]   # NAG-linked ASN, detected from 1CX8 in v1

MAXASA = {'A':129,'R':274,'N':195,'D':193,'C':167,'Q':225,'E':223,'G':104,'H':224,'I':197,
          'L':201,'K':236,'M':224,'F':240,'P':159,'S':155,'T':172,'W':285,'Y':263,'V':174}
T3 = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E','GLY':'G',
      'HIS':'H','ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S',
      'THR':'T','TRP':'W','TYR':'Y','VAL':'V'}

# complexes to map: file -> (label, exclude?)
COMPLEXES = [
    ("3s9l", "transferrin", True),
    ("1de4", "HFE", True),
    ("6wrv", "designed_3DS18", False),
    ("6wrw", "designed_2DS25.5", False),
    ("6wrx", "designed_2DS25.1", False),
    ("3kas", "MachupoGP1", False),
    ("9gh7", "bicyclic_peptide", False),
]


def uniprot_seq(path):
    txt = open(path).read()
    block = txt.split("\nSQ   ", 1)[1].split("\n", 1)[1]
    return "".join(block.split("//")[0].split())


def conservation_map(human, mouse):
    n, k = len(human), len(mouse)
    GAP = -2
    prev = [GAP * j for j in range(k + 1)]
    ptr = []
    for i in range(1, n + 1):
        cur = [GAP * i] + [0] * k
        row = bytearray(k + 1)
        hi = human[i - 1]
        for j in range(1, k + 1):
            d = prev[j - 1] + (1 if hi == mouse[j - 1] else -1)
            u = prev[j] + GAP
            l = cur[j - 1] + GAP
            if d >= u and d >= l:
                cur[j] = d; row[j] = 0
            elif u >= l:
                cur[j] = u; row[j] = 1
            else:
                cur[j] = l; row[j] = 2
        ptr.append(row); prev = cur
    i, j, cons = n, k, {}
    while i > 0 and j > 0:
        t = ptr[i - 1][j]
        if t == 0:
            cons[i] = 1 if human[i - 1] == mouse[j - 1] else 0
            i -= 1; j -= 1
        elif t == 1:
            cons[i] = 0; i -= 1
        else:
            j -= 1
    while i > 0:
        cons[i] = 0; i -= 1
    return cons


def model_of(code):
    return MMCIFParser(QUIET=True).get_structure(code, os.path.join(DATA, code + ".cif"))[0]


def aa_res(chain):
    return [r for r in chain if is_aa(r, standard=True)]


def chain_seq(rl):
    return "".join(T3.get(r.get_resname(), 'X') for r in rl)


def heavy_atoms(res_list):
    pts, owner = [], []
    for r in res_list:
        for a in r:
            if a.element != 'H':
                pts.append(a.coord); owner.append(r.id[1])
    return np.array(pts), np.array(owner)


def contacts(target_res, partner_pts, cutoff=CONTACT_A):
    pts, owner = heavy_atoms(target_res)
    hits = set()
    for s in range(0, len(pts), 1500):
        sl = slice(s, s + 1500)
        d = np.linalg.norm(pts[sl][:, None, :] - partner_pts[None, :, :], axis=-1)
        close = (d < cutoff).any(axis=1)
        hits.update(int(x) for x in owner[sl][close])
    return hits


def footprint(code, human_seq):
    """TfR1 residues contacting any non-TfR1 protein chain, in TfR1 numbering."""
    m = model_of(code)
    tfr, partner = [], []
    for c in m:
        rl = aa_res(c)
        if len(rl) < 8:
            continue
        seq = chain_seq(rl)
        probe = seq[len(seq) // 2: len(seq) // 2 + 25]
        (tfr if probe and probe in human_seq else partner).append(rl)
    if not tfr or not partner:
        return set(), 0, 0
    ppts = np.vstack([heavy_atoms(rl)[0] for rl in partner])
    hits = set()
    for rl in tfr:
        hits |= contacts(rl, ppts)
    return hits, len(tfr), len(partner)


def dimer_sasa(code, chain_a="A", chain_b="B"):
    m = model_of(code)
    sel = aa_res(m[chain_a]) + aa_res(m[chain_b])
    cls = freesasa.Classifier()
    radii, coords, keys = [], [], []
    for r in sel:
        for a in r:
            if a.element == 'H':
                continue
            try:
                rad = cls.radius(r.get_resname(), a.get_name())
            except Exception:
                rad = 1.8
            radii.append(rad if rad > 0 else 1.8)
            coords.extend(float(x) for x in a.coord)
            keys.append((r.get_parent().id, r.id[1]))
    res = freesasa.calcCoord(coords, radii)
    per = defaultdict(float)
    for i, k in enumerate(keys):
        per[k] += res.atomArea(i)
    out, cb = {}, {}
    for r in aa_res(m[chain_a]):
        aa = T3.get(r.get_resname())
        if not aa:
            continue
        s = per.get((chain_a, r.id[1]), 0.0)
        out[r.id[1]] = (aa, s, s / MAXASA[aa])
        at = r['CB'] if 'CB' in r else (r['CA'] if 'CA' in r else None)
        if at is not None:
            cb[r.id[1]] = at.coord
    return out, cb


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--window", default="0-0")
    ap.add_argument("--min-ident", type=float, default=0.85)
    ap.add_argument("--min-hyd", type=int, default=3)
    ap.add_argument("--min-exposed", type=float, default=0.25)
    ap.add_argument("--radius", type=float, default=10.0)
    ap.add_argument("--top", type=int, default=8)
    ap.add_argument("--out", default=os.path.join(RESULTS, "epitope_scan2.json"))
    a = ap.parse_args()
    lo, hi = (int(x) for x in a.window.split("-"))

    human = uniprot_seq(os.path.join(DATA, "P02786.txt"))
    mouse = uniprot_seq(os.path.join(DATA, "Q62351.txt"))
    cons = conservation_map(human, mouse)

    sasa, cb = dimer_sasa("11qc")
    modeled = sorted(sasa)

    fps, excluded = {}, set()
    for code, label, do_exclude in COMPLEXES:
        hits, nt, npart = footprint(code, human)
        fps[label] = dict(pdb=code.upper(), residues=sorted(hits), n=len(hits),
                          tfr_chains=nt, partner_chains=npart, excluded=do_exclude)
        if do_exclude:
            excluded |= hits

    anchor = cb.get(modeled[-1])
    tag_near = {n for n, c in cb.items() if anchor is not None
                and float(np.linalg.norm(c - anchor)) < TAG_EXCLUSION_A}
    hard = excluded | set(GLYCAN_SITES) | tag_near

    allowed = {n for n in modeled
               if CONSTRUCT[0] <= n <= CONSTRUCT[1] and sasa[n][2] >= 0.15 and n not in hard}
    centres = [n for n in modeled
               if CONSTRUCT[0] <= n <= CONSTRUCT[1] and sasa[n][2] >= a.min_exposed
               and n not in hard and (not hi or lo <= n <= hi)]

    known = {lab: set(d["residues"]) for lab, d in fps.items() if not d["excluded"]}
    patches = []
    for c0 in centres:
        if c0 not in cb:
            continue
        members = [n for n in allowed if n in cb
                   and float(np.linalg.norm(cb[n] - cb[c0])) <= a.radius]
        if len(members) < 6:
            continue
        seq = "".join(sasa[n][0] for n in sorted(members))
        hyd = sum(1 for ch in seq if ch in HYDROPHOBIC)
        chg = sum(1 for ch in seq if ch in CHARGED)
        ident = float(np.mean([cons.get(n, 0) for n in members]))
        if ident < a.min_ident or hyd < a.min_hyd:
            continue
        expo = float(np.mean([sasa[n][2] for n in members]))
        his = [n for n in members if sasa[n][0] == 'H']
        mset = set(members)
        overlaps = {lab: len(mset & rs) for lab, rs in known.items() if mset & rs}
        d_tf = min((float(np.linalg.norm(cb[n] - cb[c0]))
                    for n in fps["transferrin"]["residues"] if n in cb), default=float('inf'))
        d_hfe = min((float(np.linalg.norm(cb[n] - cb[c0]))
                     for n in fps["HFE"]["residues"] if n in cb), default=float('inf'))
        proven = sum(v for k, v in overlaps.items() if k.startswith("designed"))
        score = (hyd * 1.0 + expo * 4.0 - chg * 0.6 + len(members) * 0.15
                 + (1.5 if his else 0.0) + (ident - a.min_ident) * 20.0)
        patches.append(dict(centre=c0, n=len(members), residues=sorted(members), seq=seq,
                            identity=round(ident, 3), hydrophobic=hyd, charged=chg,
                            mean_rel_sasa=round(expo, 3), native_his=his,
                            dist_to_transferrin=round(d_tf, 1), dist_to_HFE=round(d_hfe, 1),
                            known_site_overlap=overlaps, designed_binder_overlap=proven,
                            score=round(score, 2)))

    patches.sort(key=lambda p: -p["score"])
    chosen, used = [], set()
    for p in patches:
        if len(set(p["residues"]) & used) > len(p["residues"]) * 0.4:
            continue
        chosen.append(p); used.update(p["residues"])
        if len(chosen) >= a.top:
            break

    out = dict(meta=dict(template="11QC (2.4 A cryo-EM, chains A+B, residues %d-%d)"
                         % (modeled[0], modeled[-1]),
                         construct=CONSTRUCT, n_modeled=len(modeled),
                         glycan_sites=GLYCAN_SITES, tag_excluded=sorted(tag_near),
                         hard_excluded=len(hard), allowed=len(allowed), centres=len(centres),
                         thresholds=vars(a), known_sites=fps),
               patches=chosen)
    os.makedirs(RESULTS, exist_ok=True)
    json.dump(out, open(a.out, "w"), indent=1)

    print("template 11QC: chain A modeled %d-%d (%d residues), dimer context chain B"
          % (modeled[0], modeled[-1], len(modeled)))
    print("\nknown binding sites on human TfR1 (this run, %.0f A heavy-atom contacts):" % CONTACT_A)
    for lab, d in fps.items():
        print("  %-18s %s  %3d residues  %s" % (lab, d["pdb"], d["n"],
              "EXCLUDED" if d["excluded"] else "annotated"))
        print("      %s" % ",".join(str(x) for x in d["residues"][:28]) +
              (" ..." if d["n"] > 28 else ""))
    print("\nhard-excluded residues: %d | allowed in patches: %d | candidate centres: %d"
          % (len(hard), len(allowed), len(centres)))
    print("\n%-7s %-4s %-6s %-6s %-4s %-4s %-8s %-7s %-7s %-10s %s"
          % ("centre", "n", "score", "ident", "hyd", "chg", "relSASA", "dTf", "dHFE", "nativeHis", "overlaps"))
    for p in chosen:
        print("%-7d %-4d %-6.2f %-6.2f %-4d %-4d %-8.3f %-7.1f %-7.1f %-10s %s"
              % (p["centre"], p["n"], p["score"], p["identity"], p["hydrophobic"], p["charged"],
                 p["mean_rel_sasa"], p["dist_to_transferrin"], p["dist_to_HFE"],
                 ",".join(str(h) for h in p["native_his"]) or "-",
                 p["known_site_overlap"] or "none"))
    print("\nresidue lists:")
    for p in chosen:
        print("  centre %-4d %s" % (p["centre"], ",".join(str(x) for x in p["residues"])))
        print("             seq %s" % p["seq"])


if __name__ == "__main__":
    main()
