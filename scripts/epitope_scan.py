#!/usr/bin/env python
"""Epitope scan for the Challenge 4 target.

Intersects three masks on the TfR1 ectodomain and emits ranked candidate patches:
  (a) human/mouse conservation  (local Needleman-Wunsch, UniProt P02786 vs Q62351)
  (b) solvent exposure in the physiological homodimer (so the dimer interface is excluded)
  (c) exclusions: outside the assay construct, the transferrin footprint (from the
      TfR1-transferrin complex), N-glycosylation sites, and the C-terminal tag vicinity.

Outputs JSON + a readable table. No network access; structures and UniProt records
are read from ./pdb and ./ .
"""
import json, os, sys, math
from collections import defaultdict
import numpy as np
import freesasa
from Bio.PDB import MMCIFParser
from Bio.PDB.Polypeptide import is_aa

HERE = os.path.dirname(os.path.abspath(__file__))
PDB = os.path.join(HERE, "pdb")
CONSTRUCT = (121, 760)          # assay construct, per the official challenge page
# optional: restrict candidate patch CENTRES to a domain window, e.g. the apical domain
WINDOW = tuple(int(x) for x in os.environ.get("WINDOW", "0-0").split("-"))
MIN_IDENT = float(os.environ.get("MIN_IDENT", "0.85"))
MIN_HYD = int(os.environ.get("MIN_HYD", "3"))
TOP_N = int(os.environ.get("TOP_N", "8"))
TAG_EXCLUSION_A = 15.0          # keep this far from the C-terminal Twin-Strep anchor
TF_EXCLUSION_A = 5.0            # heavy-atom contact cutoff for the transferrin footprint
EXPOSED_REL_SASA = 0.25         # relative SASA threshold for "surface"
PATCH_RADIUS_A = 10.0           # CB-CB radius of a patch
HYDROPHOBIC = set("AVLIMFWYC")
CHARGED = set("DEKR")

MAXASA = {  # Tien et al. 2013 theoretical maxima (A^2)
 'A':129,'R':274,'N':195,'D':193,'C':167,'Q':225,'E':223,'G':104,'H':224,'I':197,
 'L':201,'K':236,'M':224,'F':240,'P':159,'S':155,'T':172,'W':285,'Y':263,'V':174}
THREE2ONE = {
 'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E','GLY':'G','HIS':'H',
 'ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S','THR':'T','TRP':'W',
 'TYR':'Y','VAL':'V'}


# ---------- conservation ----------
def uniprot_seq(path):
    txt = open(path).read()
    block = txt.split("\nSQ   ", 1)[1].split("\n", 1)[1]
    return "".join(block.split("//")[0].split())


def conservation_map(human, mouse):
    """Return {human_resnum: 1 if identical to the aligned mouse residue else 0}."""
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
    i, j = n, k
    cons = {}
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


# ---------- structure helpers ----------
def load(path, name):
    return MMCIFParser(QUIET=True).get_structure(name, path)[0]


def residues(chain):
    return [r for r in chain if is_aa(r, standard=True)]


def atoms_np(res_list):
    pts, owner = [], []
    for r in res_list:
        for a in r:
            if a.element != 'H':
                pts.append(a.coord); owner.append(r.id[1])
    return np.array(pts), np.array(owner)


def pick_dimer_partner(model):
    """TfR1 is a homodimer; find chain A's largest-interface partner."""
    chains = {c.id: residues(c) for c in model if len(residues(c)) > 100}
    if 'A' not in chains:
        raise SystemExit("no chain A with protein residues")
    a_pts, _ = atoms_np(chains['A'])
    best, best_n = None, -1
    for cid, rl in chains.items():
        if cid == 'A':
            continue
        pts, _ = atoms_np(rl)
        if not len(pts):
            continue
        # count A atoms within 5 A of this chain
        n = 0
        for chunk in np.array_split(a_pts, max(1, len(a_pts) // 2000)):
            d = np.linalg.norm(chunk[:, None, :] - pts[None, :, :], axis=-1)
            n += int((d < 5.0).any(axis=1).sum())
        if n > best_n:
            best, best_n = cid, n
    return best, best_n, list(chains.keys())


def sasa_for_chain(model, chain_id, partner_id):
    """Per-residue relative SASA of chain_id computed in the chain+partner assembly."""
    sel = []
    for c in model:
        if c.id in (chain_id, partner_id):
            sel.extend(residues(c))
    radii, coords, keys = [], [], []
    cls = freesasa.Classifier()
    for r in sel:
        for a in r:
            if a.element == 'H':
                continue
            rn = r.get_resname()
            try:
                rad = cls.radius(rn, a.get_name())
            except Exception:
                rad = 1.8
            if rad <= 0:
                rad = 1.8
            radii.append(rad)
            coords.extend([float(x) for x in a.coord])
            keys.append((r.get_parent().id, r.id[1]))
    res = freesasa.calcCoord(coords, radii)
    per_res = defaultdict(float)
    for i, k in enumerate(keys):
        per_res[k] += res.atomArea(i)
    out = {}
    for r in residues(model[chain_id]):
        aa = THREE2ONE.get(r.get_resname())
        if not aa:
            continue
        abs_sasa = per_res.get((chain_id, r.id[1]), 0.0)
        out[r.id[1]] = (aa, abs_sasa, abs_sasa / MAXASA[aa])
    return out


def glycosylation_sites(model, chain_id):
    """ASN residues with a sugar (NAG/NDG/BMA/MAN) within 2.2 A of ND2."""
    sugars = []
    for c in model:
        for r in c:
            if r.get_resname().strip() in ("NAG", "NDG", "BMA", "MAN", "FUC"):
                sugars.extend([a.coord for a in r])
    sugars = np.array(sugars) if sugars else np.zeros((0, 3))
    sites = []
    if not len(sugars):
        return sites
    for r in residues(model[chain_id]):
        if r.get_resname() != 'ASN' or 'ND2' not in r:
            continue
        d = np.linalg.norm(sugars - r['ND2'].coord, axis=-1).min()
        if d < 2.2:
            sites.append(r.id[1])
    return sites


def tf_footprint(path):
    """TfR1 residues contacting transferrin in the complex structure."""
    model = load(path, "complex")
    chain_lens = {c.id: len(residues(c)) for c in model if len(residues(c)) > 50}
    # TfR1 chains are the ~639-residue ones; transferrin lobes are ~329/345
    tfr = [cid for cid, n in chain_lens.items() if n > 500]
    tf = [cid for cid, n in chain_lens.items() if 200 < n <= 500]
    if not tfr or not tf:
        return set(), chain_lens, tfr, tf
    tf_pts = np.vstack([atoms_np(residues(model[c]))[0] for c in tf])
    hits = set()
    for c in tfr:
        pts, owner = atoms_np(residues(model[c]))
        for chunk_i in range(0, len(pts), 2000):
            sl = slice(chunk_i, chunk_i + 2000)
            d = np.linalg.norm(pts[sl][:, None, :] - tf_pts[None, :, :], axis=-1)
            close = (d < TF_EXCLUSION_A).any(axis=1)
            hits.update(int(x) for x in owner[sl][close])
    return hits, chain_lens, tfr, tf


def main():
    human = uniprot_seq(os.path.join(HERE, "P02786.txt"))
    mouse = uniprot_seq(os.path.join(HERE, "Q62351.txt"))
    cons = conservation_map(human, mouse)

    model = load(os.path.join(PDB, "1cx8.cif"), "apo")
    partner, n_contacts, all_chains = pick_dimer_partner(model)
    sasa = sasa_for_chain(model, 'A', partner)
    glyc = glycosylation_sites(model, 'A')
    tf_hits, chain_lens, tfr_chains, tf_chains = tf_footprint(os.path.join(PDB, "1suv.cif"))

    modeled = sorted(sasa.keys())
    last_modeled = modeled[-1]
    cb = {}
    for r in residues(model['A']):
        at = r['CB'] if 'CB' in r else (r['CA'] if 'CA' in r else None)
        if at is not None:
            cb[r.id[1]] = at.coord
    anchor = cb.get(last_modeled)
    tag_near = {n for n, c in cb.items() if anchor is not None and np.linalg.norm(c - anchor) < TAG_EXCLUSION_A}

    # candidate residues
    cand = []
    for n in modeled:
        aa, abs_s, rel = sasa[n]
        if not (CONSTRUCT[0] <= n <= CONSTRUCT[1]):
            continue
        if rel < EXPOSED_REL_SASA:
            continue
        if n in tf_hits or n in glyc or n in tag_near:
            continue
        if WINDOW[1] and not (WINDOW[0] <= n <= WINDOW[1]):
            continue
        cand.append(n)
    cand_set = set(cand)

    # patch membership: any residue with some exposure that is not itself excluded
    allowed = {n for n in modeled
               if CONSTRUCT[0] <= n <= CONSTRUCT[1]
               and sasa[n][2] >= 0.15
               and n not in tf_hits and n not in glyc and n not in tag_near}

    # grow patches around each candidate
    patches = []
    for c0 in cand:
        if c0 not in cb:
            continue
        members = [n for n in allowed if n in cb and np.linalg.norm(cb[n] - cb[c0]) <= PATCH_RADIUS_A]
        if len(members) < 6:
            continue
        seq = "".join(sasa[n][0] for n in sorted(members))
        hyd = sum(1 for ch in seq if ch in HYDROPHOBIC)
        chg = sum(1 for ch in seq if ch in CHARGED)
        identity = float(np.mean([cons.get(n, 0) for n in members]))
        if identity < MIN_IDENT or hyd < MIN_HYD:
            continue
        exposure = float(np.mean([sasa[n][2] for n in members]))
        # distance of the patch centre to the nearest transferrin-footprint residue
        d_tf = min((float(np.linalg.norm(cb[n] - cb[c0])) for n in tf_hits if n in cb), default=float('inf'))
        d_tag = float(np.linalg.norm(cb[c0] - anchor)) if anchor is not None else float('inf')
        his = [n for n in members if sasa[n][0] == 'H']
        # score: hydrophobic content and exposure help; heavy charge hurts (RFdiffusion guidance)
        score = (hyd * 1.0 + exposure * 4.0 - chg * 0.6 + len(members) * 0.15
                 + (1.5 if his else 0.0) + (identity - MIN_IDENT) * 20.0)
        patches.append(dict(center=c0, n=len(members), residues=sorted(members), seq=seq,
                            identity=round(identity, 3),
                            hydrophobic=hyd, charged=chg, mean_rel_sasa=round(exposure, 3),
                            native_his=his, dist_to_tf_footprint=round(d_tf, 1),
                            dist_to_tag_anchor=round(d_tag, 1), score=round(score, 2)))

    patches.sort(key=lambda p: -p['score'])
    chosen, used = [], set()
    for p in patches:
        if len(set(p['residues']) & used) > len(p['residues']) * 0.4:
            continue
        chosen.append(p); used.update(p['residues'])
        if len(chosen) >= TOP_N:
            break

    meta = dict(
        construct=CONSTRUCT, dimer_partner_chain=partner, dimer_contact_atoms=n_contacts,
        chains_in_apo=all_chains, modeled_range=[modeled[0], last_modeled],
        n_modeled=len(modeled), glycosylation_sites=glyc,
        tf_footprint_residues=sorted(tf_hits), tf_complex_chain_lengths=chain_lens,
        tfr_chains=tfr_chains, tf_chains=tf_chains,
        tag_excluded_residues=sorted(tag_near), n_candidates=len(cand),
        thresholds=dict(exposed_rel_sasa=EXPOSED_REL_SASA, patch_radius_a=PATCH_RADIUS_A,
                        tf_cutoff_a=TF_EXCLUSION_A, tag_exclusion_a=TAG_EXCLUSION_A),
    )
    out = dict(meta=meta, patches=chosen)
    with open(os.path.join(HERE, "epitope_scan.json"), "w") as fh:
        json.dump(out, fh, indent=1)

    print("apo chains: %s | dimer partner of A: %s (%d contact atoms)" % (all_chains, partner, n_contacts))
    print("modeled in chain A: %d residues, %d-%d" % (len(modeled), modeled[0], last_modeled))
    print("glycosylation sites (NAG-linked ASN): %s" % glyc)
    print("transferrin footprint: %d TfR1 residues (chains %s vs %s)" % (len(tf_hits), tfr_chains, tf_chains))
    print("tag-vicinity exclusion (<%.0f A of residue %d): %d residues" % (TAG_EXCLUSION_A, last_modeled, len(tag_near)))
    print("candidates after all three masks: %d\n" % len(cand))
    print("%-7s %-4s %-6s %-6s %-4s %-4s %-8s %-8s %-8s %s" % ("center", "n", "score", "ident", "hyd", "chg", "relSASA", "dTf(A)", "dTag(A)", "nativeHis"))
    for p in chosen:
        print("%-7d %-4d %-6.2f %-6.2f %-4d %-4d %-8.3f %-8.1f %-8.1f %s" % (
            p['center'], p['n'], p['score'], p['identity'], p['hydrophobic'], p['charged'],
            p['mean_rel_sasa'], p['dist_to_tf_footprint'], p['dist_to_tag_anchor'], p['native_his'] or '-'))
    print("\nTop patch residue lists:")
    for p in chosen[:5]:
        print("  center %d: %s" % (p['center'], ",".join(str(x) for x in p['residues'])))
        print("            seq %s" % p['seq'])


if __name__ == "__main__":
    main()
