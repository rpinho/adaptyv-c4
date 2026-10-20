#!/usr/bin/env python
"""Trim the target to an epitope neighbourhood and emit a generator-ready PDB.

The generators' runtime scales with target size, so the target is cut down to the
epitope plus its surroundings. Residues are selected by distance, then merged into
contiguous intervals (short fragments dropped, ends padded) so the trimmed target has
as few chain breaks as possible. Both chains of the physiological dimer are considered,
so a binder is never designed into space that the partner subunit occupies.

Usage:
  python trim_target.py --set A_apical_derisked --radius 12 --out ../targets/E2.pdb
  python trim_target.py --residues 688,690,695,704,705,706,711,714,730,733,737 --out ../targets/E1.pdb
"""
import argparse, json, os
import numpy as np
from Bio.PDB import MMCIFParser, PDBIO, Select
from Bio.PDB.Polypeptide import is_aa

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(HERE, "data")
TEMPLATE = "11qc"
MIN_FRAGMENT = 4
PAD = 2


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--set", dest="setname")
    ap.add_argument("--residues")
    ap.add_argument("--radius", type=float, default=12.0)
    ap.add_argument("--chains", default="A,B")
    ap.add_argument("--epitope-chain", default="A")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    if a.residues:
        hot = [int(x) for x in a.residues.split(",")]
        label = "custom"
    else:
        sets = json.load(open(os.path.join(HERE, "results", "hotspot_sets.json")))
        if a.setname not in sets:
            raise SystemExit("unknown set %r; have %s" % (a.setname, list(sets)))
        hot = sets[a.setname]["residues"]
        label = a.setname

    model = MMCIFParser(QUIET=True).get_structure(TEMPLATE,
                        os.path.join(DATA, TEMPLATE + ".cif"))[0]
    chains = a.chains.split(",")
    res_by_chain = {c: [r for r in model[c] if is_aa(r, standard=True)] for c in chains}

    hot_atoms = np.array([at.coord for r in res_by_chain[a.epitope_chain]
                          if r.id[1] in hot for at in r if at.element != 'H'])
    if not len(hot_atoms):
        raise SystemExit("no atoms found for the requested epitope residues")

    keep = {}
    for c, rl in res_by_chain.items():
        sel = []
        for r in rl:
            pts = np.array([at.coord for at in r if at.element != 'H'])
            if not len(pts):
                continue
            d = np.linalg.norm(pts[:, None, :] - hot_atoms[None, :, :], axis=-1).min()
            if d <= a.radius:
                sel.append(r.id[1])
        keep[c] = sorted(sel)

    # merge into intervals, drop short fragments, pad the ends
    merged = {}
    for c, nums in keep.items():
        if not nums:
            merged[c] = []
            continue
        ivs, start, prev = [], nums[0], nums[0]
        for n in nums[1:]:
            if n - prev <= 2:
                prev = n
            else:
                ivs.append((start, prev)); start = prev = n
        ivs.append((start, prev))
        present = {r.id[1] for r in res_by_chain[c]}
        out = []
        for lo, hi in ivs:
            if hi - lo + 1 < MIN_FRAGMENT:
                continue
            lo2, hi2 = lo - PAD, hi + PAD
            out.append([n for n in range(lo2, hi2 + 1) if n in present])
        merged[c] = sorted({n for frag in out for n in frag})

    class Keep(Select):
        def accept_residue(self, residue):
            c = residue.get_parent().id
            return c in merged and residue.id[1] in merged[c] and is_aa(residue, standard=True)

        def accept_atom(self, atom):
            return atom.element != 'H'

    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    io = PDBIO(); io.set_structure(model); io.save(a.out, select=Keep())

    # distance from the epitope to the partner chain, to show whether it is even relevant
    partner = [c for c in chains if c != a.epitope_chain]
    d_partner = float("inf")
    for c in partner:
        pts = np.array([at.coord for r in res_by_chain[c] for at in r if at.element != 'H'])
        if len(pts):
            d_partner = min(d_partner,
                            float(np.linalg.norm(hot_atoms[:, None, :] - pts[None, :, :], axis=-1).min()))

    print("epitope set: %s (%d residues: %s)" % (label, len(hot), ",".join(map(str, hot))))
    print("radius %.1f A, template %s" % (a.radius, TEMPLATE.upper()))
    for c in chains:
        n = len(merged.get(c, []))
        if n:
            frags, cur = [], []
            for x in merged[c]:
                if cur and x - cur[-1] > 1:
                    frags.append((cur[0], cur[-1])); cur = []
                cur.append(x)
            frags.append((cur[0], cur[-1]))
            print("  chain %s: %3d residues in %d fragment(s): %s"
                  % (c, n, len(frags), ", ".join("%d-%d" % f for f in frags)))
        else:
            print("  chain %s: none" % c)
    print("closest approach of the epitope to chain %s: %.1f A" % (",".join(partner), d_partner))
    print("wrote %s" % a.out)
    print("BindCraft hotspots: %s" % ",".join("A%d" % r for r in hot))


if __name__ == "__main__":
    main()
