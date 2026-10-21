#!/usr/bin/env python
"""Acid/base census of a candidate epitope - the term the linkage budget forgot.

Measured 2026-10-20 on a real BindCraft2 design (results/c3_ranked): two designed
binder histidines contributed -0.91 to the proton linkage L, but three TARGET
glutamates had their pKa RAISED on binding (one by +2.99, from 4.92 to 7.91) and
contributed +1.91, flipping L positive. A carboxylate whose pKa rises into the
6-8 window accepts a proton on binding, so binding gets STRONGER at low pH -
the opposite of what Challenge 4's gate asks for.

The target's own side chains cannot be mutated, so the only lever is epitope
choice. This counts, per epitope, the titratable groups a binder would bury:
Asp/Glu (push L the wrong way when buried and desolvated) against His (the
intended lever) and Lys/Arg/Tyr.

Usage:
  python epitope_charge.py targets/E1_his707_tight.pdb --hotspots A695,A704,A705,A706,A707,A711,A714
  python epitope_charge.py targets/ --hotspots-from campaign_hotspots.json
"""
import argparse, json, os, sys
from collections import Counter

ACID = {"ASP", "GLU"}
BASE = {"LYS", "ARG"}
SWITCH = {"HIS"}
OTHER_TITRATABLE = {"TYR", "CYS"}


def residues(path):
    """[(chain, resnum, resname, [(x,y,z) sidechain-ish atoms])] from a PDB."""
    out, cur, key = [], None, None
    for line in open(path):
        if not line.startswith("ATOM"):
            continue
        ch, num, name = line[21], int(line[22:26]), line[17:20].strip()
        k = (ch, num)
        if k != key:
            cur = (ch, num, name, [])
            out.append(cur)
            key = k
        atom = line[12:16].strip()
        if atom not in ("N", "C", "O"):           # keep CA and side chain
            cur[3].append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
    return out


def d2(a, b):
    return (a[0]-b[0])**2 + (a[1]-b[1])**2 + (a[2]-b[2])**2


def census(path, hotspots, radius):
    res = residues(path)
    hs = set()
    for h in hotspots:
        h = h.strip()
        if not h:
            continue
        hs.add((h[0], int(h[1:])))
    anchors = [r for r in res if (r[0], r[1]) in hs]
    missing = hs - {(r[0], r[1]) for r in anchors}
    pts = [p for r in anchors for p in r[3]]
    if not pts:
        raise SystemExit("%s: none of the hotspots %s are in the file" % (path, sorted(hs)))
    near, r2 = [], radius * radius
    for r in res:
        if (r[0], r[1]) in hs:
            near.append(r)
            continue
        if any(d2(p, q) <= r2 for p in r[3] for q in pts):
            near.append(r)
    c = Counter(r[2] for r in near)
    acids = [r for r in near if r[2] in ACID]
    his = [r for r in near if r[2] in SWITCH]
    bases = [r for r in near if r[2] in BASE]
    return dict(path=os.path.basename(path), n_patch=len(near),
                n_acid=len(acids), n_his=len(his), n_base=len(bases),
                n_tyr_cys=sum(c[x] for x in OTHER_TITRATABLE),
                acid_list=["%s%d %s" % (r[0], r[1], r[2]) for r in acids],
                his_list=["%s%d" % (r[0], r[1]) for r in his],
                base_list=["%s%d %s" % (r[0], r[1], r[2]) for r in bases],
                missing=sorted("%s%d" % m for m in missing),
                net_formal=len(bases) - len(acids))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("target", help="a target PDB, or a directory of them")
    ap.add_argument("--hotspots", help="comma-separated, e.g. A695,A704,A707")
    ap.add_argument("--hotspots-from", help="JSON {basename: 'A695,A704,...'}")
    ap.add_argument("--radius", type=float, default=10.0,
                    help="a binder buries roughly this far from the hotspots (default 10 A)")
    ap.add_argument("--json")
    a = ap.parse_args()

    jobs = []
    if os.path.isdir(a.target):
        table = json.load(open(a.hotspots_from)) if a.hotspots_from else {}
        for f in sorted(os.listdir(a.target)):
            if f.endswith(".pdb") and (f in table or a.hotspots):
                jobs.append((os.path.join(a.target, f),
                             (table.get(f) or a.hotspots).split(",")))
    else:
        if not a.hotspots:
            raise SystemExit("--hotspots is required for a single target")
        jobs.append((a.target, a.hotspots.split(",")))
    if not jobs:
        raise SystemExit("nothing to do: no hotspots given for any target")

    print("Within %.0f A of the hotspots. Asp/Glu push the proton linkage L POSITIVE when\n"
          "buried (wrong way for the pH gate); His is the intended lever.\n" % a.radius)
    print("%-30s %5s %5s %5s %5s %6s" % ("target", "patch", "acid", "His", "base", "net"))
    rows = []
    for path, hs in jobs:
        r = census(path, hs, a.radius)
        rows.append(r)
        print("%-30s %5d %5d %5d %5d %+6d" % (r["path"][:30], r["n_patch"], r["n_acid"],
                                              r["n_his"], r["n_base"], r["net_formal"]))
    print()
    for r in rows:
        print("%s" % r["path"])
        print("  acids buried : %s" % (", ".join(r["acid_list"]) or "none"))
        print("  His in patch : %s" % (", ".join(r["his_list"]) or "none"))
        print("  bases        : %s" % (", ".join(r["base_list"]) or "none"))
        if r["missing"]:
            print("  NOT IN FILE  : %s" % ", ".join(r["missing"]))
        ratio = r["n_acid"] / max(1, r["n_his"])
        print("  acid : His = %.1f  %s" % (ratio,
              "acid-dominated, the His term will be swamped" if ratio >= 3 else
              "workable" if ratio <= 1.5 else "marginal"))
    if a.json:
        json.dump(rows, open(a.json, "w"), indent=1)
        print("\nwrote %s" % a.json)


if __name__ == "__main__":
    main()
