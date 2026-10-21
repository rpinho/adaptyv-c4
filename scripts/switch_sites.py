#!/usr/bin/env python
"""Rank binder positions by TWO-STATE switchability, not by burial.

Why this exists (measured 2026-10-20, results/c3_ranked):
a real design's two histidines were depressed by 2.67 and 4.33 pKa units on binding -
enormous shifts - and contributed only -0.52 and -0.39 to the proton linkage L, because
their FREE-state pKa was 6.43 and 6.20. Meanwhile three target glutamates whose pKa rose
into the pH 6.0-7.4 window contributed +1.91 and flipped the sign of the whole design.

The lesson: a group contributes to L only in proportion to how much of the 6.0-7.4 window
its titration crosses. A large dpKa that starts and ends below the window is worth nothing
(measured: A67 ASP, dpKa -1.57, contribution -0.000; A88 TYR, dpKa +3.84, contribution
+0.000). The leverage is in the FREE-state pKa:

    binder His free pKa 6.0, driven to 2.0  ->  -0.28   (what burial alone buys)
    binder His free pKa 7.4, driven to 3.5  ->  -1.12   (4x more, same burial)

So a switch position must be elevated-pKa when free AND depressed when bound. Ranking by
interface burial - what install_switch.py did - only optimises the second half.

A histidine's pKa is raised by nearby carboxylates and by solvent exposure, and depressed
by desolvation and by nearby cations. This scores both states from coordinates alone:

    free state (binder alone)   : exposed + acidic neighbours  -> high pKa
    bound state (complex)       : buried  + basic  neighbours  -> low  pKa

No mutation modelling and no PyRosetta, so it runs on CPU while the GPUs design. Its output
is the prioritised position list to hand to install_switch.py, which then does the
expensive Rosetta work on a short list instead of on every interface position.

Usage:
  python switch_sites.py design.pdb --binder-chain B --target-chain A
  python switch_sites.py dir_of_pdbs/ --csv sites.csv --top 8
"""
import argparse, csv, math, os

ACID = {"ASP", "GLU"}
BASE = {"LYS", "ARG"}
ANCHOR_AA = {"TRP", "PHE", "TYR", "LEU", "ILE", "MET", "VAL"}
SKIP_AA = {"GLY", "PRO", "CYS"}
PH_LOW, PH_HIGH = 6.0, 7.4
# a bare histidine in water; propka3's model value
HIS_MODEL_PKA = 6.5
# propka-style empirical weights, used here only to RANK positions, never as a pKa claim
W_ACID_FREE = 0.45     # per carboxylate within 8 A, raises the free-state pKa
W_BASE_BOUND = -0.40   # per cation within 8 A in the complex, lowers the bound pKa
W_DESOLV = -2.60       # per unit of fractional burial gained on binding


def atoms(path):
    out = []
    with open(path) as fh:
        for line in fh:
            if line.startswith("ATOM"):
                out.append((line[21], int(line[22:26]), line[17:20].strip(),
                            line[12:16].strip(),
                            float(line[30:38]), float(line[38:46]), float(line[46:54])))
    return out


def by_residue(ats):
    res = {}
    for ch, num, name, atom, x, y, z in ats:
        res.setdefault((ch, num), [name, []])[1].append((atom, x, y, z))
    return res


def d2(a, b):
    return (a[0]-b[0])**2 + (a[1]-b[1])**2 + (a[2]-b[2])**2


def sidechain_centre(entry):
    name, al = entry
    pts = [(x, y, z) for atom, x, y, z in al if atom not in ("N", "C", "O")]
    if not pts:
        pts = [(x, y, z) for atom, x, y, z in al]
    n = len(pts)
    return (sum(p[0] for p in pts)/n, sum(p[1] for p in pts)/n, sum(p[2] for p in pts)/n)


def neighbour_count(centre, res, exclude, radius=8.0, kinds=None):
    r2 = radius * radius
    n = 0
    for key, entry in res.items():
        if key in exclude:
            continue
        if kinds and entry[0] not in kinds:
            continue
        if any(d2(centre, (x, y, z)) <= r2 for atom, x, y, z in entry[1]):
            n += 1
    return n


def burial(centre, res, exclude, radius=10.0):
    """Heavy atoms within 10 A, normalised - a cheap, monotonic proxy for burial."""
    r2 = radius * radius
    n = 0
    for key, entry in res.items():
        if key in exclude:
            continue
        for atom, x, y, z in entry[1]:
            if d2(centre, (x, y, z)) <= r2:
                n += 1
    return min(1.0, n / 180.0)


def term(pb, pf, ph):
    return math.log10(1 + 10 ** (pb - ph)) - math.log10(1 + 10 ** (pf - ph))


def contribution(pf, pb):
    return term(pb, pf, PH_LOW) - term(pb, pf, PH_HIGH)


def score_design(path, binder_chain, target_chain):
    ats = atoms(path)
    res = by_residue(ats)
    binder = {k: v for k, v in res.items() if k[0] in binder_chain}
    target = {k: v for k, v in res.items() if k[0] in target_chain}
    if not binder or not target:
        raise SystemExit("%s: chains %s / %s not both present"
                         % (path, binder_chain, target_chain))

    rows = []
    for key, entry in sorted(binder.items(), key=lambda kv: kv[0][1]):
        name = entry[0]
        c = sidechain_centre(entry)
        # is this position at the interface at all?
        if not any(d2(c, (x, y, z)) <= 100.0 for e in target.values()
                   for atom, x, y, z in e[1]):
            continue

        b_free = burial(c, binder, {key})
        b_bound = burial(c, res, {key})
        d_burial = b_bound - b_free

        acids_free = neighbour_count(c, binder, {key}, 8.0, ACID)
        bases_bound = neighbour_count(c, res, {key}, 8.0, BASE)
        bases_target = neighbour_count(c, target, set(), 8.0, BASE)

        # a ranking estimate, NOT a pKa prediction: free state raised by its own
        # carboxylates and by exposure, bound state lowered by desolvation and cations
        pka_free_est = HIS_MODEL_PKA + W_ACID_FREE * acids_free + W_DESOLV * b_free * 0.35
        pka_bound_est = pka_free_est + W_DESOLV * max(0.0, d_burial) + W_BASE_BOUND * bases_target
        est = contribution(pka_free_est, pka_bound_est)

        flag = ""
        if name in SKIP_AA:
            flag = "skip (backbone/cys role)"
        elif name in ANCHOR_AA and b_bound > 0.55:
            flag = "anchor - mutating costs binding"
        rows.append(dict(design=os.path.basename(path), chain=key[0], resnum=key[1],
                         resname=name,
                         burial_free=round(b_free, 3), burial_bound=round(b_bound, 3),
                         d_burial=round(d_burial, 3),
                         acids_near_free=acids_free, bases_near_target=bases_target,
                         pka_free_est=round(pka_free_est, 2),
                         pka_bound_est=round(pka_bound_est, 2),
                         est_contribution=round(est, 3), note=flag))
    rows.sort(key=lambda r: (r["note"] != "", r["est_contribution"]))
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("path", help="a complex PDB, or a directory of them")
    ap.add_argument("--binder-chain", default="B")
    ap.add_argument("--target-chain", default="A")
    ap.add_argument("--top", type=int, default=10)
    ap.add_argument("--csv")
    a = ap.parse_args()

    pdbs = ([os.path.join(a.path, f) for f in sorted(os.listdir(a.path)) if f.endswith(".pdb")]
            if os.path.isdir(a.path) else [a.path])
    if not pdbs:
        raise SystemExit("no PDB files at %s" % a.path)

    print("Two-state switch ranking. A position is good when it is EXPOSED with acidic")
    print("neighbours in the free binder (high pKa) and BURIED against the target's cations")
    print("when bound (low pKa). est_contribution is a RANKING, not a pKa prediction -")
    print("propka3 via ph_linkage.py is what decides a design.\n")

    allrows = []
    for p in pdbs:
        rows = score_design(p, a.binder_chain, a.target_chain)
        allrows += rows
        usable = [r for r in rows if not r["note"]]
        print("%s  (%d interface positions, %d usable)" % (os.path.basename(p), len(rows), len(usable)))
        print("  %-6s %-5s %6s %6s %7s %6s %6s %8s %8s %8s"
              % ("pos", "aa", "b_free", "b_bnd", "dburial", "acids", "bases",
                 "pKa_f", "pKa_b", "est_L"))
        for r in usable[:a.top]:
            print("  %s%-5d %-5s %6.2f %6.2f %7.2f %6d %6d %8.2f %8.2f %+8.3f"
                  % (r["chain"], r["resnum"], r["resname"], r["burial_free"],
                     r["burial_bound"], r["d_burial"], r["acids_near_free"],
                     r["bases_near_target"], r["pka_free_est"], r["pka_bound_est"],
                     r["est_contribution"]))
        if usable:
            best = usable[:2]
            print("  best pair, if both installed: est L = %+.2f  %s"
                  % (sum(r["est_contribution"] for r in best),
                     "" if sum(r["est_contribution"] for r in best) <= -2
                     else "(short of the -2 gate on its own; needs a third site or a"
                          " charge-compensated interface)"))
        print()

    if a.csv:
        with open(a.csv, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(allrows[0].keys()))
            w.writeheader()
            w.writerows(allrows)
        print("wrote %s (%d rows)" % (a.csv, len(allrows)))


if __name__ == "__main__":
    main()
