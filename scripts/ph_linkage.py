#!/usr/bin/env python
"""Proton linkage score for a designed binder-target complex.

This is the entry's design metric. The pH dependence of binding comes entirely from the
difference in proton affinity between the bound and the free states (Wyman linkage), so
for every titratable group i with pKa_f free and pKa_b in the complex:

    K_a(pH) = K_a(ref) * PROD_i (1 + 10^(pKa_i,b - pH)) / (1 + 10^(pKa_i,f - pH))

and therefore the quantity the challenge's first criterion cares about is

    L = log10 K_a(pH_low) - log10 K_a(pH_high)
      = SUM_i [ log10((1+10^(pKa_b-pH_low)) /(1+10^(pKa_f-pH_low)))
               - log10((1+10^(pKa_b-pH_high))/(1+10^(pKa_f-pH_high))) ]

L is in log10 units of affinity. A binder that holds at pH 7.4 and lets go at pH 6.0
needs L <= -2 (a hundredfold weaker binding at the lower pH). L < 0 requires the
protonatable groups to be LESS protonated in the complex than when free, i.e. pKa
depressed on binding - which is what burying a neutral histidine against a
desolvated, apolar or positively polarised interface does.

pKa values come from PROPKA3 run three times: on the complex, on the isolated binder and
on the isolated target. Residues are matched by (chain, number). Groups whose pKa barely
moves contribute nothing, so the sum is dominated by real interface switches.

Usage:
  python ph_linkage.py complex.pdb --binder-chain B --target-chain A
  python ph_linkage.py dir_of_pdbs/ --binder-chain B --target-chain A --csv out.csv
"""
import argparse, csv, math, os, re, shutil, subprocess, sys, tempfile

TITRATABLE = {"HIS", "ASP", "GLU", "LYS", "TYR", "CYS", "ARG", "N+", "C-"}
PKA_LINE = re.compile(r"^\s*([A-Z0-9+-]{1,3})\s+(\d+)\s+([A-Z])\s+(-?\d+\.\d+)")


def propka_exe():
    """propka3 next to the running interpreter, so the env need not be activated."""
    cand = os.path.join(os.path.dirname(os.path.abspath(sys.executable)), "propka3")
    return cand if os.path.exists(cand) else "propka3"


def run_propka(pdb_path, workdir):
    """Run propka3 and return {(chain, resnum, resname): pKa}."""
    dst = os.path.join(workdir, os.path.basename(pdb_path))
    if os.path.abspath(dst) != os.path.abspath(pdb_path):
        shutil.copy(pdb_path, dst)
    try:
        subprocess.run([propka_exe(), os.path.basename(dst)], cwd=workdir,
                       check=True, capture_output=True, timeout=1800)
    except subprocess.CalledProcessError as e:
        sys.stderr.write("propka3 failed on %s: %s\n"
                         % (pdb_path, e.stderr.decode()[-400:]))
        return None
    pka_file = os.path.splitext(dst)[0] + ".pka"
    if not os.path.exists(pka_file):
        return None
    out, in_summary = {}, False
    for line in open(pka_file):
        if line.startswith("SUMMARY OF THIS PREDICTION"):
            in_summary = True
            continue
        if in_summary:
            if line.startswith("---") or not line.strip():
                if out:
                    break
                continue
            m = PKA_LINE.match(line)
            if m:
                resname, num, chain, pka = m.group(1), int(m.group(2)), m.group(3), float(m.group(4))
                out[(chain, num, resname)] = pka
    return out


def split_chains(pdb_path, keep, out_path):
    with open(pdb_path) as fh, open(out_path, "w") as out:
        for line in fh:
            if line.startswith(("ATOM", "HETATM")):
                if line[21] in keep:
                    out.write(line)
            elif line.startswith("TER"):
                out.write(line)
        out.write("END\n")


def term(pka_b, pka_f, ph):
    return (math.log10(1 + 10 ** (pka_b - ph)) - math.log10(1 + 10 ** (pka_f - ph)))


def linkage(bound, free, ph_low, ph_high, min_shift):
    """Return (L, per-residue contributions) for groups present in both states."""
    rows, total = [], 0.0
    for key, pb in bound.items():
        if key not in free:
            continue
        pf = free[key]
        if abs(pb - pf) < min_shift:
            continue
        contrib = term(pb, pf, ph_low) - term(pb, pf, ph_high)
        total += contrib
        rows.append(dict(chain=key[0], resnum=key[1], resname=key[2],
                         pka_free=round(pf, 2), pka_bound=round(pb, 2),
                         dpka=round(pb - pf, 2), contribution=round(contrib, 3)))
    rows.sort(key=lambda r: r["contribution"])
    return total, rows


def score_one(pdb, binder_chain, target_chain, ph_low, ph_high, min_shift, verbose=True):
    with tempfile.TemporaryDirectory() as td:
        cx = run_propka(pdb, td)
        if not cx:
            return None
        bp = os.path.join(td, "binder.pdb")
        tp = os.path.join(td, "target.pdb")
        split_chains(pdb, set(binder_chain), bp)
        split_chains(pdb, set(target_chain), tp)
        fb = run_propka(bp, td) or {}
        ft = run_propka(tp, td) or {}
        free = dict(fb); free.update(ft)
        L, rows = linkage(cx, free, ph_low, ph_high, min_shift)
        binder_his = [r for r in rows if r["resname"] == "HIS" and r["chain"] in binder_chain]
        target_his = [r for r in rows if r["resname"] == "HIS" and r["chain"] in target_chain]
        res = dict(design=os.path.basename(pdb), L=round(L, 3),
                   n_groups=len(rows), n_binder_his=len(binder_his),
                   n_target_his=len(target_his),
                   switch_ok=bool(L <= -2.0), rows=rows)
        if verbose:
            print("%-44s L=%+.2f  groups=%d  binderHis=%d  targetHis=%d  %s"
                  % (res["design"], L, len(rows), len(binder_his), len(target_his),
                     "PASS" if res["switch_ok"] else ""))
            show = rows[:6] + ([None] + rows[-6:] if len(rows) > 12 else rows[6:])
            for r in show:
                if r is None:
                    print("    ... the other end of the sum (these push L the WRONG way) ...")
                    continue
                print("    %s%-4d %-4s pKa %5.2f free -> %5.2f bound (d%+.2f)  contributes %+.3f"
                      % (r["chain"], r["resnum"], r["resname"], r["pka_free"],
                         r["pka_bound"], r["dpka"], r["contribution"]))
        return res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("path", help="a complex PDB, or a directory of them")
    ap.add_argument("--binder-chain", default="B")
    ap.add_argument("--target-chain", default="A")
    ap.add_argument("--ph-low", type=float, default=6.0)
    ap.add_argument("--ph-high", type=float, default=7.4)
    ap.add_argument("--min-shift", type=float, default=0.1)
    ap.add_argument("--csv")
    a = ap.parse_args()

    pdbs = []
    if os.path.isdir(a.path):
        pdbs = [os.path.join(a.path, f) for f in sorted(os.listdir(a.path))
                if f.endswith(".pdb")]
    else:
        pdbs = [a.path]
    if not pdbs:
        raise SystemExit("no PDB files found at %s" % a.path)

    print("proton linkage L = log10 Ka(%.1f) - log10 Ka(%.1f); the pH gate needs L <= -2\n"
          % (a.ph_low, a.ph_high))
    results = [r for r in (score_one(p, a.binder_chain, a.target_chain,
                                     a.ph_low, a.ph_high, a.min_shift) for p in pdbs) if r]
    if not results:
        raise SystemExit("nothing scored")
    results.sort(key=lambda r: r["L"])
    print("\nranked (most pH-selective first):")
    for r in results:
        print("  %-44s L=%+.2f %s" % (r["design"], r["L"], "PASS" if r["switch_ok"] else ""))
    passing = [r for r in results if r["switch_ok"]]
    print("\n%d of %d designs clear L <= -2" % (len(passing), len(results)))

    if a.csv:
        with open(a.csv, "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["design", "L", "n_groups", "n_binder_his", "n_target_his", "switch_ok"])
            for r in results:
                w.writerow([r["design"], r["L"], r["n_groups"], r["n_binder_his"],
                            r["n_target_his"], r["switch_ok"]])
        print("wrote %s" % a.csv)
        rp = os.path.splitext(a.csv)[0] + "_residues.csv"
        with open(rp, "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["design", "chain", "resnum", "resname", "pka_free",
                        "pka_bound", "dpka", "contribution"])
            for r in results:
                for row in r["rows"]:
                    w.writerow([r["design"], row["chain"], row["resnum"], row["resname"],
                                row["pka_free"], row["pka_bound"], row["dpka"],
                                row["contribution"]])
        print("wrote %s" % rp)


if __name__ == "__main__":
    main()
