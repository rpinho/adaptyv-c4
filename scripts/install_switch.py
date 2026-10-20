#!/usr/bin/env python
"""Install pH switches into an accepted binder, as an optimisation rather than an edit.

For one designed complex (target chain A, binder chain B) this:
  1. finds the binder's interface positions,
  2. ranks them by how buried they are at the interface - a histidine that gets buried
     against a desolvated interface is the one whose pKa drops on binding, which is what
     makes the switch - skipping positions that carry the interface's hydrophobic core,
  3. builds, for each selected position, a HIS mutant (the switch) and a GLN mutant (the
     isosteric, non-titratable matched control twin),
  4. repacks and minimises locally, then scores each variant's binding energy with
     Rosetta's InterfaceAnalyzer so designs that lose their pH 7.4 binding are dropped,
  5. writes every surviving variant as a complex PDB for ph_linkage.py to score.

The pH half of the loop is deliberately NOT Rosetta's: it comes from PROPKA linkage in
ph_linkage.py, so the two estimators that decide a design are independent of each other.

Run inside the BindCraft conda env (it provides PyRosetta).

Usage:
  python install_switch.py design.pdb --out-dir variants/ [--max-positions 8]
"""
import argparse, json, os, sys

import pyrosetta
from pyrosetta import pose_from_pdb
from pyrosetta.rosetta.core.scoring import ScoreType
from pyrosetta.rosetta.protocols.analysis import InterfaceAnalyzerMover
from pyrosetta.rosetta.core.select.residue_selector import (
    ChainSelector, InterGroupInterfaceByVectorSelector)

# positions whose side chain is the interface's hydrophobic anchor are left alone:
# mutating them away usually costs more binding than the switch is worth
ANCHOR_AA = set("WFYLIMV")
SKIP_AA = set("GPC")          # glycine/proline backbone roles, cysteine pairing


def init():
    pyrosetta.init(" ".join([
        "-mute all",
        "-ex1", "-ex2aro",
        "-use_input_sc",
        "-ignore_unrecognized_res",
        "-ignore_zero_occupancy false",
    ]))


def interface_positions(pose, target_chain="A", binder_chain="B"):
    sel = InterGroupInterfaceByVectorSelector()
    sel.group1_selector(ChainSelector(target_chain))
    sel.group2_selector(ChainSelector(binder_chain))
    sel.cb_dist_cut(11.0)
    mask = sel.apply(pose)
    info = pose.pdb_info()
    out = []
    for i in range(1, pose.total_residue() + 1):
        if mask[i] and info.chain(i) == binder_chain:
            out.append(i)
    return out


def burial(pose, resi):
    """Side-chain neighbour count: a cheap, robust proxy for interface burial."""
    xyz = pose.residue(resi).nbr_atom_xyz()
    n = 0
    for j in range(1, pose.total_residue() + 1):
        if j == resi:
            continue
        if xyz.distance(pose.residue(j).nbr_atom_xyz()) < 10.0:
            n += 1
    return n


def dG_separated(pose, target_chain="A", binder_chain="B"):
    iam = InterfaceAnalyzerMover()
    iam.set_interface("%s_%s" % (target_chain, binder_chain))
    iam.set_compute_packstat(False)
    iam.set_pack_separated(True)
    iam.apply(pose)
    data = iam.get_all_data()
    return float(iam.get_separated_interface_energy()), float(data.sasa[1] if data.sasa else 0.0)


def mutate_and_relax(pose, resi, aa, scorefxn, radius=8.0):
    from pyrosetta.toolbox import mutate_residue
    work = pose.clone()
    mutate_residue(work, resi, aa, pack_radius=radius, pack_scorefxn=scorefxn)
    mm = pyrosetta.rosetta.core.kinematics.MoveMap()
    mm.set_bb(False)
    mm.set_chi(False)
    info = work.pdb_info()
    centre = work.residue(resi).nbr_atom_xyz()
    for j in range(1, work.total_residue() + 1):
        if centre.distance(work.residue(j).nbr_atom_xyz()) < radius:
            mm.set_chi(j, True)
    minmover = pyrosetta.rosetta.protocols.minimization_packing.MinMover(
        mm, scorefxn, "lbfgs_armijo_nonmonotone", 0.001, True)
    minmover.apply(work)
    return work


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("pdb")
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--target-chain", default="A")
    ap.add_argument("--binder-chain", default="B")
    ap.add_argument("--max-positions", type=int, default=8)
    ap.add_argument("--max-ddg", type=float, default=3.0,
                    help="discard a variant that loses more than this much binding energy")
    a = ap.parse_args()

    init()
    scorefxn = pyrosetta.get_fa_scorefxn()
    pose = pose_from_pdb(a.pdb)
    info = pose.pdb_info()
    base_dg, _ = dG_separated(pose, a.target_chain, a.binder_chain)
    print("%s: interface dG_separated = %.2f" % (os.path.basename(a.pdb), base_dg))

    positions = interface_positions(pose, a.target_chain, a.binder_chain)
    scored = []
    for i in positions:
        aa = pose.residue(i).name1()
        if aa in SKIP_AA:
            continue
        scored.append((burial(pose, i), i, aa))
    scored.sort(reverse=True)

    chosen = [(i, aa) for b, i, aa in scored if aa not in ANCHOR_AA][:a.max_positions]
    anchors = [(i, aa) for b, i, aa in scored if aa in ANCHOR_AA]
    print("interface positions on chain %s: %d | candidates: %s | anchors left alone: %s"
          % (a.binder_chain, len(positions),
             ",".join("%s%d" % (aa, info.number(i)) for i, aa in chosen) or "none",
             ",".join("%s%d" % (aa, info.number(i)) for i, aa in anchors[:6]) or "none"))

    os.makedirs(a.out_dir, exist_ok=True)
    stem = os.path.splitext(os.path.basename(a.pdb))[0]
    rows = []
    for resi, wt in chosen:
        num = info.number(resi)
        for aa, kind in (("H", "switch"), ("Q", "control")):
            if wt == aa:
                continue
            try:
                var = mutate_and_relax(pose, resi, aa, scorefxn)
                dg, _ = dG_separated(var, a.target_chain, a.binder_chain)
            except Exception as e:
                print("  %s%d%s FAILED: %s" % (wt, num, aa, type(e).__name__))
                continue
            ddg = dg - base_dg
            name = "%s_%s%d%s_%s" % (stem, wt, num, aa, kind)
            keep = ddg <= a.max_ddg
            if keep:
                out = os.path.join(a.out_dir, name + ".pdb")
                var.dump_pdb(out)
            rows.append(dict(design=stem, position=num, wt=wt, mut=aa, kind=kind,
                             dG=round(dg, 2), ddG=round(ddg, 2), kept=bool(keep),
                             pdb=name + ".pdb" if keep else None))
            print("  %-22s dG %8.2f  ddG %+6.2f  %s"
                  % (name.replace(stem + "_", ""), dg, ddg, "" if keep else "DISCARDED"))

    manifest = os.path.join(a.out_dir, stem + "_variants.json")
    json.dump(dict(design=stem, source_pdb=os.path.abspath(a.pdb),
                   base_dG=round(base_dg, 2), variants=rows), open(manifest, "w"), indent=1)
    kept = [r for r in rows if r["kept"]]
    print("kept %d of %d variants (ddG <= %.1f); manifest %s"
          % (len(kept), len(rows), a.max_ddg, manifest))
    print("next: python ph_linkage.py %s --binder-chain %s --target-chain %s"
          % (a.out_dir, a.binder_chain, a.target_chain))


if __name__ == "__main__":
    main()
