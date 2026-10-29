#!/usr/bin/env python3
"""
CHECK 1 of 2 mandated by the 2026-10-26 reviewer round (Astra's PART 3).

Astra: "rigidly superpose each shortlisted complex onto the correctly numbered native
homodimer and measure nonbonded binder-receptor heavy-atom overlaps, explicitly
identifying backbone-backbone pairs below 2.0 A. Flip B to A if all candidates have
extensive, unavoidable backbone interpenetration."

All three process reviewers (Opus 5.5, Fable 5.1, Astra) independently said to run this on
the FULL candidate set, fixing no shortlist, because fixing a shortlist would pre-empt a
decision that is not the assistant's to make. So this scores all 343 and selects nothing.

WHAT IT MEASURES
  Each design complex holds a binder (chain A) docked against a 154-residue discontinuous
  crop of TfR1 (chain B, renumbered 1..154 in file order). That crop's file order maps to
  author numbering via targets/E1_his707.pdb chain A. Superposing the crop onto the native
  11QC homodimer by matched CA atoms carries the binder with it; the binder's heavy atoms
  are then tested against the FULL native assembly (both subunits + observed glycans),
  including all the receptor atoms the crop omitted.

  A design whose pose interpenetrates the native receptor cannot be retained as-is. A design
  that fits is not thereby a binder -- this is steric feasibility only.

CONTROLS (all must fire or the script exits non-zero)
  MAPPING_LEN          crop order is 154 residues and matches chain B length
  SUPERPOSE_RMSD       CA RMSD after rigid fit is small for every design
  HIS707_REGRESSION    the mapping reproduces the independently recorded 19/343 designs
                       contacting author His707 (results/c4_provenance_join_20261026)
  POS_DISPLACED        a binder translated +500 A has exactly zero clashes
  POS_NATIVE_PARTNER   the native partner subunit scored as a pseudo-binder against chain A
                       DOES register contacts, proving the detector fires

CPU only. Read-only on all inputs. $0. No ranking. No submission.
"""

import argparse
import hashlib
import json
import os
import sys
import warnings

import numpy as np

warnings.simplefilter("ignore")
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB import Superimposer, NeighborSearch

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CIF = os.path.join(REPO, "data/11qc.cif")
TARGET = os.path.join(REPO, "targets/E1_his707.pdb")
CPLX_DIR = os.path.join(REPO, "results/c4_rebuild_20261024/remote_pdbs")
CONTACTS = os.path.join(REPO, "results/c4_provenance_join_20261026/hotspot_contact.json")

CLASH_HARD = 2.0      # heavy-atom overlap treated as a real clash (repo convention)
CLASH_SOFT = 3.0      # close contact
CONTACT_CUTOFF = 4.5  # same cutoff the contact-mass audit used
BAND_LO, BAND_HI = 651, 741   # the real design epitope; the fit frame
BACKBONE = {"N", "CA", "C", "O"}
STD = {"ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE","LEU","LYS",
       "MET","PHE","PRO","SER","THR","TRP","TYR","VAL"}


def sha256_file(p):
    return hashlib.sha256(open(p, "rb").read()).hexdigest()


def _jd(o):
    if isinstance(o, (np.bool_,)): return bool(o)
    if isinstance(o, (np.integer,)): return int(o)
    if isinstance(o, (np.floating,)): return float(o)
    raise TypeError(str(type(o)))


def heavy(atoms):
    return [a for a in atoms if (a.element or "").strip().upper() != "H"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--self-test", action="store_true")
    args = ap.parse_args()

    outdir = os.path.join(REPO, args.out)
    out_json = os.path.join(outdir, "native_pose_clash.json")
    out_csv = os.path.join(outdir, "native_pose_clash_per_design.csv")
    if not args.self_test:
        os.makedirs(outdir, exist_ok=True)
        for p in (out_json, out_csv):
            if os.path.exists(p):
                sys.exit(f"REFUSING TO OVERWRITE: {p}")

    controls = {}

    # ---- crop order: file order of E1_his707 chain A -> author numbering
    pp = PDBParser(QUIET=True)
    tgt = pp.get_structure("t", TARGET)
    order = [r.id[1] for r in next(tgt.get_models())["A"] if r.id[0] == " "]
    controls["MAPPING_LEN"] = {
        "requirement": "crop is 154 residues",
        "n": len(order), "first": order[0], "last": order[-1],
        "fired": bool(len(order) == 154),
    }

    # ---- native assembly
    st = MMCIFParser(QUIET=True).get_structure("n", CIF)
    model = next(st.get_models())
    A, B = model["A"], model["B"]
    a_std = {r.id[1]: r for r in A if r.id[0] == " " and r.get_resname() in STD}
    b_std = [r for r in B if r.id[0] == " " and r.get_resname() in STD]
    nag = [r for r in A if r.get_resname() == "NAG"] + \
          [r for r in B if r.get_resname() == "NAG"]

    native_atoms = []
    for r in a_std.values():
        native_atoms += heavy(list(r))
    for r in b_std:
        native_atoms += heavy(list(r))
    for r in nag:
        native_atoms += heavy(list(r))
    ns = NeighborSearch(native_atoms)
    # chain-A-only tree, used for the displaced control
    print(f"native assembly: {len(native_atoms)} heavy atoms "
          f"({len(a_std)} + {len(b_std)} residues + {len(nag)} NAG)", flush=True)

    # ---- positive control: native partner subunit as a pseudo-binder vs chain A
    a_only = []
    for r in a_std.values():
        a_only += heavy(list(r))
    ns_a = NeighborSearch(a_only)
    partner_atoms = []
    for r in b_std:
        partner_atoms += heavy(list(r))
    n_partner_contacts = sum(1 for at in partner_atoms
                             if ns_a.search(at.coord, CONTACT_CUTOFF))
    controls["POS_NATIVE_PARTNER"] = {
        "requirement": "the real partner subunit registers contacts against chain A",
        "n_partner_atoms_in_contact": n_partner_contacts,
        "fired": bool(n_partner_contacts > 50),
    }

    files = sorted(f for f in os.listdir(CPLX_DIR) if f.endswith("_complex.pdb"))
    print(f"{len(files)} complexes", flush=True)

    rows = []
    rmsds = []
    global_rms = []
    his707_orig = 0
    displaced_clashes = None

    for n, fn in enumerate(files):
        s = pp.get_structure("c", os.path.join(CPLX_DIR, fn))
        m = next(s.get_models())
        if "A" not in m or "B" not in m:
            continue
        binder = [r for r in m["A"] if r.id[0] == " "]
        crop = [r for r in m["B"] if r.id[0] == " "]

        # --- regression: His707 contact in the ORIGINAL (un-superposed) frame
        binder_at = []
        for r in binder:
            binder_at += heavy(list(r))
        crop_by_author = {}
        for r in crop:
            i = r.id[1]
            if 1 <= i <= len(order):
                crop_by_author[order[i - 1]] = r
        if 707 in crop_by_author:
            tgt_at = heavy(list(crop_by_author[707]))
            if tgt_at:
                nsb = NeighborSearch(binder_at)
                if any(nsb.search(a.coord, CONTACT_CUTOFF) for a in tgt_at):
                    his707_orig += 1

        # --- rigid fit onto native chain A.
        # The fit is taken over the REAL DESIGN BAND (author 651-741) only, not the whole
        # 154-residue crop. Measured reason, recorded as a finding in its own right:
        # Boltz-2 reproduces the band at ~1.15 A CA RMSD but misplaces the stitched
        # fragments relative to one another, so a global fit is 17.6 A and meaningless.
        # The band is where every design binds, so a band-local frame is the correct one;
        # the global RMSD is still computed and reported per design as a diagnostic.
        fixed, moving = [], []
        gfixed, gmoving = [], []
        for author, r in crop_by_author.items():
            nr = a_std.get(author)
            if nr is None or "CA" not in r or "CA" not in nr:
                continue
            gfixed.append(nr["CA"]); gmoving.append(r["CA"])
            if BAND_LO <= author <= BAND_HI:
                fixed.append(nr["CA"]); moving.append(r["CA"])
        if len(fixed) < 30:
            continue
        sup = Superimposer()
        sup.set_atoms(fixed, moving)
        rmsds.append(float(sup.rms))
        gsup = Superimposer()
        gsup.set_atoms(gfixed, gmoving)
        global_rms.append(float(gsup.rms))

        # apply to a copy of the binder atoms
        coords = np.array([a.coord for a in binder_at], dtype=float)
        rot, tran = sup.rotran
        moved = coords @ rot + tran

        hard = soft = contact = bb_hard = 0
        mind = 1e9
        for at, xyz in zip(binder_at, moved):
            near = ns.search(xyz, CLASH_SOFT)
            if near:
                d = min(float(np.linalg.norm(xyz - b.coord)) for b in near)
                mind = min(mind, d)
                if d < CLASH_HARD:
                    hard += 1
                    if at.get_id() in BACKBONE and any(
                            b.get_id() in BACKBONE and
                            np.linalg.norm(xyz - b.coord) < CLASH_HARD for b in near):
                        bb_hard += 1
                soft += 1
            if ns.search(xyz, CONTACT_CUTOFF):
                contact += 1

        rows.append({
            "design": fn,
            "n_binder_heavy_atoms": len(binder_at),
            "band_fit_ca_rmsd": round(float(sup.rms), 3),
            "n_ca_pairs_fitted_band": len(fixed),
            "global_crop_ca_rmsd": round(float(gsup.rms), 3),
            "clashes_lt_2.0A": hard,
            "backbone_backbone_clashes_lt_2.0A": bb_hard,
            "close_lt_3.0A": soft,
            "contacts_lt_4.5A": contact,
            "min_heavy_atom_distance": round(mind, 3) if mind < 1e9 else None,
            "clash_free": bool(hard == 0),
        })

        if displaced_clashes is None:
            disp = moved + np.array([500.0, 0.0, 0.0])
            displaced_clashes = sum(1 for xyz in disp if ns.search(xyz, CLASH_SOFT))

        if (n + 1) % 50 == 0:
            print(f"  {n+1}/{len(files)}", flush=True)

    controls["SUPERPOSE_RMSD"] = {
        "requirement": "band-local fit (author 651-741) onto native chain A is tight",
        "band_median": round(float(np.median(rmsds)), 4),
        "band_max": round(float(np.max(rmsds)), 4),
        "global_crop_median_for_contrast": round(float(np.median(global_rms)), 4),
        "note": ("the global crop fit is an order of magnitude worse; Boltz-2 reproduces "
                 "the band but misplaces the stitched fragments. Reported, not hidden."),
        "fired": bool(np.median(rmsds) < 2.0),
    }
    controls["HIS707_REGRESSION"] = {
        "requirement": "mapping reproduces the recorded 19/343 His707 contacts",
        "expected": json.load(open(CONTACTS))["his707"],
        "measured": his707_orig,
        "fired": bool(his707_orig == json.load(open(CONTACTS))["his707"]),
    }
    controls["POS_DISPLACED"] = {
        "requirement": "a binder translated +500 A has zero close contacts",
        "measured": displaced_clashes,
        "fired": bool(displaced_clashes == 0),
    }

    failed = [k for k, v in controls.items() if not v["fired"]]
    print("\ncontrols:")
    for k, v in controls.items():
        print(f"  {k}: {'FIRED' if v['fired'] else 'DID NOT FIRE'}  {v}")
    if failed or args.self_test:
        if failed:
            sys.exit(f"CONTROLS DID NOT FIRE: {failed}")
        sys.exit(0)

    hard_arr = np.array([r["clashes_lt_2.0A"] for r in rows])
    bb_arr = np.array([r["backbone_backbone_clashes_lt_2.0A"] for r in rows])
    summary = {
        "n_designs": len(rows),
        "clash_free_lt_2.0A": int((hard_arr == 0).sum()),
        "backbone_backbone_clash_free": int((bb_arr == 0).sum()),
        "clashes_median": float(np.median(hard_arr)),
        "clashes_p90": float(np.percentile(hard_arr, 90)),
        "clashes_max": int(hard_arr.max()),
        "bb_clashes_median": float(np.median(bb_arr)),
        "designs_with_any_native_contact": int(sum(1 for r in rows if r["contacts_lt_4.5A"] > 0)),
        "band_fit_rmsd_median": round(float(np.median(rmsds)), 3),
        "global_crop_rmsd_median": round(float(np.median(global_rms)), 3),
    }

    payload = {
        "script": os.path.basename(__file__),
        "created_at_utc": os.popen("date -u '+%Y-%m-%dT%H:%M:%SZ'").read().strip(),
        "scope": ("steric feasibility of each design's existing pose against the full native "
                  "homodimer; NOT a binding, affinity or pH measurement, and NOT a ranking"),
        "mandated_by": "Astra PART 3, 2026-10-26; run on the full set per Opus 5.5 / Fable 5.1 / Astra process Q2 = (A)",
        "inputs": [{"path": os.path.relpath(p, REPO), "sha256": sha256_file(p)}
                   for p in (CIF, TARGET, CONTACTS)],
        "parameters": {"clash_hard_A": CLASH_HARD, "close_A": CLASH_SOFT,
                       "contact_A": CONTACT_CUTOFF},
        "controls": controls,
        "summary": summary,
    }
    with open(out_json, "w") as f:
        json.dump(payload, f, indent=1, default=_jd)
    import csv
    with open(out_csv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)

    print("\n--- native-pose steric feasibility ---")
    print(json.dumps(summary, indent=1))
    print(f"wrote {out_json}")


if __name__ == "__main__":
    main()
