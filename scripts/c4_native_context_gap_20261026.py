#!/usr/bin/env python3
"""
C4 decision evidence: what the official 640-aa tag-free homodimer changes
versus the 195-aa apical surrogate everything on disk was built against.

Everything in this repo models human TfR1 as a contiguous 195-residue apical
domain (P02786 author 189-383). The challenge target is P02786-1 residues
121-760, 640 aa per subunit, a tag-free HOMODIMER. This script measures, rather
than asserts, what the surrogate cannot see.

Method: all four contexts are built from the SAME experimental coordinates
(11QC chain A), so the only variable is context, not conformation:

  CTX1  apical only        residues 189-383 of chain A, alone        <- the surrogate
  CTX2  full monomer       residues 122-760 of chain A               <- adds intra-subunit domains
  CTX3  homodimer          chain A + chain B                         <- adds the partner subunit
  CTX4  homodimer+glycan   chain A + chain B + all 6 NAG             <- adds observed glycans

For every apical position: solvent accessibility in each context, distance to the
nearest partner-subunit atom, distance to the nearest glycan atom, and distance to
the membrane-proximal end of the construct (residue 122).

Then, for each of the 342 designs, the fraction of its measured human contact
footprint that is still accessible in the native context.

CONTROLS (must fire):
  BURIAL_MONOTONIC   accessibility must never increase as context is added
  CTX1_MATCHES_SEQ   CTX1 must contain exactly 195 residues, 189-383
  GLYCAN_PRESENT     6 NAG residues must be found
  SELF_CONSISTENT    CTX4 recomputed twice must be identical

CPU only. Read-only. No ranking. No GPU.
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
from Bio.PDB.SASA import ShrakeRupley
from Bio.PDB.Structure import Structure
from Bio.PDB.Model import Model
from Bio.PDB.Chain import Chain

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CIF = os.path.join(REPO, "data/11qc.cif")
STRAT = os.path.join(REPO, "results/c4_epitope_strategy_20261026/epitope_strategy.json")

APICAL_LO, APICAL_HI = 189, 383
MEMBRANE_PROXIMAL_RES = 122  # N-terminal end of the crystallised ectodomain
EXPOSED_REL = 0.20           # relative SASA above which a position counts as exposed
OCCLUDED_DROP = 0.50         # >=50% relative loss counts as materially occluded

MAX_ASA = {
    "ALA": 129.0, "ARG": 274.0, "ASN": 195.0, "ASP": 193.0, "CYS": 167.0,
    "GLN": 225.0, "GLU": 223.0, "GLY": 104.0, "HIS": 224.0, "ILE": 197.0,
    "LEU": 201.0, "LYS": 236.0, "MET": 224.0, "PHE": 240.0, "PRO": 159.0,
    "SER": 155.0, "THR": 172.0, "TRP": 285.0, "TYR": 263.0, "VAL": 174.0,
}


def sha256_file(p):
    return hashlib.sha256(open(p, "rb").read()).hexdigest()


def build(groups):
    """Assemble a detached Structure from lists of residues; returns (struct, index map)."""
    st = Structure("ctx")
    md = Model(0)
    st.add(md)
    for i, res in enumerate(groups):
        ch = Chain(chr(ord("A") + i))
        md.add(ch)
        for r in res:
            ch.add(r.copy())
    return st, md


def sasa_of(groups):
    """Shrake-Rupley over the assembled body; returns {(grp, resid): total_sasa}."""
    st, md = build(groups)
    ShrakeRupley(n_points=100).compute(st, level="A")
    out = {}
    for i, ch in enumerate(md):
        for r in ch:
            tot = sum(a.sasa for a in r if a.element != "H")
            out[(i, r.id)] = tot
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    outdir = os.path.join(REPO, args.out)
    os.makedirs(outdir, exist_ok=True)
    out_json = os.path.join(outdir, "native_context_gap.json")
    out_csv = os.path.join(outdir, "native_context_per_position.csv")
    for p in (out_json, out_csv):
        if os.path.exists(p):
            sys.exit(f"REFUSING TO OVERWRITE: {p}")

    st = MMCIFParser(QUIET=True).get_structure("11qc", CIF)
    model = next(st.get_models())
    A, B = model["A"], model["B"]

    a_std = [r for r in A if r.id[0] == " " and r.get_resname() in MAX_ASA]
    b_std = [r for r in B if r.id[0] == " " and r.get_resname() in MAX_ASA]
    nag_a = [r for r in A if r.get_resname() == "NAG"]
    nag_b = [r for r in B if r.get_resname() == "NAG"]
    nag = nag_a + nag_b  # NAGs reuse residue ids across chains; keep them in separate groups
    apical = [r for r in a_std if APICAL_LO <= r.id[1] <= APICAL_HI]

    controls = {
        "CTX1_MATCHES_SEQ": {
            "n_apical": len(apical),
            "lo": min(r.id[1] for r in apical),
            "hi": max(r.id[1] for r in apical),
            "fired": len(apical) == 195 and min(r.id[1] for r in apical) == APICAL_LO
                     and max(r.id[1] for r in apical) == APICAL_HI,
        },
        "GLYCAN_PRESENT": {"n_nag": len(nag), "fired": len(nag) == 6},
    }
    if not (controls["CTX1_MATCHES_SEQ"]["fired"] and controls["GLYCAN_PRESENT"]["fired"]):
        print(json.dumps(controls, indent=2))
        sys.exit("INPUT CONTROLS DID NOT FIRE")

    print("computing CTX1 apical only ...", flush=True)
    s1 = sasa_of([apical])
    print("computing CTX2 full monomer ...", flush=True)
    s2 = sasa_of([a_std])
    print("computing CTX3 homodimer ...", flush=True)
    s3 = sasa_of([a_std, b_std])
    print("computing CTX4 homodimer + glycans ...", flush=True)
    s4 = sasa_of([a_std, b_std, nag_a, nag_b])
    print("recomputing CTX4 for the self-control ...", flush=True)
    s4b = sasa_of([a_std, b_std, nag_a, nag_b])

    # geometric distances, measured once
    b_coords = np.array([a.get_coord() for r in b_std for a in r if a.element != "H"])
    nag_coords = np.array([a.get_coord() for r in nag for a in r])
    mem = [r for r in a_std if r.id[1] == MEMBRANE_PROXIMAL_RES]
    mem_coords = np.array([a.get_coord() for a in mem[0] if a.element != "H"]) if mem else None

    per_pos = []
    monotonic_violations = []
    self_violations = 0
    for r in apical:
        rid = r.id
        mx = MAX_ASA[r.get_resname()]
        v1 = s1[(0, rid)]
        v2 = s2[(0, rid)]
        v3 = s3[(0, rid)]
        v4 = s4[(0, rid)]
        if abs(v4 - s4b[(0, rid)]) > 1e-9:
            self_violations += 1
        # context may only ever remove accessibility
        if not (v1 + 1e-6 >= v2 >= v3 - 1e-6 and v3 + 1e-6 >= v4):
            monotonic_violations.append(
                {"res": rid[1], "ctx1": round(v1, 2), "ctx2": round(v2, 2),
                 "ctx3": round(v3, 2), "ctx4": round(v4, 2)})

        heavy = np.array([a.get_coord() for a in r if a.element != "H"])
        d_partner = float(np.sqrt(((heavy[:, None, :] - b_coords[None, :, :]) ** 2).sum(-1)).min())
        d_glycan = float(np.sqrt(((heavy[:, None, :] - nag_coords[None, :, :]) ** 2).sum(-1)).min())
        d_mem = (float(np.sqrt(((heavy[:, None, :] - mem_coords[None, :, :]) ** 2).sum(-1)).min())
                 if mem_coords is not None else None)

        rel1, rel4 = v1 / mx, v4 / mx
        per_pos.append({
            "human_auth": rid[1],
            "resname": r.get_resname(),
            "sasa_ctx1_apical_only": round(v1, 2),
            "sasa_ctx2_full_monomer": round(v2, 2),
            "sasa_ctx3_homodimer": round(v3, 2),
            "sasa_ctx4_homodimer_glycan": round(v4, 2),
            "rel_ctx1": round(rel1, 4),
            "rel_ctx4": round(rel4, 4),
            "abs_loss_A2": round(v1 - v4, 2),
            "frac_loss": round((v1 - v4) / v1, 4) if v1 > 1.0 else None,
            "loss_to_intra_monomer_A2": round(v1 - v2, 2),
            "loss_to_partner_subunit_A2": round(v2 - v3, 2),
            "loss_to_glycan_A2": round(v3 - v4, 2),
            "exposed_in_surrogate": bool(rel1 >= EXPOSED_REL),
            "exposed_in_native": bool(rel4 >= EXPOSED_REL),
            "materially_occluded": bool(v1 > 1.0 and (v1 - v4) / v1 >= OCCLUDED_DROP),
            "dist_to_partner_subunit_A": round(d_partner, 2),
            "dist_to_nearest_glycan_A": round(d_glycan, 2),
            "dist_to_membrane_proximal_122_A": round(d_mem, 2) if d_mem is not None else None,
        })

    controls["BURIAL_MONOTONIC"] = {
        "n_violations": len(monotonic_violations),
        "examples": monotonic_violations[:5],
        "fired": len(monotonic_violations) == 0,
    }
    controls["SELF_CONSISTENT"] = {"n_violations": self_violations, "fired": self_violations == 0}
    controls["ALL_FIRED"] = all(v["fired"] for k, v in controls.items() if isinstance(v, dict))
    if not controls["ALL_FIRED"]:
        print(json.dumps(controls, indent=2))
        sys.exit("CONTROLS DID NOT ALL FIRE - refusing to emit measurements")
    print("all controls fired", flush=True)

    byres = {p["human_auth"]: p for p in per_pos}

    # ---- summary over the apical surface ----
    exposed_surr = [p for p in per_pos if p["exposed_in_surrogate"]]
    lost = [p for p in exposed_surr if not p["exposed_in_native"]]
    occl = [p for p in per_pos if p["materially_occluded"]]
    near_partner = [p for p in per_pos if p["dist_to_partner_subunit_A"] <= 8.0]
    near_glycan = [p for p in per_pos if p["dist_to_nearest_glycan_A"] <= 8.0]

    summary = {
        "apical_positions": len(per_pos),
        "total_sasa_ctx1_A2": round(sum(p["sasa_ctx1_apical_only"] for p in per_pos), 1),
        "total_sasa_ctx2_A2": round(sum(p["sasa_ctx2_full_monomer"] for p in per_pos), 1),
        "total_sasa_ctx3_A2": round(sum(p["sasa_ctx3_homodimer"] for p in per_pos), 1),
        "total_sasa_ctx4_A2": round(sum(p["sasa_ctx4_homodimer_glycan"] for p in per_pos), 1),
        "apical_area_lost_to_intra_monomer_A2": round(sum(p["loss_to_intra_monomer_A2"] for p in per_pos), 1),
        "apical_area_lost_to_partner_subunit_A2": round(sum(p["loss_to_partner_subunit_A2"] for p in per_pos), 1),
        "apical_area_lost_to_glycan_A2": round(sum(p["loss_to_glycan_A2"] for p in per_pos), 1),
        "n_exposed_in_surrogate": len(exposed_surr),
        "n_exposed_in_native": sum(1 for p in per_pos if p["exposed_in_native"]),
        "n_exposed_in_surrogate_but_not_native": len(lost),
        "positions_exposed_in_surrogate_but_not_native": [p["human_auth"] for p in lost],
        "n_materially_occluded_ge50pct": len(occl),
        "positions_materially_occluded": [p["human_auth"] for p in occl],
        "n_within_8A_of_partner_subunit": len(near_partner),
        "positions_within_8A_of_partner_subunit": [p["human_auth"] for p in near_partner],
        "n_within_8A_of_glycan": len(near_glycan),
        "positions_within_8A_of_glycan": [p["human_auth"] for p in near_glycan],
        "min_dist_apical_to_partner_A": round(min(p["dist_to_partner_subunit_A"] for p in per_pos), 2),
        "min_dist_apical_to_glycan_A": round(min(p["dist_to_nearest_glycan_A"] for p in per_pos), 2),
        "min_dist_apical_to_membrane_proximal_A": round(
            min(p["dist_to_membrane_proximal_122_A"] for p in per_pos), 2),
    }

    # ---- per-design footprint exposure in the native context ----
    designs = []
    if os.path.exists(STRAT):
        sd = json.load(open(STRAT))
        for r in sd["designs"]:
            fp = [int(x) for x in r["footprint_human_auth"].split(";") if x]
            fp = [a for a in fp if a in byres]
            if not fp:
                continue
            n_occ = sum(1 for a in fp if byres[a]["materially_occluded"])
            n_lost = sum(1 for a in fp if byres[a]["exposed_in_surrogate"]
                         and not byres[a]["exposed_in_native"])
            n_near_p = sum(1 for a in fp if byres[a]["dist_to_partner_subunit_A"] <= 8.0)
            n_near_g = sum(1 for a in fp if byres[a]["dist_to_nearest_glycan_A"] <= 8.0)
            a1 = sum(byres[a]["sasa_ctx1_apical_only"] for a in fp)
            a4 = sum(byres[a]["sasa_ctx4_homodimer_glycan"] for a in fp)
            designs.append({
                "name": r["name"],
                "sequence_sha256": r["sequence_sha256"],
                "in_v4": r["in_v4"],
                "L6_apical": r["L6_apical"],
                "iptm_apical": r["iptm_apical"],
                "n_footprint": len(fp),
                "n_footprint_materially_occluded": n_occ,
                "frac_footprint_occluded": round(n_occ / len(fp), 4),
                "n_footprint_exposed_surrogate_not_native": n_lost,
                "n_footprint_within_8A_partner": n_near_p,
                "n_footprint_within_8A_glycan": n_near_g,
                "footprint_sasa_surrogate_A2": round(a1, 1),
                "footprint_sasa_native_A2": round(a4, 1),
                "frac_footprint_area_retained_native": round(a4 / a1, 4) if a1 > 1 else None,
            })

    fa = np.array([d["frac_footprint_area_retained_native"] for d in designs
                   if d["frac_footprint_area_retained_native"] is not None])
    fo = np.array([d["frac_footprint_occluded"] for d in designs])
    design_summary = {
        "n_designs": len(designs),
        "frac_footprint_area_retained_native": {
            "median": round(float(np.median(fa)), 4),
            "p10": round(float(np.percentile(fa, 10)), 4),
            "p90": round(float(np.percentile(fa, 90)), 4),
            "min": round(float(fa.min()), 4),
            "max": round(float(fa.max()), 4),
        },
        "frac_footprint_materially_occluded": {
            "median": round(float(np.median(fo)), 4),
            "max": round(float(fo.max()), 4),
            "n_designs_with_any_occluded_contact": int((fo > 0).sum()),
            "n_designs_with_ge25pct_occluded": int((fo >= 0.25).sum()),
        },
    }

    payload = {
        "created_at_utc": os.popen("date -u +%Y-%m-%dT%H:%M:%SZ").read().strip(),
        "purpose": "Decision evidence: 640-aa tag-free homodimer versus 195-aa apical surrogate.",
        "target_note": "11QC is P02786 122-760 per subunit (639 aa); the official target is 121-760 "
                       "(640 aa). Residue 121 is absent locally and is not measured here.",
        "method": {
            "contexts": {
                "CTX1": "apical 189-383 of 11QC chain A, alone (the surrogate)",
                "CTX2": "full 11QC chain A monomer, 122-760",
                "CTX3": "11QC chain A + chain B homodimer",
                "CTX4": "11QC chain A + chain B + all 6 observed NAG",
            },
            "note": "All contexts use the same experimental coordinates, so the only variable is "
                    "context, not conformation. This deliberately does NOT model the predicted "
                    "binder poses in the native assembly.",
            "sasa": "Bio.PDB.SASA.ShrakeRupley n_points=100, heavy atoms only",
            "exposed_rel_sasa": EXPOSED_REL,
            "materially_occluded_frac_loss": OCCLUDED_DROP,
        },
        "inputs": {"cif": "data/11qc.cif", "cif_sha256": sha256_file(CIF),
                   "strategy_json_sha256": sha256_file(STRAT) if os.path.exists(STRAT) else None},
        "controls": controls,
        "summary": summary,
        "design_summary": design_summary,
        "per_position": per_pos,
        "designs": designs,
    }
    json.dump(payload, open(out_json, "w"), indent=1)

    import csv as _csv
    with open(out_csv, "w", newline="") as fh:
        w = _csv.DictWriter(fh, fieldnames=list(per_pos[0].keys()))
        w.writeheader()
        for p in per_pos:
            w.writerow(p)

    print(f"\nwrote {out_json}\nwrote {out_csv}")
    print(json.dumps(summary, indent=1)[:2600])
    print(json.dumps(design_summary, indent=1))


if __name__ == "__main__":
    main()
