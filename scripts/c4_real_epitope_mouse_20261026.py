#!/usr/bin/env python3
"""
C4: mouse conservation AND native accessibility of the REAL design epitope.

`results/c4_native_context_20261026/FINDING_WRONG_DOMAIN.md` established that the 342
`tfr1_c4_acid_site1/2_*` designs were generated against human TfR1 author residues
~651-741 -- the helical / dimerisation domain, as a TWO-SUBUNIT groove around His707 --
and that every evaluation from v3 onward scored them against apical 189-383. That finding
closes with:

    "Mouse conservation of the real epitope, author 651-741 across both subunits,
     has not been measured by anyone."

This measures it, and the companion question the same finding leaves open: whether that
composite groove is solvent-accessible in the native homodimer at all, or whether it was
an artifact of the 154-residue crop the designs were generated against.

Both answers bear directly on the open construct-scope decision request
(`results/decision_requests/20261026-c4-construct-scope.md`). C4 requires mouse
cross-reactivity. If the real epitope is poorly conserved in mouse, no amount of
re-evaluation against the correct construct can make this design pool satisfy that
objective, and options B/C buy nothing. If it is well conserved, the pool has upside that
only a correct-construct measurement can reach.

Method
  1. Mouse conservation: full-length global BLOSUM62 alignment of human P02786 against
     mouse Q62351 (same parameters as results/apical_context/: open -10, extend -0.5),
     retaining EVERY equally-optimal alignment. A human position counts as conserved only
     if it is identical in ALL optima; "identical in at least one optimum" is reported
     alongside, so the ambiguity is visible rather than hidden by a tie-break.
  2. Native accessibility: the same four nested contexts as the apical measurement, built
     from identical 11QC coordinates so context is the only variable --
       CTX1 band alone / CTX2 full 639-aa monomer / CTX3 homodimer / CTX4 + 6 NAG.
  3. Per design: each design's own contact footprint (already mapped to author numbering
     in results/c4_provenance_join_20261026/hotspot_contact.json), scored for mouse
     identity over the whole footprint and over its natively-exposed part only.

CONTROLS (every one must fire or this exits non-zero)
  ALIGN_SELF          P02786 aligned to itself is 100% identical at every position
  NUMBERING_11QC      11QC chain A author numbering == P02786 sequence numbering
  HOTSPOT_IDENTITY    the decisive 7/7 config test: 707 HIS, 659 ASP, 682/694/714/728 GLU
  APICAL_REGRESSION   apical 189-383 reproduces the independently measured 146/195
  NEG_WINDOW          a deliberately shifted window must NOT give the same identity
  BURIAL_MONOTONIC    accessibility never increases as context is added
  GLYCAN_PRESENT      6 NAG residues found

CPU only. Read-only on every input. No ranking. No GPU. No submission.
"""

import argparse
import hashlib
import json
import os
import sys
import warnings
from collections import Counter

import numpy as np

warnings.simplefilter("ignore")
from Bio import Align
from Bio.Align import substitution_matrices
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.SASA import ShrakeRupley
from Bio.PDB.Structure import Structure
from Bio.PDB.Model import Model
from Bio.PDB.Chain import Chain

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CIF = os.path.join(REPO, "data/11qc.cif")
HUMAN_TXT = os.path.join(REPO, "data/P02786.txt")
MOUSE_TXT = os.path.join(REPO, "data/Q62351.txt")
CONTACTS = os.path.join(REPO, "results/c4_provenance_join_20261026/hotspot_contact.json")

# The real design epitope, author numbering, from the two RFD4 contigs translated out of
# label numbering (author = label + 121) in FINDING_WRONG_DOMAIN.md:
#   site1  A555-620 -> author 676-741 ; site2  A530-600 -> author 651-721
#   both   B555-575 -> author 676-696  (partner subunit)
REAL_LO, REAL_HI = 651, 741
PARTNER_LO, PARTNER_HI = 676, 696
APICAL_LO, APICAL_HI = 189, 383
MONOMER_LO, MONOMER_HI = 122, 760

# Named hotspots from the generator configs.
CARBOXYLATE_HOTSPOTS = {659: "ASP", 682: "GLU", 694: "GLU", 714: "GLU", 728: "GLU"}
GROOVE_HOTSPOTS = [685, 701, 702, 707, 709]
ANCHOR_HIS = 707

EXPOSED_REL = 0.20    # relative SASA above which a position counts as exposed
OCCLUDED_DROP = 0.50  # >=50% relative loss counts as materially occluded

THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
    "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V",
}
MAX_ASA = {
    "ALA": 129.0, "ARG": 274.0, "ASN": 195.0, "ASP": 193.0, "CYS": 167.0,
    "GLN": 225.0, "GLU": 223.0, "GLY": 104.0, "HIS": 224.0, "ILE": 197.0,
    "LEU": 201.0, "LYS": 236.0, "MET": 224.0, "PHE": 240.0, "PRO": 159.0,
    "SER": 155.0, "THR": 172.0, "TRP": 285.0, "TYR": 263.0, "VAL": 174.0,
}
MAX_ALIGNMENTS = 64  # enumeration cap; recorded so a tie explosion cannot pass silently


def sha256_file(p):
    return hashlib.sha256(open(p, "rb").read()).hexdigest()


def _json_default(o):
    if isinstance(o, (np.bool_,)):
        return bool(o)
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating,)):
        return float(o)
    raise TypeError(f"not JSON serialisable: {type(o)}")


def read_uniprot_seq(path):
    """Pull the SQ block out of a UniProt flat file."""
    seq, inside = [], False
    for line in open(path):
        if line.startswith("SQ "):
            inside = True
            continue
        if inside:
            if line.startswith("//"):
                break
            seq.append(line.strip().replace(" ", ""))
    return "".join(seq)


def aligner():
    al = Align.PairwiseAligner()
    al.substitution_matrix = substitution_matrices.load("BLOSUM62")
    al.open_gap_score = -10
    al.extend_gap_score = -0.5
    al.mode = "global"
    return al


def map_human_to_mouse(human, mouse):
    """
    Return (per_position, meta).

    per_position[human_pos_1based] = {
        "human_aa", "mouse_aa_set", "identical_all", "identical_any", "gapped_any"
    }
    computed over EVERY equally-optimal alignment, so a position whose mouse partner is
    ambiguous between optima is visible instead of silently tie-broken.
    """
    al = aligner()
    alns = al.align(human, mouse)
    best = alns.score
    n_opt = len(alns)
    enumerated = min(n_opt, MAX_ALIGNMENTS)

    partners = {i: set() for i in range(len(human))}
    for k in range(enumerated):
        a = alns[k]
        assert abs(a.score - best) < 1e-9, "non-optimal alignment enumerated"
        # a.aligned gives paired blocks of (start, end) index ranges
        hblocks, mblocks = a.aligned
        paired = set()
        for (hs, he), (ms, me) in zip(hblocks, mblocks):
            for off in range(he - hs):
                partners[hs + off].add(mouse[ms + off])
                paired.add(hs + off)
        for i in range(len(human)):
            if i not in paired:
                partners[i].add("-")

    per_position = {}
    for i, h in enumerate(human):
        s = partners[i]
        per_position[i + 1] = {
            "human_aa": h,
            "mouse_aa_set": sorted(s),
            "identical_all": s == {h},
            "identical_any": h in s,
            "gapped_any": "-" in s,
        }
    meta = {
        "best_score": float(best),
        "equally_best_alignments": int(n_opt),
        "enumerated": int(enumerated),
        "enumeration_capped": bool(n_opt > MAX_ALIGNMENTS),
        "open_gap_score": -10,
        "extend_gap_score": -0.5,
        "matrix": "BLOSUM62",
        "mode": "global",
    }
    return per_position, meta


def band_identity(per_position, lo, hi):
    """Conservation summary over an inclusive author-numbering band."""
    pos = [p for p in range(lo, hi + 1) if p in per_position]
    ident_all = [p for p in pos if per_position[p]["identical_all"]]
    ident_any = [p for p in pos if per_position[p]["identical_any"]]
    gapped = [p for p in pos if per_position[p]["gapped_any"]]
    return {
        "lo": lo, "hi": hi, "n": len(pos),
        "identical_all_optima": len(ident_all),
        "identical_any_optimum": len(ident_any),
        "gapped_in_some_optimum": len(gapped),
        "identity_fraction": (len(ident_all) / len(pos)) if pos else None,
        "substituted_positions": [p for p in pos if not per_position[p]["identical_any"]],
    }


def build(groups):
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
    st, md = build(groups)
    ShrakeRupley(n_points=100).compute(st, level="A")
    out = {}
    for i, ch in enumerate(md):
        for r in ch:
            out[(i, r.id)] = sum(a.sasa for a in r if a.element != "H")
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--self-test", action="store_true",
                    help="run controls only, write nothing")
    args = ap.parse_args()

    outdir = os.path.join(REPO, args.out)
    out_json = os.path.join(outdir, "real_epitope_mouse.json")
    out_csv = os.path.join(outdir, "real_epitope_per_position.csv")
    out_design_csv = os.path.join(outdir, "real_epitope_per_design.csv")
    if not args.self_test:
        os.makedirs(outdir, exist_ok=True)
        for p in (out_json, out_csv, out_design_csv):
            if os.path.exists(p):
                sys.exit(f"REFUSING TO OVERWRITE: {p}")

    controls = {}

    # ---------------------------------------------------------------- sequences
    human = read_uniprot_seq(HUMAN_TXT)
    mouse = read_uniprot_seq(MOUSE_TXT)
    print(f"human P02786 {len(human)} aa; mouse Q62351 {len(mouse)} aa", flush=True)

    # CONTROL: self-alignment must be perfect at every position.
    self_map, self_meta = map_human_to_mouse(human, human)
    n_self_ident = sum(1 for v in self_map.values() if v["identical_all"])
    controls["ALIGN_SELF"] = {
        "requirement": "P02786 vs P02786 is identical at all positions",
        "n_positions": len(human),
        "n_identical": n_self_ident,
        "fired": bool(n_self_ident == len(human)),
    }

    print("aligning human vs mouse (all optima) ...", flush=True)
    per_position, aln_meta = map_human_to_mouse(human, mouse)
    print(f"  {aln_meta['equally_best_alignments']} equally-best alignments, "
          f"score {aln_meta['best_score']}", flush=True)

    # ---------------------------------------------------------------- structure
    st = MMCIFParser(QUIET=True).get_structure("11qc", CIF)
    model = next(st.get_models())
    A, B = model["A"], model["B"]
    a_std = [r for r in A if r.id[0] == " " and r.get_resname() in MAX_ASA]
    b_std = [r for r in B if r.id[0] == " " and r.get_resname() in MAX_ASA]
    nag_a = [r for r in A if r.get_resname() == "NAG"]
    nag_b = [r for r in B if r.get_resname() == "NAG"]
    resname_a = {r.id[1]: r.get_resname() for r in a_std}

    controls["GLYCAN_PRESENT"] = {
        "requirement": "6 NAG residues in 11QC",
        "n_nag": len(nag_a) + len(nag_b),
        "fired": bool(len(nag_a) + len(nag_b) == 6),
    }

    # CONTROL: 11QC author numbering must BE UniProt numbering. Everything in this script
    # that maps a structural residue onto the alignment depends on it.
    mism = [(p, rn, human[p - 1]) for p, rn in sorted(resname_a.items())
            if 1 <= p <= len(human) and THREE_TO_ONE.get(rn) != human[p - 1]]
    controls["NUMBERING_11QC"] = {
        "requirement": "11QC chain A author residue identity == P02786 at the same number",
        "n_checked": len(resname_a),
        "n_mismatched": len(mism),
        "examples": mism[:5],
        "fired": bool(len(mism) == 0),
    }

    # CONTROL: the decisive config test from FINDING_WRONG_DOMAIN.md, re-run here.
    hs_checks = {}
    for pos, want in CARBOXYLATE_HOTSPOTS.items():
        hs_checks[pos] = {"want": want, "got": resname_a.get(pos),
                          "ok": resname_a.get(pos) == want}
    hs_checks[ANCHOR_HIS] = {"want": "HIS", "got": resname_a.get(ANCHOR_HIS),
                             "ok": resname_a.get(ANCHOR_HIS) == "HIS"}
    controls["HOTSPOT_IDENTITY"] = {
        "requirement": "707 HIS, 659 ASP, 682/694/714/728 GLU under author numbering",
        "checks": hs_checks,
        "n_ok": sum(1 for v in hs_checks.values() if v["ok"]),
        "n_total": len(hs_checks),
        "fired": bool(all(v["ok"] for v in hs_checks.values())),
    }

    # ---------------------------------------------------------------- conservation
    bands = {
        "real_epitope_651_741": band_identity(per_position, REAL_LO, REAL_HI),
        "partner_band_676_696": band_identity(per_position, PARTNER_LO, PARTNER_HI),
        "apical_189_383": band_identity(per_position, APICAL_LO, APICAL_HI),
        "official_target_121_760": band_identity(per_position, 121, 760),
    }

    # CONTROL: regression against the independently measured apical figure (146/195).
    controls["APICAL_REGRESSION"] = {
        "requirement": "apical 189-383 identical count == 146 (measured in a prior cycle)",
        "expected": 146,
        "measured": bands["apical_189_383"]["identical_all_optima"],
        "fired": bool(bands["apical_189_383"]["identical_all_optima"] == 146),
    }

    # CONTROL: a shifted window must not return the same identity, or the mapping is not
    # actually position-sensitive and every band number above would be meaningless.
    neg = band_identity(per_position, REAL_LO - 200, REAL_HI - 200)
    controls["NEG_WINDOW"] = {
        "requirement": "a 200-residue-shifted window differs from the real band",
        "real_identity": bands["real_epitope_651_741"]["identity_fraction"],
        "shifted_identity": neg["identity_fraction"],
        "fired": bool(neg["identity_fraction"] is not None
                      and abs(neg["identity_fraction"]
                              - bands["real_epitope_651_741"]["identity_fraction"]) > 1e-9),
    }

    # ---------------------------------------------------------------- accessibility
    band_res = [r for r in a_std if REAL_LO <= r.id[1] <= REAL_HI]
    mono = [r for r in a_std if MONOMER_LO <= r.id[1] <= MONOMER_HI]
    print(f"real-epitope band: {len(band_res)} residues resolved in 11QC chain A",
          flush=True)

    print("CTX1 band alone ...", flush=True)
    s1 = sasa_of([band_res])
    print("CTX2 full monomer ...", flush=True)
    s2 = sasa_of([mono])
    print("CTX3 homodimer ...", flush=True)
    s3 = sasa_of([a_std, b_std])
    print("CTX4 homodimer + glycans ...", flush=True)
    s4 = sasa_of([a_std, b_std, nag_a, nag_b])

    rows, viol = [], []
    for r in band_res:
        rid, rn = r.id, r.get_resname()
        pos = rid[1]
        mx = MAX_ASA[rn]
        v1, v2 = s1[(0, rid)], s2[(0, rid)]
        v3, v4 = s3[(0, rid)], s4[(0, rid)]
        # Monotonicity: adding atoms can only remove accessible surface. A small positive
        # tolerance absorbs Shrake-Rupley point-sampling noise.
        for a, b in ((v1, v2), (v2, v3), (v3, v4)):
            if b - a > 0.5:
                viol.append({"pos": pos, "from": round(a, 3), "to": round(b, 3)})
        pp = per_position.get(pos, {})
        rows.append({
            "author_pos": pos,
            "resname": rn,
            "human_aa": THREE_TO_ONE.get(rn, "X"),
            "mouse_aa_set": "|".join(pp.get("mouse_aa_set", [])),
            "identical_all_optima": bool(pp.get("identical_all", False)),
            "sasa_ctx1_band_alone": round(v1, 2),
            "sasa_ctx2_monomer": round(v2, 2),
            "sasa_ctx3_dimer": round(v3, 2),
            "sasa_ctx4_dimer_glycan": round(v4, 2),
            "rel_sasa_ctx1": round(v1 / mx, 4),
            "rel_sasa_ctx4": round(v4 / mx, 4),
            "exposed_ctx1": bool(v1 / mx >= EXPOSED_REL),
            "exposed_ctx4": bool(v4 / mx >= EXPOSED_REL),
            "frac_area_lost": round(1.0 - (v4 / v1), 4) if v1 > 1e-9 else None,
            "materially_occluded": bool(v1 > 1e-9 and (1.0 - v4 / v1) >= OCCLUDED_DROP),
            "is_named_hotspot": bool(pos in CARBOXYLATE_HOTSPOTS or pos in GROOVE_HOTSPOTS),
        })
    controls["BURIAL_MONOTONIC"] = {
        "requirement": "accessibility never increases as context is added",
        "n_violations": len(viol),
        "examples": viol[:5],
        "fired": bool(len(viol) == 0),
    }

    tot = {k: sum(r[k] for r in rows) for k in
           ("sasa_ctx1_band_alone", "sasa_ctx2_monomer",
            "sasa_ctx3_dimer", "sasa_ctx4_dimer_glycan")}
    accessibility = {
        "n_band_residues_resolved": len(rows),
        "band_sasa_A2": {
            "ctx1_band_alone": round(tot["sasa_ctx1_band_alone"], 1),
            "ctx2_monomer": round(tot["sasa_ctx2_monomer"], 1),
            "ctx3_dimer": round(tot["sasa_ctx3_dimer"], 1),
            "ctx4_dimer_glycan": round(tot["sasa_ctx4_dimer_glycan"], 1),
        },
        "frac_retained_ctx4_vs_ctx1": round(
            tot["sasa_ctx4_dimer_glycan"] / tot["sasa_ctx1_band_alone"], 4),
        "lost_to_own_subunit_A2": round(
            tot["sasa_ctx1_band_alone"] - tot["sasa_ctx2_monomer"], 1),
        "lost_to_partner_subunit_A2": round(
            tot["sasa_ctx2_monomer"] - tot["sasa_ctx3_dimer"], 1),
        "lost_to_glycan_A2": round(tot["sasa_ctx3_dimer"] - tot["sasa_ctx4_dimer_glycan"], 1),
        "n_exposed_ctx1": sum(1 for r in rows if r["exposed_ctx1"]),
        "n_exposed_ctx4": sum(1 for r in rows if r["exposed_ctx4"]),
        "n_exposed_ctx1_not_ctx4": sum(1 for r in rows
                                       if r["exposed_ctx1"] and not r["exposed_ctx4"]),
        "n_materially_occluded": sum(1 for r in rows if r["materially_occluded"]),
    }

    # Conservation restricted to what is actually presented in the native assembly --
    # the apical lesson was that burial inflates apparent conservation.
    exposed_pos = [r["author_pos"] for r in rows if r["exposed_ctx4"]]
    buried_pos = [r["author_pos"] for r in rows if not r["exposed_ctx4"]]

    def ident_frac(positions):
        if not positions:
            return None
        n = sum(1 for p in positions if per_position.get(p, {}).get("identical_all"))
        return {"n": len(positions), "n_identical": n, "identity_fraction": n / len(positions)}

    conservation_by_exposure = {
        "natively_exposed_ctx4": ident_frac(exposed_pos),
        "natively_buried_ctx4": ident_frac(buried_pos),
        "whole_band": ident_frac([r["author_pos"] for r in rows]),
    }

    hotspot_detail = []
    for pos in sorted(set(list(CARBOXYLATE_HOTSPOTS) + GROOVE_HOTSPOTS)):
        pp = per_position.get(pos, {})
        row = next((r for r in rows if r["author_pos"] == pos), None)
        hotspot_detail.append({
            "author_pos": pos,
            "human": pp.get("human_aa"),
            "mouse_aa_set": pp.get("mouse_aa_set"),
            "conserved": bool(pp.get("identical_all", False)),
            "rel_sasa_ctx4": row["rel_sasa_ctx4"] if row else None,
            "exposed_in_native": row["exposed_ctx4"] if row else None,
            "role": ("carboxylate_hotspot" if pos in CARBOXYLATE_HOTSPOTS
                     else "groove_hotspot"),
        })

    # ---------------------------------------------------------------- per design
    cj = json.load(open(CONTACTS))
    exposed_set = set(exposed_pos)
    design_rows = []
    for rec in cj["records"]:
        contacted = rec["contacted_author_residues"]
        in_band = [p for p in contacted if REAL_LO <= p <= REAL_HI]
        in_band_exposed = [p for p in in_band if p in exposed_set]
        cons_band = [p for p in in_band
                     if per_position.get(p, {}).get("identical_all")]
        cons_exposed = [p for p in in_band_exposed
                        if per_position.get(p, {}).get("identical_all")]
        design_rows.append({
            "design": rec["design"],
            "n_contacted_total": rec["n_contacted"],
            "n_in_real_band": len(in_band),
            "n_in_real_band_natively_exposed": len(in_band_exposed),
            "mouse_identity_real_band": round(len(cons_band) / len(in_band), 4) if in_band else None,
            "mouse_identity_real_band_exposed": (round(len(cons_exposed) / len(in_band_exposed), 4)
                                                 if in_band_exposed else None),
            "contacts_His707": rec["contacts_His707"],
            "contacts_any_named_hotspot": rec["contacts_any_named_hotspot"],
        })

    def summarise(key):
        vals = [r[key] for r in design_rows if r[key] is not None]
        if not vals:
            return None
        a = np.array(vals, dtype=float)
        return {
            "n_designs_with_value": len(vals),
            "median": round(float(np.median(a)), 4),
            "p10": round(float(np.percentile(a, 10)), 4),
            "p90": round(float(np.percentile(a, 90)), 4),
            "min": round(float(a.min()), 4),
            "max": round(float(a.max()), 4),
        }

    per_design_summary = {
        "n_designs": len(design_rows),
        "designs_with_any_real_band_contact": sum(1 for r in design_rows if r["n_in_real_band"] > 0),
        "designs_with_any_natively_exposed_real_band_contact":
            sum(1 for r in design_rows if r["n_in_real_band_natively_exposed"] > 0),
        "n_in_real_band": summarise("n_in_real_band"),
        "n_in_real_band_natively_exposed": summarise("n_in_real_band_natively_exposed"),
        "mouse_identity_real_band": summarise("mouse_identity_real_band"),
        "mouse_identity_real_band_exposed": summarise("mouse_identity_real_band_exposed"),
    }

    # ---------------------------------------------------------------- controls gate
    failed = [k for k, v in controls.items() if not v["fired"]]
    print("\ncontrols:")
    for k, v in controls.items():
        print(f"  {k}: {'FIRED' if v['fired'] else 'DID NOT FIRE'}")
    if args.self_test:
        print(json.dumps(controls, indent=2, default=_json_default))
        sys.exit(0 if not failed else f"CONTROLS DID NOT FIRE: {failed}")
    if failed:
        print(json.dumps(controls, indent=2, default=_json_default))
        sys.exit(f"CONTROLS DID NOT FIRE: {failed}")

    payload = {
        "script": os.path.basename(__file__),
        "created_at_utc": os.popen("date -u '+%Y-%m-%dT%H:%M:%SZ'").read().strip(),
        "scope": ("mouse conservation and native accessibility of the REAL design epitope "
                  "(author 651-741 plus the 676-696 partner band); not a binding, affinity "
                  "or pH measurement and not a ranking"),
        "inputs": [{"path": os.path.relpath(p, REPO), "sha256": sha256_file(p)}
                   for p in (HUMAN_TXT, MOUSE_TXT, CIF, CONTACTS)],
        "alignment": aln_meta,
        "controls": controls,
        "conservation_by_band": bands,
        "conservation_by_exposure": conservation_by_exposure,
        "native_accessibility": accessibility,
        "named_hotspots": hotspot_detail,
        "per_design_summary": per_design_summary,
    }
    with open(out_json, "w") as f:
        json.dump(payload, f, indent=1, default=_json_default)

    import csv
    with open(out_csv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    with open(out_design_csv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(design_rows[0].keys()))
        w.writeheader()
        w.writerows(design_rows)

    print("\n--- real epitope, author 651-741 ---")
    b = bands["real_epitope_651_741"]
    print(f"mouse identity: {b['identical_all_optima']}/{b['n']} = "
          f"{b['identity_fraction']:.4f}")
    print(f"apical 189-383 for contrast: "
          f"{bands['apical_189_383']['identity_fraction']:.4f}")
    print(f"native retention of band surface (CTX4/CTX1): "
          f"{accessibility['frac_retained_ctx4_vs_ctx1']:.4f}")
    print(f"exposed alone {accessibility['n_exposed_ctx1']} -> "
          f"native {accessibility['n_exposed_ctx4']}")
    ne = conservation_by_exposure["natively_exposed_ctx4"]
    if ne:
        print(f"identity among natively exposed: {ne['n_identical']}/{ne['n']} = "
              f"{ne['identity_fraction']:.4f}")
    print(f"\nwrote {out_json}")


if __name__ == "__main__":
    main()
