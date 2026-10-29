#!/usr/bin/env python3
"""
C4 independent check: foldability and interface chemistry.

PURPOSE
-------
The v4 ranking (results/c4_ranked_v4_apical_L6_20261025/) selected designs on
  (a) apical iPTM >= 0.60, (b) apical L6 > 0, ranked by apical L6 descending,
  (c) a 2.0 A heavy-atom clash gate,
and then reported those same quantities. This script computes signals that
entered NONE of that selection, so the v4 top rows can be tested out-of-sample:

  Family 1  binder monomer core quality   (is the binder a foldable protein at all)
  Family 2  interface chemistry           (does the interface bury apolar surface)
  Family 3  sequence liabilities          (composition / low-complexity, pose-free)

None of these is a variant of L6, of iPTM, or of the clash screen. Family 1 and 3
do not look at the interface at all; Family 2 measures a different physical
property (buried apolar area) than an interface confidence score.

Net charge is reported as a COMPOSITION DESCRIPTOR only. It is deliberately NOT
framed as a pH-mechanism signal and is not combined into any score, so that this
check stays clean of the pH axis that L6 already occupies.

CONTROLS (all must fire; see --self-test)
-----------------------------------------
  POS_HETERO   1suv chain C (HFE) vs chain A (TfR1): real experimental interface
  POS_HOMO     11qc chain A vs chain B: the official TfR1 homodimer interface
  NEG_APART    binder translated +500 A: interface area must collapse to 0
  SELF         SASA recomputed twice on one design must be bit-identical
  NAT_SEQ      natural human TfR1 189-383 sequence as the composition baseline

Outputs are written once; the script refuses to overwrite an existing file.
CPU only. No network, no GPU, no model calls.

Usage:
  python3 scripts/c4_independent_foldability_20261026.py --out results/c4_independent_check_20261026
  python3 scripts/c4_independent_foldability_20261026.py --self-test
"""

import argparse
import hashlib
import json
import math
import os
import sys
import warnings
from collections import Counter

import numpy as np

warnings.simplefilter("ignore")

from Bio.PDB import PDBParser
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.SASA import ShrakeRupley

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

POOL_DIR = os.path.join(REPO, "results/c4_apical_pdb_audit_20261025/all343")
LINKAGE = os.path.join(REPO, "results/c4_linkage_apical_20261025/linkage_apical.json")
V4_JSON = os.path.join(
    REPO,
    "results/c4_ranked_v4_apical_L6_20261025/c4_acid_ranked_v4_apical_L6_20261025.json",
)
CIF_1SUV = os.path.join(REPO, "data/1suv.cif")
CIF_11QC = os.path.join(REPO, "data/11qc.cif")

# Tien et al. 2013, theoretical maximum ASA (A^2), used for relative SASA.
MAX_ASA = {
    "ALA": 129.0, "ARG": 274.0, "ASN": 195.0, "ASP": 193.0, "CYS": 167.0,
    "GLN": 225.0, "GLU": 223.0, "GLY": 104.0, "HIS": 224.0, "ILE": 197.0,
    "LEU": 201.0, "LYS": 236.0, "MET": 224.0, "PHE": 240.0, "PRO": 159.0,
    "SER": 155.0, "THR": 172.0, "TRP": 285.0, "TYR": 263.0, "VAL": 174.0,
}

THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
    "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V",
}

APOLAR_AA = set("AVLIMFWCGP")
CHARGED_AA = set("DEKR")
ACIDIC_AA = set("DE")

BURIED_REL_SASA = 0.10  # residue is "buried" below this relative SASA

# Henderson-Hasselbalch pKa set (composition descriptor only).
PKA = {"D": 3.65, "E": 4.25, "H": 6.00, "C": 8.18, "Y": 10.07, "K": 10.53, "R": 12.48}
POS_RES = set("HKR")
NEG_RES = set("DECY")


def _json_default(o):
    """numpy scalars are not JSON-serializable; coerce them explicitly."""
    if isinstance(o, (np.bool_,)):
        return bool(o)
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating,)):
        return float(o)
    raise TypeError(f"Object of type {o.__class__.__name__} is not JSON serializable")


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def sha256_text(s):
    return hashlib.sha256(s.encode("utf-8")).hexdigest()


def atom_is_apolar(atom):
    """Standard apolar partition: carbon and sulfur."""
    el = (atom.element or "").strip().upper()
    if not el:
        el = atom.get_id()[0].upper()
    return el in ("C", "S")


def load_chain_atoms(structure, chain_id):
    """Return list of standard-residue atoms for one chain, first model."""
    model = next(structure.get_models())
    chain = model[chain_id]
    residues = [r for r in chain if r.id[0] == " " and r.get_resname() in MAX_ASA]
    atoms = []
    for r in residues:
        for a in r:
            if a.element != "H":
                atoms.append(a)
    return residues, atoms


class SASACalc:
    """Shrake-Rupley over an explicit atom list, via a detached container."""

    def __init__(self, n_points=100):
        self.sr = ShrakeRupley(n_points=n_points)

    def run(self, entity):
        self.sr.compute(entity, level="A")


def compute_sasa_for(entities):
    """
    Compute SASA in-place for a list of Bio.PDB entities treated as ONE body.
    Uses a synthetic Structure so Shrake-Rupley sees exactly these atoms.
    """
    from Bio.PDB.Structure import Structure
    from Bio.PDB.Model import Model
    from Bio.PDB.Chain import Chain

    st = Structure("tmp")
    md = Model(0)
    st.add(md)
    for i, residues in enumerate(entities):
        ch = Chain(chr(ord("A") + i))
        md.add(ch)
        for r in residues:
            ch.add(r.copy())
    sr = ShrakeRupley(n_points=100)
    sr.compute(st, level="A")
    # map back by (chain_index, residue_id, atom_name)
    out = {}
    for i, ch in enumerate(md):
        for r in ch:
            for a in r:
                out[(i, r.id, a.get_id())] = a.sasa
    return out


def residue_sasa(residues, sasa_map, chain_index):
    """Per-residue total SASA and relative SASA from an atom-level map."""
    per = []
    for r in residues:
        tot = 0.0
        apolar = 0.0
        for a in r:
            if a.element == "H":
                continue
            v = sasa_map.get((chain_index, r.id, a.get_id()))
            if v is None:
                continue
            tot += v
            if atom_is_apolar(a):
                apolar += v
        rel = tot / MAX_ASA[r.get_resname()]
        per.append({"resname": r.get_resname(), "sasa": tot, "apolar": apolar, "rel": rel})
    return per


def analyse_complex(binder_res, target_res):
    """
    Core measurement. Returns Family-1 and Family-2 metrics.
    Three SASA passes: binder alone, target alone, complex.
    """
    sasa_binder_alone = compute_sasa_for([binder_res])
    sasa_target_alone = compute_sasa_for([target_res])
    sasa_complex = compute_sasa_for([binder_res, target_res])

    per_b_alone = residue_sasa(binder_res, sasa_binder_alone, 0)
    per_b_cplx = residue_sasa(binder_res, sasa_complex, 0)
    per_t_alone = residue_sasa(target_res, sasa_target_alone, 0)
    per_t_cplx = residue_sasa(target_res, sasa_complex, 1)

    # ---- Family 1: binder monomer core quality (isolated binder only) ----
    buried = [p for p in per_b_alone if p["rel"] < BURIED_REL_SASA]
    n_buried = len(buried)
    buried_aa = [THREE_TO_ONE[p["resname"]] for p in buried]
    n_buried_charged = sum(1 for c in buried_aa if c in CHARGED_AA)
    n_buried_acidic = sum(1 for c in buried_aa if c in ACIDIC_AA)
    n_buried_apolar = sum(1 for c in buried_aa if c in APOLAR_AA)

    total_b_alone = sum(p["sasa"] for p in per_b_alone)
    apolar_b_alone = sum(p["apolar"] for p in per_b_alone)

    fam1 = {
        "binder_n_res": len(binder_res),
        "binder_sasa_alone_A2": round(total_b_alone, 1),
        "binder_apolar_sasa_frac": round(apolar_b_alone / total_b_alone, 4) if total_b_alone else None,
        "n_buried": n_buried,
        "buried_frac_of_binder": round(n_buried / len(binder_res), 4) if binder_res else None,
        "n_buried_apolar": n_buried_apolar,
        "n_buried_charged": n_buried_charged,
        "n_buried_acidic": n_buried_acidic,
        "core_apolar_frac": round(n_buried_apolar / n_buried, 4) if n_buried else None,
        "core_charged_frac": round(n_buried_charged / n_buried, 4) if n_buried else None,
    }

    # ---- Family 2: interface chemistry ----
    d_b_tot = sum(a["sasa"] - c["sasa"] for a, c in zip(per_b_alone, per_b_cplx))
    d_b_apo = sum(a["apolar"] - c["apolar"] for a, c in zip(per_b_alone, per_b_cplx))
    d_t_tot = sum(a["sasa"] - c["sasa"] for a, c in zip(per_t_alone, per_t_cplx))
    d_t_apo = sum(a["apolar"] - c["apolar"] for a, c in zip(per_t_alone, per_t_cplx))

    bsa_total = (d_b_tot + d_t_tot) / 2.0
    bsa_apolar = (d_b_apo + d_t_apo) / 2.0

    n_iface_res = sum(
        1 for a, c in zip(per_b_alone, per_b_cplx) if (a["sasa"] - c["sasa"]) > 1.0
    )

    fam2 = {
        "bsa_total_A2": round(bsa_total, 1),
        "bsa_apolar_A2": round(bsa_apolar, 1),
        "bsa_polar_A2": round(bsa_total - bsa_apolar, 1),
        "iface_apolar_frac": round(bsa_apolar / bsa_total, 4) if bsa_total > 1.0 else None,
        "binder_dsasa_A2": round(d_b_tot, 1),
        "n_iface_res_binder": n_iface_res,
    }
    return fam1, fam2


def seq_metrics(seq):
    """Family 3: pose-free sequence liabilities."""
    n = len(seq)
    cnt = Counter(seq)
    # Shannon entropy over observed amino acids, in bits
    ent = -sum((c / n) * math.log2(c / n) for c in cnt.values())
    # longest homopolymer run
    best, cur, best_c = 1, 1, seq[0]
    for i in range(1, n):
        if seq[i] == seq[i - 1]:
            cur += 1
            if cur > best:
                best, best_c = cur, seq[i]
        else:
            cur = 1
    # net charge, Henderson-Hasselbalch (COMPOSITION DESCRIPTOR ONLY)
    def net_charge(ph):
        q = 0.0
        for c in seq:
            if c in POS_RES:
                q += 1.0 / (1.0 + 10 ** (ph - PKA[c]))
            elif c in NEG_RES:
                q -= 1.0 / (1.0 + 10 ** (PKA[c] - ph))
        q += 1.0 / (1.0 + 10 ** (ph - 8.0))    # N-terminus
        q -= 1.0 / (1.0 + 10 ** (3.65 - ph))   # C-terminus
        return q

    return {
        "length": n,
        "entropy_bits": round(ent, 4),
        "n_distinct_aa": len(cnt),
        "max_homopolymer_run": best,
        "max_homopolymer_aa": best_c,
        "frac_acidic_DE": round(sum(cnt[c] for c in "DE") / n, 4),
        "frac_basic_KR": round(sum(cnt[c] for c in "KR") / n, 4),
        "frac_his": round(cnt["H"] / n, 4),
        "frac_apolar": round(sum(cnt[c] for c in APOLAR_AA) / n, 4),
        "net_charge_pH7.4_composition_only": round(net_charge(7.4), 2),
        "net_charge_pH6.0_composition_only": round(net_charge(6.0), 2),
    }


def chain_sequence(residues):
    return "".join(THREE_TO_ONE[r.get_resname()] for r in residues)


# ---------------------------------------------------------------- controls


def run_controls():
    """Every control must fire. Returns (all_ok, report_dict)."""
    rep = {}
    ok = True

    # POS_HETERO: real experimental TfR1-HFE interface (1suv A vs C)
    st = MMCIFParser(QUIET=True).get_structure("1suv", CIF_1SUV)
    a_res, _ = load_chain_atoms(st, "A")
    c_res, _ = load_chain_atoms(st, "C")
    f1, f2 = analyse_complex(c_res, a_res)
    rep["POS_HETERO_1suv_C_vs_A"] = {
        "bsa_total_A2": f2["bsa_total_A2"],
        "iface_apolar_frac": f2["iface_apolar_frac"],
        "fired": bool(f2["bsa_total_A2"] > 300 and f2["iface_apolar_frac"] is not None),
    }
    ok &= rep["POS_HETERO_1suv_C_vs_A"]["fired"]

    # POS_HOMO: official TfR1 homodimer interface (11qc A vs B)
    st2 = MMCIFParser(QUIET=True).get_structure("11qc", CIF_11QC)
    qa, _ = load_chain_atoms(st2, "A")
    qb, _ = load_chain_atoms(st2, "B")
    f1h, f2h = analyse_complex(qa, qb)
    rep["POS_HOMO_11qc_A_vs_B"] = {
        "bsa_total_A2": f2h["bsa_total_A2"],
        "iface_apolar_frac": f2h["iface_apolar_frac"],
        "chain_A_n_res": len(qa),
        "fired": bool(f2h["bsa_total_A2"] > 300),
    }
    ok &= rep["POS_HOMO_11qc_A_vs_B"]["fired"]
    rep["natural_core_reference_11qc_A"] = {
        "core_apolar_frac": f1h["core_apolar_frac"],
        "core_charged_frac": f1h["core_charged_frac"],
        "buried_frac_of_chain": f1h["buried_frac_of_binder"],
    }

    # NEG_APART: translate the binder far away; interface must collapse
    pdb = os.path.join(POOL_DIR, "tfr1_c4_acid_site1_123_complex.pdb")
    sd = PDBParser(QUIET=True).get_structure("d", pdb)
    b_res, _ = load_chain_atoms(sd, "A")
    t_res, _ = load_chain_atoms(sd, "B")
    _, f2_near = analyse_complex(b_res, t_res)

    sd2 = PDBParser(QUIET=True).get_structure("d", pdb)
    b2, _ = load_chain_atoms(sd2, "A")
    t2, _ = load_chain_atoms(sd2, "B")
    for r in b2:
        for a in r:
            a.set_coord(a.get_coord() + np.array([500.0, 0.0, 0.0]))
    _, f2_far = analyse_complex(b2, t2)
    rep["NEG_APART_translated_500A"] = {
        "bsa_total_A2_contact": f2_near["bsa_total_A2"],
        "bsa_total_A2_apart": f2_far["bsa_total_A2"],
        "fired": bool(f2_near["bsa_total_A2"] > 100 and abs(f2_far["bsa_total_A2"]) < 1.0),
    }
    ok &= rep["NEG_APART_translated_500A"]["fired"]

    # SELF: recomputation must be identical
    sd3 = PDBParser(QUIET=True).get_structure("d", pdb)
    b3, _ = load_chain_atoms(sd3, "A")
    t3, _ = load_chain_atoms(sd3, "B")
    _, f2_again = analyse_complex(b3, t3)
    rep["SELF_recompute_identical"] = {
        "first": f2_near["bsa_total_A2"],
        "second": f2_again["bsa_total_A2"],
        "fired": bool(f2_near["bsa_total_A2"] == f2_again["bsa_total_A2"]),
    }
    ok &= rep["SELF_recompute_identical"]["fired"]

    # NAT_SEQ: natural human TfR1 apical span as the composition baseline
    nat = chain_sequence(t_res)
    rep["NAT_SEQ_human_apical_189_383"] = seq_metrics(nat)
    rep["NAT_SEQ_human_apical_189_383"]["n_res"] = len(nat)
    rep["NAT_SEQ_human_apical_189_383"]["fired"] = bool(len(nat) == 195)

    ok &= rep["NAT_SEQ_human_apical_189_383"]["fired"]

    rep["ALL_CONTROLS_FIRED"] = bool(ok)
    return bool(ok), rep


# ---------------------------------------------------------------- main


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", help="output directory (must not already contain outputs)")
    ap.add_argument("--self-test", action="store_true", help="run controls only")
    ap.add_argument("--limit", type=int, default=0, help="debug: only N designs")
    args = ap.parse_args()

    if args.self_test:
        ok, rep = run_controls()
        print(json.dumps(rep, indent=2, default=_json_default))
        print("\nALL CONTROLS FIRED:", ok)
        sys.exit(0 if ok else 1)

    if not args.out:
        ap.error("--out is required (or use --self-test)")

    os.makedirs(args.out, exist_ok=True)
    out_json = os.path.join(args.out, "independent_metrics.json")
    out_csv = os.path.join(args.out, "independent_metrics.csv")
    for p in (out_json, out_csv):
        if os.path.exists(p):
            sys.exit(f"REFUSING TO OVERWRITE: {p}")

    print("running controls first ...", flush=True)
    ok, controls = run_controls()
    if not ok:
        sys.exit("CONTROLS DID NOT ALL FIRE - refusing to emit measurements")
    print("controls OK", flush=True)

    link = json.load(open(LINKAGE))
    recs = {r["name"]: r for r in link["records"]}
    v4 = json.load(open(V4_JSON))
    v4_names = [r["name"] for r in v4["chosen_by_apical_L6"]]
    v4_rank = {r["name"]: r["rank"] for r in v4["chosen_by_apical_L6"]}

    names = sorted(recs)
    if args.limit:
        # keep all v4 rows plus a sample, so the comparison stays meaningful
        rest = [n for n in names if n not in v4_rank][: args.limit]
        names = sorted(set(v4_names) | set(rest))

    rows = []
    for i, name in enumerate(names, 1):
        pdb = os.path.join(POOL_DIR, f"{name}_complex.pdb")
        if not os.path.exists(pdb):
            continue
        st = PDBParser(QUIET=True).get_structure(name, pdb)
        b_res, _ = load_chain_atoms(st, "A")
        t_res, _ = load_chain_atoms(st, "B")
        if not b_res or not t_res:
            continue
        seq = chain_sequence(b_res)
        f1, f2 = analyse_complex(b_res, t_res)
        f3 = seq_metrics(seq)
        row = {
            "name": name,
            "sequence_sha256": sha256_text(seq),
            "pdb_sha256": sha256_file(pdb),
            "in_v4": name in v4_rank,
            "v4_rank": v4_rank.get(name),
            "L6_apical": recs[name]["L6"],
            "iptm_apical": round(recs[name]["iptm"], 4),
            "receptor_ca": recs[name]["receptor_ca"],
        }
        row.update(f1)
        row.update(f2)
        row.update(f3)
        rows.append(row)
        if i % 25 == 0:
            print(f"  {i}/{len(names)}", flush=True)

    payload = {
        "created_at_utc": os.popen("date -u +%Y-%m-%dT%H:%M:%SZ").read().strip(),
        "purpose": "Independent out-of-sample check on the v4 C4 ranking. "
                   "Signals here entered none of v4's selection (iPTM, apical L6, 2.0 A clash gate).",
        "method": {
            "sasa": "Bio.PDB.SASA.ShrakeRupley, n_points=100, heavy atoms only",
            "buried_rel_sasa_threshold": BURIED_REL_SASA,
            "max_asa_reference": "Tien et al. 2013 theoretical",
            "apolar_atom_partition": "element C or S",
            "bsa": "((dSASA_binder) + (dSASA_target)) / 2, binder/target alone vs complex",
            "net_charge": "Henderson-Hasselbalch; COMPOSITION DESCRIPTOR ONLY, not a pH-mechanism signal",
        },
        "inputs": {
            "pool_dir": os.path.relpath(POOL_DIR, REPO),
            "linkage_json": os.path.relpath(LINKAGE, REPO),
            "linkage_sha256": sha256_file(LINKAGE),
            "v4_json": os.path.relpath(V4_JSON, REPO),
            "v4_sha256": sha256_file(V4_JSON),
            "cif_1suv_sha256": sha256_file(CIF_1SUV),
            "cif_11qc_sha256": sha256_file(CIF_11QC),
        },
        "controls": controls,
        "n_designs": len(rows),
        "n_in_v4": sum(1 for r in rows if r["in_v4"]),
        "designs": rows,
    }

    with open(out_json, "w") as fh:
        json.dump(payload, fh, indent=1, default=_json_default)

    if rows:
        import csv
        keys = list(rows[0].keys())
        with open(out_csv, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=keys)
            w.writeheader()
            for r in rows:
                w.writerow(r)

    print(f"\nwrote {out_json}  ({len(rows)} designs, {payload['n_in_v4']} in v4)")
    print(f"wrote {out_csv}")


if __name__ == "__main__":
    main()
