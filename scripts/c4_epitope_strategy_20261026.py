#!/usr/bin/env python3
"""
C4 strategy question: can the current design approach cross-react with mouse TfR1?

The challenge ranks mouse cross-reactivity as its own objective. The prior cycle
measured per-design mouse pose compatibility for the 20 v4 seats. This script asks
the strategic question instead, across the whole 342-design pool:

  Where on the human apical surface do these designs actually bind, how mouse-
  conserved is that surface, and does a contactable mouse-conserved patch exist
  that the campaign is not using?

That distinguishes two very different answers:
  (i) the surface itself is substituted everywhere a binder could sit
      -> the apical domain cannot support cross-reactivity, change target domain
  (ii) the campaign converged on a substituted patch while a conserved patch exists
      -> the strategy is not capped; the epitope choice is, and it is fixable

Numbering: chain B residue i (1..195) maps to human P02786 author position 188+i.
Verified against the human_aa column of the conservation alignment before use;
the script exits if that control does not fire.

Contacts: any binder heavy atom within 4.5 A of any target heavy atom (same
convention as results/c4_mouse_compatibility_20261026).

CPU only. Read-only on all inputs. No ranking is produced.
"""

import argparse
import ast
import csv
import hashlib
import json
import os
import sys
import warnings
from collections import Counter, defaultdict

import numpy as np

warnings.simplefilter("ignore")
from Bio.PDB import PDBParser

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
POOL_DIR = os.path.join(REPO, "results/c4_apical_pdb_audit_20261025/all343")
ALIGN_CSV = os.path.join(REPO, "results/mouse_epitope_conservation/apical_alignment.csv")
LINKAGE = os.path.join(REPO, "results/c4_linkage_apical_20261025/linkage_apical.json")
V4_JSON = os.path.join(
    REPO,
    "results/c4_ranked_v4_apical_L6_20261025/c4_acid_ranked_v4_apical_L6_20261025.json",
)

CONTACT_CUTOFF = 4.5
APICAL_START = 189  # human author numbering of chain B residue 1

THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
    "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V",
}
CHARGE = {"D": -1, "E": -1, "K": +1, "R": +1, "H": 0}


def sha256_file(p):
    return hashlib.sha256(open(p, "rb").read()).hexdigest()


def load_alignment():
    """human_auth -> dict(human_aa, mouse_aa_options, identical, ambiguous)"""
    out = {}
    with open(ALIGN_CSV) as fh:
        for r in csv.DictReader(fh):
            auth = int(r["human_auth"])
            opts = ast.literal_eval(r["mouse_aa_options"])
            out[auth] = {
                "human_aa": r["human_aa"],
                "mouse_aa_options": [o for o in opts],
                "identical": r["identical"] == "True",
                "ambiguous": r["mapping_ambiguous"] == "True",
            }
    return out


def chain_residues(struct, cid):
    model = next(struct.get_models())
    return [r for r in model[cid] if r.id[0] == " " and r.get_resname() in THREE_TO_ONE]


def contacts(binder_res, target_res, cutoff=CONTACT_CUTOFF):
    """Return set of target residue numbers contacted by any binder heavy atom."""
    bcoords = []
    for r in binder_res:
        for a in r:
            if a.element != "H":
                bcoords.append(a.get_coord())
    if not bcoords:
        return set()
    B = np.asarray(bcoords, float)
    hit = set()
    c2 = cutoff * cutoff
    for r in target_res:
        for a in r:
            if a.element == "H":
                continue
            d2 = ((B - a.get_coord()) ** 2).sum(axis=1)
            if (d2 <= c2).any():
                hit.add(r.id[1])
                break
    return hit


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    outdir = os.path.join(REPO, args.out)
    os.makedirs(outdir, exist_ok=True)
    out_json = os.path.join(outdir, "epitope_strategy.json")
    out_csv = os.path.join(outdir, "epitope_strategy.csv")
    for p in (out_json, out_csv):
        if os.path.exists(p):
            sys.exit(f"REFUSING TO OVERWRITE: {p}")

    align = load_alignment()
    link = json.load(open(LINKAGE))
    recs = {r["name"]: r for r in link["records"]}
    v4n = {r["name"] for r in json.load(open(V4_JSON))["chosen_by_apical_L6"]}

    # ---- CONTROL: numbering map must reproduce the target sequence ----
    probe = os.path.join(POOL_DIR, "tfr1_c4_acid_site1_123_complex.pdb")
    st = PDBParser(QUIET=True).get_structure("p", probe)
    t_res = chain_residues(st, "B")
    mism = []
    for r in t_res:
        auth = APICAL_START + r.id[1] - 1
        aa = THREE_TO_ONE[r.get_resname()]
        if auth not in align:
            mism.append((r.id[1], auth, aa, "absent_from_alignment"))
        elif align[auth]["human_aa"] != aa:
            mism.append((r.id[1], auth, aa, align[auth]["human_aa"]))
    control = {
        "n_target_res": len(t_res),
        "n_mismatch": len(mism),
        "first_mismatches": mism[:5],
        "fired": len(t_res) == 195 and len(mism) == 0,
    }
    if not control["fired"]:
        print(json.dumps(control, indent=2))
        sys.exit("NUMBERING CONTROL DID NOT FIRE - refusing to emit measurements")
    print(f"numbering control OK: 195/195 target residues match human 189-383")

    # ---- per-design footprints ----
    rows = []
    position_hits = Counter()
    position_hits_v4 = Counter()
    names = sorted(recs)
    for i, name in enumerate(names, 1):
        pdb = os.path.join(POOL_DIR, f"{name}_complex.pdb")
        if not os.path.exists(pdb):
            continue
        s = PDBParser(QUIET=True).get_structure(name, pdb)
        b = chain_residues(s, "A")
        t = chain_residues(s, "B")
        if not b or not t:
            continue
        hits = contacts(b, t)
        auths = sorted(APICAL_START + h - 1 for h in hits)
        n_ident = sum(1 for a in auths if align[a]["identical"])
        n_amb = sum(1 for a in auths if align[a]["ambiguous"])
        # charge reversal: human charged -> mouse opposite sign, any optimum
        n_rev = 0
        n_sub = 0
        for a in auths:
            e = align[a]
            if e["identical"]:
                continue
            n_sub += 1
            hq = CHARGE.get(e["human_aa"], 0)
            for m in e["mouse_aa_options"]:
                mq = CHARGE.get(m, 0)
                if hq != 0 and mq != 0 and hq != mq:
                    n_rev += 1
                    break
        seq = "".join(THREE_TO_ONE[r.get_resname()] for r in b)
        row = {
            "name": name,
            "site": "site1" if "_site1_" in name else ("site2" if "_site2_" in name else "other"),
            "sequence_sha256": hashlib.sha256(seq.encode()).hexdigest(),
            "in_v4": name in v4n,
            "L6_apical": recs[name]["L6"],
            "iptm_apical": round(recs[name]["iptm"], 4),
            "n_contact_residues": len(auths),
            "n_identical_in_mouse": n_ident,
            "n_substituted": n_sub,
            "n_charge_reversal": n_rev,
            "n_ambiguous": n_amb,
            "frac_identical": round(n_ident / len(auths), 4) if auths else None,
            "footprint_human_auth": ";".join(str(a) for a in auths),
        }
        rows.append(row)
        for a in auths:
            position_hits[a] += 1
            if name in v4n:
                position_hits_v4[a] += 1
        if i % 50 == 0:
            print(f"  {i}/{len(names)}", flush=True)

    # ---- per-position summary: contactability x mouse identity ----
    per_pos = []
    for auth in sorted(align):
        if auth < APICAL_START or auth > APICAL_START + 194:
            continue
        e = align[auth]
        per_pos.append({
            "human_auth": auth,
            "human_aa": e["human_aa"],
            "mouse_aa_options": "/".join(e["mouse_aa_options"]) if e["mouse_aa_options"] else "gap",
            "identical": e["identical"],
            "ambiguous": e["ambiguous"],
            "n_designs_contacting": position_hits.get(auth, 0),
            "n_v4_contacting": position_hits_v4.get(auth, 0),
            "frac_pool_contacting": round(position_hits.get(auth, 0) / len(rows), 4) if rows else None,
        })

    n_ident_pos = sum(1 for p in per_pos if p["identical"])
    contacted = [p for p in per_pos if p["n_designs_contacting"] > 0]
    hot = [p for p in per_pos if p["frac_pool_contacting"] and p["frac_pool_contacting"] >= 0.10]

    fi = np.array([r["frac_identical"] for r in rows if r["frac_identical"] is not None])

    summary = {
        "n_designs": len(rows),
        "apical_positions": len(per_pos),
        "apical_positions_identical_in_mouse": n_ident_pos,
        "apical_identity_fraction": round(n_ident_pos / len(per_pos), 4),
        "positions_contacted_by_at_least_one_design": len(contacted),
        "positions_contacted_identical": sum(1 for p in contacted if p["identical"]),
        "positions_contacted_identity_fraction": round(
            sum(1 for p in contacted if p["identical"]) / len(contacted), 4) if contacted else None,
        "hot_positions_ge10pct_of_pool": len(hot),
        "hot_positions_identical": sum(1 for p in hot if p["identical"]),
        "hot_positions_identity_fraction": round(
            sum(1 for p in hot if p["identical"]) / len(hot), 4) if hot else None,
        "per_design_frac_identical": {
            "median": round(float(np.median(fi)), 4),
            "p10": round(float(np.percentile(fi, 10)), 4),
            "p90": round(float(np.percentile(fi, 90)), 4),
            "min": round(float(fi.min()), 4),
            "max": round(float(fi.max()), 4),
            "n_designs_fully_conserved_footprint": int((fi == 1.0).sum()),
            "n_designs_ge_0.90": int((fi >= 0.90).sum()),
            "n_designs_ge_0.80": int((fi >= 0.80).sum()),
        },
        "by_site": {},
    }
    for site in ("site1", "site2", "other"):
        sub = [r for r in rows if r["site"] == site]
        if not sub:
            continue
        f = np.array([r["frac_identical"] for r in sub if r["frac_identical"] is not None])
        summary["by_site"][site] = {
            "n": len(sub),
            "median_frac_identical": round(float(np.median(f)), 4),
            "max_frac_identical": round(float(f.max()), 4),
        }

    payload = {
        "created_at_utc": os.popen("date -u +%Y-%m-%dT%H:%M:%SZ").read().strip(),
        "question": "Can a fixed scaffold family on the human apical patch cross-react with mouse TfR1?",
        "method": {
            "contacts": f"any binder heavy atom within {CONTACT_CUTOFF} A of any target heavy atom",
            "numbering": "chain B residue i -> human P02786 author 188+i; verified against alignment",
            "mouse_identity": "all optimal BLOSUM62 alignments retained; identical only if identical in every optimum",
            "charge_reversal": "human charged residue -> opposite-sign mouse residue in any optimum",
        },
        "inputs": {
            "pool_dir": os.path.relpath(POOL_DIR, REPO),
            "alignment_csv": os.path.relpath(ALIGN_CSV, REPO),
            "alignment_sha256": sha256_file(ALIGN_CSV),
            "linkage_sha256": sha256_file(LINKAGE),
            "v4_sha256": sha256_file(V4_JSON),
        },
        "numbering_control": control,
        "summary": summary,
        "per_position": per_pos,
        "designs": rows,
    }
    json.dump(payload, open(out_json, "w"), indent=1)

    with open(out_csv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        for r in rows:
            w.writerow(r)

    print(f"\nwrote {out_json}\nwrote {out_csv}")
    print(json.dumps(summary, indent=1))


if __name__ == "__main__":
    main()
