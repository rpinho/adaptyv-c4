#!/usr/bin/env python3
"""C4 rule-v4 pre-registration freeze: candidate-source manifest + model inventory.

Runs BEFORE any scoring, gate or ranking. Records SHA256 of every input the rule
v4 ranking will read, plus what models actually exist for each design. Writes
results/c4_ranking_v4_acid_off_20261026/manifest_frozen.json.

$0, CPU only, read-only with respect to every input.
"""
import hashlib
import json
import os
import re
import pathlib
import subprocess

ROOT = pathlib.Path(__file__).resolve().parents[1]
OUT = ROOT / "results/c4_ranking_v4_acid_off_20261026/manifest_frozen.json"


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def rel(p):
    return str(pathlib.Path(p).resolve().relative_to(ROOT))


INPUTS = [
    "results/decisions/20261026-c4-selection-rule-v4.md",
    "results/decisions/20261026-c4-selection-rule-v3.md",
    "results/c4_linkage_apical_20261025/linkage_apical.json",
    "results/c4_acid_combined.csv",
    "results/c4_acid_r2_binders.csv",
    "results/c4_apical_pdb_audit_20261025/structural_recheck.json",
    "results/c4_mouse_compatibility_20261026/per_design_assessment.csv",
    "results/c4_novelty_screen_20261025/summary.json",
    "results/c4_novelty_screen_20261025/membership.json",
    "results/c4_provenance_20261025/summary.json",
    "results/c4_acid_off_scored.json",
    "results/c4_independent_check_20261026/independent_metrics.json",
    "results/c4_ranked_v4_apical_L6_20261025/c4_acid_ranked_v4_apical_L6_20261025.json",
    "data/P02786.txt",
    "data/Q62351.txt",
    "scripts/c4_v4_freeze_20261026.py",
]

# Verified to exist and match on 2026-10-26. v1 of this manifest used wrong paths and was
# superseded before any scoring ran (manifest_frozen_v1_SUPERSEDED.json).
FROZEN_GUARD = {
    "results/c3_submission_v8_20261025.csv": "631ed0bf",
    "results/c3_submission_v7_20261025.csv": "1d99f843",
    "results/c3_submission_v6_20261023.csv": "8e9b2d74",
    "results/submission_v9.csv": "31d98a41",
}

APICAL_DIR = ROOT / "results/c4_apical_pdb_audit_20261025/all343"


def seq_of_pdb_chain(path, chain):
    """Chain sequence from CA records, in file order."""
    out = []
    seen = set()
    for line in open(path, errors="replace"):
        if line.startswith("ATOM") and line[12:16].strip() == "CA" and line[21] == chain:
            key = line[22:27]
            if key in seen:
                continue
            seen.add(key)
            out.append(line[17:20].strip())
    return "".join({"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
                    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
                    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
                    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V"}.get(a, "X")
                   for a in out)


def main():
    if OUT.exists():
        raise SystemExit(f"REFUSING TO OVERWRITE frozen manifest: {OUT}")

    link = json.load(open(ROOT / "results/c4_linkage_apical_20261025/linkage_apical.json"))
    recs = link["records"]

    # binder sequences keyed by sha256(canonical uppercase sequence)
    seq_by_name = {}
    for src in ("results/c4_acid_combined.csv", "results/c4_acid_r2_binders.csv"):
        for line in open(ROOT / src).read().splitlines()[1:]:
            if not line.strip():
                continue
            name, seq = line.split(",")[0], line.split(",")[1]
            seq_by_name.setdefault(name, "".join(seq.split()).upper())

    inv = []
    for r in recs:
        name = r["name"]
        pdb = APICAL_DIR / f"{name}_complex.pdb"
        ent = {
            "name": name,
            "sequence_sha256": hashlib.sha256(seq_by_name.get(name, "").encode()).hexdigest()
            if seq_by_name.get(name) else None,
            "length": len(seq_by_name.get(name, "")),
            "pdb": rel(pdb) if pdb.exists() else None,
            "pdb_sha256": sha(pdb) if pdb.exists() else None,
            "models_available": [],
            "iptm_models": 1 if r.get("iptm") is not None else 0,
            "l6_models": 1 if r.get("L6") is not None else 0,
        }
        if pdb.exists():
            ent["models_available"].append({
                "model": "boltz2_apical_195",
                "source": rel(pdb),
                "sha256": ent["pdb_sha256"],
                "chain_A_residues": len(set(l[22:27] for l in open(pdb, errors="replace")
                                             if l.startswith("ATOM") and l[21] == "A")),
                "chain_B_residues": len(set(l[22:27] for l in open(pdb, errors="replace")
                                             if l.startswith("ATOM") and l[21] == "B")),
                "has_bfactor_plddt": True,
            })
            got = seq_of_pdb_chain(pdb, "A")
            ent["pdb_chainA_matches_csv"] = (got == seq_by_name.get(name, ""))
        inv.append(ent)

    target = None
    lines = (ROOT / "data/P02786.txt").read_text().splitlines()
    inseq = False
    buf = []
    for l in lines:
        if l.startswith("SQ   "):
            inseq = True
            continue
        if inseq:
            if l.startswith("//"):
                break
            buf.append("".join(c for c in l if c.isupper()))
    full = "".join(buf)
    target = full[188:383]

    payload = {
        "created_at_utc": subprocess.check_output(["date", "-u", "+%Y-%m-%dT%H:%M:%SZ"]).decode().strip(),
        "purpose": "Pre-registration freeze for C4 selection rule v4 (acid-OFF). Written BEFORE scoring.",
        "rule": {
            "file": "results/decisions/20261026-c4-selection-rule-v4.md",
            "sha256": sha(ROOT / "results/decisions/20261026-c4-selection-rule-v4.md"),
            "parent_v3_sha256": sha(ROOT / "results/decisions/20261026-c4-selection-rule-v3.md"),
            "eligibility_direction": "controlled human apical L6 < 0 (ACID-OFF)",
            "score": "S = 30P + 25M + 20D + 15C + 10N",
            "seats": 20,
            "greedy": "by S, recompute D and N per seat",
        },
        "candidate_source_manifest": {
            "pool_definition": "342 designs with a recomputed apical L6 record "
                               "(results/c4_linkage_apical_20261025/linkage_apical.json records)",
            "n_pool": len(recs),
            "sequence_source": ["results/c4_acid_combined.csv", "results/c4_acid_r2_binders.csv"],
            "target": {
                "accession": "P02786",
                "residues": "189-383 (1-based, inclusive)",
                "length": len(target),
                "sha256": hashlib.sha256(target.encode()).hexdigest(),
            },
            "inputs": [{"path": p, "sha256": sha(ROOT / p)} for p in INPUTS if (ROOT / p).exists()],
            "missing_inputs": [p for p in INPUTS if not (ROOT / p).exists()],
        },
        "available_model_inventory": {
            "complex_models_per_design": 1,
            "model_names": ["boltz2_apical_195"],
            "note": "Every design has exactly ONE original human apical complex prediction "
                    "(Boltz-2, chain B = 195-residue P02786 189-383). No design has a second "
                    "original complex model, so 'all available models' == one model and "
                    "'minimum over available models' is not a min-5.",
            "n_designs_with_1_model": sum(1 for e in inv if len(e["models_available"]) == 1),
            "n_designs_with_0_models": sum(1 for e in inv if not e["models_available"]),
            "designs_by_model_count": {"1": sum(1 for e in inv if len(e["models_available"]) == 1),
                                       "0": sum(1 for e in inv if not e["models_available"])},
            "mouse_models": {
                "audit": "results/c4_mouse_compatibility_20261026/per_design_assessment.csv",
                "n_designs_audited": 20,
                "coverage": "the 20 incumbent acid-ON seats ONLY; 0 of the v4 acid-OFF eligible set",
            },
            "novelty_screen": {
                "file": "results/c4_novelty_screen_20261025/summary.json",
                "n_designs_screened": 35,
                "control_status": "strict search gave zero hits in all four DBs; validated by a "
                                  "positive control in the same run (real swissprot sequences "
                                  "self-hit at fident 1.000)",
            },
            "second_ph_estimator": {
                "exists_on_disk": True,
                "file": "results/c4_acid_off_scored.json",
                "implementation": "scripts/score_acid_off.py (His-distance/burial switch score)",
                "controls_recorded": False,
                "consequence_for_P": "P=0.5 requires a CONTROLLED second estimator. This file has "
                                     "no recorded controls, so it cannot by itself grant P=0.5.",
            },
            "gpu": "none used; no folding engine in this stack accepts pH as input",
        },
        "design_inventory": inv,
        "frozen_files_guard": {},
    }

    for p, pref in FROZEN_GUARD.items():
        full_p = ROOT / p
        if full_p.exists():
            payload["frozen_files_guard"][p] = {"sha256": sha(full_p), "prefix_expected": pref,
                                                "prefix_matches": sha(full_p).startswith(pref)}
        else:
            payload["frozen_files_guard"][p] = {"sha256": None, "note": "not present at this path"}

    OUT.write_text(json.dumps(payload, indent=1))
    print("wrote", rel(OUT))
    print("pool", len(recs), "| 1-model", payload["available_model_inventory"]["n_designs_with_1_model"],
          "| 0-model", payload["available_model_inventory"]["n_designs_with_0_models"])
    print("target sha", payload["candidate_source_manifest"]["target"]["sha256"])
    print("manifest sha", sha(OUT))
    for p, v in payload["frozen_files_guard"].items():
        print("  ", p, v)


if __name__ == "__main__":
    main()
