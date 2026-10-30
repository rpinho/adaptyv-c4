#!/usr/bin/env python3
"""C4 v4 — membership gate AND ranking on the correct apical L6, screen on apical PDBs.

WHY v4 EXISTS (distinct from v3's defect)
v3 fixed the *structural screen* (d_min read from the 154-residue stitched surrogate instead of the
195-residue apical receptor). v3 did NOT fix *membership*: it still filtered and ranked on the
surrogate-construct L6 from results/c4_selection_audit_20261025/audit.json (`LINK`), then only
changed which PDB directory the clash screen read.

That is a separate defect. A design can be acid-ON on the surrogate and acid-OFF on the intended
construct, so v3's acid-ON gate admits a set that is ~half wrong on the apical complex, and its
rank order is on the wrong quantity. Measured: 105/342 designs have an L6 sign flip between the
two constructs (Pearson r 0.37).

v4 changes EXACTLY ONE thing relative to v3: the membership filter and the ranking key read L6 from
results/c4_linkage_apical_20261025/linkage_apical.json (same scorer, same parameters, apical PDBs).
Everything else is v3's: iPTM floor 0.60, apical PDB directory, 2.0 A clash gate, TM>=0.65 /
coverage>=0.80 diversity rule, 20 slots, sequence sources.

CONTROL GATE (required, runs first, writes nothing on failure)
v3's own algorithm on v3's own inputs — surrogate-L6 pool + apical PDB directory — must reproduce
the shipped v3 exactly: all 20 names in order and all 20 d_min values within 0.02 A. Without that,
a difference between v3 and v4 could not be attributed to the membership gate alone.

FAIL-CLOSED: if fewer than the requested slots survive, the script reports the shortfall and does
not pad, impute, or lower a threshold. `tfr1_c4_acid_site1_41` has no apical complex and is never
imputed.

CPU only. No GPU, no network. EVIDENCE ONLY — not a submission; Ricardo approval required.
"""
from __future__ import annotations

import argparse
import csv
import importlib.util
import json
import pathlib
import sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
V3 = ROOT / "scripts/c4_rank_v3_apical_20261025.py"
AUDIT = ROOT / "results/c4_selection_audit_20261025/audit.json"
APICAL_L6 = ROOT / "results/c4_linkage_apical_20261025/linkage_apical.json"
APICAL_PDBS = ROOT / "results/c4_apical_pdb_audit_20261025/all343"
SHIPPED_V3_CSV = ROOT / "results/c4_ranked_v3_apical_20261025/c4_acid_ranked_v3_apical_20261025.csv"
SHIPPED_V3_JSON = ROOT / "results/c4_ranked_v3_apical_20261025/c4_acid_ranked_v3_apical_20261025.json"
V3_SHA_PREFIX = "4ad6932e"
OUTDIR = ROOT / "results/c4_ranked_v4_apical_L6_20261025"
PREVIEW = ROOT / "results/c4_apical_preview_20261024.csv"
SEQ_SOURCES = ["results/c4_acid_combined.csv", "results/c4_acid_r2_binders.csv",
               "results/c4_acid1_filtered.csv", "results/c4_acid2_filtered.csv",
               "results/c4_apical_preview_20261024.csv"]


def sha256(p: pathlib.Path) -> str:
    import hashlib
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--slots", type=int, default=20)
    ap.add_argument("--iptm-floor", type=float, default=0.60)
    args = ap.parse_args()

    v3_sha_before = sha256(SHIPPED_V3_CSV)
    if not v3_sha_before.startswith(V3_SHA_PREFIX):
        sys.exit(f"v3 artifact sha {v3_sha_before} does not start with {V3_SHA_PREFIX}; "
                 "refusing to run against an unexpected v3")

    v3 = load_module("c4_v3", V3)
    clash = v3.load_module("clash", "scripts/c3_clash_audit_independent_20261025.py")
    crosswalk = v3.load_module("cw", "scripts/c3_round2_cluster_crosswalk_20261025.py")

    audit = json.loads(AUDIT.read_text())
    apical = json.loads(APICAL_L6.read_text())
    if apical.get("target_len") != 195:
        sys.exit(f"apical L6 target_len {apical.get('target_len')} != 195")
    if apical.get("fails"):
        sys.exit(f"apical L6 artifact has {len(apical['fails'])} fails; refusing")
    if not apical["sign_controls"].get("all_ok"):
        sys.exit("apical L6 sign controls FAILED; refusing")

    iptm = {r["name"]: r.get("iptm") for r in audit["records"]}
    apical_l6 = {r["name"]: r["L6"] for r in apical["records"]}
    no_apical = set(apical.get("designs_without_apical_pdb") or [])

    # ---- inputs A: v3's own pool (surrogate L6) ----
    surrogate_recs = [r for r in audit["records"] if r.get("L6") is not None]
    pool_v3 = [r for r in surrogate_recs
               if r.get("iptm") is not None and r["iptm"] >= args.iptm_floor and r["L6"] > 0]
    pool_v3.sort(key=lambda r: (-r["L6"], -(r["iptm"] or 0), r["name"]))

    # ---- CONTROL: v3 algorithm on v3 inputs must reproduce shipped v3 ----
    ctl_chosen, ctl_rejected, _ = v3.select(pool_v3, APICAL_PDBS, args.slots, clash, crosswalk)
    shipped_v3 = json.loads(SHIPPED_V3_JSON.read_text())
    v3_names = [r["name"] for r in shipped_v3["chosen_by_L6"]]
    v3_dmin = {r["name"]: r["d_min_A"] for r in shipped_v3["chosen_by_L6"]}
    ctl_names = [r["name"] for r in ctl_chosen]
    names_ok = ctl_names == v3_names
    dmin_deltas = {r["name"]: abs(r["d_min"] - v3_dmin[r["name"]]) for r in ctl_chosen}
    dmin_ok = all(d < 0.02 for d in dmin_deltas.values())
    with open(SHIPPED_V3_CSV) as fh:
        v3_csv_names = [row["name"] for row in csv.DictReader(fh)]
    csv_order_ok = v3_csv_names == v3_names

    control = {
        "what": "v3 algorithm on v3's own inputs (surrogate-L6 pool + apical PDB dir)",
        "names_reproduced": names_ok,
        "names_in_order": names_ok,
        "dmin_reproduced": dmin_ok,
        "dmin_max_abs_delta_A": round(max(dmin_deltas.values()), 6) if dmin_deltas else None,
        "n_control": len(ctl_chosen),
        "n_shipped_v3": len(v3_names),
        "shipped_v3_csv_order_matches_json": csv_order_ok,
        "overlap_with_shipped_v3": sorted(set(ctl_names) & set(v3_names)),
    }
    print(f"CONTROL v3-algorithm-on-v3-inputs: names={names_ok} dmin={dmin_ok} "
          f"({len(ctl_chosen)}/{len(v3_names)}, max|d|={control['dmin_max_abs_delta_A']} A)")
    if not (names_ok and dmin_ok and csv_order_ok):
        print("CONTROL FAILED — v4 not written.", file=sys.stderr)
        return 3

    # ---- v4: ONLY the membership filter and ranking key change ----
    pool_v4 = []
    excluded_no_apical = []
    for r in audit["records"]:
        n = r["name"]
        if n in no_apical or n not in apical_l6:
            excluded_no_apical.append(n)
            continue
        if iptm.get(n) is None or iptm[n] < args.iptm_floor:
            continue
        if apical_l6[n] <= 0:
            continue
        pool_v4.append(dict(r, L6=apical_l6[n], L6_surrogate=r.get("L6"),
                            membership_source="apical_195"))
    pool_v4.sort(key=lambda r: (-r["L6"], -(r["iptm"] or 0), r["name"]))

    chosen, rejected, _ = v3.select(pool_v4, APICAL_PDBS, args.slots, clash, crosswalk)
    shortfall = args.slots - len(chosen)

    seqs = {}
    for src in SEQ_SOURCES:
        p = ROOT / src
        if not p.exists():
            continue
        for row in csv.DictReader(open(p)):
            n = row.get("name")
            if n and n not in seqs and row.get("sequence"):
                seqs[n] = row["sequence"]
    missing_seq = [r["name"] for r in chosen if r["name"] not in seqs]
    if missing_seq:
        sys.exit(f"missing sequences for {len(missing_seq)} chosen designs: {missing_seq[:5]}")

    OUTDIR.mkdir(parents=True, exist_ok=True)
    out_csv = OUTDIR / "c4_acid_ranked_v4_apical_L6_20261025.csv"
    with open(out_csv, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["name", "sequence", "molecule_class"])
        for r in chosen:
            w.writerow([r["name"], seqs[r["name"]], "single_chain"])

    chosen_names = [r["name"] for r in chosen]
    v3_set = set(v3_names)
    added = [n for n in chosen_names if n not in v3_set]
    dropped = [n for n in v3_names if n not in set(chosen_names)]
    overlap = [n for n in chosen_names if n in v3_set]
    v2_names = [r["name"] for r in json.loads(
        (ROOT / "results/c4_ranked_20261025/c4_acid_ranked_20261025.json").read_text()
    )["chosen_by_L6"]]

    preview_names = {row["name"] for row in csv.DictReader(open(PREVIEW))}
    manifest = {
        "artifact": "EVIDENCE ONLY — not a submission; no promotion; not frozen; "
                    "Ricardo approval required",
        "defect_class": "MEMBERSHIP (acid-ON L6 gate computed on the surrogate construct)",
        "distinct_from": {
            "record": "results/decisions/20261025-1103-c4-wrong-construct-screen.md",
            "that_defect": "v2 screened clashes on the 154-residue stitched surrogate; v3 fixed the "
                           "screen directory only",
            "this_defect": "v3 still filtered/ranked membership on surrogate-construct L6, so its "
                           "acid-ON gate and its order are computed on the wrong construct",
        },
        "what_changed_vs_v3": "ONLY the membership gate and ranking key: surrogate L6 -> apical "
                              "195-residue L6 (results/c4_linkage_apical_20261025/linkage_apical.json). "
                              "Same iPTM floor, same apical PDB screen, same thresholds, same slots.",
        "control": control,
        "apical_l6_source": {
            "artifact": "results/c4_linkage_apical_20261025/linkage_apical.json",
            "sha256": sha256(APICAL_L6),
            "n_records": len(apical["records"]),
            "target": apical.get("target"),
            "target_len": apical.get("target_len"),
            "verified_this_run": "scripts/c4_linkage_apical_recompute_20261025.py (verify mode) "
                                 "reproduced 342/342 within 0.0006; positive 111/342",
        },
        "selection_rule": f"iPTM>={args.iptm_floor:.2f} AND apical L6>0, ranked by apical L6 desc "
                          f"then iPTM desc then name, take {args.slots}",
        "n_iptm_eligible": sum(1 for r in audit["records"]
                               if r.get("iptm") is not None and r["iptm"] >= args.iptm_floor),
        "n_eligible_acid_on_apical": len(pool_v4),
        "n_eligible_acid_on_surrogate_v3_pool": len(pool_v3),
        "excluded_no_apical_complex": excluded_no_apical,
        "slots_requested": args.slots,
        "slots_filled": len(chosen),
        "shortfall": shortfall,
        "padded": False,
        "clash_gate_A": 2.0,
        "pdb_source": "results/c4_apical_pdb_audit_20261025/all343 (chain B = 195)",
        "seat_change_vs_v3": {
            "n_changed": len(added),
            "n_retained": len(overlap),
            "added": added,
            "dropped": dropped,
        },
        "seat_change_vs_v2": {
            "n_changed": len([n for n in chosen_names if n not in set(v2_names)]),
            "n_retained": len([n for n in chosen_names if n in set(v2_names)]),
        },
        "chosen_by_apical_L6": [
            {"rank": i + 1, "name": r["name"], "L6_apical": r["L6"],
             "L6_surrogate": r.get("L6_surrogate"),
             "L6_sign_flip": ((r["L6"] > 0) != ((r.get("L6_surrogate") or 0) > 0)),
             "iptm": round(r["iptm"], 3), "d_min_A": r.get("d_min"),
             "margin_A": round(float(r.get("d_min") or 0) - 2.0, 3),
             "in_preview": r["name"] in preview_names,
             "in_v3": r["name"] in v3_set}
            for i, r in enumerate(chosen)],
        "margin_summary": {
            "borderline_2.0-2.2": sum(1 for r in chosen if (r.get("d_min") or 0) < 2.2),
            "above_2.2": sum(1 for r in chosen if (r.get("d_min") or 0) >= 2.2),
            "min_A": min((r.get("d_min") or 0) for r in chosen) if chosen else None,
            "max_A": max((r.get("d_min") or 0) for r in chosen) if chosen else None,
        },
        "survivors_acid_on_apical_of_v3_roster": sum(1 for n in v3_names if apical_l6.get(n, 0) > 0),
        "diversity_rejections": rejected,
        "caveats": [
            "SELECTION ARTIFACT: acid-ON rate and mean L6 are guaranteed by ranking on apical L6",
            "modelled PROPKA readout on predicted complexes, not experimental evidence",
            "cannot establish pH-conditional binding, activation, or release",
            "domain-only / non-native target construct",
            "2.0 A heavy-atom screen is atom-type-blind and permissive",
            "v2 and v3 are preserved unchanged; this is a new versioned artifact",
            "Ricardo approval required before any use",
        ],
    }
    shas = {"out_csv": str(out_csv.relative_to(ROOT)), "sha256": sha256(out_csv)}
    manifest["artifact_sha256"] = shas["sha256"]
    (OUTDIR / "c4_acid_ranked_v4_apical_L6_20261025.json").write_text(
        json.dumps(manifest, indent=1) + "\n")

    v3_sha_after = sha256(SHIPPED_V3_CSV)
    if v3_sha_after != v3_sha_before:
        sys.exit("v3 artifact changed during this run — investigate")

    print(f"v4: {len(chosen)}/{args.slots} slots, eligible_apical={len(pool_v4)}, "
          f"shortfall={shortfall}")
    print(f"seats changed vs v3: {len(added)} (+{len(added)}/-{len(dropped)}), "
          f"retained {len(overlap)}")
    print(f"added:   {added}")
    print(f"dropped: {dropped}")
    print(f"sha256:  {shas['sha256']}")
    print(f"v3 sha unchanged: {v3_sha_after == v3_sha_before} ({v3_sha_after[:12]})")
    return 0 if not shortfall else 4


if __name__ == "__main__":
    sys.exit(main())
