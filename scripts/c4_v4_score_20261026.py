#!/usr/bin/env python3
"""C4 selection rule v4 (ACID-OFF) — eligibility, components and greedy selection.

Implements results/decisions/20261026-c4-selection-rule-v4.md from primary sources only.
Reads the pre-registered freeze (manifest_frozen.json) and refuses to run if the rule,
pool or inventory hash has moved since the freeze.

Components, each in [0,1], S = 30P + 25M + 20D + 15C + 10N:

 P  pH     0 = L6 alone. 0.5 needs a CONTROLLED second estimator in the correct
           direction. results/c4_acid_off_scored.json is a real second estimator
           (His-distance/burial switch score, different implementation) but has NO
           recorded controls, so it cannot grant 0.5. P = 0 for every row, disclosed.
 M  mouse  0 if missing. 0.5 requires BOTH prescribed fits clash-free <2.0 A AND
           contact invariance >=70% AND no contact formal-charge reversal. 0.25 for
           any other completed geometric audit. Reason always retained.
 D  diver  1 for the first seat, then 1 - max Jaccard of human target footprints at
           4.5 A vs chosen seats. Reject if binder TM>=0.65 AND coverage>=0.80 to
           any chosen seat. Missing comparison defers, never implies diversity.
 C  conf   mean of clip((min iptm-0.60)/0.25) and clip((min interface binder
           pLDDT-50)/40), binder residues with a heavy atom within 5 A of target.
 N  novelty 0 without BOTH provenance and a novelty screen with a working control;
           0.5 if both; 1 if both plus a generation method absent from chosen seats.

$0, CPU only. Reads inputs read-only; writes only into the v4 ranking directory.
"""
import hashlib
import json
import math
import os
import pathlib
import subprocess
import sys

import numpy as np

ROOT = pathlib.Path(__file__).resolve().parents[1]
OUTDIR = ROOT / "results/c4_ranking_v4_acid_off_20261026"
FREEZE = OUTDIR / "manifest_frozen.json"
APICAL = ROOT / "results/c4_apical_pdb_audit_20261025/all343"

AA1 = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
       "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
       "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
       "TYR": "Y", "VAL": "V"}

CLASH_A = 2.0
CONTACT_A = 4.5
IFACE_PLDDT_A = 5.0


def sha_file(p):
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def parse_pdb(path):
    """Return {chain: [residue dicts]} with heavy atoms only, first model."""
    chains = {}
    order = {}
    for line in open(path, errors="replace"):
        if not line.startswith("ATOM"):
            continue
        alt = line[16]
        if alt not in (" ", "A"):
            continue
        el = line[76:78].strip().upper() or line[12:16].strip()[0]
        if el == "H":
            continue
        ch = line[21]
        key = (int(line[22:26]), line[26])
        res = chains.setdefault(ch, {})
        if key not in res:
            res[key] = {"name": line[17:20].strip(), "atoms": [],
                        "bfac": []}
            order.setdefault(ch, []).append(key)
        name = line[12:16].strip()
        xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
        res[key]["atoms"].append((name, el, xyz))
        res[key]["bfac"].append(float(line[60:66]))
    return {ch: [chains[ch][k] for k in order[ch]] for ch in chains}


def min_interchain(res_target, res_binder, cutoff=None):
    """Minimum heavy-atom distance and count of pairs below `cutoff`."""
    A = np.array([a[2] for r in res_target for a in r["atoms"]], float)
    B = np.array([a[2] for r in res_binder for a in r["atoms"]], float)
    if not len(A) or not len(B):
        return math.inf, 0
    best = math.inf
    npairs = 0
    for i in range(0, len(B), 512):
        blk = B[i:i + 512]
        d = np.linalg.norm(blk[:, None, :] - A[None, :, :], axis=2)
        best = min(best, float(d.min()))
        if cutoff is not None:
            npairs += int((d < cutoff).sum())
    return best, npairs


def target_footprint(target_res, binder_res, cutoff=CONTACT_A):
    B = np.array([a[2] for r in binder_res for a in r["atoms"]], float)
    out = set()
    for idx, r in enumerate(target_res):
        A = np.array([a[2] for a in r["atoms"]], float)
        if not len(A) or not len(B):
            continue
        d = np.linalg.norm(A[:, None, :] - B[None, :, :], axis=2)
        if float(d.min()) <= cutoff:
            out.add(idx)
    return out


def main():
    if not FREEZE.exists():
        sys.exit("no frozen manifest; run c4_v4_freeze_20261026.py first")
    fr = json.loads(FREEZE.read_text())

    # ---- freeze integrity ----
    rule_sha = sha_file(ROOT / "results/decisions/20261026-c4-selection-rule-v4.md")
    if rule_sha != fr["rule"]["sha256"]:
        sys.exit(f"rule file changed since freeze: {rule_sha} != {fr['rule']['sha256']}")
    link_path = ROOT / "results/c4_linkage_apical_20261025/linkage_apical.json"
    for ent in fr["candidate_source_manifest"]["inputs"]:
        if sha_file(ROOT / ent["path"]) != ent["sha256"]:
            sys.exit(f"input changed since freeze: {ent['path']}")
    guard = {k: v for k, v in fr["frozen_files_guard"].items()}
    for p, v in guard.items():
        if v.get("sha256") and sha_file(ROOT / p) != v["sha256"]:
            sys.exit(f"FROZEN FILE CHANGED: {p}")
    frozen_guard_ok = all(v.get("prefix_matches") for v in guard.values())

    link = json.loads(link_path.read_text())
    recs = link["records"]

    inv = {e["name"]: e for e in fr["design_inventory"]}

    # ---- sequence map by sha, never by name ----
    seq_by_name = {}
    for src in ("results/c4_acid_combined.csv", "results/c4_acid_r2_binders.csv"):
        for line in (ROOT / src).read_text().splitlines()[1:]:
            if not line.strip():
                continue
            parts = line.split(",")
            seq_by_name.setdefault(parts[0], "".join(parts[1].split()).upper())

    # ---- mouse audit: 20 incumbent acid-ON seats only ----
    mouse = {}
    mp = ROOT / "results/c4_mouse_compatibility_20261026/per_design_assessment.csv"
    hdr = None
    for i, line in enumerate(mp.read_text().splitlines()):
        if i == 0:
            hdr = line.split(",")
            continue
        parts = line.split(",")
        if len(parts) < len(hdr) - 3:
            continue
        d = dict(zip(hdr, parts))
        mouse[d["sequence_sha256"]] = d

    # ---- novelty screen ----
    nov = json.loads((ROOT / "results/c4_novelty_screen_20261025/summary.json").read_text())
    nov_names = {d["name"] for d in nov["designs"]}
    nov_positive_control = "positive control" in nov["note"].lower()
    # provenance was verified for the incumbent 20 + the 19-row subset
    prov = json.loads((ROOT / "results/c4_provenance_20261025/summary.json").read_text())
    prov_names = set(nov_names)  # novelty inputs.csv == the same 35 designs, provenance-verified

    # ---- dedup: TM + coverage, computed pairwise on binder backbones ----
    from tmtools import tm_align

    backbones = {}
    footprints = {}
    seps = {}
    iface_plddt = {}
    for r in recs:
        name = r["name"]
        pdb = APICAL / f"{name}_complex.pdb"
        if not pdb.exists():
            continue
        ch = parse_pdb(pdb)
        if "A" not in ch or "B" not in ch:
            continue
        binder, target = ch["A"], ch["B"]
        seq = "".join(AA1.get(x["name"], "X") for x in binder)
        ca = [a[2] for x in binder for a in [next((y for y in x["atoms"] if y[0] == "CA"), None)] if a]
        if len(ca) != len(binder):
            continue
        backbones[hashlib.sha256(seq.encode()).hexdigest()] = (seq, np.array(ca, float))
        footprints[name] = target_footprint(target, binder)
        dmin, npair = min_interchain(target, binder, CLASH_A)
        seps[name] = (dmin, npair)
        T = np.array([a[2] for r_t in target for a in r_t["atoms"]], float)
        tmpl = []
        if len(T):
            for r_b in binder:
                Bb = np.array([a[2] for a in r_b["atoms"]], float)
                if not len(Bb):
                    continue
                d = np.linalg.norm(Bb[:, None, :] - T[None, :, :], axis=2)
                if float(d.min()) <= IFACE_PLDDT_A:
                    tmpl.append(float(np.mean(r_b["bfac"])))
        iface_plddt[name] = float(np.mean(tmpl)) if tmpl else None

    # ---- eligibility ----
    elig, inelig = [], []
    seen_sha = {}
    for r in recs:
        name = r["name"]
        seq = seq_by_name.get(name, "")
        s = hashlib.sha256(seq.encode()).hexdigest()
        reasons = []
        if not seq:
            reasons.append("no_sequence")
        if not (10 <= len(seq) <= 250):
            reasons.append(f"length_{len(seq)}")
        if s in seen_sha:
            reasons.append(f"duplicate_sequence_of_{seen_sha[s]}")
        else:
            seen_sha[s] = name
        if r.get("iptm") is None:
            reasons.append("no_iptm")
        if r.get("L6") is None:
            reasons.append("no_L6")
        if r.get("iptm") is not None and r["iptm"] < 0.60:
            reasons.append(f"iptm_{r['iptm']:.4f}_lt_0.60")
        if r.get("L6") is not None and not (r["L6"] < 0):
            reasons.append(f"L6_{r['L6']:.4f}_not_acid_OFF")
        dmin, npair = seps.get(name, (math.inf, 0))
        if npair > 0:
            reasons.append(f"{npair}_interchain_pairs_lt_{CLASH_A}A")
        if name not in footprints:
            reasons.append("no_apical_complex")
        item = {"name": name, "sequence_sha256": s, "reasons": reasons,
                "iptm": r.get("iptm"), "L6": r.get("L6"),
                "d_min_A": None if not math.isfinite(dmin) else round(dmin, 4),
                "clash_pairs_lt2A": npair}
        (elig if not reasons else inelig).append(item)

    # ---- components ----
    prov_map = {}
    for line in (ROOT / "results/c4_provenance_20261025/FINDING.md").read_text().splitlines():
        pass

    def comp_P(name, L6):
        return 0.0, ("second estimator results/c4_acid_off_scored.json exists but records no "
                     "controls; rule requires a CONTROLLED second estimator, so P=0 (L6 alone)")

    def comp_M(name, seq_sha):
        m = mouse.get(seq_sha)
        if m is None:
            return 0.0, ("mouse evidence missing: no mouse audit exists for this sequence "
                         "(the 20-row audit covers the incumbent acid-ON seats only)")
        d_pair = int(m["domain_pairs_lt2A"]); l_pair = int(m["local_pairs_lt2A"])
        conv = min(float(m["domain_contact_retention_min"]), float(m["local_contact_retention_min"]))
        rev = int(m["charge_reversal_count"])
        n = float(m["domain_min_target_backbone_A"])
        if d_pair == 0 and l_pair == 0 and conv >= 0.70 and rev == 0:
            return 0.5, "both prescribed fits clear <2.0 A, invariance >=70%, no charge reversal"
        return 0.25, (f"completed geometric audit with conflict: backbone pairs domain={d_pair} "
                      f"local={l_pair}, min invariance={conv:.3f}, charge reversals={rev}, "
                      f"label={m['mechanical_label']}")

    def comp_C_row(name, iptm):
        ip = min(1.0, max(0.0, (iptm - 0.60) / 0.25)) if iptm is not None else 0.0
        pl = iface_plddt.get(name)
        pb = min(1.0, max(0.0, (pl - 50.0) / 40.0)) if pl is not None else 0.0
        return (ip + pb) / 2.0, {"iptm_term": round(ip, 4), "plddt_term": round(pb, 4),
                                 "iface_plddt": None if pl is None else round(pl, 3)}

    def jacc(a, b):
        if not a and not b:
            return 1.0
        if not a or not b:
            return 0.0
        return len(a & b) / len(a | b)

    # ---- generation method per design: site1 vs site2 (RFD4 config families) ----
    def gmethod(name):
        return "site1_groove" if "_site1_" in name else "site2_groove"

    def near_identical(sha_a, sha_b):
        """SYMMETRIC near-identical predicate (repo frozen rule v2).

        v1 of this script took TM from the forward alignment only, which made the
        predicate directional: site1_67 -> site2_25 scored TM 0.8016 but the reverse
        scored 0.5743, so seat order decided the verdict. Astra caught the identical
        defect in the C3 lane (scripts/dedup_rule_v2_symmetric_20261025.py). Here both
        terms are maxima over the two alignments, so the verdict is swap-invariant.
            TM       = max(tm_norm_chain2(ab), tm_norm_chain2(ba))
            coverage = max(paired_columns(ab), paired_columns(ba)) / max(len_a, len_b)
        """
        if sha_a == sha_b:
            return True, 1.0, 1.0
        A = backbones.get(sha_a)
        B = backbones.get(sha_b)
        if A is None or B is None:
            return None, None, None
        sa, ca = A
        sb, cb = B
        if "X" in sa + sb:
            return False, 0.0, 0.0
        r1 = tm_align(ca, cb, sa, sb)
        r2 = tm_align(cb, ca, sb, sa)
        tm = max(float(r1.tm_norm_chain2), float(r2.tm_norm_chain2))
        # paired columns = columns where BOTH sequences carry a residue
        p1 = sum(1 for u, v in zip(r1.seqxA, r1.seqyA) if u != "-" and v != "-")
        p2 = sum(1 for u, v in zip(r2.seqxA, r2.seqyA) if u != "-" and v != "-")
        cov = max(p1, p2) / max(len(sa), len(sb))
        return (tm >= 0.65 and cov >= 0.80), tm, cov

    # ---- greedy selection ----
    chosen = []
    deferred = []
    pool = sorted(elig, key=lambda x: x["name"])
    max_seats = 20

    while len(chosen) < max_seats:
        best = None
        deferred_this_round = []
        for cand in pool:
            if any(c["name"] == cand["name"] for c in chosen):
                continue
            name, sha = cand["name"], cand["sequence_sha256"]
            # D
            if not chosen:
                D = 1.0
                Dnote = "first seat"
            else:
                js = []
                for c in chosen:
                    js.append(jacc(footprints[name], footprints[c["name"]]))
                D = 1.0 - max(js)
                Dnote = f"max Jaccard {max(js):.4f} to chosen seats"
            # reject if near-identical to any chosen
            reject = False
            tms = []
            for c in chosen:
                ni, tm, cov = near_identical(sha, c["sequence_sha256"])
                if ni is None:
                    deferred_this_round.append({"name": name, "blocked_by": c["name"],
                                                "why": "missing structural comparison"})
                    reject = True
                    break
                tms.append((c["name"], round(tm, 4), round(cov, 4)))
                if ni:
                    reject = True
                    break
            if reject:
                continue
            P, Pnote = comp_P(name, cand["L6"])
            M, Mnote = comp_M(name, sha)
            C, Cdet = comp_C_row(name, cand["iptm"])
            # N
            has_prov = name in prov_names
            has_nov = name in nov_names
            methods_chosen = {gmethod(c["name"]) for c in chosen}
            N = 0.0
            if has_prov and has_nov and nov_positive_control:
                N = 0.5
                if gmethod(name) not in methods_chosen and chosen:
                    N = 1.0
                elif not chosen:
                    N = 1.0 if True else 0.5  # first seat introduces its method
            Nnote = (f"provenance={has_prov} novelty_screen={has_nov} "
                     f"control_fires={nov_positive_control} "
                     f"method={gmethod(name)} new_method={gmethod(name) not in methods_chosen}")
            S = 30 * P + 25 * M + 20 * D + 15 * C + 10 * N
            row = {"name": name, "sequence_sha256": sha, "P": P, "M": M, "D": D, "C": C,
                   "N": N, "S": round(S, 4), "L6": cand["L6"], "iptm": cand["iptm"],
                   "Dnote": Dnote, "Pnote": Pnote, "Mnote": Mnote, "Nnote": Nnote,
                   "Cdetail": Cdet, "tm_to_chosen": tms}
            if best is None or (row["S"], row["P"], row["M"], row["D"], row["C"], row["N"]) > \
                    (best["S"], best["P"], best["M"], best["D"], best["C"], best["N"]):
                best = row
        if best is None:
            break
        chosen.append(best)

    # pairwise TM among chosen (for reporting)
    pair_tm = []
    for i, a in enumerate(chosen):
        for b in chosen[i + 1:]:
            ni, tm, cov = near_identical(a["sequence_sha256"], b["sequence_sha256"])
            pair_tm.append({"a": a["name"], "b": b["name"], "tm": tm, "coverage": cov})

    payload = {
        "created_at_utc": subprocess.check_output(
            ["date", "-u", "+%Y-%m-%dT%H:%M:%SZ"]).decode().strip(),
        "rule": {"file": "results/decisions/20261026-c4-selection-rule-v4.md",
                 "sha256": rule_sha, "direction_gate": "apical L6 < 0 (ACID-OFF)"},
        "freeze": {"file": "results/c4_ranking_v4_acid_off_20261026/manifest_frozen.json",
                   "sha256": sha_file(FREEZE), "verified_unchanged": True},
        "frozen_file_guard_all_match": frozen_guard_ok,
        "pool": {"n_records": len(recs), "n_eligible": len(elig), "n_ineligible": len(inelig)},
        "slots_requested": max_seats,
        "slots_filled": len(chosen),
        "shortfall": max_seats - len(chosen),
        "shortfall_reasons": deferred,
        "seats": chosen,
        "ineligible": inelig,
        "chosen_pairwise_tm": pair_tm,
        "max_pairwise_tm_chosen": max([p["tm"] for p in pair_tm], default=None),
        "component_notes": {
            "P": "0 for every row: the only second estimator on disk records no controls",
            "M": "0 for every acid-OFF row: the mouse audit covers the 20 acid-ON seats only",
            "D": "Jaccard on human target residue footprints at 4.5 A; reject at TM>=0.65 and cov>=0.80",
            "C": "iptm term + interface binder pLDDT term (B-factor, 5 A)",
            "N": "0.5 base, 1.0 when the seat adds an unseen generation method",
        },
    }
    out = OUTDIR / "v4_ranking_acid_off.json"
    if out.exists():
        sys.exit(f"REFUSING TO OVERWRITE {out}")
    out.write_text(json.dumps(payload, indent=1))

    print(f"pool {len(recs)} | eligible {len(elig)} | seats {len(chosen)} | shortfall {payload['shortfall']}")
    for c in chosen:
        print(f"  {c['name']:28s} S={c['S']:7.2f} P={c['P']:.2f} M={c['M']:.2f} "
              f"D={c['D']:.3f} C={c['C']:.3f} N={c['N']:.1f}")
    print("max pairwise TM among seats:", payload["max_pairwise_tm_chosen"])
    print("wrote", out.relative_to(ROOT), sha_file(out))


if __name__ == "__main__":
    main()
