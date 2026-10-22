#!/usr/bin/env python
"""Rank and shortlist scored candidates into a submission-ready ordered list.

This is the selection policy, kept deliberately separate from scripts/build_submission.py
so that policy can change late without touching the format emitter. It answers three
questions the emitter must never contain:

  WHICH designs qualify   -- campaign gates are reported, and enforced only on request,
                             because the organizers state that no single in silico metric
                             decides and do not impose our internal thresholds.
  IN WHAT ORDER           -- a transparent weighted composite over normalised confidence
                             metrics, with every component dumped for audit. Row order is
                             the submitted rank, so this ordering is a claim.
  HOW MANY OF A KIND      -- diversity is a documented hard filter on their side, and our
                             pools are severely clustered: a single scaffold family can
                             fill all 20 slots. This is the one axis where the script is
                             STRICTER than the raw score, by design.

Honesty rule: a gate that cannot be evaluated from the input is reported as NOT
EVALUABLE, never as passed. The three mandatory termini gates
(N/C_Terminus_Away_Cosine, Termini_Away_Cosine) are absent from the scored candidate
CSVs, so designs drawn from those files are NOT certified against the Twin-Strep tag
accessibility constraint. That absence is printed loudly rather than rounded up.

  python scripts/rank_candidates.py results/v3_all/\\!_Refolded.csv \\
      --gates results/modal_native_replicates_v14/replicate_2/mirrored/campaign.json \\
      -o shortlist.csv
"""
import argparse
import csv
import json
import os
import re
import sys
from datetime import datetime, timezone

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from novelty_audit import smith_waterman  # reuse the audited local aligner

# Mandatory campaign gates routinely absent from scored candidate CSVs. Their absence is
# a reporting obligation, not a pass.
TERMINI_GATES = ("N_Terminus_Away_Cosine", "C_Terminus_Away_Cosine", "Termini_Away_Cosine")

# Visible, auditable ranking composite. Each component is min-max normalised across the
# candidate pool, then weighted; "higher False" means smaller raw values are better.
RANK_WEIGHTS = {
    "i_pTM":                         (0.35, True),
    "Unbound_Binder_pLDDT":          (0.20, True),
    "Hotspot_Contact_Fraction":      (0.15, True),
    "i_pAE":                         (0.15, False),
    "Interface_BuriedArea_Fraction": (0.10, True),
    "Off_Epitope_Contact_Fraction":  (0.05, False),
}


def fnum(row, key):
    try:
        return float(row[key])
    except (KeyError, TypeError, ValueError):
        return None


def load_candidates(paths):
    """Read BindCraft-style or Proteinbase-style scored CSVs into uniform records."""
    out = []
    for path in paths:
        with open(path) as fh:
            rows = list(csv.DictReader(fh))
        if not rows:
            print("note: %s has no rows, skipped" % path)
            continue
        cols = rows[0].keys()
        scol = next((c for c in ("Binder_Sequence", "sequence", "Sequence") if c in cols), None)
        ncol = next((c for c in ("design", "name", "Design") if c in cols), None)
        if scol is None:
            sys.exit("error: %s has no sequence column" % path)
        for i, r in enumerate(rows):
            seq = (r.get(scol) or "").strip()
            if seq:
                out.append({"source_file": os.path.relpath(path),
                            "source_id": (r.get(ncol) or "row%d" % i).strip(),
                            "sequence": seq, "metrics": r})
    return out


def apply_gates(rec, gates):
    passed, failed, uncheckable = [], [], []
    for key, spec in sorted(gates.items()):
        val = fnum(rec["metrics"], key)
        if val is None:
            uncheckable.append(key)
            continue
        thr = float(spec["threshold"])
        good = val >= thr if spec.get("higher") else val <= thr
        (passed if good else failed).append(key)
    return passed, failed, uncheckable


def rank_scores(records):
    pools = {}
    for key in RANK_WEIGHTS:
        vals = [v for v in (fnum(r["metrics"], key) for r in records) if v is not None]
        pools[key] = (min(vals), max(vals)) if vals else None
    for rec in records:
        total, parts = 0.0, {}
        for key, (weight, higher) in RANK_WEIGHTS.items():
            span, val = pools[key], fnum(rec["metrics"], key)
            if span is None or val is None:
                parts[key] = None
                continue
            lo, hi = span
            norm = 0.5 if hi == lo else (val - lo) / (hi - lo)
            if not higher:
                norm = 1.0 - norm
            parts[key] = round(norm * weight, 6)
            total += norm * weight
        rec["rank_score"] = round(total, 6)
        rec["rank_parts"] = parts


def identity(a, b, min_aln):
    """Best local identity fraction, ignoring alignments too short to mean anything."""
    _, ids, aln = smith_waterman(a, b)
    return 0.0 if aln < min_aln else ids / float(aln)


def cluster(records, threshold, min_aln):
    """Greedy single-linkage clustering by local identity; best-ranked design seeds each."""
    clusters = []
    for rec in records:
        for idx, members in enumerate(clusters):
            if any(identity(rec["sequence"], m["sequence"], min_aln) >= threshold
                   for m in members):
                members.append(rec)
                rec["cluster"] = idx
                break
        else:
            rec["cluster"] = len(clusters)
            clusters.append([rec])
    return clusters


def slugify(text, fallback):
    s = re.sub(r"[^A-Za-z0-9_-]+", "_", (text or "").strip()).strip("_-")
    return s[:64] or fallback


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("candidates", nargs="+")
    ap.add_argument("-o", "--out", required=True, help="ranked shortlist CSV to write")
    ap.add_argument("--max-rows", type=int, default=20, help="shortlist length (Track 3 cap)")
    ap.add_argument("--molecule-class", default="protein")
    ap.add_argument("--gates", default=None, help="campaign.json whose filters classify rows")
    ap.add_argument("--require-gates", action="store_true",
                    help="drop rows failing any evaluable gate (default: keep and report)")
    ap.add_argument("--max-per-cluster", type=int, default=3,
                    help="diversity cap per near-identical family (0 = unlimited)")
    ap.add_argument("--cluster-identity", type=float, default=0.60)
    ap.add_argument("--min-aln", type=int, default=20)
    ap.add_argument("--name-prefix", default="c4_tfr1")
    ap.add_argument("--keep-source-names", action="store_true")
    a = ap.parse_args()

    records = load_candidates(a.candidates)
    print("loaded %d candidate rows from %d file(s)" % (len(records), len(a.candidates)))
    if not records:
        sys.exit("error: no candidates loaded")

    seen, unique = {}, []
    for rec in records:
        if rec["sequence"] in seen:
            continue
        seen[rec["sequence"]] = rec
        unique.append(rec)
    if len(unique) != len(records):
        print("deduplicated %d identical sequence(s)" % (len(records) - len(unique)))

    gates, uncheckable_global, mandatory_missing = {}, [], []
    if a.gates:
        gates = json.load(open(a.gates)).get("filters", {})
        print("loaded %d campaign gates from %s" % (len(gates), a.gates))
    for rec in unique:
        rec["gates_passed"], rec["gates_failed"], rec["gates_uncheckable"] = (
            apply_gates(rec, gates) if gates else ([], [], []))
    if gates:
        uncheckable_global = sorted({g for r in unique for g in r["gates_uncheckable"]})
        mandatory_missing = [g for g in uncheckable_global
                            if gates.get(g, {}).get("mandatory") or g in TERMINI_GATES]
        if mandatory_missing:
            print("")
            print("!! NOT EVALUABLE from this input, and MANDATORY in the campaign filters:")
            for g in mandatory_missing:
                print("     %s" % g)
            print("   Rows from this input are NOT certified against those gates. The termini")
            print("   gates encode the Twin-Strep tag accessibility constraint; absence of the")
            print("   column is not a pass.")
            print("")
        if a.require_gates:
            before = len(unique)
            unique = [r for r in unique if not r["gates_failed"]]
            print("--require-gates: %d of %d rows pass every evaluable gate"
                  % (len(unique), before))
        else:
            clean = sum(1 for r in unique if not r["gates_failed"])
            print("gate report: %d of %d rows pass every evaluable gate (not enforced)"
                  % (clean, len(unique)))

    if not unique:
        sys.exit("error: no rows survived filtering")

    rank_scores(unique)
    unique.sort(key=lambda r: (-r["rank_score"], r["source_id"]))
    clusters = cluster(unique, a.cluster_identity, a.min_aln)
    print("ranked %d rows into %d sequence famil%s at >=%.0f%% local identity"
          % (len(unique), len(clusters), "y" if len(clusters) == 1 else "ies",
             a.cluster_identity * 100))

    selected, per_cluster, deferred = [], {}, []
    for rec in unique:
        if len(selected) >= a.max_rows:
            deferred.append(rec)
            continue
        used = per_cluster.get(rec["cluster"], 0)
        if a.max_per_cluster and used >= a.max_per_cluster:
            rec["deferred_reason"] = ("diversity cap: family %d already contributes %d row(s)"
                                      % (rec["cluster"], used))
            deferred.append(rec)
            continue
        per_cluster[rec["cluster"]] = used + 1
        selected.append(rec)

    short = min(a.max_rows, len(unique)) - len(selected)
    if short > 0:
        print("")
        print("!! the diversity cap held back %d row(s) of a %d-row package: this pool is too"
              % (short, a.max_rows))
        print("   clustered to fill the track legitimately. %d famil%s for %d slots."
              % (len(clusters), "y" if len(clusters) == 1 else "ies", a.max_rows))
        print("")

    rows = []
    for i, rec in enumerate(selected, 1):
        name = (slugify(rec["source_id"], "design_%02d" % i) if a.keep_source_names
                else "%s_%02d" % (a.name_prefix, i))
        rec["submission_name"], rec["submission_rank"] = name, i
        rows.append({"name": name, "sequence": rec["sequence"],
                     "molecule_class": a.molecule_class,
                     "rank": i, "cluster": rec["cluster"],
                     "rank_score": rec["rank_score"],
                     "source_id": rec["source_id"]})

    outdir = os.path.dirname(os.path.abspath(a.out))
    if outdir:
        os.makedirs(outdir, exist_ok=True)
    cols = ["name", "sequence", "molecule_class", "rank", "cluster", "rank_score", "source_id"]
    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, lineterminator="\n")
        w.writeheader()
        w.writerows(rows)

    prov = os.path.splitext(a.out)[0] + ".provenance.json"
    with open(prov, "w") as fh:
        json.dump({
            "generated_utc": datetime.now(timezone.utc).isoformat(),
            "generator": "scripts/rank_candidates.py",
            "inputs": [os.path.relpath(p) for p in a.candidates],
            "gates_file": a.gates,
            "gates_enforced": bool(a.gates and a.require_gates),
            "gates_not_evaluable": uncheckable_global,
            "mandatory_gates_not_evaluable": mandatory_missing,
            "ranking_weights": {k: {"weight": v[0], "higher_is_better": v[1]}
                                for k, v in RANK_WEIGHTS.items()},
            "diversity": {"cluster_identity": a.cluster_identity,
                          "max_per_cluster": a.max_per_cluster,
                          "min_alignment": a.min_aln,
                          "families_found": len(clusters),
                          "slots_held_back_by_cap": max(0, short)},
            "counts": {"loaded": len(records), "unique": len(unique),
                       "shortlisted": len(rows), "deferred": len(deferred)},
            "shortlist": [{
                "rank": r["submission_rank"], "name": r["submission_name"],
                "sequence": r["sequence"], "length": len(r["sequence"]),
                "source_file": r["source_file"], "source_id": r["source_id"],
                "cluster": r["cluster"], "rank_score": r["rank_score"],
                "rank_components": r["rank_parts"],
                "gates_passed": r["gates_passed"], "gates_failed": r["gates_failed"],
                "gates_not_evaluable": r["gates_uncheckable"],
                "metrics": {k: r["metrics"].get(k) for k in RANK_WEIGHTS},
            } for r in selected],
            "deferred": [{"source_id": r["source_id"], "cluster": r.get("cluster"),
                          "rank_score": r.get("rank_score"),
                          "reason": r.get("deferred_reason", "below the row cap")}
                         for r in deferred[:50]],
        }, fh, indent=2)

    print("wrote %s  (%d row%s)" % (a.out, len(rows), "" if len(rows) == 1 else "s"))
    print("wrote %s" % prov)
    print("")
    print("rank  name                   len  family  score  i_pTM  pLDDT  gates_failed")
    for r in selected:
        m = r["metrics"]
        print("%4d  %-20s  %4d  %6d  %5.3f  %5s  %5s  %d" % (
            r["submission_rank"], r["submission_name"], len(r["sequence"]), r["cluster"],
            r["rank_score"], (m.get("i_pTM") or "-")[:5],
            (m.get("Unbound_Binder_pLDDT") or "-")[:5], len(r["gates_failed"])))
    print("")
    print("next: python scripts/build_submission.py %s -o submission/round4.csv" % a.out)


if __name__ == "__main__":
    main()
