#!/usr/bin/env python
"""Emit the Proteinbase submission CSV from an already-ranked shortlist.

This script does ONE job: take a ranked list of (name, sequence, molecule_class) and
write a file that conforms to the competition's format contract -- or refuse to write
anything at all.

It deliberately contains NO qualification logic. No campaign gates, no confidence
ranking, no diversity clustering. Those decisions live upstream in
scripts/rank_candidates.py, for a concrete reason: the native-context clash gate that
rejected 66 candidates is itself unverified, and if the 6W3H control falsifies it the
eligible pool widens again. Selection policy must be free to change late; the emitter
must not. Keeping them apart means a reversed gate decision never touches this file.

What it DOES enforce is the format contract itself -- column names, the length window,
the residue alphabet, the row cap, name/sequence uniqueness. That is not qualification;
that is the difference between a submission and a disqualified upload, and refusing to
emit a non-conformant file is precisely this script's reason to exist.

Input row order IS the submitted rank. FAQ 5 asks for the CSV "ordered by how you would
rank your molecules (top row higher)", so order is meaning, not incidental. This script
preserves input order exactly and never re-sorts.

  python scripts/rank_candidates.py results/v3_all/\\!_Refolded.csv -o shortlist.csv
  python scripts/build_submission.py shortlist.csv -o submission/round4.csv
  python scripts/novelty_audit.py submission/round4.csv --max-identity 0.30
"""
import argparse
import csv
import json
import os
import re
import sys
from datetime import datetime, timezone

# --------------------------------------------------------------------------------------
# OFFICIAL FORMAT CONTRACT
#
# Confidence tags are load-bearing. The second-opinion pass that produced this spec read
# the competition pages through a delegating subagent that capped quotes at 125 chars,
# never retrieved the selection page's prompt body, and got only navigation from the
# submit page. So:
#
#   [DOC] reported as page content (still second-hand -- re-fetch before freezing)
#   [INF] inference, defensible but not stated anywhere
#   [UNK] not present in anything read; MUST NOT be encoded as a rule
#
# Re-verify [DOC] directly, and resolve every [UNK] from the AUTHENTICATED submit page
# once Ricardo's Proteinbase account exists. Four unknowns collapse in that one visit:
# upload mechanics, file naming, whether a submission can be replaced, and the
# Twin-Strep convention.
# --------------------------------------------------------------------------------------
REQUIRED_COLUMNS = ("name", "sequence", "molecule_class")  # [DOC] FAQ 5: "at least" these
EXTRA_COLUMNS_ALLOWED = True        # [DOC] "at least" implies extras are permitted
HEADER_ROW_REQUIRED = True          # [INF] columns referred to by name; near-certain
LENGTH_MIN, LENGTH_MAX = 10, 250    # [DOC] challenge page hard requirement
TRACK_ROW_CAP = {"1": 40, "2": 20, "3": 20}  # [DOC] 40 for Track 1, 20 for Tracks 2/3
ROW_ORDER_IS_RANK = True            # [DOC] FAQ 5: ordered by rank, top row higher

# [DOC] with a paraphrase caveat -- medium-high confidence, not verbatim-certain.
# One named constant so re-verification is a one-line edit.
MOLECULE_CLASSES = ("protein", "nanobody", "scfv", "fab_kappa", "fab_lambda")
MOLECULE_CLASS_DEFAULT = "protein"  # our Track 3 designs are de novo single chains

# [DOC] Fab rows carry VH and VL joined by ":". Any other class must not contain one --
# a stray colon is a malformed row, never something to silently strip.
CHAIN_SEPARATOR = ":"
MULTI_CHAIN_CLASSES = ("fab_kappa", "fab_lambda")

AA20 = set("ACDEFGHIKLMNPQRSTVWY")   # [INF] nothing documented on X/lowercase/gaps

# [UNK] -- entirely absent from every page read. Self-imposed so names stay portable and
# collision-free; declared as OUR convention in the methods write-up, not as their rule.
NAME_PATTERN = re.compile(r"^[A-Za-z0-9_-]{1,64}$")

# [UNK] and the highest-consequence unknown in the format. Designs are captured through a
# C-terminal Twin-Strep tag that "must not be part of the binding" [DOC], but whether the
# tag belongs IN the submitted sequence is stated nowhere. Default excluded ([INF]: it is
# Adaptyv's expression-time capture handle and would eat ~30 residues of the length
# budget), recorded explicitly in the manifest, and confirmable by one flag rather than a
# repo-wide find-and-replace at the deadline.
TWIN_STREP_DEFAULT = "excluded"


def validate_row(idx, row, track_cap, seen_names, seen_seqs, strict_class):
    """Return a list of format violations for one row. Empty list means conformant."""
    bad = []
    name = (row.get("name") or "").strip()
    seq = (row.get("sequence") or "").strip()
    cls = (row.get("molecule_class") or "").strip()
    where = "row %d" % idx

    if not name:
        bad.append("%s: empty name" % where)
    elif not NAME_PATTERN.match(name):
        bad.append("%s: name %r violates our [A-Za-z0-9_-]{1,64} convention" % (where, name))
    elif name in seen_names:
        bad.append("%s: duplicate name %r (also row %d)" % (where, name, seen_names[name]))

    if not cls:
        bad.append("%s: empty molecule_class" % where)
    elif cls not in MOLECULE_CLASSES:
        msg = "%s: molecule_class %r not in %s" % (where, cls, "/".join(MOLECULE_CLASSES))
        bad.append(msg) if strict_class else print("warning: %s (spelling is [DOC] but "
                                                  "paraphrased -- verify)" % msg)

    if not seq:
        bad.append("%s: empty sequence" % where)
        return bad

    if CHAIN_SEPARATOR in seq:
        if cls in MULTI_CHAIN_CLASSES:
            chains = seq.split(CHAIN_SEPARATOR)
            if len(chains) != 2 or not all(chains):
                bad.append("%s: %s expects exactly VH%sVL" % (where, cls, CHAIN_SEPARATOR))
            seq_parts = chains
        else:
            bad.append("%s: %r contains %r, which only %s may use -- a colon here is a "
                       "malformed row, not something to strip"
                       % (where, name, CHAIN_SEPARATOR, "/".join(MULTI_CHAIN_CLASSES)))
            seq_parts = [seq.replace(CHAIN_SEPARATOR, "")]
    else:
        seq_parts = [seq]

    for part in seq_parts:
        if not (LENGTH_MIN <= len(part) <= LENGTH_MAX):
            bad.append("%s: length %d outside the documented %d-%d window"
                       % (where, len(part), LENGTH_MIN, LENGTH_MAX))
        illegal = sorted(set(part) - AA20)
        if illegal:
            lower = [c for c in illegal if c.islower()]
            hint = " (lowercase not accepted)" if lower else ""
            bad.append("%s: non-standard residues %s%s" % (where, ",".join(illegal), hint))

    if seq in seen_seqs:
        bad.append("%s: duplicate sequence, identical to row %d -- uniqueness is a "
                   "documented hard filter" % (where, seen_seqs[seq]))
    return bad


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("shortlist", help="ranked CSV with at least name,sequence,molecule_class")
    ap.add_argument("-o", "--out", required=True, help="submission CSV to write")
    ap.add_argument("--track", default="3", choices=sorted(TRACK_ROW_CAP))
    ap.add_argument("--molecule-class", default=None,
                    help="fill/override molecule_class for every row (default: take the "
                         "column, or %r when absent)" % MOLECULE_CLASS_DEFAULT)
    ap.add_argument("--twin-strep", default=TWIN_STREP_DEFAULT,
                    choices=("excluded", "included"),
                    help="whether submitted sequences already carry the C-terminal "
                         "Twin-Strep tag. UNVERIFIED against the rules; recorded in the "
                         "manifest so it is one flag to change, not a find-and-replace")
    ap.add_argument("--extra-columns", default="",
                    help="comma-separated extra columns to carry through (allowed, and an "
                         "asset: shared info reaches the selector for Tracks 2/3)")
    ap.add_argument("--strict-class", action="store_true",
                    help="treat an unrecognised molecule_class as fatal rather than a warning")
    ap.add_argument("--format-test", action="store_true",
                    help="mark the output as a FORMAT TEST, not a submission candidate")
    a = ap.parse_args()

    cap = TRACK_ROW_CAP[a.track]

    with open(a.shortlist) as fh:
        reader = csv.DictReader(fh)
        incoming = list(reader)
        present = list(reader.fieldnames or [])

    if not incoming:
        sys.exit("error: %s has no rows" % a.shortlist)

    missing_cols = [c for c in REQUIRED_COLUMNS
                    if c not in present and not (c == "molecule_class" and a.molecule_class)]
    if missing_cols:
        sys.exit("error: %s is missing required column(s): %s\n  present: %s"
                 % (a.shortlist, ", ".join(missing_cols), ", ".join(present)))

    extras = [c.strip() for c in a.extra_columns.split(",") if c.strip()]
    unknown_extras = [c for c in extras if c not in present]
    if unknown_extras:
        sys.exit("error: --extra-columns names absent column(s): %s" % ", ".join(unknown_extras))

    # Build rows in input order. Order is rank; never re-sort.
    rows, problems, seen_names, seen_seqs = [], [], {}, {}
    for i, src in enumerate(incoming, 1):
        # Never normalise the sequence. Silently upcasing lowercase input would hide the
        # very violation the validator exists to catch; the emitter passes bytes through
        # and refuses bad ones instead of repairing them.
        row = {
            "name": (src.get("name") or "").strip(),
            "sequence": (src.get("sequence") or "").strip(),
            "molecule_class": (a.molecule_class
                               or (src.get("molecule_class") or "").strip()
                               or MOLECULE_CLASS_DEFAULT),
        }
        raw = row
        bad = validate_row(i, raw, cap, seen_names, seen_seqs, a.strict_class)
        problems.extend(bad)
        if raw["name"]:
            seen_names.setdefault(raw["name"], i)
        if raw["sequence"]:
            seen_seqs.setdefault(raw["sequence"], i)
        for c in extras:
            row[c] = src.get(c, "")
        rows.append(row)

    if len(rows) > cap:
        problems.insert(0, "row count %d exceeds the documented Track %s cap of %d"
                        % (len(rows), a.track, cap))

    if problems:
        print("")
        for p in problems:
            print("INVALID: %s" % p)
        print("")
        sys.exit("refusing to write a non-conformant submission file (%d problem(s))"
                 % len(problems))

    fieldnames = list(REQUIRED_COLUMNS) + extras
    outdir = os.path.dirname(os.path.abspath(a.out))
    if outdir:
        os.makedirs(outdir, exist_ok=True)
    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames, lineterminator="\n")
        w.writeheader()
        for row in rows:
            w.writerow({k: row.get(k, "") for k in fieldnames})

    manifest = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "generator": "scripts/build_submission.py",
        "is_format_test": bool(a.format_test),
        "source_shortlist": os.path.relpath(a.shortlist),
        "track": a.track,
        "row_cap": cap,
        "rows": len(rows),
        "row_order_is_rank": ROW_ORDER_IS_RANK,
        "columns": fieldnames,
        "twin_strep_tag": a.twin_strep,
        "twin_strep_status": "UNVERIFIED against the official rules -- confirm on the "
                             "authenticated submit page before submission",
        "unverified_spec_items": [
            "molecule_class spellings are paraphrased from a capped fetch",
            "header-row requirement is inferred, not documented",
            "residue alphabet restrictions are inferred",
            "name charset/length/uniqueness is our own convention",
            "whether an uploaded submission can be edited or replaced is unknown "
            "-- assume the first upload is final",
        ],
        "designs": [{"rank": i, "name": r["name"], "length": len(r["sequence"]),
                     "molecule_class": r["molecule_class"]}
                    for i, r in enumerate(rows, 1)],
    }
    man_path = os.path.splitext(a.out)[0] + ".manifest.json"
    with open(man_path, "w") as fh:
        json.dump(manifest, fh, indent=2)

    label = "FORMAT TEST -- NOT A SUBMISSION" if a.format_test else "submission"
    print("wrote %s  (%s, %d/%d rows, Track %s)" % (a.out, label, len(rows), cap, a.track))
    print("wrote %s" % man_path)
    print("")
    print("Twin-Strep tag: %s  [UNVERIFIED -- confirm in the portal]" % a.twin_strep)
    if a.format_test:
        print("")
        print("This file proves the FORMAT only. Its sequences are fixtures, not candidates.")
    print("")
    print("rank  name                              len  class")
    for i, r in enumerate(rows, 1):
        print("%4d  %-32s %4d  %s" % (i, r["name"], len(r["sequence"]), r["molecule_class"]))
    print("")
    print("next: python scripts/novelty_audit.py %s --max-identity 0.30" % a.out)


if __name__ == "__main__":
    main()
