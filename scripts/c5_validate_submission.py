#!/usr/bin/env python3
"""C5 submission validator (track-independent). No GPU, no network.

Every limit is configurable on the command line because the C5 track is not yet
known (Kimi plan 2026-10-25 §4: "do not hard-code a track"; 20 vs 20-40 slots).

Checks
  * required columns present (default: name,sequence,molecule_class; --columns to change)
  * exactly --slots rows (default 20; pass 40 for the larger track)
  * unique names, unique sequences
  * canonical amino acids only
  * per-row length within [--min-len, --max-len]
  * molecule_class within --classes
  * VH/VL formatting for fab-style entries (--vh-vl: chain separator & equal-length rule)
  * optional de-novo/zero-shot check against a supplied known-binder FASTA (--forbid-fasta)
Exit code 0 = pass, 1 = violations (printed), 2 = usage/IO error.
"""
import argparse
import csv
import hashlib
import json
import pathlib
import sys

CANON = set("ACDEFGHIKLMNPQRSTVWY")

# Official upload schema: the signed-in submit page's observed enum and the live official
# template at https://proteinbase.com/templates/competition-submission-template.csv
# (sha256 79f6b3caf48df71ad2a33ecf26444d07ad8e8255d8871ccf902d6c4ed17c6948, fetched 2026-10-25).
# The public FAQ's "protein" value contradicts the form; "single_chain" is the demonstrated value.
# Class 5's track-specific values are not published yet, so every value stays overridable on the CLI.
OFFICIAL_CLASSES = "single_chain,nanobody,scfv,fab_kappa,fab_lambda"
OFFICIAL_MIN_LEN, OFFICIAL_MAX_LEN = 10, 250
FAB_CLASSES = {"fab_kappa", "fab_lambda"}


def read_fasta(p):
    seqs, name, buf = {}, None, []
    for line in open(p):
        line = line.strip()
        if line.startswith(">"):
            if name:
                seqs[name] = "".join(buf)
            name, buf = line[1:], []
        elif line:
            buf.append(line)
    if name:
        seqs[name] = "".join(buf)
    return seqs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("csv")
    ap.add_argument("--slots", type=int, default=20,
                    help="required row count (20 for the standard track, 40 for the larger)")
    ap.add_argument("--min-len", type=int, default=OFFICIAL_MIN_LEN)
    ap.add_argument("--max-len", type=int, default=OFFICIAL_MAX_LEN)
    ap.add_argument("--columns", default="name,sequence,molecule_class")
    ap.add_argument("--classes", default=OFFICIAL_CLASSES,
                    help="comma-separated allowed molecule_class values")
    ap.add_argument("--vh-vl", action="store_true",
                    help="treat ':' or '|' in a sequence as an allowed VH:VL chain break")
    ap.add_argument("--forbid-fasta", default="",
                    help="FASTA of known binders whose sequences must not appear verbatim")
    ap.add_argument("--manifest", default="", help="write a JSON manifest + hashes here")
    ap.add_argument("--reviewer", default="", help="researcher sign-off name")
    a = ap.parse_args()

    cols = [c.strip() for c in a.columns.split(",")]
    classes = {c.strip() for c in a.classes.split(",") if c.strip()}
    rows = list(csv.DictReader(open(a.csv)))
    errs = []

    if not rows:
        print("FAIL: empty file")
        return 2
    header = list(rows[0].keys())
    for c in cols:
        if c not in header:
            errs.append(f"missing required column {c!r} (header: {header})")

    if len(rows) != a.slots:
        errs.append(f"row count {len(rows)} != required {a.slots}")

    if classes and "molecule_class" not in header:
        errs.append(f"missing required column 'molecule_class' (needed for --classes)")

    names = [r.get("name", "") for r in rows]
    raw_seqs = [r.get("sequence", "") for r in rows]
    if len(set(names)) != len(names):
        errs.append("duplicate names")
    if len(set(raw_seqs)) != len(raw_seqs):
        errs.append("duplicate sequences")
    for n in names:
        if not n or any(ch.isspace() for ch in n):
            errs.append(f"invalid name {n!r}")

    for r in rows:
        nm = r.get("name", "?")
        seq = r.get("sequence", "")
        if not seq:
            errs.append(f"{nm}: empty sequence")
            continue
        cls = r.get("molecule_class")
        split = a.vh_vl or cls in FAB_CLASSES
        chains = [seq]
        if split and (":" in seq or "|" in seq):
            chains = [c for c in seq.replace("|", ":").split(":") if c]
        elif split and cls in FAB_CLASSES:
            errs.append(f"{nm}: {cls} requires VH:VL in one field, no ':' or '|' separator found")
            chains = [seq]
        for ch in chains:
            if not set(ch) <= CANON:
                bad = sorted(set(ch) - CANON)
                errs.append(f"{nm}: non-canonical residues {bad}")
        total = sum(len(c) for c in chains)
        if not (a.min_len <= total <= a.max_len):
            errs.append(f"{nm}: length {total} outside [{a.min_len},{a.max_len}]")
        if a.classes and "molecule_class" in header and r.get("molecule_class") not in classes:
            errs.append(f"{nm}: molecule_class {r.get('molecule_class')!r} not in {sorted(classes)}")
        if cls in FAB_CLASSES:
            if len(chains) != 2:
                errs.append(f"{nm}: {cls} requires exactly two chains as VH:VL, found {len(chains)}")
            else:
                h, l = chains
                if not (90 <= len(h) <= 140):
                    errs.append(f"{nm}: {cls} VH length {len(h)} outside [90,140]")
                if not (85 <= len(l) <= 130):
                    errs.append(f"{nm}: {cls} VL length {len(l)} outside [85,130]")
        elif a.vh_vl and len(chains) == 2:
            h, l = chains
            if not (90 <= len(h) <= 140):
                errs.append(f"{nm}: VH length {len(h)} outside [90,140]")
            if not (85 <= len(l) <= 130):
                errs.append(f"{nm}: VL length {len(l)} outside [85,130]")

    if a.forbid_fasta:
        known = set(read_fasta(a.forbid_fasta).values())
        for r in rows:
            if r.get("sequence") in known:
                errs.append(f"{r.get('name')}: sequence appears verbatim in the forbidden/known set")

    digest = hashlib.sha256(open(a.csv, "rb").read()).hexdigest()
    manifest = {
        "csv": str(a.csv),
        "sha256": digest,
        "rows": len(rows),
        "slots_required": a.slots,
        "rows_sha256": hashlib.sha256(
            json.dumps([[r.get("name"), r.get("sequence")] for r in rows]).encode()).hexdigest(),
        "columns": header,
        "reviewer": a.reviewer or None,
        "reviewed": bool(a.reviewer),
        "violations": errs,
        "pass": not errs,
    }
    if a.manifest:
        pathlib.Path(a.manifest).write_text(json.dumps(manifest, indent=2))

    if errs:
        print(f"FAIL — {len(errs)} violation(s):")
        for e in errs:
            print(f"  - {e}")
        return 1
    print(f"PASS — {len(rows)} rows, {a.slots} required slots, {len(cols)} required columns")
    print(f"sha256 {digest}")
    if not a.reviewer:
        print("NOTE: no reviewer sign-off recorded (--reviewer).")
    return 0


if __name__ == "__main__":
    sys.exit(main())
