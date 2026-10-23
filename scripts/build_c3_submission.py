"""Build C3 (GLP1R) submission CSV from ranked designs.

Reads the BindCraft ranked CSV + ESMFold refold results, applies quality gates,
and emits a Proteinbase-ready submission CSV.

Usage:
    python3 -I scripts/build_c3_submission.py [--dry-run]
"""
import csv
import hashlib
import json
import re
import sys
from pathlib import Path

AA = set("ACDEFGHIKLMNPQRSTVWY")

RANKED_CSV = Path("results/pods/c3_glp1r/3_Ranked/!_Ranked.csv")
REFOLD_CSV = Path("results/c3_refold_results.csv")
OUTPUT_DIR = Path("results")

DRY_RUN = "--dry-run" in sys.argv


def load_ranked():
    designs = []
    with open(RANKED_CSV) as f:
        for row in csv.DictReader(f):
            designs.append({
                "name": row["design"],
                "sequence": row["Binder_Sequence"],
                "molecule_class": "single_chain",
                "rank": int(row["rank"]),
                "length": int(row["length"]),
                "iPTM": float(row["i_pTM"]),
                "pLDDT": float(row["pLDDT"]),
                "hotspot_contact": float(row["Hotspot_Contact_Fraction"]),
                "clashes": int(float(row["All_Atom_Clashes"])),
            })
    return designs


def load_refold():
    if not REFOLD_CSV.exists():
        return {}
    refold = {}
    with open(REFOLD_CSV) as f:
        for row in csv.DictReader(f):
            refold[row["name"]] = {
                "refold_plddt": float(row["mean_plddt"]),
                "refold_status": row["status"],
            }
    return refold


def validate(designs):
    names = set()
    sequences = set()
    for d in designs:
        name, seq, kind = d["name"], d["sequence"], d["molecule_class"]
        if not isinstance(name, str) or not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]{0,127}', name):
            raise ValueError(f"Invalid name: {name}")
        if name in names or seq in sequences:
            raise ValueError(f"Duplicate: {name}")
        if not isinstance(seq, str) or not 10 <= len(seq) <= 250 or set(seq) - AA:
            raise ValueError(f"Bad sequence for {name}: len={len(seq)}")
        if kind != "single_chain":
            raise ValueError(f"Bad molecule_class for {name}")
        names.add(name)
        sequences.add(seq)


def main():
    designs = load_ranked()
    refold = load_refold()

    for d in designs:
        if d["name"] in refold:
            d.update(refold[d["name"]])

    # Quality gates
    kept = []
    dropped = []
    for d in designs:
        reasons = []
        if d.get("refold_status") == "FAIL":
            reasons.append(f"refold FAIL (pLDDT {d.get('refold_plddt', '?')})")
        if d["clashes"] > 10:
            reasons.append(f"{d['clashes']} all-atom clashes")
        if reasons:
            dropped.append((d, reasons))
        else:
            kept.append(d)

    print(f"Designs: {len(designs)} total, {len(kept)} kept, {len(dropped)} dropped")
    for d, reasons in dropped:
        print(f"  DROPPED {d['name']}: {', '.join(reasons)}")

    validate(kept)

    tag = "dry_run" if DRY_RUN else "final"
    out_csv = OUTPUT_DIR / f"c3_submission_{tag}_20261023.csv"
    if out_csv.exists():
        raise FileExistsError(f"{out_csv} exists")

    rows = [{"name": d["name"], "sequence": d["sequence"], "molecule_class": d["molecule_class"]}
            for d in kept]

    with out_csv.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["name", "sequence", "molecule_class"])
        writer.writeheader()
        writer.writerows(rows)

    audit = {
        "status": f"c3_{tag}",
        "row_count": len(rows),
        "designs_dropped": len(dropped),
        "source": str(RANKED_CSV),
        "refold_source": str(REFOLD_CSV) if REFOLD_CSV.exists() else "NOT_AVAILABLE",
        "csv_sha256": hashlib.sha256(out_csv.read_bytes()).hexdigest(),
    }
    audit_path = out_csv.with_suffix(".audit.json")
    audit_path.write_text(json.dumps(audit, indent=2) + "\n")

    print(f"\nWrote {out_csv} ({len(rows)} designs)")
    print(f"Audit: {audit_path}")

    for d in kept:
        refold_info = f", refold pLDDT={d.get('refold_plddt', 'N/A')}" if 'refold_plddt' in d else ""
        print(f"  {d['rank']:2d}. {d['name']} (len={d['length']}, iPTM={d['iPTM']}, pLDDT={d['pLDDT']}{refold_info})")


if __name__ == "__main__":
    main()
