#!/usr/bin/env python3
"""Merge pod-scored C3 AF2M results with the 17 Modal-scored designs.

Reads the Modal table (results/c3_screen_af2m_scores_20261024.json) and every
/tmp/podpull/p*/out/*/scores.json, and writes the union to
results/c3_screen_af2m_scores_combined_20261025.json.

No interpretation: this only unions identical-schema rows and reports which
designs (if any) are still missing from the 44-design input.
"""
import csv, glob, json, pathlib, sys

ROOT = pathlib.Path(".")
MODAL = ROOT / "results/c3_screen_af2m_scores_20261024.json"
INPUT = ROOT / "results/c3_screen_input_20261024.csv"
OUT = ROOT / "results/c3_screen_af2m_scores_combined_20261025.json"
PULL = "/tmp/podpull"


def valid(d):
    return isinstance(d.get("iptm"), (int, float)) and isinstance(d.get("ptm"), (int, float))


def main():
    rows = {}
    for r in json.load(open(MODAL)):
        if valid(r):
            rows[r["name"]] = r
    modal_n = len(rows)
    pod_rows = 0
    for f in sorted(glob.glob(f"{PULL}/p*/out/*/scores.json")):
        try:
            d = json.load(open(f))
        except Exception as ex:
            print("skip unreadable", f, ex)
            continue
        if not valid(d):
            continue
        if d["name"] in rows:
            continue
        d.setdefault("run_host", "runpod")
        d.setdefault("gpu", "A40")
        rows[d["name"]] = d
        pod_rows += 1
    ordered = [rows[k] for k in sorted(rows)]
    with open(OUT, "w") as f:
        json.dump(ordered, f, indent=2)
    names = {r["name"] for r in csv.DictReader(open(INPUT))}
    missing = sorted(names - set(rows))
    print(f"modal={modal_n} pod_new={pod_rows} total={len(rows)}/{len(names)}")
    print(f"missing {len(missing)}: {' '.join(missing) if missing else '(none)'}")
    print(f"wrote {OUT}")
    return 0 if not missing else 1


if __name__ == "__main__":
    sys.exit(main())
