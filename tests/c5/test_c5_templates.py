#!/usr/bin/env python3
"""Fail-closed checks on the C5 templates.

Purpose: the templates must not be mistakable for completed deliverables. A blank template is
an INCOMPLETE state, and the stop/go checklist must contain the kill conditions that matter.
No GPU, no network.
"""
import pathlib, sys

ROOT = pathlib.Path(__file__).resolve().parents[2]
T = ROOT / "c5_prep" / "templates"
res = []


def check(name, cond, detail=""):
    res.append(bool(cond))
    print(("PASS " if cond else "FAIL ") + name + ("" if cond else f"  {detail}"))


required = ["provenance_methodology.md", "pilot_budget.md", "stop_go_checklist.md", "constraints_table.md"]
for f in required:
    p = T / f
    check(f"template exists: {f}", p.is_file())
    if not p.is_file():
        continue
    s = p.read_text()
    check(f"{f} is declared a TEMPLATE", "TEMPLATE" in s.upper())
    # a template may legitimately contain the word "completed"; it must not claim to BE one
    check(f"{f} does not claim completion",
          "THIS FILE IS COMPLETE" not in s.upper())
    check(f"{f} is not empty", len(s.strip()) > 400)

prov = (T / "provenance_methodology.md").read_text()
for k in ["method", "version", "sha256", "seeds", "de novo", "limitations", "sign-off"]:
    check(f"provenance covers {k!r}", k.lower() in prov.lower())

bud = (T / "pilot_budget.md").read_text()
check("budget states the $0 default", "$0" in bud)
check("budget forbids new Modal spend without approval", "no new modal spend" in bud.lower())
check("budget carries measured unit costs", "0.017/design" in bud)
check("budget requires a ceiling before the first pod",
      "approved **before** the first pod" in bud.lower())

sg = (T / "stop_go_checklist.md").read_text().lower()
for k in ["explicitly approved", "track confirmed", "budget approved", "hash", "nothing uploaded",
          "no measurable progress", "overwritten"]:
    check(f"stop/go covers {k!r}", k in sg)

cons = (T / "constraints_table.md").read_text().lower()
for k in ["selection prompt", "snapshot", "slot limit", "de novo", "decline"]:
    check(f"constraints covers {k!r}", k in cons)
check("constraints refuses to infer from C3/C4", "inferred from c3/c4" in cons)

nc = (ROOT / "results" / "c5_naming_convention.md").read_text().lower()
for k in ["never edited", "overwritten", "sha256", "one line", "without ricardo"]:
    check(f"naming convention covers {k!r}", k in nc)

print(f"\n{sum(res)}/{len(res)} passed")
sys.exit(0 if all(res) else 1)
