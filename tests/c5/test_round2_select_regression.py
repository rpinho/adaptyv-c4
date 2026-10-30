#!/usr/bin/env python3
"""Regression: the refactored selector library API must reproduce the FROZEN selection.

The 24-design selection is upstream of every round-2 screen result. If `select()` ever changes
which 24 designs come out (or their order), the screen evidence no longer describes what ran.
This recomputes it from the same inputs and demands an exact name-for-name, order-for-order match.
"""
import hashlib
import json
import pathlib
import sys

ROOT = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
import c3_round2_select_20261025 as S  # noqa: E402
from c3_soluble_filters_20261025 import failures  # noqa: E402

R2 = ROOT / "results/c3_round2_111_20261025"
FROZEN = R2 / "selected24_filtered_20261025.json"
FROZEN_SHA = "ee2d4f70af5e38e2ad13575a2b64624948f277e5fca2e4c1a9a9d28ec4926416"

res = []


def check(name, cond, extra=""):
    res.append(bool(cond))
    print(("PASS " if cond else "FAIL ") + name + ("" if cond else f"  {extra}"))


refold = json.loads((R2 / "refold_checks_20261025.json").read_text())
clusters = json.loads((R2 / "clusters_20261025.json").read_text())["cluster_of"]
seqs = json.loads((R2 / "mpnn_sequences_20261025.json").read_text())
recs = []
for s in seqs:
    if s["sample"] == 0:
        continue
    name = f"{s['backbone']}_s{s['sample']}"
    r = refold.get(name)
    if not r:
        continue
    recs.append({"name": name, "seq_id": name, "backbone": s["backbone"], "sequence": s["seq"],
                 "refold_ok": r["pass"], "rmsd": r["ca_rmsd_vs_designed_backbone"],
                 "mean_plddt": r["refold_mean_plddt"],
                 "frac70": r["refold_frac_residues_plddt_ge70"], "length": s["len"]})
recs = [r for r in recs if not failures(r["sequence"])]

chosen, rep = S.select(recs, clusters)
frozen = json.loads(FROZEN.read_text())
fn = [x["name"] for x in frozen]
cn = [x["name"] for x in chosen]

check("frozen selection file is byte-unchanged",
      hashlib.sha256(FROZEN.read_bytes()).hexdigest() == FROZEN_SHA)
check("selector reproduces the frozen set exactly", set(fn) == set(cn), set(fn) ^ set(cn))
check("selector reproduces the frozen ORDER exactly", fn == cn)
check("24 selected, 24 clusters, <=2 per cluster",
      rep["n_selected"] == 24 and rep["n_clusters"] >= 8 and rep["max_per_cluster"] <= 2, rep)
check("every selected design passes the eligibility gate", all(S.eligible(r) for r in chosen))
check("chosen records are ranked non-decreasing",
      all(S.rank_key(a) <= S.rank_key(b) for a, b in zip(chosen, chosen[1:])))

n, ok = len(res), sum(res)
print(f"\n{ok}/{n} selector regression checks passed")
sys.exit(0 if ok == n else 1)
