#!/usr/bin/env python3
"""Fail-closed tests for the C3 round-2 72->24 selection. CPU only, no GPU."""
import sys, pathlib
ROOT = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
import c3_round2_select_20261025 as S

RESULTS = []


def make(n_clusters=13, per_cluster=4, **over):
    recs, cluster_of = [], {}
    for c in range(n_clusters):
        for j in range(per_cluster):
            bb = f"bb{c}"
            cluster_of[bb] = f"cl{c}"
            rec = {"seq_id": f"s{c}_{j}", "backbone": bb,
                   "sequence": f"SEQ{c}_{j}", "refold_ok": True,
                   "mean_plddt": 90.0 - j, "frac70": 0.95, "rmsd": 0.5 + j}
            rec.update(over)
            recs.append(rec)
    # 2 sequences per backbone
    extra = []
    for r in recs:
        e = dict(r)
        e["seq_id"] = r["seq_id"] + "b"
        e["sequence"] = r["sequence"] + "B"
        e["rmsd"] = r["rmsd"] + 0.01
        extra.append(e)
    return recs + extra, cluster_of


def check(name, fn):
    try:
        fn(); RESULTS.append((name, True, ""))
    except AssertionError as e:
        RESULTS.append((name, False, f"assertion: {e}"))
    except Exception as e:
        RESULTS.append((name, False, f"{type(e).__name__}: {e}"))


def expect_stop(fn, label):
    try:
        fn()
    except S.SelectionError:
        return
    raise AssertionError(f"{label} did not stop")


def t_selects_exactly_24():
    recs, cl = make()
    chosen, rep = S.select(recs, cl)
    assert len(chosen) == 24, len(chosen)
    assert rep["n_selected"] == 24


def t_never_more_than_two_per_cluster():
    recs, cl = make()
    chosen, _ = S.select(recs, cl)
    counts = {}
    for r in chosen:
        counts[r["cluster"]] = counts.get(r["cluster"], 0) + 1
    assert max(counts.values()) <= 2, counts


def t_disjoint_sequences():
    recs, cl = make()
    chosen, _ = S.select(recs, cl)
    assert len({r["sequence"] for r in chosen}) == 24


def t_rank_order_respected():
    recs, cl = make()
    chosen, _ = S.select(recs, cl)
    # pass 1 takes the best of each cluster (13); pass 2 fills to 24
    first_pass = [r for r in chosen if r["rmsd"] < 1.0]
    assert len(first_pass) >= 13, [r["rmsd"] for r in chosen]
    assert all(S.rank_key(a) <= S.rank_key(b) for a, b in zip(chosen, chosen[1:])), "not ranked"


def t_ineligible_excluded():
    recs, cl = make()
    recs[0]["mean_plddt"] = 70.0      # below 80
    recs[1]["frac70"] = 0.5           # below 0.90
    recs[2]["rmsd"] = 5.0             # above 3.0
    recs[3]["refold_ok"] = False
    chosen, _ = S.select(recs, cl)
    ids = {r["seq_id"] for r in chosen}
    assert not ids & {"s0_0", "s0_1", "s0_2", "s0_3"}, ids
    assert all(S.eligible(r) for r in chosen)


def t_stops_when_cannot_reach_24():
    recs, cl = make(n_clusters=4, per_cluster=4)   # only 8 per backbone-group
    expect_stop(lambda: S.select(recs, cl), "insufficient pool")


def t_stops_on_too_few_clusters():
    # plenty of sequences but they all come from 6 clusters -> <8
    recs, cl = make(n_clusters=6, per_cluster=8)
    expect_stop(lambda: S.select(recs, cl), "too few clusters")


def t_stops_on_duplicate_ids():
    recs, cl = make()
    recs[1]["seq_id"] = recs[0]["seq_id"]
    expect_stop(lambda: S.select(recs, cl), "duplicate ids")


def t_stops_on_unknown_backbone():
    recs, cl = make()
    recs[0]["backbone"] = "nope"
    expect_stop(lambda: S.select(recs, cl), "unknown backbone")


def t_stops_on_no_eligible():
    recs, cl = make()
    for r in recs:
        r["refold_ok"] = False
    expect_stop(lambda: S.select(recs, cl), "no eligible")


def t_threshold_boundaries_inclusive():
    recs, cl = make(n_clusters=13, per_cluster=4)
    r = recs[0]
    r.update(mean_plddt=80.0, frac70=0.90, rmsd=3.0)
    assert S.eligible(r), "boundary values must pass (>= / <=)"


TESTS = [(k[2:], v) for k, v in sorted(globals().items()) if k.startswith("t_") and callable(v)]
for n, f in TESTS:
    check(n, f)
ok = sum(1 for _, p, _ in RESULTS if p)
for n, p, e in RESULTS:
    print(f"{'PASS' if p else 'FAIL'} {n}{'' if p else '  <- ' + e}")
print(f"\n{ok}/{len(RESULTS)} passed")
sys.exit(0 if ok == len(RESULTS) else 1)
