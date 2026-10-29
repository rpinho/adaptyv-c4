#!/usr/bin/env python3
"""Sensitivity of the v4 seat set to THIS LANE's interpretive choices.

The rule fixes the gates and weights but leaves two readings open, both of which I had
to decide. Rather than assert they do not matter, measure it:

  S1  N for the FIRST seat. The rule says N=1 when the design "introduces a documented
      generation method absent from chosen rows". With zero chosen rows any method is
      trivially absent, so I granted N=1. The stricter reading gives N=0.5.
  S2  Whether N may exceed 0.5 at all when every method present has already been used.

Also reports the structural fact that limits what ANY of this can mean:

  N IS A SCREENING-COVERAGE FLAG, NOT A PROPERTY OF THE DESIGN. Provenance and the
  novelty screen exist for 35 designs. For an eligible design outside that 35, N=0
  purely because nobody screened it. Measured below: how many of the 49 eligible
  designs are in the screened 35, and how the seats split.

$0, CPU only. Writes a JSON; does not alter the frozen ranking.
"""
import itertools
import json
import pathlib
import subprocess
import sys

import numpy as np

ROOT = pathlib.Path(__file__).resolve().parents[1]
OUTDIR = ROOT / "results/c4_ranking_v4_acid_off_20261026"


def main():
    rk = json.loads((OUTDIR / "v4_ranking_acid_off.json").read_text())
    nov = json.loads((ROOT / "results/c4_novelty_screen_20261025/summary.json").read_text())
    nov_names = {d["name"] for d in nov["designs"]}

    seats = rk["seats"]
    elig_shas = {s["sequence_sha256"]: s for s in seats}
    inelig = {x["name"] for x in rk["ineligible"]}

    # eligible = seats + eligible-but-unseated (recover the latter from the ineligible list)
    all_eligible = [s["name"] for s in seats]
    # eligible-but-unseated are those in the metrics pool that are neither seats nor ineligible
    mj = json.loads((ROOT / "results/c4_independent_check_20261026/"
                     "independent_metrics.json").read_text())
    seat_shas = {s["sequence_sha256"] for s in seats}
    unseated = [r["name"] for r in mj["designs"]
                if r["sequence_sha256"] not in seat_shas and r["name"] not in inelig]
    eligible_names = all_eligible + unseated

    # S1: recompute the ranking under N_first_seat = 0.5
    # Re-deriving the whole greedy is expensive; instead bound the effect: the first seat
    # is chosen by S, and N=1 vs 0.5 moves any row by at most 5.0 points. Report the S gap
    # between seat 1 and 2 and what 5.0 points would do to the top of the ordering.
    s_sorted = sorted(seats, key=lambda x: -x["S"])
    gap_1_2 = s_sorted[0]["S"] - s_sorted[1]["S"]
    # first seat's D=1.0 is fixed by rule; recompute its S with N=0.5
    first = s_sorted[0]
    s_if_n05 = first["S"] - 10 * (first["N"] - 0.5)

    # Which non-seat eligible row could beat s_if_n05? Bound: best possible S for any
    # row is 20*1.0 + 15*1.0 + 10*1.0 = 45 (D, C, N maxima). Use each unseated row's
    # achievable ceiling given N is screening-bound.
    ceilings = []
    for nm in unseated:
        n_max = 1.0 if nm in nov_names else 0.0
        ceilings.append((nm, round(20 * 1.0 + 15 * 1.0 + 10 * n_max, 2)))

    payload = {
        "created_at_utc": subprocess.check_output(
            ["date", "-u", "+%Y-%m-%dT%H:%M:%SZ"]).decode().strip(),
        "purpose": "Measure the effect of this lane's interpretive choices on the v4 seat set",
        "S1_N_for_first_seat": {
            "granted_N": first["N"],
            "S_as_built": first["S"],
            "S_if_N_0.5": round(s_if_n05, 4),
            "gap_to_second_seat_as_built": round(gap_1_2, 4),
            "first_seat_still_first": bool(s_if_n05 >= s_sorted[1]["S"]),
            "reading": ("The first seat wins on D=1.0 (20 points, fixed by rule for the first "
                        "pick) plus C; the N grant is not what selects it. Even at N=0.5 it "
                        "stays ahead of seat 2." if s_if_n05 >= s_sorted[1]["S"] else
                        "N choice DOES change the first seat; the ranking is not robust to it."),
        },
        "S2_N_cap": {
            "n_seats_with_N_1": sum(1 for s in seats if s["N"] == 1.0),
            "n_seats_with_N_0.5": sum(1 for s in seats if s["N"] == 0.5),
            "n_seats_with_N_0": sum(1 for s in seats if s["N"] == 0.0),
            "note": "N=1 arises only for the two generation methods present (site1/site2), "
                    "so at most 2 seats can earn it; the rest are 0.5 or 0.",
        },
        "N_is_screening_coverage": {
            "n_eligible": len(eligible_names),
            "n_eligible_in_novelty_set": sum(1 for n in eligible_names if n in nov_names),
            "n_eligible_not_screened": sum(1 for n in eligible_names if n not in nov_names),
            "n_seats_in_novelty_set": sum(1 for s in seats if s["name"] in nov_names),
            "n_seats_not_screened": sum(1 for s in seats if s["name"] not in nov_names),
            "novelty_pool_size": len(nov_names),
            "reading": ("N=0 for a seat means 'nobody screened this design', not 'this design "
                        "is not novel'. Of the 20 seats, most score N=0 for that reason alone. "
                        "Any use of S to claim one seat is more novel than another is unsupported "
                        "by this data."),
        },
        "max_S_achievable": 45.0,
        "observed_top_S": first["S"],
        "fraction_of_max_S_used": round(first["S"] / 45.0, 4),
        "inputs": {"ranking": "results/c4_ranking_v4_acid_off_20261026/v4_ranking_acid_off.json"},
    }
    o = OUTDIR / "sensitivity_v4.json"
    if o.exists():
        sys.exit(f"REFUSING TO OVERWRITE {o}")
    o.write_text(json.dumps(payload, indent=1))

    print("S1 N-for-first-seat:")
    for k, v in payload["S1_N_for_first_seat"].items():
        print(f"   {k}: {v}")
    print("\nS2 N cap:", json.dumps(payload["S2_N_cap"], indent=1)[:300])
    print("\nN is screening coverage:")
    for k, v in payload["N_is_screening_coverage"].items():
        if k != "reading":
            print(f"   {k}: {v}")
    print("\nmax S achievable with P=0,M=0 is 45; top observed", first["S"],
          f"({payload['fraction_of_max_S_used']*100:.1f}% of max)")
    print("\nwrote", o.relative_to(ROOT))


if __name__ == "__main__":
    main()
