"""Build C4 v7: fix c4_v_33_4x clustering + replace near-pI designs.

Changes from v6:
1. c4_v_33_4x → c4_v_33_4x_v2: spread out mutations (D30N, D43N, E48Q, E50Q)
   per Fable review — original clustered 4 mutations at positions 42-48 next to
   H37/H38, lowering His pKa and shifting pH-switch wrong direction.

2. Drop his015_55 (pI 5.56, very near pH 5.5 measurement) — aggregation risk.
3. Drop his015_190 (pI 5.59, very near pH 5.5) — same issue.
4. Add 2 pool replacements with pI well above 5.7.

Usage:
    python3 -I scripts/build_c4_v7.py
"""
import csv
import math
from pathlib import Path

RESULTS = Path(__file__).resolve().parent.parent / "results"

PK = {
    'D': 3.65, 'E': 4.25, 'C': 8.18, 'Y': 10.07,
    'H': 6.00, 'K': 10.53, 'R': 12.48,
    'N_term': 8.0, 'C_term': 3.1,
}

DROP = {"c4_v_33_4x", "tfr1_e1_his015_55", "tfr1_e1_his015_190"}

# Fable-recommended spread-out variant of c4_v_33_4x
# Parent: FLHHVHELFLSFIHKVREVHPHLPFHVIDDLVTIVEHHAWEDDETAMETEVSTLLHEEHLSHETQSKLLSVSTDYYHDLRTLLTS
# Mutations: D30N, D43N, E48Q, E50Q (spans 30-50 instead of 42-48)
C4_V_33_4X_V2_SEQ = "FLHHVHELFLSFIHKVREVHPHLPFHVIDNLVTIVEHHAWEDNETAMQTQVSTLLHEEHLSHETQSKLLSVSTDYYHDLRTLLTS"

ALREADY_IN_V6 = {
    "tfr1_his200m_20261023__tfr1_apical_his_148",
    "c4_v_186_max",
    "tfr1_e1_his015_115", "tfr1_e1_his015_148", "tfr1_e1_his015_57",
    "tfr1_e1_his015_26", "tfr1_e1_his015_61", "tfr1_e1_his015_55",
    "tfr1_e1_his015_14", "tfr1_e1_his015_190",
    "c4_v_185_4x", "c4_v_62_max",
    "tfr1_e1_his015_140", "tfr1_e1_his015_20", "tfr1_e1_his015_174",
    "tfr1_e1_plain_98",
    "c4_v_58_4x", "tfr1_e1_his015_131", "c4_v_33_4x",
    "tfr1_e1_his015_145",
}


def approx_pi(seq):
    pos_pks = [PK['N_term']]
    neg_pks = [PK['C_term']]
    for aa in seq:
        if aa == 'K': pos_pks.append(PK['K'])
        elif aa == 'R': pos_pks.append(PK['R'])
        elif aa == 'H': pos_pks.append(PK['H'])
        elif aa == 'D': neg_pks.append(PK['D'])
        elif aa == 'E': neg_pks.append(PK['E'])
        elif aa == 'C': neg_pks.append(PK['C'])
        elif aa == 'Y': neg_pks.append(PK['Y'])

    lo, hi = 0.0, 14.0
    for _ in range(200):
        mid = (lo + hi) / 2
        q = sum(1 / (1 + 10 ** (mid - pk)) for pk in pos_pks)
        q -= sum(1 / (1 + 10 ** (pk - mid)) for pk in neg_pks)
        if q > 0:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


def net_charge(seq):
    return seq.count("K") + seq.count("R") - seq.count("D") - seq.count("E")


def charge_at_ph(seq, ph):
    q = 0.0
    q += 1 / (1 + 10 ** (ph - PK['N_term']))
    q -= 1 / (1 + 10 ** (PK['C_term'] - ph))
    for aa in seq:
        if aa == 'K': q += 1 / (1 + 10 ** (ph - PK['K']))
        elif aa == 'R': q += 1 / (1 + 10 ** (ph - PK['R']))
        elif aa == 'H': q += 1 / (1 + 10 ** (ph - PK['H']))
        elif aa == 'D': q -= 1 / (1 + 10 ** (PK['D'] - ph))
        elif aa == 'E': q -= 1 / (1 + 10 ** (PK['E'] - ph))
        elif aa == 'C': q -= 1 / (1 + 10 ** (PK['C'] - ph))
        elif aa == 'Y': q -= 1 / (1 + 10 ** (PK['Y'] - ph))
    return q


def has_poly_aa(seq, n=4):
    for i in range(len(seq) - n + 1):
        if len(set(seq[i:i+n])) == 1:
            return True
    return False


def has_exact_repeat(seq, min_len=10):
    for length in range(min_len, len(seq) // 2 + 1):
        seen = set()
        for i in range(len(seq) - length + 1):
            sub = seq[i:i+length]
            if sub in seen:
                return True
            seen.add(sub)
    return False


def main():
    # 1. Verify the v2 variant
    v2 = C4_V_33_4X_V2_SEQ
    assert len(v2) == 85
    his_frac = v2.count("H") / len(v2)
    nc = net_charge(v2)
    pi = approx_pi(v2)
    q55 = charge_at_ph(v2, 5.5)
    q74 = charge_at_ph(v2, 7.4)
    dq = q55 - q74
    print(f"c4_v_33_4x_v2: len={len(v2)}, His={v2.count('H')} ({his_frac*100:.1f}%), "
          f"charge={nc}, pI={pi:.2f}, q@5.5={q55:.1f}, q@7.4={q74:.1f}, Δq={dq:.1f}")

    # 2. Find pool replacements for his015_55 and his015_190
    pool_names_in_v6 = {n.replace("tfr1_e1_", "") for n in ALREADY_IN_V6
                        if n.startswith("tfr1_e1_")}
    pool_names_in_v6.add("plain_98")  # also from the pool

    candidates = []
    with open(RESULTS / "e1_his015_sequences.csv") as f:
        for r in csv.DictReader(f):
            if r["chain"] != "C":
                continue
            name = r["design"]
            short = name.replace("tfr1_e1_", "")
            if short in pool_names_in_v6:
                continue
            seq = r["sequence"]
            if len(seq) < 50 or len(seq) > 120:
                continue
            hf = float(r["his_fraction"])
            if hf < 0.06 or hf > 0.16:
                continue
            if "C" in seq:
                continue
            if has_poly_aa(seq, 4):
                continue
            if has_exact_repeat(seq, 10):
                continue

            pi_val = approx_pi(seq)
            if abs(pi_val - 5.5) < 0.3 or abs(pi_val - 7.4) < 0.3:
                continue

            nc_val = net_charge(seq)
            if abs(nc_val) > 11:
                continue

            q55 = charge_at_ph(seq, 5.5)
            q74 = charge_at_ph(seq, 7.4)
            dq_val = q55 - q74

            if dq_val < 8:
                continue

            candidates.append({
                "name": name, "sequence": seq, "length": len(seq),
                "his_count": int(r["his_count"]), "his_frac": hf,
                "net_charge": nc_val, "pI": round(pi_val, 2),
                "q_5.5": round(q55, 1), "q_7.4": round(q74, 1),
                "delta_q": round(dq_val, 1),
            })

    # Rank by Δq (40%) + pI distance from danger zones (30%) + His midchain (30%)
    for c in candidates:
        dq_score = min(c["delta_q"] / 15, 1.0)
        pi_dist = min(abs(c["pI"] - 5.5), abs(c["pI"] - 7.4))
        pi_score = min(pi_dist / 1.5, 1.0)
        c["score"] = 0.40 * dq_score + 0.30 * pi_score + 0.30 * (c["his_frac"] / 0.16)

    candidates.sort(key=lambda x: -x["score"])

    print(f"\nPool candidates (pI-safe, Δq≥8, His 6-16%): {len(candidates)}")
    print(f"{'Name':30s} Len  His%   pI   Δq   Charge  Score")
    for c in candidates[:10]:
        print(f"  {c['name']:28s} {c['length']:3d}  {c['his_frac']*100:4.1f}  "
              f"{c['pI']:5.2f}  {c['delta_q']:+5.1f}  {c['net_charge']:+4d}    {c['score']:.3f}")

    # Pick top 2 diverse candidates
    selected = []
    for c in candidates:
        diverse = True
        for s in selected:
            n = min(len(c["sequence"]), len(s["sequence"]))
            if n > 0:
                d = sum(1 for i in range(n) if c["sequence"][i] != s["sequence"][i]) / max(len(c["sequence"]), len(s["sequence"]))
                if d < 0.40:
                    diverse = False
                    break
        if diverse:
            selected.append(c)
            print(f"\n  SELECTED: {c['name']} (len={c['length']}, His={c['his_frac']*100:.1f}%, "
                  f"pI={c['pI']}, Δq={c['delta_q']}, charge={c['net_charge']})")
            if len(selected) == 2:
                break

    # 3. Build the ESMFold candidates CSV
    esmfold_rows = [
        {"name": "c4_v_33_4x_v2", "sequence": C4_V_33_4X_V2_SEQ},
    ]
    for s in selected:
        esmfold_rows.append({"name": s["name"], "sequence": s["sequence"]})

    out = RESULTS / "c4_v7_candidates.csv"
    with open(out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["name", "sequence"])
        w.writeheader()
        w.writerows(esmfold_rows)
    print(f"\nSaved {len(esmfold_rows)} candidates to {out}")

    # 4. Build v7 submission
    rows = []
    with open(RESULTS / "submission_v6_20261023.csv") as f:
        for r in csv.DictReader(f):
            if r["name"] in DROP:
                print(f"\n  DROPPED: {r['name']}")
                continue
            rows.append(r)
    print(f"\nKept from v6: {len(rows)}")

    rows.append({"name": "c4_v_33_4x_v2", "sequence": C4_V_33_4X_V2_SEQ, "molecule_class": "single_chain"})
    print(f"  ADDED: c4_v_33_4x_v2 (spread-out mutations)")

    for s in selected:
        rows.append({"name": s["name"], "sequence": s["sequence"], "molecule_class": "single_chain"})
        print(f"  ADDED: {s['name']} (pI={s['pI']}, Δq={s['delta_q']})")

    assert len(rows) == 20, f"Expected 20 designs, got {len(rows)}"

    # Validate
    for r in rows:
        seq = r["sequence"]
        assert 10 <= len(seq) <= 250, f"{r['name']}: length {len(seq)} out of range"
        assert set(seq) <= set("ACDEFGHIKLMNPQRSTVWY"), f"{r['name']}: invalid chars"

    names = [r["name"] for r in rows]
    assert len(names) == len(set(names)), "Duplicate names"

    out_sub = RESULTS / "submission_v7_20261023.csv"
    with open(out_sub, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["name", "sequence", "molecule_class"])
        w.writeheader()
        w.writerows(rows)
    print(f"\nSaved v7 to {out_sub}")

    # Summary
    print("\nV7 design summary:")
    for i, r in enumerate(rows):
        seq = r["sequence"]
        nc_val = net_charge(seq)
        pi_val = approx_pi(seq)
        hf = seq.count("H") / len(seq) * 100
        q55 = charge_at_ph(seq, 5.5)
        q74 = charge_at_ph(seq, 7.4)
        dq_val = q55 - q74
        tag = ""
        if r["name"] == "c4_v_33_4x_v2":
            tag = " [NEW: spread-out mutations]"
        elif r["name"] in [s["name"] for s in selected]:
            tag = " [NEW: pI-safe replacement]"
        print(f"  {i+1:2d}. {r['name']:<48s} len={len(seq):3d} charge={nc_val:+4d} "
              f"pI={pi_val:5.2f} His={hf:4.1f}% Δq={dq_val:+5.1f}{tag}")


if __name__ == "__main__":
    main()
