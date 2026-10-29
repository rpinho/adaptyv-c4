"""Select 4 C4 replacements for designs flagged by Fable review:
- his015_197 (39-res tandem repeat)
- plain_170 (25-res tandem repeat)
- apical_his_70 (VHH, 1 His, can't pH-switch)
- his015_177 (pI at pH 7.4, precipitation risk, only 4 His)

Selection criteria per Fable recommendations:
- His 8-15% (not too low, not extreme)
- pI NOT near 5.5 or 7.4 (±0.3)
- No exact internal repeats ≥10 residues
- No poly-AA runs ≥4
- Prefer mid-chain His clustering over terminal
- Large charge swing between pH 5.5 and 7.4
"""
import csv
import math
from pathlib import Path

RESULTS = Path(__file__).resolve().parent.parent / "results"

# Current v4 names (excluding the 4 we're dropping)
KEEP_NAMES = {
    "tfr1_his200m_20261023__tfr1_apical_his_148",
    "c4_v_186_max",
    "tfr1_e1_his015_115",
    "tfr1_e1_his015_148",
    "tfr1_e1_his015_57",
    "tfr1_e1_his015_26",
    "tfr1_e1_his015_61",
    "tfr1_e1_his015_55",
    "tfr1_e1_his015_14",
    "tfr1_e1_his015_190",
    "c4_v_185_4x",
    "c4_v_62_max",
    "tfr1_e1_his015_140",
    "tfr1_e1_his015_20",
    "tfr1_e1_his015_174",
    "tfr1_e1_plain_98",
}

DROP_NAMES = {
    "tfr1_e1_his015_197",
    "tfr1_e1_plain_170",
    "tfr1_his030_200m_20261023__tfr1_apical_his_70",
    "tfr1_e1_his015_177",
}

# Already selected (to avoid duplicates)
ALREADY_SELECTED = set()


def longest_exact_repeat(seq):
    n = len(seq)
    best = 0
    for length in range(10, n // 2 + 1):
        for i in range(n - 2 * length + 1):
            sub = seq[i:i+length]
            if sub in seq[i+length:]:
                best = max(best, length)
    return best


def poly_run(seq, threshold=4):
    max_run = 1
    cur = 1
    for i in range(1, len(seq)):
        if seq[i] == seq[i-1]:
            cur += 1
            max_run = max(max_run, cur)
        else:
            cur = 1
    return max_run >= threshold


def his_midchain_score(seq):
    """Score how 'mid-chain' the His residues are (0-1, higher = more mid-chain)."""
    positions = [i for i, aa in enumerate(seq) if aa == "H"]
    if not positions:
        return 0
    n = len(seq)
    scores = []
    for p in positions:
        rel = p / (n - 1) if n > 1 else 0.5
        scores.append(1 - abs(rel - 0.5) * 2)  # 1 at center, 0 at termini
    return sum(scores) / len(scores)


def approx_pI(seq):
    """Rough pI estimate using Henderson-Hasselbalch for K/R/H/D/E."""
    pKa = {'D': 3.9, 'E': 4.1, 'H': 6.04, 'K': 10.5, 'R': 12.5,
            'Nterm': 8.0, 'Cterm': 3.1}
    for pH in [x * 0.1 for x in range(10, 140)]:
        charge = 0
        charge += 1 / (1 + 10 ** (pH - pKa['Nterm']))  # N-term
        charge -= 1 / (1 + 10 ** (pKa['Cterm'] - pH))  # C-term
        for aa in seq:
            if aa == 'K':
                charge += 1 / (1 + 10 ** (pH - pKa['K']))
            elif aa == 'R':
                charge += 1 / (1 + 10 ** (pH - pKa['R']))
            elif aa == 'H':
                charge += 1 / (1 + 10 ** (pH - pKa['H']))
            elif aa == 'D':
                charge -= 1 / (1 + 10 ** (pKa['D'] - pH))
            elif aa == 'E':
                charge -= 1 / (1 + 10 ** (pKa['E'] - pH))
        if abs(charge) < 0.05:
            return pH
    return 7.0


def charge_at_pH(seq, pH):
    pKa = {'D': 3.9, 'E': 4.1, 'H': 6.04, 'K': 10.5, 'R': 12.5,
            'Nterm': 8.0, 'Cterm': 3.1}
    charge = 0
    charge += 1 / (1 + 10 ** (pH - pKa['Nterm']))
    charge -= 1 / (1 + 10 ** (pKa['Cterm'] - pH))
    for aa in seq:
        if aa == 'K':
            charge += 1 / (1 + 10 ** (pH - pKa['K']))
        elif aa == 'R':
            charge += 1 / (1 + 10 ** (pH - pKa['R']))
        elif aa == 'H':
            charge += 1 / (1 + 10 ** (pH - pKa['H']))
        elif aa == 'D':
            charge -= 1 / (1 + 10 ** (pKa['D'] - pH))
        elif aa == 'E':
            charge -= 1 / (1 + 10 ** (pKa['E'] - pH))
    return charge


def hamming_dist(a, b):
    n = min(len(a), len(b))
    if n == 0:
        return 1.0
    return sum(1 for i in range(n) if a[i] != b[i]) / max(len(a), len(b))


def main():
    # Load existing sequences we're keeping
    existing_seqs = []
    with open(RESULTS / "submission_v4_20261023.csv") as f:
        for r in csv.DictReader(f):
            if r["name"] not in DROP_NAMES:
                existing_seqs.append(r["sequence"])
    print(f"Keeping {len(existing_seqs)} designs")

    # Load full pool
    candidates = []
    pool_names_seen = set()
    for csv_name in ["e1_his015_sequences.csv", "e1_plain_sequences.csv"]:
        path = RESULTS / csv_name
        if not path.exists():
            print(f"Skipping {csv_name} (not found)")
            continue
        with open(path) as f:
            for r in csv.DictReader(f):
                if r["chain"] != "C":
                    continue
                name = r["design"]
                if name in pool_names_seen:
                    continue
                pool_names_seen.add(name)
                seq = r["sequence"]

                # Skip already in submission
                full_name = name
                if full_name in KEEP_NAMES or full_name in DROP_NAMES:
                    continue

                if len(seq) < 60 or len(seq) > 120:
                    continue
                if "C" in seq:
                    continue

                his_pct = seq.count("H") / len(seq)
                if his_pct < 0.06 or his_pct > 0.16:
                    continue

                if poly_run(seq):
                    continue

                rep_len = longest_exact_repeat(seq)
                if rep_len >= 10:
                    continue

                pI = approx_pI(seq)
                if abs(pI - 5.5) < 0.3 or abs(pI - 7.4) < 0.3:
                    continue

                q55 = charge_at_pH(seq, 5.5)
                q74 = charge_at_pH(seq, 7.4)
                delta_q = q55 - q74

                nc = seq.count("K") + seq.count("R") - seq.count("D") - seq.count("E")
                if abs(nc) > 12:
                    continue

                midchain = his_midchain_score(seq)

                candidates.append({
                    "name": name, "sequence": seq, "length": len(seq),
                    "his_pct": round(his_pct * 100, 1),
                    "net_charge": nc, "pI": round(pI, 2),
                    "charge_55": round(q55, 1), "charge_74": round(q74, 1),
                    "delta_q": round(delta_q, 1),
                    "midchain_his": round(midchain, 3),
                    "repeat_len": rep_len,
                })

    print(f"Candidates after filters: {len(candidates)}")

    # Score: delta_q (40%) + midchain_his (25%) + pI distance from both measurement points (20%) + length pref (15%)
    for c in candidates:
        dq_score = min(c["delta_q"] / 12, 1.0)
        mid_score = c["midchain_his"]
        pI_dist = min(abs(c["pI"] - 5.5), abs(c["pI"] - 7.4))
        pI_score = min(pI_dist / 1.0, 1.0)
        len_score = 1 - abs(c["length"] - 85) / 35
        c["score"] = 0.40 * dq_score + 0.25 * mid_score + 0.20 * pI_score + 0.15 * max(len_score, 0)

    candidates.sort(key=lambda x: -x["score"])

    # Diversity-filtered selection
    selected = []
    for c in candidates:
        diverse_existing = all(hamming_dist(c["sequence"], s) > 0.30 for s in existing_seqs)
        diverse_sel = all(hamming_dist(c["sequence"], s["sequence"]) > 0.30 for s in selected)
        if diverse_existing and diverse_sel:
            selected.append(c)
            print(f"  SELECTED: {c['name']} len={c['length']} His={c['his_pct']}% "
                  f"charge={c['net_charge']} pI={c['pI']} Δq={c['delta_q']} "
                  f"midchain={c['midchain_his']:.2f} score={c['score']:.3f}")
            if len(selected) == 6:  # get 6, ESMFold top 4
                break

    print(f"\nSelected {len(selected)} candidates for ESMFold validation")

    out = RESULTS / "c4_v5_candidates.csv"
    with open(out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["name", "sequence", "length", "his_pct",
                                           "net_charge", "pI", "charge_55", "charge_74",
                                           "delta_q", "midchain_his", "repeat_len", "score"])
        w.writeheader()
        w.writerows(selected)
    print(f"Saved to {out}")


if __name__ == "__main__":
    main()
