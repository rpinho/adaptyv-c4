"""Select aromatic-rich C3 RFD4 candidates to replace redundant l100 variants.

Fable review: "I would trade two of those [l100 variants] for fresh backbones"
and flagged missing Trp/Tyr aromatic complementarity for the W33/W39/Y69 groove.

Filters: same as select_c3_replacements.py + WY% > 8%
"""
import csv
from pathlib import Path

RESULTS = Path(__file__).resolve().parent.parent / "results"

KD_SCALE = {
    'A': 1.8, 'R': -4.5, 'N': -3.5, 'D': -3.5, 'C': 2.5,
    'Q': -3.5, 'E': -3.5, 'G': -0.4, 'H': -3.2, 'I': 4.5,
    'L': 3.8, 'K': -3.9, 'M': 1.9, 'F': 2.8, 'P': -1.6,
    'S': -0.8, 'T': -0.7, 'W': -0.9, 'Y': -1.3, 'V': 4.2
}

ALREADY_IN = {
    "c3_glp1r_ecd_22", "c3_glp1r_ecd_172", "c3_glp1r_ecd_171",
}


def gravy(seq):
    return sum(KD_SCALE.get(aa, 0) for aa in seq) / len(seq)


def max_hydro_window(seq, w=19):
    if len(seq) < w:
        return gravy(seq)
    vals = [KD_SCALE.get(aa, 0) for aa in seq]
    return max(sum(vals[i:i+w]) / w for i in range(len(seq) - w + 1))


def net_charge(seq):
    return seq.count("K") + seq.count("R") - seq.count("D") - seq.count("E")


def hamming_dist(a, b):
    n = min(len(a), len(b))
    if n == 0:
        return 1.0
    return sum(1 for i in range(n) if a[i] != b[i]) / max(len(a), len(b))


def main():
    # Load existing sequences (keeping all except the 2 l100 single variants)
    existing_seqs = []
    with open(RESULTS / "c3_submission_v4_20261023.csv") as f:
        for r in csv.DictReader(f):
            if r["name"] not in ("c3_v_l100_E93Q", "c3_v_l100_E79Q"):
                existing_seqs.append(r["sequence"])
    print(f"Existing (keeping): {len(existing_seqs)}")

    candidates = []
    with open(RESULTS / "c3_rfd4_sequences.csv") as f:
        for r in csv.DictReader(f):
            if r["chain"] != "A":
                continue
            seq = r["sequence"]
            name = r["design"]
            if name in ALREADY_IN:
                continue
            if len(seq) < 45 or len(seq) > 100:
                continue
            if "C" in seq:
                continue
            g = gravy(seq)
            mhw = max_hydro_window(seq)
            if g > 0.3 or mhw > 1.6:
                continue
            nc = net_charge(seq)
            if abs(nc) > 9:
                continue
            wy = (seq.count("W") + seq.count("Y")) / len(seq)
            if wy < 0.08:
                continue
            charged = sum(1 for aa in seq if aa in "KRDEH") / len(seq)
            if charged < 0.15:
                continue
            candidates.append({
                "name": name, "sequence": seq, "length": len(seq),
                "w_count": seq.count("W"), "y_count": seq.count("Y"),
                "wy_pct": round(wy * 100, 1),
                "net_charge": nc, "gravy": round(g, 2),
                "charged_frac": round(charged * 100, 1),
            })

    print(f"WY-rich candidates: {len(candidates)}")

    # Score: WY% (50%) + charge neutrality (25%) + GRAVY toward -0.5 (25%)
    for c in candidates:
        wy_score = min(c["wy_pct"] / 12, 1.0)
        charge_score = 1 - abs(c["net_charge"]) / 10
        gravy_score = 1 - abs(c["gravy"] - (-0.5)) / 1.0
        c["score"] = 0.50 * wy_score + 0.25 * charge_score + 0.25 * gravy_score

    candidates.sort(key=lambda x: -x["score"])

    selected = []
    for c in candidates:
        diverse = all(hamming_dist(c["sequence"], s) > 0.40 for s in existing_seqs)
        diverse_sel = all(hamming_dist(c["sequence"], s["sequence"]) > 0.40 for s in selected)
        if diverse and diverse_sel:
            selected.append(c)
            print(f"  SELECTED: {c['name']} len={c['length']} W={c['w_count']} Y={c['y_count']} "
                  f"WY={c['wy_pct']}% charge={c['net_charge']} GRAVY={c['gravy']}")
            if len(selected) == 4:
                break

    out = RESULTS / "c3_aromatic_candidates.csv"
    with open(out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["name", "sequence", "length", "w_count", "y_count",
                                           "wy_pct", "net_charge", "gravy", "charged_frac", "score"])
        w.writeheader()
        w.writerows(selected)
    print(f"Saved {len(selected)} to {out}")


if __name__ == "__main__":
    main()
