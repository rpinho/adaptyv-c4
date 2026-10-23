"""Select 2 C3 RFD4 replacements for ecd_96 (TM-like) and fff0d0 (MARGINAL refold).

Filters:
- Chain A, 45-100 residues, no Cys
- NOT TM-like: max 19-mer Kyte-Doolittle window < 1.6, GRAVY < 0.3
- Charged fraction > 15%
- |net charge| < 10
- Hamming diversity > 40% against existing 18 designs

Ranks by: charge neutrality (30%) + length near 70 (20%) + charged fraction (20%) + GRAVY toward -0.5 (30%)
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

EXISTING_NAMES = {
    "c3_glp1r_ecd_22",  # keeping this one
}
DROP_NAMES = {
    "c3_glp1r_ecd_96",  # TM-like, dropping
}


def gravy(seq):
    return sum(KD_SCALE.get(aa, 0) for aa in seq) / len(seq)


def max_hydro_window(seq, w=19):
    if len(seq) < w:
        return gravy(seq)
    vals = [KD_SCALE.get(aa, 0) for aa in seq]
    best = -999
    for i in range(len(seq) - w + 1):
        s = sum(vals[i:i+w]) / w
        if s > best:
            best = s
    return best


def net_charge(seq):
    return seq.count("K") + seq.count("R") - seq.count("D") - seq.count("E")


def charged_frac(seq):
    return sum(1 for aa in seq if aa in "KRDEH") / len(seq)


def hamming_dist(a, b):
    n = min(len(a), len(b))
    if n == 0:
        return 1.0
    return sum(1 for i in range(n) if a[i] != b[i]) / max(len(a), len(b))


def main():
    # Load existing submission sequences (minus the 2 we're dropping)
    existing_seqs = []
    with open(RESULTS / "c3_submission_v3_20261023.csv") as f:
        for r in csv.DictReader(f):
            if r["name"] not in ("c3_glp1r_ecd_96", "c3_bc_l52_fff0d0"):
                existing_seqs.append(r["sequence"])
    print(f"Existing designs (keeping): {len(existing_seqs)}")

    # Load full RFD4 pool
    candidates = []
    with open(RESULTS / "c3_rfd4_sequences.csv") as f:
        for r in csv.DictReader(f):
            if r["chain"] != "A":
                continue
            seq = r["sequence"]
            name = r["design"]
            if name in EXISTING_NAMES or name in DROP_NAMES:
                continue
            if len(seq) < 45 or len(seq) > 100:
                continue
            if "C" in seq:
                continue
            g = gravy(seq)
            mhw = max_hydro_window(seq)
            if g > 0.3:
                continue
            if mhw > 1.6:
                continue
            cf = charged_frac(seq)
            if cf < 0.15:
                continue
            nc = net_charge(seq)
            if abs(nc) > 9:
                continue
            candidates.append({
                "name": name, "sequence": seq, "length": len(seq),
                "gravy": g, "max_hydro_window": mhw, "charged_frac": cf,
                "net_charge": nc,
            })

    print(f"Candidates after filters: {len(candidates)}")

    # Score
    for c in candidates:
        charge_score = 1 - abs(c["net_charge"]) / 10
        length_score = 1 - abs(c["length"] - 70) / 30
        cf_score = min(c["charged_frac"] / 0.3, 1.0)
        gravy_score = 1 - abs(c["gravy"] - (-0.5)) / 1.0
        c["score"] = 0.3 * charge_score + 0.2 * length_score + 0.2 * cf_score + 0.3 * gravy_score

    candidates.sort(key=lambda x: -x["score"])

    # Diversity filter
    selected = []
    for c in candidates:
        diverse = all(hamming_dist(c["sequence"], s) > 0.40 for s in existing_seqs)
        diverse_sel = all(hamming_dist(c["sequence"], s["sequence"]) > 0.40 for s in selected)
        if diverse and diverse_sel:
            selected.append(c)
            print(f"  SELECTED: {c['name']} len={c['length']} charge={c['net_charge']} "
                  f"GRAVY={c['gravy']:.2f} MHW={c['max_hydro_window']:.2f} "
                  f"charged={c['charged_frac']:.1%} score={c['score']:.3f}")
            if len(selected) == 4:  # get top 4 for ESMFold, pick best 2
                break

    print(f"\nSelected {len(selected)} candidates for ESMFold validation")

    # Save
    out = RESULTS / "c3_rfd4_v4_candidates.csv"
    with open(out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["name", "sequence", "length", "net_charge", "gravy",
                                           "max_hydro_window", "charged_frac", "score"])
        w.writeheader()
        w.writerows(selected)
    print(f"Saved to {out}")


if __name__ == "__main__":
    main()
