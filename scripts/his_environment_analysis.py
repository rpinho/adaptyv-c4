"""Analyze His residue local sequence environment for C4 v8 designs.

For each His, reports:
- Position (N-term, mid, C-term)
- Local hydrophobicity (±3 residue window)
- Flanking charged/polar residues
- His-His pairs (within 8 residues, potential cooperative switching)
- His-D/E pairs (within 4 residues, intramolecular pairing risk)

Usage:
    python3 -I scripts/his_environment_analysis.py
"""
import csv
import math
from pathlib import Path

RESULTS = Path(__file__).resolve().parent.parent / "results"
CSV = RESULTS / "submission_v8_20261023.csv"

HYDROPHOBIC = set("AILMFWVP")
CHARGED_NEG = set("DE")
CHARGED_POS = set("KRH")
POLAR = set("STNQYC")

KD = {
    "I": 4.5, "V": 4.2, "L": 3.8, "F": 2.8, "C": 2.5,
    "M": 1.9, "A": 1.8, "G": -0.4, "T": -0.7, "S": -0.8,
    "W": -0.9, "Y": -1.3, "P": -1.6, "H": -3.2, "D": -3.5,
    "E": -3.5, "N": -3.5, "Q": -3.5, "K": -3.9, "R": -4.5,
}


def charge_at_ph(seq, ph):
    charge = 0
    for aa in seq:
        if aa == "K":
            charge += 1 / (1 + 10**(ph - 10.5))
        elif aa == "R":
            charge += 1 / (1 + 10**(ph - 12.4))
        elif aa == "H":
            charge += 1 / (1 + 10**(ph - 6.0))
        elif aa == "D":
            charge -= 1 / (1 + 10**(3.9 - ph))
        elif aa == "E":
            charge -= 1 / (1 + 10**(4.1 - ph))
    charge += 1 / (1 + 10**(ph - 8.0))  # N-term
    charge -= 1 / (1 + 10**(3.1 - ph))  # C-term
    return charge


def analyze_his_environment(name, seq):
    n = len(seq)
    his_positions = [i for i, aa in enumerate(seq) if aa == "H"]

    if not his_positions:
        return {"name": name, "his_count": 0}

    results = []
    for pos in his_positions:
        # Position classification
        frac = pos / (n - 1) if n > 1 else 0.5
        if frac < 0.2:
            region = "N-term"
        elif frac > 0.8:
            region = "C-term"
        else:
            region = "mid"

        # Local hydrophobicity (±3 window)
        window_start = max(0, pos - 3)
        window_end = min(n, pos + 4)
        window = seq[window_start:window_end]
        local_kd = sum(KD.get(aa, 0) for aa in window) / len(window)

        # Flanking residues
        left = seq[pos-1] if pos > 0 else "-"
        right = seq[pos+1] if pos < n-1 else "-"

        # Count neighbors in ±4 window
        near_start = max(0, pos - 4)
        near_end = min(n, pos + 5)
        neighbors = seq[near_start:pos] + seq[pos+1:near_end]
        near_de = sum(1 for aa in neighbors if aa in CHARGED_NEG)
        near_kr = sum(1 for aa in neighbors if aa in CHARGED_POS)
        near_hydrophobic = sum(1 for aa in neighbors if aa in HYDROPHOBIC)

        # His-His pairs within 8 residues
        his_neighbors = [p for p in his_positions if p != pos and abs(p - pos) <= 8]

        # His-D/E pairs within 4 residues (intramolecular pairing risk)
        de_close = []
        for j in range(max(0, pos-4), min(n, pos+5)):
            if j != pos and seq[j] in CHARGED_NEG:
                de_close.append((j, seq[j]))

        # Environment classification
        if local_kd > 0.5:
            env = "BURIED"
        elif near_de >= 2:
            env = "ACIDIC"
        elif near_hydrophobic >= 4:
            env = "HYDROPHOBIC"
        else:
            env = "SURFACE"

        results.append({
            "pos": pos + 1,  # 1-indexed
            "region": region,
            "left": left,
            "right": right,
            "local_kd": round(local_kd, 2),
            "env": env,
            "near_de": near_de,
            "near_kr": near_kr,
            "his_neighbors": [p + 1 for p in his_neighbors],
            "de_close": [(j + 1, aa) for j, aa in de_close],
        })

    # Aggregate metrics
    surface_his = sum(1 for r in results if r["env"] == "SURFACE")
    buried_his = sum(1 for r in results if r["env"] == "BURIED")
    acidic_his = sum(1 for r in results if r["env"] == "ACIDIC")
    mid_his = sum(1 for r in results if r["region"] == "mid")
    his_pairs = sum(1 for r in results if r["his_neighbors"]) // 2
    de_paired = sum(1 for r in results if r["de_close"])

    # Cooperative switching estimate: His pairs within 8 res can cooperate
    cooperative_score = his_pairs * 2  # each pair contributes ~2x to pH switch

    dq = charge_at_ph(seq, 5.5) - charge_at_ph(seq, 7.4)

    return {
        "name": name,
        "length": n,
        "his_count": len(his_positions),
        "his_pct": round(len(his_positions) / n * 100, 1),
        "surface_his": surface_his,
        "buried_his": buried_his,
        "acidic_his": acidic_his,
        "mid_his": mid_his,
        "mid_pct": round(mid_his / len(his_positions) * 100, 1) if his_positions else 0,
        "his_pairs": his_pairs,
        "de_paired": de_paired,
        "cooperative_score": cooperative_score,
        "dq": round(dq, 1),
        "details": results,
    }


def main():
    designs = []
    with open(CSV) as f:
        for r in csv.DictReader(f):
            designs.append((r["name"], r["sequence"]))

    print(f"His environment analysis for {len(designs)} C4 v8 designs")
    print("=" * 120)

    all_results = []
    for name, seq in designs:
        result = analyze_his_environment(name, seq)
        all_results.append(result)

    # Summary table
    print(f"\n{'Design':<30s} {'Len':>4s} {'His':>4s} {'%':>5s} {'Surf':>5s} {'Bur':>4s} "
          f"{'Acid':>5s} {'Mid':>4s} {'Mid%':>5s} {'Pairs':>5s} {'DE-p':>5s} {'Coop':>5s} {'Δq':>5s}")
    print("-" * 120)

    for r in sorted(all_results, key=lambda x: -x.get("cooperative_score", 0)):
        if r["his_count"] == 0:
            continue
        print(f"{r['name']:<30s} {r['length']:>4d} {r['his_count']:>4d} {r['his_pct']:>5.1f} "
              f"{r['surface_his']:>5d} {r['buried_his']:>4d} {r['acidic_his']:>5d} "
              f"{r['mid_his']:>4d} {r['mid_pct']:>5.1f} {r['his_pairs']:>5d} "
              f"{r['de_paired']:>5d} {r['cooperative_score']:>5d} {r['dq']:>5.1f}")

    # Detailed per-His report for top candidates
    print("\n\nDETAILED His ENVIRONMENT — top 5 by cooperative score:")
    print("=" * 120)

    top5 = sorted(all_results, key=lambda x: -x.get("cooperative_score", 0))[:5]
    for r in top5:
        print(f"\n{r['name']} (len={r['length']}, {r['his_count']} His, Δq={r['dq']}):")
        for h in r.get("details", []):
            de_str = ",".join(f"{aa}{p}" for p, aa in h["de_close"]) if h["de_close"] else "none"
            hh_str = ",".join(str(p) for p in h["his_neighbors"]) if h["his_neighbors"] else "none"
            print(f"  H{h['pos']:3d} [{h['region']:>6s}] {h['left']}-H-{h['right']} "
                  f"KD={h['local_kd']:+.2f} env={h['env']:<10s} "
                  f"near_DE={h['near_de']} near_KR={h['near_kr']} "
                  f"His_nbrs=[{hh_str}] DE_close=[{de_str}]")

    # Risk assessment
    print("\n\nRISK FLAGS:")
    print("-" * 80)
    for r in all_results:
        flags = []
        if r.get("buried_his", 0) > 0:
            flags.append(f"{r['buried_his']} buried His (may not contact target)")
        if r.get("acidic_his", 0) > r.get("his_count", 0) * 0.5:
            flags.append(f">{50}% His in acidic environment (intramolecular D/E pairing)")
        if r.get("de_paired", 0) > r.get("his_count", 0) * 0.7:
            flags.append(f"{r['de_paired']}/{r['his_count']} His near D/E (intramolecular pairing risk)")
        if r.get("mid_pct", 0) < 40:
            flags.append(f"only {r['mid_pct']}% mid-chain His (terminal His less likely at interface)")
        if flags:
            print(f"  {r['name']}: {'; '.join(flags)}")


if __name__ == "__main__":
    main()
