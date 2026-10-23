"""Select E1 binder candidates from RFD4 generation runs.

Reads the extracted sequence CSVs, applies filters (His cap, length),
and ranks for submission alongside the E2 clash-filtered survivors.
"""
import csv
import sys
from collections import Counter
from pathlib import Path

RESULTS = Path(__file__).resolve().parent.parent / "results"
HIS_CAP = 0.15
MIN_LEN = 60
MAX_LEN = 250

def hamming_distance(s1, s2):
    """Normalised Hamming distance between two sequences (0=identical, 1=totally different)."""
    if len(s1) != len(s2):
        return 1.0
    if len(s1) == 0:
        return 0.0
    return sum(a != b for a, b in zip(s1, s2)) / len(s1)

def load_binders(csv_path, source_tag):
    """Load chain-C binder rows from a sequence CSV."""
    rows = []
    with open(csv_path) as f:
        for r in csv.DictReader(f):
            if r["chain"] != "C":
                continue
            seq = r["sequence"]
            his_frac = float(r["his_fraction"])
            n_res = int(r["n_residues"])
            rows.append({
                "name": r["design"],
                "sequence": seq,
                "n_residues": n_res,
                "his_count": int(r["his_count"]),
                "his_fraction": his_frac,
                "source": source_tag,
            })
    return rows

def has_poly_his(seq, run_len=5):
    """Check for poly-His runs of length >= run_len."""
    return "H" * run_len in seq

def longest_repeat(seq, aa="H"):
    """Length of longest homopolymeric run of the given amino acid."""
    max_run = 0
    current = 0
    for c in seq:
        if c == aa:
            current += 1
            max_run = max(max_run, current)
        else:
            current = 0
    return max_run

def main():
    plain = load_binders(RESULTS / "e1_plain_sequences.csv", "e1_plain")
    his015 = load_binders(RESULTS / "e1_his015_sequences.csv", "e1_his015")
    all_binders = plain + his015

    print(f"Total E1 binders: {len(all_binders)} ({len(plain)} plain, {len(his015)} his015)")

    # Filter: His cap
    passed_his = [r for r in all_binders if r["his_fraction"] <= HIS_CAP]
    print(f"After His cap ({HIS_CAP*100:.0f}%): {len(passed_his)}")

    # Filter: no poly-His runs >= 4
    passed_poly = [r for r in passed_his if longest_repeat(r["sequence"], "H") < 4]
    print(f"After poly-His filter (max run < 4): {len(passed_poly)}")

    # Filter: length
    passed_len = [r for r in passed_poly if MIN_LEN <= r["n_residues"] <= MAX_LEN]
    print(f"After length filter ({MIN_LEN}-{MAX_LEN}): {len(passed_len)}")

    # No cysteines (disulfide bonds complicate expression)
    passed_cys = [r for r in passed_len if "C" not in r["sequence"]]
    print(f"After Cys filter (0 Cys): {len(passed_cys)}")

    candidates = passed_cys

    # Score: prefer designs with some His (for pH switching) but not too much
    # Also prefer moderate length (70-110 residues) for expression
    for r in candidates:
        his_score = r["his_fraction"] / HIS_CAP  # 0-1, higher = more His
        len_penalty = abs(r["n_residues"] - 90) / 40  # penalty for being far from 90
        len_penalty = min(len_penalty, 1.0)
        r["score"] = his_score * 0.6 + (1 - len_penalty) * 0.4
        r["longest_his_run"] = longest_repeat(r["sequence"], "H")

    candidates.sort(key=lambda r: -r["score"])

    # Diversity filter: greedy selection, skip if too similar to already selected
    selected = []
    MIN_DIVERSITY = 0.3  # at least 30% different from all selected
    for r in candidates:
        if not selected:
            selected.append(r)
            continue
        # Check diversity against all selected of similar length
        too_similar = False
        for s in selected:
            if abs(len(r["sequence"]) - len(s["sequence"])) <= 5:
                d = hamming_distance(r["sequence"][:min(len(r["sequence"]), len(s["sequence"]))],
                                     s["sequence"][:min(len(r["sequence"]), len(s["sequence"]))])
                if d < MIN_DIVERSITY:
                    too_similar = True
                    break
        if not too_similar:
            selected.append(r)
        if len(selected) >= 17:  # leave 3 slots for E2 survivors
            break

    print(f"\nSelected {len(selected)} diverse E1 candidates:")
    print(f"{'Rank':<5} {'Name':<35} {'Len':<5} {'His%':<7} {'HisRun':<7} {'Score':<7} {'Source':<12}")
    print("-" * 85)
    for i, r in enumerate(selected):
        print(f"{i+1:<5} {r['name']:<35} {r['n_residues']:<5} "
              f"{r['his_fraction']*100:5.1f}% {r['longest_his_run']:<7} "
              f"{r['score']:.3f}   {r['source']}")

    # Summary stats
    if selected:
        lengths = [r["n_residues"] for r in selected]
        his_fracs = [r["his_fraction"] for r in selected]
        sources = Counter(r["source"] for r in selected)
        print(f"\nLength: min={min(lengths)}, max={max(lengths)}, mean={sum(lengths)/len(lengths):.0f}")
        print(f"His%: min={min(his_fracs)*100:.1f}%, max={max(his_fracs)*100:.1f}%, "
              f"mean={sum(his_fracs)/len(his_fracs)*100:.1f}%")
        print(f"Sources: {dict(sources)}")

    # Write selection CSV
    out = RESULTS / "e1_selected_candidates.csv"
    with open(out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["rank", "name", "sequence", "n_residues",
                                           "his_count", "his_fraction", "longest_his_run",
                                           "score", "source"])
        w.writeheader()
        for i, r in enumerate(selected):
            w.writerow({
                "rank": i + 1,
                "name": r["name"],
                "sequence": r["sequence"],
                "n_residues": r["n_residues"],
                "his_count": r["his_count"],
                "his_fraction": f"{r['his_fraction']:.4f}",
                "longest_his_run": r["longest_his_run"],
                "score": f"{r['score']:.4f}",
                "source": r["source"],
            })
    print(f"\nSaved to {out}")

if __name__ == "__main__":
    main()
