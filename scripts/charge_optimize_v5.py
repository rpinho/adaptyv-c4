"""Generate charge-optimized variants for the 2 new C4 designs with charge -12.

Conservative approach (4 mutations each, like c4_v_185_4x):
- E→Q and D→N at positions ≥4 residues from any His
- Pick the 4 positions most central in the sequence (away from termini)

Designs:
  his015_58:  charge -12, His 13.5%, 111 res
  his015_33:  charge -12, His 15.3%, 85 res
"""
import csv
from pathlib import Path

RESULTS = Path(__file__).resolve().parent.parent / "results"

DESIGNS = {
    "tfr1_e1_his015_58": None,
    "tfr1_e1_his015_33": None,
}

MIN_HIS_DISTANCE = 4
NUM_MUTATIONS = 4


def load_sequences():
    for csv_name in ["e1_his015_sequences.csv"]:
        with open(RESULTS / csv_name) as f:
            for r in csv.DictReader(f):
                if r["chain"] == "C" and r["design"] in DESIGNS:
                    DESIGNS[r["design"]] = r["sequence"]


def find_safe_positions(seq, num_mutations):
    his_positions = {i for i, aa in enumerate(seq) if aa == "H"}
    candidates = []
    for i, aa in enumerate(seq):
        if aa not in ("E", "D"):
            continue
        if any(abs(i - h) < MIN_HIS_DISTANCE for h in his_positions):
            continue
        centrality = 1 - abs(i / (len(seq) - 1) - 0.5) * 2
        candidates.append((i, aa, centrality))

    candidates.sort(key=lambda x: -x[2])
    return candidates[:num_mutations]


def mutate(seq, positions):
    s = list(seq)
    mutations = []
    for pos, orig, _ in positions:
        new = "Q" if orig == "E" else "N"
        s[pos] = new
        mutations.append(f"{orig}{pos+1}{new}")
    return "".join(s), mutations


def net_charge(seq):
    return seq.count("K") + seq.count("R") - seq.count("D") - seq.count("E")


def main():
    load_sequences()

    variants = []
    for name, seq in DESIGNS.items():
        if seq is None:
            print(f"WARNING: {name} not found!")
            continue

        short = name.split("_")[-1]  # e.g. "58" or "33"
        positions = find_safe_positions(seq, NUM_MUTATIONS)
        mutated, mutation_list = mutate(seq, positions)

        variant_name = f"c4_v_{short}_4x"
        old_charge = net_charge(seq)
        new_charge = net_charge(mutated)

        print(f"{name} → {variant_name}")
        print(f"  Parent charge: {old_charge}, Variant charge: {new_charge}")
        print(f"  Mutations: {', '.join(mutation_list)}")
        print(f"  His preserved: {mutated.count('H')} (was {seq.count('H')})")

        variants.append({
            "name": variant_name,
            "sequence": mutated,
            "molecule_class": "single_chain",
            "parent": name,
            "mutations": "+".join(mutation_list),
            "charge_kr_de": new_charge,
        })

    out = RESULTS / "c4_v5_charge_variants.csv"
    with open(out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["name", "sequence", "molecule_class",
                                           "parent", "mutations", "charge_kr_de"])
        w.writeheader()
        w.writerows(variants)
    print(f"\nSaved {len(variants)} variants to {out}")


if __name__ == "__main__":
    main()
