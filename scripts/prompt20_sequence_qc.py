#!/usr/bin/env python3
"""Descriptive QC of the 20 supplied sequences; no affinity predictions/ranking.

Usage: python3 scripts/prompt20_sequence_qc.py --out-dir results/prompt20_sequence_qc
All residue coordinates are one-based, inclusive. Sequence separation means
abs(i-j), so 'within five' includes separations 1 through 5. Acidic-neighbor
statistics refer only to residues in the SAME sequence, not receptor contacts.
The supplied charge labels are retained as metadata, not validated net charges.
"""

import argparse
import csv
import json
import re
from pathlib import Path


# name, supplied length, supplied charge, supplied His percentage, sequence
DATA = [
    ("apical_his_148", 78, -8, 4, "SSESTVEFKGHTIYVHTFSTLEEALDFYNSLPDSYTRRLFKNSAGEYVVLYTEDSGPSYDSLLEDYAESKGWKPFEHP"),
    ("apical_his_70", 117, -1, 1, "EVQLVQSGGGLVQPGGSLRLSCAASGNGYTFTDYAMHWFRQAPGKEREWVGWIYPYNGNTNYADSVKGRFTITADTSKNTAYLQMNSLRAEDTAVYYCARGGEYSDYWGQGTQVTVS"),
    ("his015_115", 74, -8, 15, "EHHLHAILEKKGNHWHYQVFQDGELVFKGEGFGSREHALATMHHDIERRFGKILEEELLEETEDHLHFHIVAQP"),
    ("his015_148", 88, -5, 11, "TDKYYETVHQIPDEFHVTAVLILQFFLEHPDETFTLDEILEHLKKKYPHLRFHHHTVSKALEYLVEKHLLETSKDSAGHTLYKLNKDF"),
    ("his015_57", 84, -6, 12, "SLHHLISETAEEIAHLVYEHRSTDRYFSIISSVLLKHFSSLLTEDQIHSLALEIERALRHHIFDKDKIKEIVFSELHSFLSDHL"),
    ("his015_26", 82, -3, 12, "SAHHERTMAILHILQSHGPLTAFQLLEHLRKDYPHLQISDATLYHTLFHLREDGLIESTVKPHTLTTYYALTEKGREFLDSH"),
    ("his015_61", 98, -6, 12, "SLATKTRATDHLSNMLHSLASEILSKKGETAAQFFAEKVKAYVETHHHDILHISFHKHPDGSFSFHVLPHSTATDDQVHSAVDEITAYATSLHDEILS"),
    ("his015_55", 93, -8, 11, "SLTQAEYDILEAIHFLSHHNNTEFVHVSDVLKETNRLRMAKGHSPLTEDEFLTHLHTLHDKGLIEFDPSDSLLVRLTDKGYAAFHHFTGESLL"),
    ("his015_14", 99, -7, 12, "SHRIQDLVTALEALPHDDREELHSLASEYLDSLSHEDHFTRSHSSAKLEAYARHLAANYPTFKAFLEKHYPHFLSDEFVSSNTLTHFAKDVVTAHSLLH"),
    ("his015_190", 89, -7, 10, "TLHVTEPAPDEWHGAVNVYLHKGSHLLALSFFRSGKYSLSHRHGSTLDYFHGLTRAEVIHHLSEFLDISPATAEDIVDELELAADKYFA"),
    ("his015_197", 91, -11, 10, "STELHHYLDDIESSLLKAKSAETETDRNHFIDEALHKVHELREKAPEFSHYLDDIESSLLKAKSAETETDRNHFIDEALHKVHELREKVGE"),
    ("his015_140", 80, -4, 11, "SHWKHILATAKENWESGDLFSLSLKEQLYAHAYSYSVAHHLNEDQSIELLKHVDTYHKEHPSHSTTEFFDFLKSKIADLK"),
    ("his015_177", 66, 0, 6, "SMSKQVLTAIVSQSQGKDTISLKEAAHYAGHSQSTVKMVFENLGIHPDSHGKYKVSDIKDALSSLE"),
    ("his015_20", 75, -2, 8, "SPLKQRVMTLLSHAHSAHHIKLLTFIYLKYVEEYTDKEIAEHLEVSLATVKSYTKKITEHGGVWSYLDVLETEIL"),
    ("his015_174", 73, -3, 12, "SHRYIISKHGTKITVTGHHGLHVGDRVTVSKGSLHWLGTVTHVDPDTDTAEIQLDSPAPHLHAGDEVTVQVLP"),
    ("c4_v_186_max", 97, -3, 14, "HITAQQLLTYGPQKLLNLLQQQYQQGGNKNFILGRLHVIEKNLQRFKHHLEEHYEHAMEHTHDPEMKEHLTHEFEEHMDHVDHLLSRVHTLKNKVQA"),
    ("c4_v_185_4x", 85, -9, 11, "SSSLEDLHLHLEVFILAGQQLNSMLLHYLSHHYHLTEQEIKEHWKAKLNELHITDEFFAQKRQERKEHWEEYLSEVLTGLEEYFS"),
    ("c4_v_62_max", 111, -7, 14, "SSHERLHHLMEVFTHLVESAFANNQTTFKEHATEFKEHFLALSHEEQHTFLEAAHHMFPSLHATLTQALQMLQAGQSSSDVLHTVLYGVANLSHPSHHRTFHEIVNQVMNS"),
    ("plain_170", 88, -10, 11, "SSYERVLHHLVTAKSNLEKGHPDEALHFVEKAINHLPEDASEKLHEYLLTAKSNLEKGHPDEALHFVEKAINHLLSEHPELESEYESM"),
    ("plain_98", 107, -4, 10, "SSESAKARLRATVKSMLHEYGLSEHHAETAVLVFEIFVEHSNESLTAKQVLDIYNKRHPTHTLDLDTIKLHLKYLEDLGLLESHLNSEHTRLYKLSSKGLELVESHH"),
]


def runs(sequence, residue, minimum=2):
    return [
        {"start": m.start() + 1, "end": m.end(), "length": len(m.group())}
        for m in re.finditer(residue + "{" + str(minimum) + ",}", sequence)
    ]


def longest_nonoverlapping_repeats(sequence, minimum=8):
    """Return all longest exact repeat pairs with disjoint copies.

    These are maximal exact substrings, not inferred structural repeat domains.
    A mismatch can split a larger approximate repeat into shorter exact blocks.
    """
    best = minimum - 1
    matches = []
    for i in range(len(sequence)):
        for j in range(i + minimum, len(sequence)):
            size = 0
            while (i + size < j and j + size < len(sequence)
                   and sequence[i + size] == sequence[j + size]):
                size += 1
            if size >= minimum and size >= best:
                record = {
                    "length": size,
                    "copy1": [i + 1, i + size],
                    "copy2": [j + 1, j + size],
                    "sequence": sequence[i:i + size],
                }
                if size > best:
                    best, matches = size, []
                matches.append(record)
    return matches


def analyze(record, radius=5):
    name, supplied_length, supplied_charge, supplied_his_pct, sequence = record
    if not sequence or set(sequence) - set("ACDEFGHIKLMNPQRSTVWY"):
        raise ValueError(f"Invalid amino-acid sequence: {name}")
    histidines = [i for i, aa in enumerate(sequence, 1) if aa == "H"]
    acids = [i for i, aa in enumerate(sequence, 1) if aa in "DE"]
    near = sum(any(abs(h - a) <= radius for a in acids) for h in histidines)
    isolated = len(histidines) - near
    # 'Isolated' here means no D/E in the sequence window, not isolated from H.
    ratio = near / isolated if isolated else ("infinity" if near else None)
    return {
        "name": name,
        "sequence": sequence,
        "supplied_length": supplied_length,
        "length": len(sequence),
        "length_matches": supplied_length == len(sequence),
        "supplied_charge_unvalidated": supplied_charge,
        "supplied_his_percent": supplied_his_pct,
        "his_count": len(histidines),
        "his_percent": round(100 * len(histidines) / len(sequence), 3),
        "his_positions": histidines,
        "his_index_gaps": [b - a for a, b in zip(histidines, histidines[1:])],
        "consecutive_his_runs": runs(sequence, "H"),
        "acidic_residue_count": len(acids),
        "sequence_window_radius": radius,
        "his_near_sequence_acid_count": near,
        "his_without_sequence_acid_count": isolated,
        "near_to_without_ratio_descriptive_only": ratio,
        "fraction_his_near_sequence_acid": round(near / len(histidines), 4) if histidines else None,
        "gln_runs_ge_2": runs(sequence, "Q"),
        "cys_count": sequence.count("C"),
        "longest_exact_nonoverlapping_repeats": longest_nonoverlapping_repeats(sequence),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out-dir", type=Path, required=True)
    args = parser.parse_args()
    results = [analyze(record) for record in DATA]
    args.out_dir.mkdir(parents=True, exist_ok=True)
    (args.out_dir / "sequence_qc.json").write_text(
        json.dumps({
            "interpretation": "Descriptive sequence QC only; no receptor-contact, affinity, or pH-switch prediction.",
            "coordinates": "One-based inclusive; gaps are differences between consecutive His indices.",
            "results": results,
        }, indent=2, allow_nan=False) + "\n"
    )
    with (args.out_dir / "sequence_qc.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(results[0]))
        writer.writeheader()
        for result in results:
            writer.writerow({
                key: json.dumps(value) if isinstance(value, (list, dict)) else value
                for key, value in result.items()
            })
    for result in results:
        print(f"{result['name']:18} L={result['length']:3} H={result['his_count']:2} "
              f"H%={result['his_percent']:6.2f} "
              f"near/without={result['his_near_sequence_acid_count']}/"
              f"{result['his_without_sequence_acid_count']} "
              f"gaps={result['his_index_gaps']}")


if __name__ == "__main__":
    main()
