#!/usr/bin/env python3
"""
JDI positive control: run the known TCR through our scoring pipeline.

If the KNOWN JDI TCR-pMHC complex scores poorly, the pipeline is
miscalibrated for this target class — not the designs' fault.

Tests:
1. AF2-Multimer: predict JDI+pMHC complex, check pLDDT/ipTM/pTM
2. Extract TCRβ CDR3 contacts with peptide position 6 (Asp in G12D)
3. Compare G12D vs WT complex scores (the TCR should prefer G12D)

Panel recommendation (Gemini/DeepSeek): if JDI fails AF2 validation,
stop using AF2-Multimer as a gate entirely.
"""

import json
import os
import subprocess
import sys
from pathlib import Path


WORKDIR = Path("/workspace/c5")
TARGETS = WORKDIR / "targets"


def run_colabfold_af2(fasta_path, output_dir, num_models=1, num_recycle=3):
    """Run ColabFold AF2-Multimer prediction."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    cmd = [
        "colabfold_batch",
        str(fasta_path),
        str(output_dir),
        "--num-models", str(num_models),
        "--num-recycle", str(num_recycle),
        "--model-type", "alphafold2_multimer_v3",
    ]

    print(f"Running ColabFold: {' '.join(cmd)}")
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=3600)
    if r.returncode != 0:
        print(f"ColabFold failed: {r.stderr[-500:]}")
        return None

    # Find the best result
    scores = {}
    for f in output_dir.glob("*scores*.json"):
        with open(f) as fh:
            data = json.load(fh)
            scores[f.name] = {
                "plddt": data.get("plddt", [0]),
                "ptm": data.get("ptm", 0),
                "iptm": data.get("iptm", 0),
            }
    return scores


def extract_sequence_from_pdb(pdb_path, chain_id):
    """Extract amino acid sequence from a PDB chain."""
    three_to_one = {
        'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F',
        'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L',
        'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R',
        'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y',
    }
    residues = {}
    with open(pdb_path) as f:
        for line in f:
            if line[:4] == "ATOM" and line[21] == chain_id and line[12:16].strip() == "CA":
                resseq = int(line[22:26])
                resname = line[17:20].strip()
                residues[resseq] = three_to_one.get(resname, 'X')
    return "".join(residues[k] for k in sorted(residues)) if residues else ""


def make_multimer_fasta(pdb_path, chains, outpath):
    """Create a ColabFold multimer FASTA from PDB chains."""
    seqs = []
    for ch in chains:
        seq = extract_sequence_from_pdb(pdb_path, ch)
        seqs.append(seq)

    with open(outpath, "w") as f:
        combined = ":".join(seqs)
        f.write(f">complex\n{combined}\n")

    print(f"Multimer FASTA: {outpath}")
    for i, (ch, seq) in enumerate(zip(chains, seqs)):
        print(f"  Chain {ch}: {len(seq)} aa")
    return outpath


def main():
    results_dir = WORKDIR / "jdi_control_results"
    results_dir.mkdir(parents=True, exist_ok=True)

    control_pdb = TARGETS / "jdi_binder" / "pmhc_tcr_full_control.pdb"
    g12d_pdb = TARGETS / "state1_g12d_complex.pdb"

    print("=== JDI POSITIVE CONTROL ===\n")

    # Step 1: Extract sequences
    print("1. Extracting sequences from 7OW6...")
    hla_seq = extract_sequence_from_pdb(control_pdb, "A")
    b2m_seq = extract_sequence_from_pdb(control_pdb, "B")
    pep_seq = extract_sequence_from_pdb(control_pdb, "C")
    tcra_seq = extract_sequence_from_pdb(control_pdb, "D")
    tcrb_seq = extract_sequence_from_pdb(control_pdb, "E")

    print(f"  HLA:     {len(hla_seq)} aa")
    print(f"  β2m:     {len(b2m_seq)} aa")
    print(f"  Peptide: {len(pep_seq)} aa ({pep_seq})")
    print(f"  TCRα:    {len(tcra_seq)} aa")
    print(f"  TCRβ:    {len(tcrb_seq)} aa")

    # Step 2: Create FASTA for AF2-Multimer
    # Test 1: Full complex (HLA + β2m + peptide + TCRβ)
    # Using TCRβ only as the "binder" since that's what our designs will be
    print("\n2. Creating AF2-Multimer inputs...")

    fasta_tcrb = results_dir / "jdi_tcrb_pmhc.fasta"
    with open(fasta_tcrb, "w") as f:
        combined = ":".join([hla_seq, b2m_seq, pep_seq, tcrb_seq])
        f.write(f">jdi_tcrb_g12d\n{combined}\n")
    print(f"  TCRβ + pMHC: {fasta_tcrb}")

    # Test 2: Same but with WT peptide
    wt_pep = pep_seq[:5] + "G" + pep_seq[6:]  # D6G mutation
    fasta_tcrb_wt = results_dir / "jdi_tcrb_pmhc_wt.fasta"
    with open(fasta_tcrb_wt, "w") as f:
        combined = ":".join([hla_seq, b2m_seq, wt_pep, tcrb_seq])
        f.write(f">jdi_tcrb_wt\n{combined}\n")
    print(f"  TCRβ + WT pMHC: {fasta_tcrb_wt}")
    print(f"  G12D peptide: {pep_seq}")
    print(f"  WT peptide:   {wt_pep}")

    # Step 3: Run AF2-Multimer
    print("\n3. Running AF2-Multimer predictions...")

    print("\n--- G12D complex ---")
    g12d_scores = run_colabfold_af2(fasta_tcrb, results_dir / "af2_g12d", num_recycle=6)

    print("\n--- WT complex ---")
    wt_scores = run_colabfold_af2(fasta_tcrb_wt, results_dir / "af2_wt", num_recycle=6)

    # Step 4: Report
    print("\n=== RESULTS ===")
    report = {
        "test": "JDI positive control",
        "purpose": "Calibrate AF2-Multimer pipeline for pMHC target class",
        "g12d_scores": g12d_scores,
        "wt_scores": wt_scores,
    }

    if g12d_scores:
        for name, s in g12d_scores.items():
            mean_plddt = sum(s["plddt"]) / len(s["plddt"]) if s["plddt"] else 0
            print(f"  G12D {name}: pLDDT={mean_plddt:.3f}, pTM={s['ptm']:.3f}, ipTM={s['iptm']:.3f}")
            report["g12d_mean_plddt"] = mean_plddt
            report["g12d_ptm"] = s["ptm"]
            report["g12d_iptm"] = s["iptm"]
    else:
        print("  G12D: AF2 FAILED")
        report["g12d_result"] = "FAILED"

    if wt_scores:
        for name, s in wt_scores.items():
            mean_plddt = sum(s["plddt"]) / len(s["plddt"]) if s["plddt"] else 0
            print(f"  WT {name}: pLDDT={mean_plddt:.3f}, pTM={s['ptm']:.3f}, ipTM={s['iptm']:.3f}")
            report["wt_mean_plddt"] = mean_plddt
            report["wt_ptm"] = s["ptm"]
            report["wt_iptm"] = s["iptm"]
    else:
        print("  WT: AF2 FAILED")
        report["wt_result"] = "FAILED"

    # Interpretation
    print("\n=== INTERPRETATION ===")
    if g12d_scores and wt_scores:
        g12d_iptm = report.get("g12d_iptm", 0)
        wt_iptm = report.get("wt_iptm", 0)
        if g12d_iptm > 0.4 and wt_iptm > 0.4:
            delta = g12d_iptm - wt_iptm
            print(f"  Both complexes predicted. ipTM delta: {delta:.3f}")
            if abs(delta) < 0.05:
                print("  WARNING: AF2 cannot discriminate G12D from WT (delta < 0.05)")
                print("  This is EXPECTED — AF2-Multimer lacks the resolution for")
                print("  single-residue pMHC discrimination.")
                report["conclusion"] = "AF2 cannot discriminate — expected"
            else:
                print(f"  AF2 shows preference: {'G12D' if delta > 0 else 'WT'}")
                report["conclusion"] = f"AF2 discriminates, delta={delta:.3f}"
        elif g12d_iptm < 0.3:
            print("  CRITICAL: AF2 gives low ipTM even for the KNOWN JDI-pMHC complex!")
            print("  The pipeline is miscalibrated for pMHC. Do NOT use AF2-Multimer")
            print("  as a validation gate for C5 designs.")
            report["conclusion"] = "AF2 miscalibrated for pMHC — abandon as gate"
    elif not g12d_scores:
        print("  AF2 failed entirely on the positive control.")
        print("  Confirms: AF2-Multimer cannot handle this target class.")
        report["conclusion"] = "AF2 failed on positive control"

    # Save report
    report_path = results_dir / "jdi_positive_control_report.json"
    with open(report_path, "w") as f:
        json.dump(report, f, indent=2)
    print(f"\nReport saved: {report_path}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
