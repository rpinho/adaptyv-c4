#!/usr/bin/env python3
"""
C5 AF2 diagnostic: JDI positive control + BindCraft sequence triage.

1. JDI positive control: run the known JDI TCR through AF2 to calibrate
   whether our validation pipeline works for pMHC targets at all.
2. BindCraft diagnostic: test hallucinated vs MPNN sequences through AF2
   monomer (fold check) and AF2 multimer (binding check).

Panel consensus (Gemini/DeepSeek): the 0/160 MPNN failure may be
AF2-Multimer's limitation on pMHC, not a design problem.
"""

import json
import os
import subprocess
import sys
import urllib.request
from pathlib import Path

WORKDIR = Path("/workspace/c5_diagnostic")

# HLA + b2m + peptide sequences from 7OW6 (extracted by bootstrap)
# These will be extracted fresh from 7OW6

# Hallucinated sequences from BindCraft v2 (from submission CSV)
HALLUCINATED_SEQS = {
    "c5_g12d_l71_s494197": "MKVDEERVRNWEPEEKEFWKRYIRVELPHDRANLAYLEDLAKQGEPNAPQMAKWMKSRIEAFVRMMEILNS",
    "c5_g12d_l74_s496495": "SPKPSEEKLEKYFDEWIQAEVVFWKSGNPETGMAAYMAVFTKMTPMERKEMKKLMTAWGMWTWWFNMVAAMVAS",
    "c5_g12d_l92_s775831": "MILEAWEEFHEAYVKMYVESWRMSKAYPPSRDEKNPYTQFHIKKMKEMIPYMKKITQLRQKEGHKPPTAEFWGWHFPYEVQPFAVKYMKEVV",
    "c5_g12d_l94_s451111": "MPWEDLTPGEKAFNTEAWRKIDAVWAVHMAVYKGMEKEYRRMTSPAEKAGMKARMVDVFEQYDRDTIAVIREMSEVGKKRGLRHPVHFSWYPTQ",
    "c5_g12d_l100_s567490": "SMSPQERTIDMMRAELEHLDEEYPELYEERRARDPWDAHYYFGFQYEHILDYIITTVFGPFSEEAVAWWVAFAAVWMSEPEAMREVVKEAIKAATEKVRS",
    "c5_g12d_l102_s622442": "MVISELKKKMRGIVKEITEAPPGKLFYAIMKGQAQMDAIIYKGMQRMVKETGMSWEEVIDATYEVMVMTYRPIAERFRAWWNSPKSWEAFKDPNMKMPVVEW",
    "c5_g12d_l109_s634072": "MASKIKMPSYEEMLKNWKEWLYYDAGYYEYQFYNMPDNPKEQKKELSKMNPQDVYYLALYFMFFLGKSPEPVIKEYMAMVEAGMSEKEMVAVLTARVVEYRAEFMEMFM",
}

MPNN_SEQS = {
    "c5_g12d_l94_s451111_mpnn6": "MDWSELTEEEKAFNEKVWKRIDEVWEEYMARYRALEEELRKLTDEAEKERMRAEMVAIFDEYDKATIEVIREASEEGKKLGLKHPVHFSWYPTP",
    "c5_g12d_l94_s451111_mpnn20": "MDWSELTEEEKEFNEQVWKKIDAVWEKYMARYEALEEELRKLTDEAEKAALRAEMVKIFDEYDKDTIAVIREASEEGKKQGLKHPVHFSWYPTP",
}


def extract_sequence_from_pdb(pdb_path, chain_id):
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


def download_pdb(pdb_id, outpath):
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    print(f"Downloading {pdb_id}...")
    urllib.request.urlretrieve(url, outpath)


def run_colabfold(fasta_path, output_dir, num_recycle=3):
    """Run ColabFold AF2 prediction."""
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    cmd = [
        "colabfold_batch",
        str(fasta_path),
        str(output_dir),
        "--num-models", "1",
        "--num-recycle", str(num_recycle),
        "--model-type", "alphafold2_multimer_v3",
    ]
    print(f"  colabfold_batch {fasta_path} -> {output_dir}")
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=3600)
    if r.returncode != 0:
        print(f"  FAILED: {r.stderr[-300:]}")
        return None

    for f in Path(output_dir).glob("*scores*.json"):
        with open(f) as fh:
            data = json.load(fh)
            plddt = data.get("plddt", [])
            mean_plddt = sum(plddt) / len(plddt) if plddt else 0
            return {
                "plddt": mean_plddt,
                "ptm": data.get("ptm", 0),
                "iptm": data.get("iptm", 0),
            }
    return None


def run_monomer_prediction(name, sequence, output_dir):
    """Run AF2 monomer prediction on a single sequence."""
    fasta = Path(output_dir) / f"{name}.fasta"
    with open(fasta, "w") as f:
        f.write(f">{name}\n{sequence}\n")

    result_dir = Path(output_dir) / name
    result_dir.mkdir(parents=True, exist_ok=True)

    cmd = [
        "colabfold_batch",
        str(fasta),
        str(result_dir),
        "--num-models", "1",
        "--num-recycle", "3",
    ]
    print(f"  Monomer: {name} ({len(sequence)} aa)")
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=1800)
    if r.returncode != 0:
        print(f"    FAILED: {r.stderr[-200:]}")
        return None

    for f in result_dir.glob("*scores*.json"):
        with open(f) as fh:
            data = json.load(fh)
            plddt = data.get("plddt", [])
            return sum(plddt) / len(plddt) if plddt else 0
    return None


def main():
    WORKDIR.mkdir(parents=True, exist_ok=True)
    targets = WORKDIR / "targets"
    targets.mkdir(exist_ok=True)

    # Download 7OW6
    pdb_7ow6 = targets / "7OW6.pdb"
    if not pdb_7ow6.exists():
        download_pdb("7OW6", pdb_7ow6)

    # Extract sequences
    hla_seq = extract_sequence_from_pdb(pdb_7ow6, "A")
    b2m_seq = extract_sequence_from_pdb(pdb_7ow6, "B")
    pep_g12d = extract_sequence_from_pdb(pdb_7ow6, "C")
    tcrb_seq = extract_sequence_from_pdb(pdb_7ow6, "E")

    pep_wt = pep_g12d[:5] + "G" + pep_g12d[6:]

    print(f"HLA: {len(hla_seq)} aa")
    print(f"b2m: {len(b2m_seq)} aa")
    print(f"Peptide G12D: {pep_g12d}")
    print(f"Peptide WT:   {pep_wt}")
    print(f"TCRb: {len(tcrb_seq)} aa")

    results = {}

    # =========================================================
    # PART 1: JDI POSITIVE CONTROL
    # =========================================================
    print("\n" + "="*60)
    print("PART 1: JDI POSITIVE CONTROL")
    print("="*60)

    jdi_dir = WORKDIR / "jdi_control"
    jdi_dir.mkdir(exist_ok=True)

    # Test 1a: TCRb + G12D pMHC
    print("\n--- JDI TCRb + G12D pMHC ---")
    fasta = jdi_dir / "jdi_tcrb_g12d.fasta"
    with open(fasta, "w") as f:
        f.write(f">jdi_tcrb_g12d\n{hla_seq}:{b2m_seq}:{pep_g12d}:{tcrb_seq}\n")

    g12d_result = run_colabfold(fasta, jdi_dir / "af2_g12d", num_recycle=6)
    results["jdi_g12d"] = g12d_result
    if g12d_result:
        print(f"  pLDDT={g12d_result['plddt']:.3f} pTM={g12d_result['ptm']:.3f} ipTM={g12d_result['iptm']:.3f}")

    # Test 1b: TCRb + WT pMHC
    print("\n--- JDI TCRb + WT pMHC ---")
    fasta_wt = jdi_dir / "jdi_tcrb_wt.fasta"
    with open(fasta_wt, "w") as f:
        f.write(f">jdi_tcrb_wt\n{hla_seq}:{b2m_seq}:{pep_wt}:{tcrb_seq}\n")

    wt_result = run_colabfold(fasta_wt, jdi_dir / "af2_wt", num_recycle=6)
    results["jdi_wt"] = wt_result
    if wt_result:
        print(f"  pLDDT={wt_result['plddt']:.3f} pTM={wt_result['ptm']:.3f} ipTM={wt_result['iptm']:.3f}")

    # Test 1c: pMHC only (no binder)
    print("\n--- pMHC only (baseline) ---")
    fasta_pmhc = jdi_dir / "pmhc_only.fasta"
    with open(fasta_pmhc, "w") as f:
        f.write(f">pmhc_only\n{hla_seq}:{b2m_seq}:{pep_g12d}\n")

    pmhc_result = run_colabfold(fasta_pmhc, jdi_dir / "af2_pmhc", num_recycle=3)
    results["pmhc_only"] = pmhc_result
    if pmhc_result:
        print(f"  pLDDT={pmhc_result['plddt']:.3f} pTM={pmhc_result['ptm']:.3f} ipTM={pmhc_result['iptm']:.3f}")

    # =========================================================
    # PART 2: BINDCRAFT DIAGNOSTIC
    # =========================================================
    print("\n" + "="*60)
    print("PART 2: BINDCRAFT SEQUENCE DIAGNOSTIC")
    print("="*60)

    bc_dir = WORKDIR / "bindcraft_diagnostic"
    bc_dir.mkdir(exist_ok=True)

    # Step 2a: Monomer predictions (does the sequence fold on its own?)
    print("\n--- Monomer fold check ---")
    monomer_dir = bc_dir / "monomer"
    monomer_dir.mkdir(exist_ok=True)

    monomer_results = {}
    for name, seq in {**HALLUCINATED_SEQS, **MPNN_SEQS}.items():
        plddt = run_monomer_prediction(name, seq, monomer_dir)
        monomer_results[name] = plddt
        if plddt is not None:
            tag = "RAW" if name in HALLUCINATED_SEQS else "MPNN"
            print(f"    [{tag}] {name}: monomer pLDDT = {plddt:.3f}")
        else:
            print(f"    {name}: FAILED")

    results["monomer"] = monomer_results

    # Step 2b: Complex predictions (best 3 hallucinated + both MPNN)
    print("\n--- Complex binding check (top 3 hallucinated + 2 MPNN) ---")
    complex_dir = bc_dir / "complex"
    complex_dir.mkdir(exist_ok=True)

    # Pick top 3 hallucinated by length diversity
    test_seqs = {
        "c5_g12d_l71_s494197": HALLUCINATED_SEQS["c5_g12d_l71_s494197"],
        "c5_g12d_l94_s451111": HALLUCINATED_SEQS["c5_g12d_l94_s451111"],
        "c5_g12d_l109_s634072": HALLUCINATED_SEQS["c5_g12d_l109_s634072"],
    }
    test_seqs.update(MPNN_SEQS)

    complex_results = {}
    for name, seq in test_seqs.items():
        print(f"\n  Testing {name} ({len(seq)} aa) vs G12D pMHC...")
        fasta = complex_dir / f"{name}_complex.fasta"
        with open(fasta, "w") as f:
            f.write(f">{name}\n{hla_seq}:{b2m_seq}:{pep_g12d}:{seq}\n")

        result = run_colabfold(fasta, complex_dir / name, num_recycle=6)
        complex_results[name] = result
        if result:
            tag = "RAW" if name in HALLUCINATED_SEQS else "MPNN"
            print(f"    [{tag}] pLDDT={result['plddt']:.3f} pTM={result['ptm']:.3f} ipTM={result['iptm']:.3f}")

    results["complex"] = complex_results

    # =========================================================
    # PART 3: REPORT
    # =========================================================
    print("\n" + "="*60)
    print("DIAGNOSTIC REPORT")
    print("="*60)

    # JDI interpretation
    print("\n--- JDI Positive Control ---")
    if results.get("jdi_g12d"):
        g = results["jdi_g12d"]
        print(f"  G12D: pLDDT={g['plddt']:.3f}, ipTM={g['iptm']:.3f}")
        if g["iptm"] < 0.3:
            print("  CRITICAL: Known binder gets low ipTM -> AF2 miscalibrated for pMHC")
            results["jdi_conclusion"] = "AF2_MISCALIBRATED"
        elif g["iptm"] < 0.5:
            print("  WARNING: Known binder gets marginal ipTM -> AF2 unreliable as gate")
            results["jdi_conclusion"] = "AF2_UNRELIABLE"
        else:
            print("  OK: Known binder gets reasonable ipTM -> pipeline may work")
            results["jdi_conclusion"] = "AF2_USABLE"

        if results.get("jdi_wt"):
            w = results["jdi_wt"]
            delta = g["iptm"] - w["iptm"]
            print(f"  WT:   pLDDT={w['plddt']:.3f}, ipTM={w['iptm']:.3f}")
            print(f"  Delta ipTM (G12D - WT): {delta:.3f}")
            if abs(delta) < 0.05:
                print("  AF2 cannot discriminate G12D from WT (as expected)")
    else:
        print("  AF2 FAILED on positive control -> do not use as gate")
        results["jdi_conclusion"] = "AF2_FAILED"

    # Monomer interpretation
    print("\n--- Monomer Fold Check ---")
    raw_pass = sum(1 for n, v in monomer_results.items()
                   if n in HALLUCINATED_SEQS and v is not None and v > 0.7)
    mpnn_pass = sum(1 for n, v in monomer_results.items()
                    if n in MPNN_SEQS and v is not None and v > 0.7)
    print(f"  Hallucinated: {raw_pass}/{len(HALLUCINATED_SEQS)} pass pLDDT > 0.7")
    print(f"  MPNN: {mpnn_pass}/{len(MPNN_SEQS)} pass pLDDT > 0.7")

    if raw_pass == 0 and mpnn_pass == 0:
        print("  -> Neither folds -> hallucination overfitting")
        results["monomer_conclusion"] = "OVERFITTING"
    elif raw_pass > 0 and mpnn_pass == 0:
        print("  -> Raw folds but MPNN breaks it -> MPNN destabilizes")
        results["monomer_conclusion"] = "MPNN_DESTABILIZES"
    elif raw_pass > 0 and mpnn_pass > 0:
        print("  -> Both fold -> sequence quality OK")
        results["monomer_conclusion"] = "SEQUENCES_OK"

    # Complex interpretation
    print("\n--- Complex Binding Check ---")
    raw_complex_pass = sum(1 for n, v in complex_results.items()
                          if n in HALLUCINATED_SEQS and v and v["iptm"] > 0.3)
    mpnn_complex_pass = sum(1 for n, v in complex_results.items()
                           if n in MPNN_SEQS and v and v["iptm"] > 0.3)
    print(f"  Hallucinated ipTM > 0.3: {raw_complex_pass}/{sum(1 for n in complex_results if n in HALLUCINATED_SEQS)}")
    print(f"  MPNN ipTM > 0.3: {mpnn_complex_pass}/{sum(1 for n in complex_results if n in MPNN_SEQS)}")

    # Save full report
    report_path = WORKDIR / "diagnostic_report.json"
    with open(report_path, "w") as f:
        json.dump(results, f, indent=2, default=str)
    print(f"\nFull report: {report_path}")

    # Flag file for event-watch
    with open(WORKDIR / "done.flag", "w") as f:
        f.write("DIAGNOSTIC_COMPLETE\n")

    return 0


if __name__ == "__main__":
    sys.exit(main())
