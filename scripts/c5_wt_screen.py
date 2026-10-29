#!/usr/bin/env python3
"""
C5 WT selectivity screen.

Takes accepted binder designs from the G12D BindCraft run and scores them
against the WT pMHC target using ColabFold AF2-Multimer. Designs that score
high against both G12D and WT fail the selectivity gate. Designs that score
high against G12D but low against WT are selectivity candidates.

Usage:
    python3 c5_wt_screen.py \
        --designs /root/c5_design/output_g12d/Accepted \
        --wt-pdb /root/c5_design/targets/state2_wt_complex.pdb \
        --pepfree-pdb /root/c5_design/targets/state3_peptide_free.pdb \
        --output /root/c5_design/wt_screen \
        --colabfold colabfold_batch
"""

import argparse
import csv
import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path


def extract_binder_sequence(pdb_path, binder_chain="D"):
    """Extract the binder chain's sequence from a PDB file."""
    three_to_one = {
        'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F',
        'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L',
        'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R',
        'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y',
    }
    residues = {}
    with open(pdb_path) as f:
        for line in f:
            if line[:4] == "ATOM" and line[21] == binder_chain and line[12:16].strip() == "CA":
                resseq = int(line[22:26])
                resname = line[17:20].strip()
                residues[resseq] = three_to_one.get(resname, 'X')
    return "".join(residues[k] for k in sorted(residues))


def extract_target_sequence(pdb_path, chains="ABC"):
    """Extract multi-chain target sequence from PDB."""
    three_to_one = {
        'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F',
        'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L',
        'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R',
        'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y',
    }
    chain_seqs = {}
    for c in chains:
        chain_seqs[c] = {}

    with open(pdb_path) as f:
        for line in f:
            if line[:4] == "ATOM" and line[21] in chains and line[12:16].strip() == "CA":
                c = line[21]
                resseq = int(line[22:26])
                resname = line[17:20].strip()
                chain_seqs[c][resseq] = three_to_one.get(resname, 'X')

    sequences = []
    for c in chains:
        if chain_seqs[c]:
            seq = "".join(chain_seqs[c][k] for k in sorted(chain_seqs[c]))
            sequences.append(seq)
    return sequences


def write_fasta(fasta_path, name, target_seqs, binder_seq):
    """Write a ColabFold-compatible FASTA with target chains + binder."""
    combined = ":".join(target_seqs + [binder_seq])
    with open(fasta_path, 'w') as f:
        f.write(f">{name}\n{combined}\n")


def parse_colabfold_scores(results_dir, name):
    """Parse ColabFold output for ipTM/pTM scores."""
    scores = []
    for f in sorted(Path(results_dir).glob(f"{name}_scores_rank_*.json")):
        with open(f) as fh:
            d = json.load(fh)
            scores.append({
                'rank': f.stem.split('rank_')[1].split('_')[0],
                'ipTM': d.get('iptm', 0),
                'pTM': d.get('ptm', 0),
                'pLDDT': d.get('mean_plddt', 0),
            })
    return scores


def main():
    parser = argparse.ArgumentParser(description='C5 WT selectivity screen')
    parser.add_argument('--designs', required=True, help='Path to Accepted designs dir')
    parser.add_argument('--wt-pdb', required=True, help='WT pMHC PDB path')
    parser.add_argument('--pepfree-pdb', required=True, help='Peptide-free HLA+b2m PDB path')
    parser.add_argument('--output', required=True, help='Output directory')
    parser.add_argument('--colabfold', default='colabfold_batch', help='ColabFold command')
    parser.add_argument('--num-models', type=int, default=3, help='AF2 models to use')
    parser.add_argument('--num-recycles', type=int, default=3, help='Recycles per model')
    args = parser.parse_args()

    designs_dir = Path(args.designs)
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Get WT target sequences
    wt_seqs = extract_target_sequence(args.wt_pdb, "ABC")
    pepfree_seqs = extract_target_sequence(args.pepfree_pdb, "AB")

    # Find all accepted design PDBs
    design_pdbs = sorted(designs_dir.glob("*.pdb"))
    print(f"Found {len(design_pdbs)} accepted designs")
    print(f"WT target: {len(wt_seqs)} chains ({'+'.join(str(len(s)) for s in wt_seqs)} residues)")
    print(f"Pepfree target: {len(pepfree_seqs)} chains ({'+'.join(str(len(s)) for s in pepfree_seqs)} residues)")

    results = []

    for pdb in design_pdbs:
        design_name = pdb.stem
        print(f"\n--- Screening {design_name} ---")

        binder_seq = extract_binder_sequence(str(pdb), binder_chain="D")
        if not binder_seq:
            # Try chain B (some BindCraft outputs use B for binder)
            binder_seq = extract_binder_sequence(str(pdb), binder_chain="B")
        if not binder_seq:
            print(f"  WARNING: no binder sequence found, skipping")
            continue

        print(f"  Binder: {len(binder_seq)} aa")

        # Create FASTAs for WT and pepfree screens
        fasta_dir = output_dir / "fastas"
        fasta_dir.mkdir(exist_ok=True)

        wt_name = f"{design_name}_vs_wt"
        pf_name = f"{design_name}_vs_pepfree"

        write_fasta(fasta_dir / f"{wt_name}.fasta", wt_name, wt_seqs, binder_seq)
        write_fasta(fasta_dir / f"{pf_name}.fasta", pf_name, pepfree_seqs, binder_seq)

        # Run ColabFold for WT
        wt_out = output_dir / "wt_results" / design_name
        wt_out.mkdir(parents=True, exist_ok=True)

        print(f"  Scoring vs WT...")
        cmd_wt = [
            args.colabfold,
            str(fasta_dir / f"{wt_name}.fasta"),
            str(wt_out),
            "--num-models", str(args.num_models),
            "--num-recycle", str(args.num_recycles),
            "--model-type", "alphafold2_multimer_v3",
        ]
        subprocess.run(cmd_wt, capture_output=True, text=True, timeout=3600)

        wt_scores = parse_colabfold_scores(wt_out, wt_name)
        best_wt_iptm = max((s['ipTM'] for s in wt_scores), default=0)

        # Run ColabFold for pepfree
        pf_out = output_dir / "pepfree_results" / design_name
        pf_out.mkdir(parents=True, exist_ok=True)

        print(f"  Scoring vs peptide-free...")
        cmd_pf = [
            args.colabfold,
            str(fasta_dir / f"{pf_name}.fasta"),
            str(pf_out),
            "--num-models", str(args.num_models),
            "--num-recycle", str(args.num_recycles),
            "--model-type", "alphafold2_multimer_v3",
        ]
        subprocess.run(cmd_pf, capture_output=True, text=True, timeout=3600)

        pf_scores = parse_colabfold_scores(pf_out, pf_name)
        best_pf_iptm = max((s['ipTM'] for s in pf_scores), default=0)

        result = {
            'design': design_name,
            'binder_length': len(binder_seq),
            'sequence': binder_seq,
            'best_wt_ipTM': round(best_wt_iptm, 4),
            'best_pepfree_ipTM': round(best_pf_iptm, 4),
            'wt_scores': wt_scores,
            'pepfree_scores': pf_scores,
        }
        results.append(result)

        print(f"  WT best ipTM: {best_wt_iptm:.3f}")
        print(f"  Pepfree best ipTM: {best_pf_iptm:.3f}")

        selectivity = "SELECTIVE" if best_wt_iptm < 0.7 and best_pf_iptm < 0.7 else "NON-SELECTIVE"
        print(f"  Verdict: {selectivity}")

    # Save results
    summary_path = output_dir / "wt_screen_summary.json"
    with open(summary_path, 'w') as f:
        json.dump(results, f, indent=2)
    print(f"\nResults saved to {summary_path}")

    # Print summary table
    print(f"\n{'='*70}")
    print(f"SELECTIVITY SCREEN SUMMARY")
    print(f"{'='*70}")
    print(f"{'Design':<30} {'Len':>4} {'WT ipTM':>8} {'PF ipTM':>8} {'Verdict':>12}")
    for r in results:
        sel = "SELECTIVE" if r['best_wt_ipTM'] < 0.7 and r['best_pepfree_ipTM'] < 0.7 else "NON-SEL"
        print(f"{r['design']:<30} {r['binder_length']:4d} {r['best_wt_ipTM']:8.3f} {r['best_pepfree_ipTM']:8.3f} {sel:>12}")


if __name__ == "__main__":
    sys.exit(main() or 0)
