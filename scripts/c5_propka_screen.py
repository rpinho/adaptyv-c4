#!/usr/bin/env python3
"""
C5 propka pKa-shift screen for binder selectivity.

The JDI TCR discriminates G12D from WT through protonation-coupled
electrostatics: TCRβ His100 sits 7.5 Å from the buried Asp6 carboxylate,
and the Asp→Gly charge change shifts its pKa. This script runs propka
on designed binder complexes to check if the binder's His/charged residues
show pKa shifts between G12D and WT contexts — the same readout that
the natural TCR uses.

Usage:
    python3 c5_propka_screen.py \
        --designs /root/c5_design/output_g12d_v2/Accepted \
        --state1 /root/c5_design/targets/state1_g12d_complex.pdb \
        --state2 /root/c5_design/targets/state2_wt_complex.pdb \
        --output /root/c5_design/propka_screen
"""

import argparse
import json
import os
import subprocess
import sys
from pathlib import Path


def run_propka(pdb_path, output_dir):
    """Run propka on a PDB file, return the .pka output path."""
    pdb_name = Path(pdb_path).stem
    pka_path = os.path.join(output_dir, f"{pdb_name}.pka")

    cmd = ["propka3", pdb_path]
    r = subprocess.run(cmd, capture_output=True, text=True, cwd=output_dir, timeout=120)

    # propka writes output to the same dir as input or cwd
    expected = os.path.join(output_dir, f"{Path(pdb_path).stem}.pka")
    if not os.path.exists(expected):
        # Try in the input dir
        expected2 = str(Path(pdb_path).with_suffix('.pka'))
        if os.path.exists(expected2):
            import shutil
            shutil.move(expected2, expected)

    return expected if os.path.exists(expected) else None


def parse_propka_pka(pka_path):
    """Parse propka .pka file for all residue pKa values."""
    residues = []
    in_summary = False

    with open(pka_path) as f:
        for line in f:
            if 'SUMMARY OF THIS PREDICTION' in line:
                in_summary = True
                continue
            if in_summary and line.strip() and not line.startswith('-') and not line.startswith('Group'):
                parts = line.split()
                if len(parts) >= 4:
                    try:
                        resname = parts[0]
                        resnum = int(parts[1])
                        chain = parts[2]
                        pka = float(parts[3])
                        residues.append({
                            'resname': resname,
                            'resnum': resnum,
                            'chain': chain,
                            'pka': pka,
                        })
                    except (ValueError, IndexError):
                        continue
    return residues


def build_complex_pdb(binder_pdb, target_pdb, output_path, binder_chain='D'):
    """Combine binder and target PDBs into one file for propka analysis."""
    with open(output_path, 'w') as out:
        # Write target chains first
        with open(target_pdb) as f:
            for line in f:
                if line.startswith('ATOM') or line.startswith('TER'):
                    out.write(line)

        # Write binder
        with open(binder_pdb) as f:
            for line in f:
                if line.startswith('ATOM'):
                    # Ensure binder chain is set to binder_chain
                    line = line[:21] + binder_chain + line[22:]
                    out.write(line)
                elif line.startswith('TER'):
                    out.write(line)

        out.write('END\n')


def main():
    parser = argparse.ArgumentParser(description='C5 propka pKa-shift selectivity screen')
    parser.add_argument('--designs', required=True, help='Accepted designs directory')
    parser.add_argument('--state1', required=True, help='G12D pMHC PDB')
    parser.add_argument('--state2', required=True, help='WT pMHC PDB')
    parser.add_argument('--output', required=True, help='Output directory')
    args = parser.parse_args()

    output = Path(args.output)
    output.mkdir(parents=True, exist_ok=True)

    designs = sorted(Path(args.designs).glob("*.pdb"))
    print(f"Found {len(designs)} accepted designs")

    results = []

    for binder_pdb in designs:
        design_name = binder_pdb.stem
        print(f"\n--- {design_name} ---")

        design_out = output / design_name
        design_out.mkdir(exist_ok=True)

        # Build complexes: binder + G12D, binder + WT
        g12d_complex = design_out / f"{design_name}_g12d_complex.pdb"
        wt_complex = design_out / f"{design_name}_wt_complex.pdb"

        build_complex_pdb(str(binder_pdb), args.state1, str(g12d_complex))
        build_complex_pdb(str(binder_pdb), args.state2, str(wt_complex))

        # Run propka on both
        g12d_pka_path = run_propka(str(g12d_complex), str(design_out))
        wt_pka_path = run_propka(str(wt_complex), str(design_out))

        if not g12d_pka_path or not wt_pka_path:
            print(f"  propka failed for {design_name}")
            continue

        g12d_pkas = parse_propka_pka(g12d_pka_path)
        wt_pkas = parse_propka_pka(wt_pka_path)

        # Find binder residues (chain D) with pKa shifts
        g12d_binder = {(r['resnum'], r['resname']): r['pka'] for r in g12d_pkas if r['chain'] == 'D'}
        wt_binder = {(r['resnum'], r['resname']): r['pka'] for r in wt_pkas if r['chain'] == 'D'}

        shifts = []
        for key in g12d_binder:
            if key in wt_binder:
                delta = g12d_binder[key] - wt_binder[key]
                if abs(delta) > 0.1:  # Significant shift
                    shifts.append({
                        'resnum': key[0],
                        'resname': key[1],
                        'pka_g12d': round(g12d_binder[key], 2),
                        'pka_wt': round(wt_binder[key], 2),
                        'delta_pka': round(delta, 2),
                    })

        # Sort by absolute shift
        shifts.sort(key=lambda x: abs(x['delta_pka']), reverse=True)

        # His residues are especially interesting (TCR uses His100)
        his_shifts = [s for s in shifts if s['resname'] == 'HIS']

        result = {
            'design': design_name,
            'total_shifts': len(shifts),
            'his_shifts': len(his_shifts),
            'max_shift': round(max(abs(s['delta_pka']) for s in shifts), 2) if shifts else 0,
            'shifts': shifts[:10],  # Top 10
        }
        results.append(result)

        print(f"  Residues with pKa shift > 0.1: {len(shifts)}")
        if his_shifts:
            for h in his_shifts:
                print(f"  HIS{h['resnum']}: pKa G12D={h['pka_g12d']}, WT={h['pka_wt']}, Δ={h['delta_pka']}")
        if shifts:
            print(f"  Max shift: {result['max_shift']} (residue {shifts[0]['resname']}{shifts[0]['resnum']})")

    # Save results
    summary_path = output / "propka_screen_summary.json"
    with open(summary_path, 'w') as f:
        json.dump(results, f, indent=2)
    print(f"\nResults saved to {summary_path}")

    # Print summary
    print(f"\n{'='*60}")
    print(f"PROPKA PKA-SHIFT SCREEN SUMMARY")
    print(f"{'='*60}")
    print(f"{'Design':<30} {'Shifts':>7} {'His':>4} {'Max Δ':>7}")
    for r in results:
        print(f"{r['design']:<30} {r['total_shifts']:7d} {r['his_shifts']:4d} {r['max_shift']:7.2f}")


if __name__ == "__main__":
    sys.exit(main() or 0)
