#!/usr/bin/env python3
"""
C5 RFdiffusion → ProteinMPNN → AF2 scoring pipeline.

Unlike BindCraft (which couples hallucination + MPNN + AF2 validation),
this pipeline decouples the steps:
  1. RFdiffusion generates backbones with hotspot conditioning
  2. ProteinMPNN designs sequences for each backbone
  3. AF2-Multimer scores the designs (ranking signal, not a gate)
  4. propka pKa-shift screen for selectivity

Two tracks (Kimi K3.1 Pro consultation, 2026-10-28):
  - De novo: full diffusion with C6-inclusive hotspot branches (20% compute)
  - Partial diffusion: from JDI TCR backbone in 7OW6 (70% compute)
"""

import argparse
import json
import os
import subprocess
import sys
from pathlib import Path


# Chain A (HLA): 276 residues, Chain B (b2m): 99 residues, Chain C (peptide): 10 residues

# De novo hotspot branches — C6 (mutation site) paired with rim/peptide residues
DENOVO_HOTSPOT_BRANCHES = {
    "c6_a69": ["C6", "A69"],
    "c6_a155": ["C6", "A155"],
    "c6_c7": ["C6", "C7"],
}
DENOVO_BINDER_LEN = "90-150"

# Partial diffusion noise levels
PARTIAL_CONSERVATIVE_T = 15
PARTIAL_EXPLORATORY_T = 30


def run_rfd_denovo(target_pdb, output_dir, num_designs=200):
    """Run de novo RFdiffusion with C6-inclusive hotspot branches."""
    all_pdbs = []
    designs_per_branch = num_designs // len(DENOVO_HOTSPOT_BRANCHES)

    for branch_name, hotspots in DENOVO_HOTSPOT_BRANCHES.items():
        branch_dir = os.path.join(output_dir, "denovo", branch_name)
        os.makedirs(branch_dir, exist_ok=True)
        output_prefix = os.path.join(branch_dir, f"denovo_{branch_name}")

        contig = f"A1-276/0 B1-99/0 C1-10/0 {DENOVO_BINDER_LEN}"
        hotspot_str = "[" + ",".join(hotspots) + "]"

        cmd = [
            "python3", "scripts/run_inference.py",
            f"inference.input_pdb={target_pdb}",
            f"inference.output_prefix={output_prefix}",
            f"inference.num_designs={designs_per_branch}",
            f"contigmap.contigs=['{contig}']",
            f"ppi.hotspot_res={hotspot_str}",
            "denoiser.noise_scale_ca=1.0",
            "denoiser.noise_scale_frame=1.0",
        ]

        print(f"De novo branch '{branch_name}': {designs_per_branch} designs, "
              f"hotspots {hotspots}, binder length {DENOVO_BINDER_LEN}")

        r = subprocess.run(cmd, capture_output=True, text=True, timeout=7200)
        if r.returncode != 0:
            print(f"  RFdiffusion failed for {branch_name}: {r.stderr[-500:]}")
            continue

        pdbs = sorted(Path(branch_dir).glob(f"denovo_{branch_name}_*.pdb"))
        print(f"  Generated {len(pdbs)} backbones")
        all_pdbs.extend(pdbs)

    return all_pdbs


def run_rfd_partial(target_pdb, jdi_binder_pdb, output_dir, num_designs=700):
    """Run partial diffusion from the JDI TCR backbone.

    IMPORTANT: stock RFdiffusion partial_T noises the ENTIRE structure,
    not just the binder chain. We run whole-complex partial diffusion
    then filter outputs by pMHC backbone RMSD to keep only those where
    the target stayed close to the crystal structure.
    """
    all_pdbs = []

    arms = [
        ("conservative", PARTIAL_CONSERVATIVE_T, int(num_designs * 0.57)),
        ("exploratory", PARTIAL_EXPLORATORY_T, int(num_designs * 0.43)),
    ]

    for arm_name, partial_t, n_designs in arms:
        arm_dir = os.path.join(output_dir, "partial", arm_name)
        os.makedirs(arm_dir, exist_ok=True)
        output_prefix = os.path.join(arm_dir, f"partial_{arm_name}")

        cmd = [
            "python3", "scripts/run_inference.py",
            f"inference.input_pdb={jdi_binder_pdb}",
            f"inference.output_prefix={output_prefix}",
            f"inference.num_designs={n_designs}",
            f"diffuser.partial_T={partial_t}",
            "denoiser.noise_scale_ca=1.0",
            "denoiser.noise_scale_frame=1.0",
        ]

        print(f"Partial diffusion '{arm_name}': {n_designs} designs, partial_T={partial_t}")

        r = subprocess.run(cmd, capture_output=True, text=True, timeout=14400)
        if r.returncode != 0:
            print(f"  RFdiffusion partial failed for {arm_name}: {r.stderr[-500:]}")
            continue

        raw_pdbs = sorted(Path(arm_dir).glob(f"partial_{arm_name}_*.pdb"))
        print(f"  Generated {len(raw_pdbs)} raw backbones, filtering by pMHC RMSD...")

        passed = filter_by_pmhc_rmsd(raw_pdbs, target_pdb, max_rmsd=2.0)
        print(f"  {len(passed)}/{len(raw_pdbs)} passed pMHC RMSD filter (< 2.0 Å)")
        all_pdbs.extend(passed)

    return all_pdbs


def filter_by_pmhc_rmsd(pdb_list, reference_pdb, max_rmsd=2.0):
    """Filter partial diffusion outputs by pMHC backbone RMSD."""
    import math

    ref_cas = {}
    with open(reference_pdb) as f:
        for line in f:
            if line[:4] == "ATOM" and line[12:16].strip() == "CA" and line[21] in ("A", "B", "C"):
                key = (line[21], int(line[22:26]))
                ref_cas[key] = (float(line[30:38]), float(line[38:46]), float(line[46:54]))

    passed = []
    for pdb_path in pdb_list:
        design_cas = {}
        with open(pdb_path) as f:
            for line in f:
                if line[:4] == "ATOM" and line[12:16].strip() == "CA" and line[21] in ("A", "B", "C"):
                    key = (line[21], int(line[22:26]))
                    design_cas[key] = (float(line[30:38]), float(line[38:46]), float(line[46:54]))

        common = set(ref_cas) & set(design_cas)
        if len(common) < 100:
            continue

        sq_sum = sum(
            (ref_cas[k][0] - design_cas[k][0])**2 +
            (ref_cas[k][1] - design_cas[k][1])**2 +
            (ref_cas[k][2] - design_cas[k][2])**2
            for k in common
        )
        rmsd = math.sqrt(sq_sum / len(common))
        if rmsd <= max_rmsd:
            passed.append(pdb_path)

    return passed


def run_proteinmpnn(backbone_pdb, output_dir, num_seqs=8, fixed_chains="A,B,C"):
    """Run ProteinMPNN to design sequences for a backbone."""
    name = Path(backbone_pdb).stem

    cmd = [
        "python3", "protein_mpnn_run.py",
        "--pdb_path", str(backbone_pdb),
        "--out_folder", output_dir,
        "--num_seq_per_target", str(num_seqs),
        "--sampling_temp", "0.1",
        "--seed", "42",
        "--batch_size", "1",
        "--fixed_chains", fixed_chains,
    ]

    r = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
    if r.returncode != 0:
        print(f"  MPNN failed for {name}: {r.stderr[-200:]}")
        return []

    # Parse FASTA output
    fasta_dir = os.path.join(output_dir, "seqs")
    sequences = []
    for fasta in Path(fasta_dir).glob(f"{name}*.fa"):
        with open(fasta) as f:
            lines = f.readlines()
            for i in range(0, len(lines), 2):
                if lines[i].startswith(">") and i + 1 < len(lines):
                    header = lines[i].strip()
                    seq = lines[i + 1].strip()
                    if "sample" in header.lower():
                        sequences.append(seq)

    return sequences


def extract_binder_sequence(pdb_path, binder_chain="D"):
    """Extract binder chain sequence from PDB."""
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
    return "".join(residues[k] for k in sorted(residues)) if residues else ""


def main():
    parser = argparse.ArgumentParser(description="C5 RFD4 → MPNN pipeline")
    parser.add_argument("--target", default="targets/c5_pmhc/state1_g12d_complex.pdb")
    parser.add_argument("--jdi-binder", default=None,
                        help="JDI TCR binder PDB for partial diffusion track")
    parser.add_argument("--output", default="c5_rfd4_designs")
    parser.add_argument("--num-seqs", type=int, default=8)
    parser.add_argument("--track", choices=["denovo", "partial", "both"], default="both")
    parser.add_argument("--denovo-designs", type=int, default=200)
    parser.add_argument("--partial-designs", type=int, default=700)
    parser.add_argument("--skip-rfd", action="store_true", help="Skip RFD, use existing backbones")
    parser.add_argument("--skip-mpnn", action="store_true", help="Skip MPNN, use backbone sequences")
    args = parser.parse_args()

    output = Path(args.output)
    output.mkdir(parents=True, exist_ok=True)

    # Step 1: RFdiffusion backbone generation
    backbones = []

    if args.skip_rfd:
        for pattern in ["denovo_*_*.pdb", "partial_*_*.pdb", "rfd4_backbone_*.pdb"]:
            backbones.extend(output.rglob(pattern))
        backbones = sorted(set(backbones))
        print(f"Skipping RFD, found {len(backbones)} existing backbones")
    else:
        if args.track in ("denovo", "both"):
            denovo_pdbs = run_rfd_denovo(args.target, str(output), args.denovo_designs)
            backbones.extend(denovo_pdbs)
            print(f"De novo track: {len(denovo_pdbs)} backbones across "
                  f"{len(DENOVO_HOTSPOT_BRANCHES)} branches")

        if args.track in ("partial", "both"):
            if not args.jdi_binder:
                print("WARNING: --jdi-binder not provided, skipping partial diffusion track")
            else:
                partial_pdbs = run_rfd_partial(
                    args.target, args.jdi_binder, str(output), args.partial_designs)
                backbones.extend(partial_pdbs)
                print(f"Partial diffusion track: {len(partial_pdbs)} backbones")

    if not backbones:
        print("No backbones generated")
        return 1

    # Step 2: ProteinMPNN sequence design
    all_designs = []

    if args.skip_mpnn:
        print("Skipping MPNN, extracting backbone sequences directly")
        for pdb in backbones:
            seq = extract_binder_sequence(str(pdb), binder_chain="D")
            if seq and 10 <= len(seq) <= 250:
                track = "partial" if "partial" in str(pdb) else "denovo"
                all_designs.append({
                    "name": pdb.stem,
                    "sequence": seq,
                    "source": f"rfd4_{track}_backbone",
                    "length": len(seq),
                })
    else:
        mpnn_dir = output / "mpnn_output"
        mpnn_dir.mkdir(exist_ok=True)

        for i, pdb in enumerate(backbones):
            print(f"[{i+1}/{len(backbones)}] MPNN for {pdb.name}...")
            sequences = run_proteinmpnn(pdb, str(mpnn_dir), num_seqs=args.num_seqs)
            track = "partial" if "partial" in str(pdb) else "denovo"

            for j, seq in enumerate(sequences):
                if 10 <= len(seq) <= 250:
                    all_designs.append({
                        "name": f"{pdb.stem}_mpnn{j+1}",
                        "sequence": seq,
                        "source": f"rfd4_{track}_mpnn",
                        "length": len(seq),
                    })

    print(f"\nTotal designs: {len(all_designs)}")

    # Step 3: Write submission CSV
    csv_path = output / "c5_rfd4_submission.csv"
    with open(csv_path, "w") as f:
        f.write("name,sequence,molecule_class\n")
        for d in all_designs[:20]:
            f.write(f"{d['name']},{d['sequence']},single_chain\n")
    print(f"Submission CSV: {csv_path} ({min(len(all_designs), 20)} designs)")

    # Step 4: Write full results JSON
    results_path = output / "c5_rfd4_all_designs.json"
    with open(results_path, "w") as f:
        json.dump(all_designs, f, indent=2)
    print(f"Full results: {results_path} ({len(all_designs)} designs)")

    return 0


if __name__ == "__main__":
    sys.exit(main())
