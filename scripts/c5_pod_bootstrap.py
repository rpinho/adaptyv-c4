#!/usr/bin/env python3
"""
C5 pod bootstrap: download 7OW6, prepare target states, extract JDI binder,
install RFdiffusion + ProteinMPNN.

Run on the pod: python3 c5_pod_bootstrap.py
"""

import os
import subprocess
import sys
import urllib.request
from pathlib import Path


WORKDIR = Path("/workspace/c5")
TARGETS = WORKDIR / "targets"
PEPTIDE_G12D = "VVVGADGVGK"
PEPTIDE_WT = "VVVGAGGVGK"


def download_pdb(pdb_id, outpath):
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    print(f"Downloading {pdb_id} from RCSB...")
    urllib.request.urlretrieve(url, outpath)
    print(f"  Saved to {outpath}")


def parse_chains(pdb_path):
    chains = {}
    with open(pdb_path) as f:
        for line in f:
            if line[:4] != "ATOM":
                continue
            ch = line[21]
            resseq = int(line[22:26])
            if ch not in chains:
                chains[ch] = {"min": resseq, "max": resseq, "ca": 0}
            chains[ch]["min"] = min(chains[ch]["min"], resseq)
            chains[ch]["max"] = max(chains[ch]["max"], resseq)
            if line[12:16].strip() == "CA":
                chains[ch]["ca"] += 1
    return chains


def extract_chains(pdb_path, keep_chains, outpath, rename_map=None):
    """Extract specific chains, optionally renaming them."""
    if rename_map is None:
        rename_map = {}
    with open(pdb_path) as f:
        lines = f.readlines()

    out = []
    prev_chain = None
    for line in lines:
        if line[:4] == "ATOM" or line[:6] == "HETATM":
            ch = line[21]
            if ch in keep_chains:
                new_ch = rename_map.get(ch, ch)
                out.append(line[:21] + new_ch + line[22:])
                if prev_chain and prev_chain != ch:
                    out.insert(-1, "TER\n")
                prev_chain = ch

    with open(outpath, "w") as f:
        f.writelines(out)
        f.write("TER\nEND\n")
    print(f"  Extracted chains {keep_chains} -> {outpath}")


def mutate_peptide_g12d_to_wt(pdb_path, outpath):
    """Mutate position 6 of chain C from ASP to GLY (G12D -> WT)."""
    with open(pdb_path) as f:
        lines = f.readlines()

    out = []
    for line in lines:
        if (line[:4] == "ATOM" and line[21] == "C" and
                int(line[22:26]) == 6):
            resname = line[17:20].strip()
            if resname == "ASP":
                atom_name = line[12:16].strip()
                if atom_name in ("N", "CA", "C", "O"):
                    out.append(line[:17] + "GLY" + line[20:])
            continue
        out.append(line)

    with open(outpath, "w") as f:
        f.writelines(out)
    print(f"  WT mutant (D6G): {outpath}")


def remove_peptide(pdb_path, outpath):
    """Remove chain C (peptide) for peptide-free state."""
    with open(pdb_path) as f:
        lines = f.readlines()

    out = [l for l in lines if not (
        (l[:4] == "ATOM" or l[:6] == "HETATM") and l[21] == "C"
    )]

    with open(outpath, "w") as f:
        f.writelines(out)
    print(f"  Peptide-free: {outpath}")


def get_contacts(pdb_path, chain1, chain2_set, cutoff=8.0):
    """Find chain1 residues with CA within cutoff of any chain2 CA."""
    import math
    target_cas = []
    binder_cas = []

    with open(pdb_path) as f:
        for line in f:
            if line[:4] != "ATOM" or line[12:16].strip() != "CA":
                continue
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            ch = line[21]
            resseq = int(line[22:26])
            if ch in chain2_set:
                target_cas.append((x, y, z))
            elif ch == chain1:
                binder_cas.append((x, y, z, resseq))

    contacts = set()
    for bx, by, bz, bresseq in binder_cas:
        for tx, ty, tz in target_cas:
            if math.sqrt((bx-tx)**2 + (by-ty)**2 + (bz-tz)**2) < cutoff:
                contacts.add(bresseq)
                break
    return sorted(contacts)


def main():
    WORKDIR.mkdir(parents=True, exist_ok=True)
    TARGETS.mkdir(parents=True, exist_ok=True)

    # Step 1: Download 7OW6
    pdb_7ow6 = TARGETS / "7OW6.pdb"
    if not pdb_7ow6.exists():
        download_pdb("7OW6", pdb_7ow6)

    chains = parse_chains(pdb_7ow6)
    print("\n7OW6 chains:")
    for ch, info in sorted(chains.items()):
        print(f"  Chain {ch}: residues {info['min']}-{info['max']} ({info['ca']} CA atoms)")

    # Step 2: Prepare the three target states
    print("\nPreparing target states...")

    # G12D complex (chains A+B+C, no TCR)
    g12d = TARGETS / "state1_g12d_complex.pdb"
    extract_chains(pdb_7ow6, {"A", "B", "C"}, g12d)

    # WT complex (mutate D6 -> G in chain C)
    wt = TARGETS / "state2_wt_complex.pdb"
    mutate_peptide_g12d_to_wt(g12d, wt)

    # Peptide-free
    pepfree = TARGETS / "state3_peptide_free.pdb"
    remove_peptide(g12d, pepfree)

    # Step 3: Extract JDI binder
    print("\nExtracting JDI binder...")

    jdi_dir = TARGETS / "jdi_binder"
    jdi_dir.mkdir(exist_ok=True)

    # Full complex with TCR for positive control
    full_complex = jdi_dir / "pmhc_tcr_full_control.pdb"
    extract_chains(pdb_7ow6, {"A", "B", "C", "D", "E"}, full_complex)

    # pMHC + TCRβ only (for partial diffusion)
    pmhc_tcrb = jdi_dir / "pmhc_tcrb_for_partial.pdb"
    extract_chains(pdb_7ow6, {"A", "B", "C", "E"}, pmhc_tcrb, rename_map={"E": "D"})

    # TCR-only
    tcr_only = jdi_dir / "tcr_only.pdb"
    extract_chains(pdb_7ow6, {"D", "E"}, tcr_only)

    # Contact analysis
    print("\nContact analysis (CA-CA):")
    pep_e = get_contacts(full_complex, "E", {"C"}, 8.0)
    pep_d = get_contacts(full_complex, "D", {"C"}, 8.0)
    pmhc_e = get_contacts(full_complex, "E", {"A", "B", "C"}, 10.0)
    pmhc_d = get_contacts(full_complex, "D", {"A", "B", "C"}, 10.0)

    print(f"  TCRβ within 8Å of peptide: {pep_e}")
    print(f"  TCRα within 8Å of peptide: {pep_d}")
    print(f"  TCRβ within 10Å of pMHC: {len(pmhc_e)} residues")
    print(f"  TCRα within 10Å of pMHC: {len(pmhc_d)} residues")

    print("\n=== Bootstrap complete ===")
    print(f"Targets: {TARGETS}")
    print(f"  G12D: {g12d}")
    print(f"  WT:   {wt}")
    print(f"  Pepfree: {pepfree}")
    print(f"  JDI control: {full_complex}")
    print(f"  JDI TCRβ partial: {pmhc_tcrb}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
