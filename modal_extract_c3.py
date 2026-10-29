"""Extract C3 GLP1R binder sequences from RFD4 generation on Modal.

The C3 CIFs are named c3_glp1r_ecd_N_generated.cif (not tfr1_c3_glp1r_*),
so the generic extract_sequences glob doesn't match. This script uses
the correct pattern.

Usage:
    modal run modal_extract_c3.py::extract
"""
import csv
import io
import os
import pathlib

import modal

app = modal.App("extract-c3-seqs")

runs_vol = modal.Volume.from_name("rfd4-runs", create_if_missing=False)

image = (
    modal.Image.debian_slim(python_version="3.10")
    .pip_install("biotite>=1.0.0")
)


@app.function(image=image, volumes={"/runs": runs_vol}, timeout=600, cpu=2.0, memory=4096)
def extract():
    import biotite.structure.io.pdbx as pdbx
    import glob

    base = "/runs/tfr1_c3_glp1r/inference_samples"
    cifs = sorted(glob.glob(f"{base}/c3_glp1r_ecd_*/c3_glp1r_ecd_*_generated.cif"))
    print(f"Found {len(cifs)} CIF files")

    AA3TO1 = {
        'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F',
        'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L',
        'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R',
        'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y',
    }

    rows = []
    for cif_path in cifs:
        name = pathlib.Path(cif_path).parent.name
        try:
            f = pdbx.CIFFile.read(cif_path)
            arr = pdbx.get_structure(f.block, model=1)
            ca = arr[arr.atom_name == "CA"]
            chains = {}
            for a in ca:
                chains.setdefault(a.chain_id, []).append(a)
            for cid in sorted(chains):
                atoms = chains[cid]
                seq = "".join(AA3TO1.get(a.res_name, "X") for a in atoms)
                rows.append({
                    "design": name,
                    "chain": cid,
                    "n_residues": len(atoms),
                    "sequence": seq,
                })
        except Exception as e:
            print(f"ERROR reading {cif_path}: {e}")

    buf = io.StringIO()
    writer = csv.DictWriter(buf, fieldnames=["design", "chain", "n_residues", "sequence"])
    writer.writeheader()
    writer.writerows(rows)
    csv_text = buf.getvalue()

    out_path = "/runs/c3_rfd4_sequences.csv"
    pathlib.Path(out_path).write_text(csv_text)
    runs_vol.commit()

    binder_rows = [r for r in rows if r["chain"] == "A"]
    print(f"\nTotal rows: {len(rows)}, binder (chain A) rows: {len(binder_rows)}")
    if binder_rows:
        lengths = [int(r["n_residues"]) for r in binder_rows]
        print(f"Binder lengths: min={min(lengths)}, max={max(lengths)}, "
              f"mean={sum(lengths)/len(lengths):.0f}")
        cys_count = sum(1 for r in binder_rows if "C" in r["sequence"])
        print(f"Binders with Cys: {cys_count}/{len(binder_rows)}")

    print(f"\nSaved to {out_path}")
    return csv_text
