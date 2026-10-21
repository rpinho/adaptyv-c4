#!/usr/bin/env python
"""Convert a BindCraft2 output mmCIF to PDB.

BindCraft2 writes a minimal _atom_site loop with no occupancy column, which
Biopython's MMCIFParser rejects outright. Everything downstream of the
generator (propka3, pdb2pqr, most scoring tools) wants PDB anyway, so this
reads the loop by its own header order - never by fixed column positions -
and emits a standard PDB.

Usage:
  python cif2pdb.py design.cif                 # -> design.pdb
  python cif2pdb.py design.cif out.pdb
  python cif2pdb.py dir_of_cifs/ --outdir pdbs/
"""
import argparse, os, sys

WANT = ("group_PDB", "type_symbol", "auth_atom_id", "label_alt_id", "auth_comp_id",
        "auth_asym_id", "auth_seq_id", "pdbx_PDB_ins_code",
        "Cartn_x", "Cartn_y", "Cartn_z", "B_iso_or_equiv", "pdbx_PDB_model_num")


def read_atom_site(path):
    """Return a list of dicts, one per atom, keyed by the _atom_site tag names."""
    cols, rows, in_loop = [], [], False
    with open(path) as fh:
        for line in fh:
            s = line.strip()
            if s.startswith("_atom_site."):
                cols.append(s.split(".", 1)[1].split()[0])
                in_loop = True
                continue
            if not in_loop:
                continue
            if s.startswith(("ATOM", "HETATM")):
                f = s.split()
                if len(f) != len(cols):
                    raise SystemExit("%s: %d fields but %d _atom_site tags"
                                     % (path, len(f), len(cols)))
                rows.append(dict(zip(cols, f)))
            elif s.startswith("#") or s.startswith("loop_") or s.startswith("_"):
                if rows:
                    break
    missing = [c for c in ("auth_atom_id", "auth_comp_id", "auth_asym_id", "auth_seq_id",
                           "Cartn_x", "Cartn_y", "Cartn_z") if c not in cols]
    if missing:
        # some writers only emit the label_* flavour
        alt = {"auth_atom_id": "label_atom_id", "auth_comp_id": "label_comp_id",
               "auth_asym_id": "label_asym_id", "auth_seq_id": "label_seq_id"}
        for m in list(missing):
            if alt.get(m) in cols:
                for r in rows:
                    r[m] = r[alt[m]]
                missing.remove(m)
    if missing:
        raise SystemExit("%s: _atom_site lacks %s" % (path, ", ".join(missing)))
    return rows


def dot(v, default=""):
    return default if v in (".", "?") else v


def atom_name_field(name, element):
    """PDB columns 13-16: a one-letter element starts in column 14."""
    if len(name) >= 4:
        return name[:4]
    if len(element) == 1 and len(name) < 4:
        return " %-3s" % name
    return "%-4s" % name


def convert(src, dst, model="1"):
    rows = read_atom_site(src)
    n = 0
    with open(dst, "w") as out:
        prev_chain = None
        for i, r in enumerate(rows, 1):
            if dot(r.get("pdbx_PDB_model_num", model), model) != model:
                continue
            chain = dot(r["auth_asym_id"], "A")[:1]
            if prev_chain is not None and chain != prev_chain:
                out.write("TER\n")
            prev_chain = chain
            elem = dot(r.get("type_symbol", ""))
            out.write("%-6s%5d %s%1s%3s %1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s\n"
                      % (r["group_PDB"], i,
                         atom_name_field(dot(r["auth_atom_id"]), elem),
                         dot(r.get("label_alt_id", "")),
                         dot(r["auth_comp_id"])[:3], chain,
                         int(float(r["auth_seq_id"])),
                         dot(r.get("pdbx_PDB_ins_code", "")),
                         float(r["Cartn_x"]), float(r["Cartn_y"]), float(r["Cartn_z"]),
                         1.00, float(dot(r.get("B_iso_or_equiv", "0"), "0")),
                         elem.upper()))
            n += 1
        out.write("TER\nEND\n")
    return n


def chain_report(path):
    rows = read_atom_site(path)
    seen = {}
    for r in rows:
        ch = dot(r["auth_asym_id"], "A")[:1]
        seen.setdefault(ch, set()).add(int(float(r["auth_seq_id"])))
    return {c: (len(v), min(v), max(v)) for c, v in seen.items()}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("src", help="a .cif or a directory of them")
    ap.add_argument("out", nargs="?", help="output .pdb (single input only)")
    ap.add_argument("--outdir")
    ap.add_argument("--model", default="1")
    a = ap.parse_args()

    if os.path.isdir(a.src):
        outdir = a.outdir or a.src
        os.makedirs(outdir, exist_ok=True)
        cifs = [f for f in sorted(os.listdir(a.src)) if f.endswith(".cif")]
        if not cifs:
            raise SystemExit("no .cif files in %s" % a.src)
        for f in cifs:
            src = os.path.join(a.src, f)
            dst = os.path.join(outdir, f[:-4] + ".pdb")
            n = convert(src, dst, a.model)
            ch = chain_report(src)
            print("%-62s %4d atoms  chains %s"
                  % (f[:62], n, ", ".join("%s:%d res %d-%d" % (c, v[0], v[1], v[2])
                                          for c, v in sorted(ch.items()))))
    else:
        dst = a.out or (a.src[:-4] + ".pdb" if a.src.endswith(".cif") else a.src + ".pdb")
        n = convert(a.src, dst, a.model)
        ch = chain_report(a.src)
        print("%s -> %s  %d atoms" % (a.src, dst, n))
        for c, v in sorted(ch.items()):
            print("  chain %s: %d residues, %d-%d" % (c, v[0], v[1], v[2]))


if __name__ == "__main__":
    main()
