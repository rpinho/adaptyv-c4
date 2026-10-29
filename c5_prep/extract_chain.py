#!/usr/bin/env python3
"""Extract one chain from an experimental PDB/mmCIF into a clean PDB (CPU only)."""
import argparse, pathlib, sys
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import template_lib as T


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("source")
    ap.add_argument("--chain", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--allow-predicted", action="store_true",
                    help="only for parser smoke tests; never for template selection")
    a = ap.parse_args()
    src = pathlib.Path(a.source)
    atoms, meta = T.parse_structure(src)
    if not a.allow_predicted:
        T.reject_if_predicted(src, meta, atoms)
    ch = T.chains(atoms)
    if a.chain not in ch:
        raise SystemExit(f"chain {a.chain!r} absent; available {sorted(ch)}")
    residues = ch[a.chain]
    out = pathlib.Path(a.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    serial = 1
    with out.open("w") as fh:
        for r in residues:
            for name in sorted(r["atoms"], key=lambda n: (n != "N", n != "CA", n != "C", n != "O", n)):
                x, y, z = r["atoms"][name]
                elem = name[0] if not name[0].isdigit() else name[1]
                aname = name if len(name) == 4 else f" {name:<3}"
                fh.write(f"ATOM  {serial:5d} {aname}{' '}{r['resn'][:3]:>3} {a.chain}"
                         f"{r['resseq']:4d}{r['icode'] or ' '}   "
                         f"{x:8.3f}{y:8.3f}{z:8.3f}{1.0:6.2f}{r['bfac'].get(name, 0.0):6.2f}"
                         f"          {elem:>2}\n")
                serial += 1
        fh.write("TER\nEND\n")
    seq = T.chain_sequence(residues)
    print(f"{src.name}:{a.chain} -> {out}  residues={len(residues)} seq_sha={T.hashlib.sha256(seq.encode()).hexdigest()[:12]}")


if __name__ == "__main__":
    main()
