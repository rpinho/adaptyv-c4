#!/usr/bin/env python
"""Limited local-reference similarity screen; not proof of de novo novelty.

Compare project sequences to explicitly available reference binders. The finite
reference list and local sequence metric cannot establish organizer novelty scores,
structural diversity or whether a design was made from scratch.

Reference binders are pulled straight out of the complex structures in data/: whichever
chain is NOT the target is a known binder, plus any sequence listed in EXTRA.

Usage:
  python novelty_audit.py designs.csv            # a Proteinbase-style CSV: name,sequence,...
  python novelty_audit.py designs.fasta --max-identity 0.30
"""
import argparse, csv, json, os, sys
from Bio.PDB import MMCIFParser
from Bio.PDB.Polypeptide import is_aa

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(HERE, "data")
T3 = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E','GLY':'G',
      'HIS':'H','ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S',
      'THR':'T','TRP':'W','TYR':'Y','VAL':'V'}

# complex structures whose non-target chains are known binders of our target
COMPLEXES = {
    "6wrv": "designed binder 3DS18",
    "6wrw": "designed binder 2DS25.5",
    "6wrx": "designed binder 2DS25.1",
    "3kas": "Machupo virus GP1",
    "9gh7": "bicyclic peptide",
    "1de4": "HFE / beta2-microglobulin",
    "3s9l": "transferrin",
}
EXTRA = {}  # name -> sequence, for references with no structure here

# a project can override the reference set with data/known_binders.json:
#   {"complexes": {"6x18": "GLP-1 peptide"}, "extra": {"exendin-4": "HGEGTFTSDLSK..."}}
_cfg = os.path.join(DATA, "known_binders.json")
if os.path.exists(_cfg):
    _d = json.load(open(_cfg))
    COMPLEXES = _d.get("complexes", COMPLEXES)
    EXTRA = _d.get("extra", EXTRA)


def chain_seqs(code):
    path = os.path.join(DATA, code + ".cif")
    if not os.path.exists(path):
        return {}
    m = MMCIFParser(QUIET=True).get_structure(code, path)[0]
    out = {}
    for c in m:
        rs = [r for r in c if is_aa(r, standard=True)]
        if len(rs) >= 8:
            out[c.id] = "".join(T3.get(r.get_resname(), "X") for r in rs)
    return out


def smith_waterman(a, b, match=2, mismatch=-1, gap=-2):
    """Return (best score, identities in the best local alignment, aligned length)."""
    n, k = len(a), len(b)
    prev = [0] * (k + 1)
    prev_id = [0] * (k + 1)
    prev_len = [0] * (k + 1)
    best = (0, 0, 0)
    for i in range(1, n + 1):
        cur = [0] * (k + 1); cur_id = [0] * (k + 1); cur_len = [0] * (k + 1)
        ai = a[i - 1]
        for j in range(1, k + 1):
            eq = ai == b[j - 1]
            d = prev[j - 1] + (match if eq else mismatch)
            u = prev[j] + gap
            l = cur[j - 1] + gap
            m = max(0, d, u, l)
            cur[j] = m
            if m == 0:
                cur_id[j] = cur_len[j] = 0
            elif m == d:
                cur_id[j] = prev_id[j - 1] + (1 if eq else 0); cur_len[j] = prev_len[j - 1] + 1
            elif m == u:
                cur_id[j] = prev_id[j]; cur_len[j] = prev_len[j] + 1
            else:
                cur_id[j] = cur_id[j - 1]; cur_len[j] = cur_len[j - 1] + 1
            if m > best[0]:
                best = (m, cur_id[j], cur_len[j])
        prev, prev_id, prev_len = cur, cur_id, cur_len
    return best


def load_designs(path):
    designs = []
    if path.endswith((".fasta", ".fa")):
        name, buf = None, []
        for line in open(path):
            line = line.strip()
            if line.startswith(">"):
                if name: designs.append((name, "".join(buf)))
                name, buf = line[1:].split()[0], []
            elif line:
                buf.append(line)
        if name: designs.append((name, "".join(buf)))
    else:
        with open(path) as fh:
            for row in csv.DictReader(fh):
                seq = (row.get("sequence") or "").strip().replace(":", "")
                if seq:
                    designs.append((row.get("name", "?"), seq))
    return designs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("designs")
    ap.add_argument("--target-chain-hint", default=None,
                    help="chain id of the TARGET in the complexes (default: the longest chain)")
    ap.add_argument("--max-identity", type=float, default=0.30,
                    help="flag any design whose best local identity to a known binder exceeds this")
    ap.add_argument("--min-aln", type=int, default=20,
                    help="ignore local alignments shorter than this; a 4-residue perfect match is noise, not homology")
    ap.add_argument("--json", help="write the full table here")
    ap.add_argument("--complex-ids",help="Explicit comma-separated subset of the configured PDB reference set")
    a = ap.parse_args()

    refs = []
    for code, label in COMPLEXES.items():
        if a.complex_ids and code not in a.complex_ids.lower().split(','):continue
        seqs = chain_seqs(code)
        if not seqs:
            print("note: %s.cif not in data/, skipped" % code)
            continue
        longest = max(seqs.values(), key=len)
        for cid, s in seqs.items():
            if s is longest or len(s) > 500:   # the receptor itself, not a binder
                continue
            if all(s != r[2] for r in refs):
                refs.append((code.upper(), "%s (chain %s)" % (label, cid), s))
    for name, s in EXTRA.items():
        refs.append(("-", name, s))
    print("reference binders: %d" % len(refs))
    for pdb, label, s in refs:
        print("  %-6s %-34s %3d aa" % (pdb, label, len(s)))

    designs = load_designs(a.designs)
    print("\ndesigns: %d (identity is judged only over alignments of at least %d residues)\n"
          % (len(designs), a.min_aln))
    rows, flagged = [], 0
    print("%-28s %-6s %-34s %s" % ("design", "ident", "closest known binder", "aln"))
    for name, seq in designs:
        # pick the alignment with the highest SCORE (a long, good alignment), then report its
        # identity; a short perfect match scores low and is correctly ignored
        best = (0.0, None, 0, 0)
        for pdb, label, ref in refs:
            score, ident, alen = smith_waterman(seq, ref)
            if alen < a.min_aln:
                continue
            if score > best[3]:
                best = (ident / alen if alen else 0.0, "%s %s" % (pdb, label), alen, score)
        rows.append(dict(design=name, length=len(seq), max_identity=round(best[0], 3),
                         closest=best[1], aligned_length=best[2], sw_score=best[3],
                         flagged=bool(best[0] > a.max_identity)))
        if best[0] > a.max_identity:
            flagged += 1
        print("%-28s %-6.2f %-34s %d aa, score %d"
              % (name[:28], best[0], (best[1] or "-")[:34], best[2], best[3]))

    print("\n%d of %d designs exceed %.0f%% local identity to a known binder"
          % (flagged, len(rows), a.max_identity * 100))
    if flagged == 0:
        print("No flags against this limited reference set; external/structural novelty remains unestablished")
    if a.json:
        json.dump(dict(threshold=a.max_identity,organizer_novelty_score=False,
                       scope='Limited local-sequence reference screen; highest-scoring alignment identity, not maximum identity over every alignment; no de novo/structural novelty proof',
                       references=[{"pdb": p, "label": l, "length": len(s)}
                                                             for p, l, s in refs], designs=rows),
                  open(a.json, "w"), indent=1)
        print("wrote %s" % a.json)


if __name__ == "__main__":
    main()
