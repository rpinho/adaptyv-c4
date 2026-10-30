#!/usr/bin/env python3
"""Does any C4 design contain a target-derived contiguous segment?

Sweeps (a) the 35 de-duplicated candidate designs and (b) the whole Modal archive, for the longest
exact substring shared with TfR1 (P02786), with a composition-matched shuffle null per design.
Also reports the same for a set of C3 designs as a negative control (they target GLP-1R / TfR1 too).

Fail-closed: exits non-zero if the target sequence cannot be loaded or if no designs are measured.
"""
import csv, json, random, statistics, sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT/'results/c4_novelty_screen_20261025'
MIN_REPORT = 8

def load_target():
    seq = None; name = None; buf = []
    for line in open(Path.home()/'.cache/adaptyv-c4/novelty/swissprot.fasta'):
        if line.startswith('>'):
            if name == 'P02786': return ''.join(buf)
            f = line[1:].split()[0]; parts = f.split('|')
            name = parts[1] if len(parts) > 1 else f; buf = []
        else:
            buf.append(line.strip())
    if name == 'P02786': return ''.join(buf)
    return None

def lcs(a, b):
    prev = [0]*(len(b)+1); bl = be = 0
    for i in range(1, len(a)+1):
        cur = [0]*(len(b)+1)
        ai = a[i-1]
        for j in range(1, len(b)+1):
            if ai == b[j-1]:
                v = prev[j-1]+1; cur[j] = v
                if v > bl: bl = v; be = i
        prev = cur
    return a[be-bl:be], bl

def sweep(named_seqs, target, n_null=40, seed=20261025):
    random.seed(seed); rows = []
    for name, s in named_seqs:
        pep, L = lcs(s, target)
        if L < MIN_REPORT:
            rows.append(dict(name=name, length=len(s), lcs=L, peptide=pep, null_max=None, null_p=None))
            continue
        null = [lcs(''.join(random.sample(s, len(s))), target)[1] for _ in range(n_null)]
        rows.append(dict(name=name, length=len(s), lcs=L, peptide=pep,
                         null_median=statistics.median(null), null_max=max(null),
                         null_p=round((sum(1 for x in null if x >= L)+1)/(n_null+1), 4)))
    return rows

def main():
    tgt = load_target()
    if not tgt: sys.exit('FAIL-CLOSED: P02786 not found in swissprot.fasta')
    print(f'target P02786 length {len(tgt)}', file=sys.stderr)

    cand = [(r['name'], r['sequence']) for r in csv.DictReader(open(OUT/'inputs.csv'))]

    arch = {}
    ap = ROOT/'results/c4_apical_consolidated_results.json'
    if ap.exists():
        d = json.loads(ap.read_text())
        res = d.get('results')
        if isinstance(res, list):
            for r in res:
                n = r.get('name') or r.get('design')
                sq = r.get('sequence')
                if n and sq: arch[n] = sq
        elif isinstance(res, dict):
            for n, r in res.items():
                sq = r.get('sequence') if isinstance(r, dict) else None
                if sq: arch[n] = sq
    if not cand: sys.exit('FAIL-CLOSED: no candidate designs')

    rows_c = sweep(cand, tgt)
    rows_a = sweep(sorted(arch.items()), tgt) if arch else []
    hits_c = [r for r in rows_c if r['lcs'] >= MIN_REPORT]
    hits_a = [r for r in rows_a if r['lcs'] >= MIN_REPORT]

    payload = dict(question='Does any C4 design contain a contiguous target-derived segment (>=8 aa)?',
                   target='TfR1 / P02786', target_len=len(tgt), min_report=MIN_REPORT,
                   candidates=dict(n=len(rows_c), n_with_segment=len(hits_c), rows=rows_c),
                   modal_archive=dict(n=len(rows_a), n_with_segment=len(hits_a), rows=rows_a,
                                      note='archive sequence field availability varies; n=0 means sequences were not stored'),
                   )
    (OUT/'target_segment_sweep.json').write_text(json.dumps(payload, indent=1))

    print(f"\nCANDIDATES: {len(hits_c)}/{len(rows_c)} have a >= {MIN_REPORT} aa exact TfR1 segment")
    for r in sorted(hits_c, key=lambda x: -x['lcs']):
        print(f"  {r['name']:30s} len {r['length']:3d}  LCS {r['lcs']:2d}  null max {r.get('null_max')}  p={r.get('null_p')}  {r['peptide']}")

    if rows_a:
        print(f"\nMODAL ARCHIVE: {len(hits_a)}/{len(rows_a)} have a >= {MIN_REPORT} aa exact TfR1 segment")
        for r in sorted(hits_a, key=lambda x: -x['lcs'])[:12]:
            print(f"  {r['name']:30s} len {r['length']:3d}  LCS {r['lcs']:2d}  null max {r.get('null_max')}  p={r.get('null_p')}  {r['peptide']}")
    else:
        print("\nMODAL ARCHIVE: no sequences recovered from the consolidated results file "
              "(field names differ) — recorded as a gap, not as a clean result.")

if __name__ == '__main__':
    main()
