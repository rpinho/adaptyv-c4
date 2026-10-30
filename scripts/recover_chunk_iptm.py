"""Recover ipTM data from chunk log files where JSON output was overwritten.

Parses the summary section at end of each log (name: ipTM=X.XXX format)
and also captures error designs.
"""
import json
import re
import sys

def parse_chunk_log(log_path, refold_json_path):
    """Parse a chunk log to extract per-design ipTM from the summary section."""
    with open(refold_json_path) as f:
        refold = json.load(f)
    refold_map = {r['name']: r for r in refold}

    with open(log_path) as f:
        log = f.read()

    results = []
    for line in log.split('\n'):
        m = re.match(r'\s+(\S+): ipTM=([\d.]+), complex_His_elevated=(\d+)', line)
        if m:
            name = m.group(1)
            iptm = float(m.group(2))
            n_elevated = int(m.group(3))

            r = refold_map.get(name, {})
            results.append({
                'name': name,
                'sequence': r.get('sequence', ''),
                'iptm': iptm,
                'n_elevated_complex': n_elevated,
            })

    # Also check for total count from "Saved N results" line
    m = re.search(r'Saved (\d+) results', log)
    total = int(m.group(1)) if m else '?'

    # Check for low-ipTM designs not in summary (ipTM<=0.5)
    m2 = re.search(r'(\d+)/(\d+) with ipTM>0\.5', log)
    if m2:
        n_good = int(m2.group(1))
        n_total = int(m2.group(2))
        n_low = n_total - n_good
    else:
        n_low = 0

    return results, total, n_low


if __name__ == '__main__':
    base = '.'
    for chunk_n in [1, 2, 4]:
        log_path = f'/tmp/chunk{chunk_n}.log'
        refold_path = f'{base}/results/c4_acid_chunk{chunk_n}_refold.json'
        out_path = f'{base}/results/c4_acid_boltz2_chunk{chunk_n}.json'

        try:
            results, total, n_low = parse_chunk_log(log_path, refold_path)
            with open(out_path, 'w') as f:
                json.dump(results, f, indent=2)
            print(f'Chunk {chunk_n}: recovered {len(results)} designs with ipTM>0.5 '
                  f'(total={total}, {n_low} low-ipTM excluded) -> {out_path}')
        except FileNotFoundError as e:
            print(f'Chunk {chunk_n}: {e}')
