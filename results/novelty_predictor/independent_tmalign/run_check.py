"""Compare standalone TM-align with preserved Foldseek scores, CPU only."""
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import time

HERE = Path(__file__).resolve().parent


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def ca_count(path):
    return sum(line.startswith('ATOM') and line[12:16].strip() == 'CA'
               for line in Path(path).read_text().splitlines())


def parse_output(text):
    scores = re.findall(r'TM-score=\s*([0-9.]+)\s*\(normalized by length of Structure_([12]):\s*L=(\d+)', text)
    if len(scores) != 2 or {s[1] for s in scores} != {'1', '2'}:
        raise ValueError('Unexpected standalone normalization schema')
    result = {f'structure_{index}': dict(tm_score=float(score), length=int(length))
              for score, index, length in scores}
    if any(not 0 <= s['tm_score'] <= 1 for s in result.values()):
        raise ValueError('Invalid TM-score')
    return result


def main():
    installation = json.loads((HERE/'INSTALLATION.json').read_text())
    binary = Path(installation['binary'])
    assert sha(binary) == installation['binary_sha256']
    inputs = json.loads((HERE/'PAIRS.json').read_text())
    pairs = inputs['pairs']
    reference_manifest = json.loads((HERE/'REFERENCE_MANIFEST.json').read_text())
    references = {r['database_key']:r for r in reference_manifest['records']}
    assert len(pairs) == 27
    for pair in pairs:
        ref = references[pair['database_key']]
        assert ref['target'] == pair['hit']['target']
        assert sha(ref['path']) == ref['sha256']
        pair['reference_path'] = ref['path']
        assert ca_count(pair['query_path']) == pair['hit']['qlen']
        assert ca_count(pair['reference_path']) == pair['hit']['tlen']
    logs = HERE/'alignments'; logs.mkdir(exist_ok=True)

    def run(item):
        i, pair = item
        cmd = [str(binary), pair['query_path'], pair['reference_path']]
        started = time.monotonic()
        completed = subprocess.run(cmd, capture_output=True, text=True, check=True, timeout=30,
            env={**os.environ, 'OMP_NUM_THREADS': '1', 'CUDA_VISIBLE_DEVICES': ''})
        log = logs/f'{i:03d}.txt'; log.write_text(completed.stdout)
        scores = parse_output(completed.stdout)
        assert scores['structure_1']['length'] == pair['hit']['qlen']
        assert scores['structure_2']['length'] == pair['hit']['tlen']
        qtm, ttm = scores['structure_1']['tm_score'], scores['structure_2']['tm_score']
        return dict(query_name=pair['query_name'], role=pair['role'], target=pair['hit']['target'],
            foldseek_query_TM=pair['hit']['qtmscore'], standalone_query_TM=qtm,
            foldseek_target_TM=pair['hit']['ttmscore'], standalone_target_TM=ttm,
            query_length=pair['hit']['qlen'], target_length=pair['hit']['tlen'],
            delta_query_TM=qtm-pair['hit']['qtmscore'],
            query_sha256=sha(pair['query_path']), reference_sha256=sha(pair['reference_path']),
            command=cmd, log_file=str(log.relative_to(HERE)), log_sha256=sha(log),
            elapsed_seconds=time.monotonic()-started)

    started = time.monotonic()
    with ThreadPoolExecutor(max_workers=4) as pool:
        results = list(pool.map(run, enumerate(pairs, 1)))
    control = next(r for r in results if r['role']=='native_positive_control')
    assert control['standalone_query_TM'] >= .98
    retired = sorted((r for r in results if r['role']=='remaining_local4'),
                     key=lambda r: (-r['standalone_query_TM'], r['target']))
    assert len(retired) == 26
    self_control = parse_output((HERE/'self_control.txt').read_text())
    assert self_control['structure_1']['tm_score'] == self_control['structure_2']['tm_score'] == 1
    old_max = max(r['foldseek_query_TM'] for r in retired)
    new_max = retired[0]['standalone_query_TM']
    category = lambda tm: 'high' if tm >= .8 else 'moderate' if tm >= .5 else 'low'
    same_class = category(old_max) == category(new_max)
    interpretation = ('repeatable_for_tested_matches' if same_class and abs(new_max-old_max)<=.01
                      else 'implementation_sensitive')
    summary = dict(status='completed', pairs=len(results), retired_reference_matches=len(retired),
        foldseek_best_query_TM=old_max, standalone_best_query_TM=new_max,
        standalone_best_target=retired[0]['target'], old_structure_class=category(old_max),
        standalone_structure_class=category(new_max), interpretation=interpretation,
        maximum_absolute_pair_delta=max(abs(r['delta_query_TM']) for r in retired),
        pairs_crossing_moderate_threshold=sum((r['foldseek_query_TM']>=.5)!=(r['standalone_query_TM']>=.5) for r in retired),
        positive_control=control, elapsed_seconds=time.monotonic()-started,
        self_control=self_control, reference_manifest_sha256=sha(HERE/'REFERENCE_MANIFEST.json'),
        binary_sha256=sha(binary), input_pairs_sha256=sha(HERE/'PAIRS.json'),
        script_sha256=sha(__file__), results=results, new_gpu_spend_usd=0,
        limitations=['Finite26 returned matches, no new database-wide discovery',
            'Same supplied coordinates, historical reference DB and whole-chain domain assumption',
            'Shared TM-align algorithm ancestry: software sensitivity check, not independent biological evidence',
            'Decoded CA coordinates rounded to PDB precision',
            'Exact organizer sequence/structure/domain policy remains unknown'])
    (HERE/'RESULTS.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps({k:v for k,v in summary.items() if k not in ('results','positive_control','limitations')}))


if __name__ == '__main__': main()
