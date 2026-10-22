"""Validate the full CPU follow-up against frozen earlier evidence."""
from collections import Counter
import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
BASE = HERE.parent / 'structural/run/predictions.json'


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def best_tm(row):
    matches = [d['best_hit']['qtmscore'] for d in
               row.get('structure_evidence', {}).get('domains', []) if d.get('best_hit')]
    return max(matches) if matches else None


def main():
    baseline = json.loads(BASE.read_text())
    path = HERE / 'full_pool_union/predictions.json'
    current = json.loads(path.read_text())
    previous = {r['name']: r for r in baseline['designs']}
    assert current['input_csv_sha256'] == baseline['input_csv_sha256']
    assert current['structure_manifest_sha256'] == baseline['structure_manifest_sha256']
    assert current['structure_database']['files'] == baseline['structure_database']['files']
    assert current['structural_search']['prefilter_modes'] == [0, 1]
    assert current['structural_search']['query_domains'] == 68
    comparisons = []
    for r in current['designs']:
        old = previous[r['name']]
        a, b = best_tm(old), best_tm(r)
        if a is not None: assert b is not None and b >= a-0.0001, r['name']
        assert not r.get('organizer_level_established', False)
        comparisons.append(dict(name=r['name'], old_local_level=old.get('local_predicted_level'),
            new_local_level=r.get('local_predicted_level'), old_best_TM=a, new_best_TM=b,
            novelty_downgraded=(old.get('local_predicted_level') is not None and
                r.get('local_predicted_level') is not None and r['local_predicted_level'] < old['local_predicted_level'])))
    retired = [r for r in comparisons if r['name'] != 'KNOWN_TFR1_POSITIVE_CONTROL']
    control = next(r for r in comparisons if r['name'] == 'KNOWN_TFR1_POSITIVE_CONTROL')
    assert control['new_best_TM'] >= 0.98 and control['new_local_level'] == 1
    result = dict(baseline_report_sha256=digest(BASE), current_report_sha256=digest(path),
        compared_retired=len(retired), before_counts=dict(Counter(str(r['old_local_level']) for r in retired)),
        after_counts=dict(Counter(str(r['new_local_level']) for r in retired)),
        downgraded=sum(r['novelty_downgraded'] for r in retired), all_previous_best_TMs_retained=True,
        elapsed_seconds=current['structural_search']['elapsed_seconds'],
        unique_query_target_matches=current['structural_search']['hit_rows'],
        raw_alignment_rows=current['structural_search']['raw_alignment_rows'],
        comparisons=comparisons, new_gpu_spend_usd=0)
    (HERE / 'FULL_POOL_COMPARISON.json').write_text(json.dumps(result, indent=2)+'\n')
    lines = ['# CPU follow-up: candidate discovery sensitivity', '',
        'The earlier single-discovery structural calls were too optimistic for some retired designs. '
        'A focused eight-query test found two of four local level4 calls drop to3, and one local3 drop to2. '
        'This is a measured change between local methods, not an organizer false-positive rate.', '',
        'The predictor now defaults to the union of Foldseek kmer0 and all-reference ungapped1 discovery, '
        'with50000 candidates per mode/query. mergedbs combines candidates for the same query database, '
        'then the existing full TM-align stage runs. Repeated target hits retain their highest query-domain TM '
        'and do not count as additional evidence. The old coordinates, domain intervals, thresholds and '
        'historical1,562,678-chain public PDB snapshot are unchanged.', '',
        f"Full follow-up:68 coordinate-bearing queries in{result['elapsed_seconds']:.3f} CPU seconds; "
        f"{result['unique_query_target_matches']} unique query/target matches from{result['raw_alignment_rows']} raw alignments. "
        'Every previous best query-normalized TM was retained or improved. All74 input sequences are accounted for; '
        'the six missing exact coordinate models remain unestablished.', '',
        '| Local level, retired designs | Before | After |', '|---|---:|---:|']
    for level in ('2', '3', '4', 'None'):
        lines.append(f"| {'Unestablished' if level == 'None' else level} | {result['before_counts'].get(level,0)} | {result['after_counts'].get(level,0)} |")
    lines += ['', f"{result['downgraded']} retired local calls were downgraded. Positive control remains local1 "
        f"with query-normalized TM{control['new_best_TM']:.4f}.", '',
        '| Changed retired design | Before TM / level | After TM / level |', '|---|---|---|']
    for r in retired:
        if r['novelty_downgraded']:
            before = 'no hit' if r['old_best_TM'] is None else f"{r['old_best_TM']:.4f}"
            lines.append(f"| {r['name']} | {before} / {r['old_local_level']} | {r['new_best_TM']:.4f} / {r['new_local_level']} |")
    lines += ['', 'No significant sequence hit is still provisionally treated as low sequence similarity. '
        'Exact organizer levels established:0. Missing patent sequences, older public snapshots, supplied '
        'complex-predicted coordinates instead of organizer monomer predictions, assumed whole-chain domains '
        'instead of the unnamed three-predictor consensus, and finite search sensitivity remain unresolved. '
        'The remaining local4 calls are not guaranteed passes. Nothing reinstates the retired pool or changes '
        'the separately vindicated native-context clash gate.', '',
        'A broad exploratory CA-spacing diagnostic found0 flags in68 structures (adjacent distances2.5–4.5 Å). '
        'This is not a peptide/sidechain/folding or functional validator. GEOMETRY.json preserves every measurement.', '',
        'Selection and CPU limits were recorded before each follow-up in PLAN.json and FULL_POOL_PLAN.json. '
        'The eight-query comparison is in COMPARISON.json; the full comparison, predictions, raw alignments, '
        'coordinate fingerprints, parameters and logs are preserved. Fifteen focused unit tests pass, including '
        'duplicate-target handling, domain coverage, normalization, sequence units and threshold equality.', '',
        'Use the same predictor/coordinate manifests as the previous structural report. Default discovery is '
        'now union, cap50000. Explicit --structure-prefilter-mode0 or1 and --structure-max-seqs remain '
        'available for sensitivity checks; higher sequence/structural similarity can lower predicted novelty. '
        'Use fresh output directories. No new inference, generator, exporter or account/submission work occurred.', '',
        'All stages ran locally,4 CPU threads,4G search memory,300seconds per stage. No GPU, Modal, '
        'Runpod or cloud CPU service was touched. New GPU spend$0.00. Recorded Runpod$56.589988 '
        '+our Modal estimate$10.718947 =~$67.31; separate RFD4 lane, storage/builds/consultations excluded.']
    (HERE / 'RESULTS.md').write_text('\n'.join(lines)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k != 'comparisons'}))


if __name__ == '__main__': main()
