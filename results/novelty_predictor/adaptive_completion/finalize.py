"""Validate completed native alignments and publish conservative mixed evidence.

No design, coordinate, scoring threshold or sequence-search result is changed.
This runs once after the detached native worker has exported its actual results.
"""
import collections
import copy
import csv
from datetime import datetime, timezone
import hashlib
import itertools
import json
from pathlib import Path
import time

from native import HERE, ROOT, BASE, n, stream_best


def completion_bounds(length, domains, best, complete):
    """Evaluate all threshold-relevant completions of unreturned pair scores.

    Zero here is an internal lower bound for missing evidence, never a reported
    measured no-hit score. A completed domain must have a real observed score.
    """
    choices = []
    matches = []
    observed = False
    for domain in domains:
        query = domain['query_id']
        hit = best.get(query)
        observed |= hit is not None
        score = hit['qtmscore'] if hit else 0.0
        if complete.get(query, False) and hit is None:
            raise ValueError('A fully scored nonempty candidate set lacks evidence')
        options = {score}
        if not complete.get(query, False):
            options |= {x for x in (n.STRUCTURE_MODERATE_TM, n.STRUCTURE_HIGH_TM, 1.0) if x >= score}
        choices.append(sorted(options))
        matches.append(dict(**n.domain_coordinates(domain), tm_score=score,
                            normalization='query_domain_length'))
    lower = n.structural_class(length, matches)
    categories = set()
    for scores in itertools.product(*choices):
        current = [dict(**n.domain_coordinates(d), tm_score=s,
                        normalization='query_domain_length') for d, s in zip(domains, scores)]
        categories.add(n.structural_class(length, current)['category'])
    return lower, sorted(categories, key=('low', 'moderate', 'high').index), observed


def refine(row, top, full_counts, expected):
    domains = row['structure_evidence']['domains']
    complete = {}
    best = {}
    details = {}
    for domain in domains:
        query = domain['query_id']
        returned = full_counts[query]
        requested = expected[query]['unique_union_targets']
        if returned > requested:
            raise ValueError('More exported rows than unique requested native pairs')
        complete[query] = requested > 0 and returned == requested
        hits = top[query]
        best[query] = hits[0] if hits else None
        domain.update(best_hit=best[query], top_hits=hits, hit_count=None,
                      hit_count_semantics='Not materialized; top ten unique observed targets retained',
                      full_union_requested_pairs=requested, full_union_exported_rows=returned,
                      full_union_complete=complete[query])
        details[query] = dict(requested_pairs=requested, exported_rows=returned,
                              complete=complete[query], observed_best_tm=(best[query] or {}).get('qtmscore'))
    lower, categories, observed = completion_bounds(row['length'], domains, best, complete)
    seq = row['sequence_evidence_class']
    if seq == 'no_significant_hit':
        seq = 'low'
    levels = sorted({n.general_level(seq, category) for category in categories})
    all_complete = all(complete.values())
    # Incomplete work can prove a fixed point through positive resemblance, but
    # cannot manufacture local3/4 solely because accepted domains cover little.
    can_state = observed and (all_complete or lower['category'] in ('moderate', 'high'))
    point = levels[0] if len(levels) == 1 and can_state else None
    coverage = sum(len(n.domain_positions(d, row['length'])) for d in domains) / row['length']
    status = 'full_declared_candidate_union' if all_complete else 'incomplete_declared_candidate_union'
    similarity = dict(lower)
    similarity.update(category=categories[0] if len(categories) == 1 else 'unestablished',
                      lower_bound_category=lower['category'], possible_categories=categories,
                      domain_segmentation_coverage=coverage, evidence_status=status,
                      tm_boundary_precision_caveat=any(
                          h and min(abs(h['qtmscore']-x) for x in (.5, .8)) <= .0001 for h in best.values()),
                      coverage_limited_low_category=(all_complete and lower['category'] == 'low' and coverage <= .7))
    row.update(structural_similarity=similarity, predicted_level=point, local_predicted_level=point,
               structure_status=status, local_level_interval_with_unfinished_alignments=levels,
               possible_levels_if_structure_policy_accepted=sorted({
                   level for category in categories for level in n.possible_levels(seq, row['sequence_database_gaps'], category)}),
               level_status='provisional_local_prediction' if point is not None else 'unestablished_incomplete_structural_evidence',
               local_gate_assessment='provisional_local_level_'+str(point) if point is not None else 'unestablished',
               organizer_level_established=False,
               prediction_confidence='low: supplied generator coordinates, finite discovery, historical reference and domain/typing/corpus gaps',
               explanation='Fixed general-scale rule on actual observed domain scores. Missing pair outputs retain completion bounds. A coverage-limited low structural category is not proof of a new fold.',
               adaptive_completion=dict(scope=status, domains=details, all_requested_pairs_exported=all_complete))
    return status


def main():
    start = time.monotonic()
    if (HERE/'predictions.json').exists():
        raise FileExistsError('The evidence report is already frozen')
    status_bytes = (HERE/'STATUS.json').read_bytes()
    status = json.loads(status_bytes)
    if status['state'] not in ('alignments_complete', 'alignments_complete_report_pending'):
        raise ValueError('Only a completed native export can be finalized')
    union_path = HERE/'CANDIDATE_UNION.json'
    union = json.loads(union_path.read_text())
    assert n.digest(union_path) == status['candidate_union_sha256']
    assert n.digest(BASE) == union['base_report_sha256']
    assert n.digest(union['query_domains']) == union['query_domains_sha256']
    assert n.digest(union['seed_tsv']) == union['seed_tsv_sha256']
    full = Path(status['raw_hits'])
    assert n.digest(full) == status['raw_hits_sha256']
    checkpoint = HERE/'ALIGNMENT_CHECKPOINT.json'
    if checkpoint.exists():
        assert checkpoint.read_bytes() == status_bytes
    else:
        checkpoint.write_bytes(status_bytes)
    base = json.loads(BASE.read_text())
    designs = n.load_designs(BASE.parents[2]/'analysis_inputs.csv')
    assert n.digest(BASE.parents[2]/'analysis_inputs.csv') == union['input_csv_sha256']
    assert [r['name'] for r in designs] == [r['name'] for r in base['designs']]
    domains = json.loads(Path(union['query_domains']).read_text())
    previous = BASE.parent/'structure/hits.tsv'
    # The bounded parent-context report's receipt identifies its actual export.
    old_search = base['structural_search']
    assert n.digest(previous) == old_search['raw_hits_sha256']
    top, counts = stream_best([previous, full], domains)
    full_counts = counts[str(full)]
    expected = union['per_domain']
    assert set(full_counts) <= set(expected)
    report = copy.deepcopy(base)
    scopes = {}
    preserved = []
    for original, row, design in zip(base['designs'], report['designs'], designs):
        assert row['sequence_sha256'] == hashlib.sha256(design['sequence'].encode()).hexdigest()
        ds = row.get('structure_evidence', {}).get('domains', [])
        if original.get('local_predicted_level') is not None:
            assert row == original
            preserved.append(row['name'])
            scope = 'point_proven_by_existing_positive_evidence_best_TM_not_exhaustive'
        elif not ds:
            assert row['predicted_level'] is None
            scope = 'no_accepted_domain_unestablished'
        else:
            assert all(d['query_id'] in expected for d in ds)
            scope = refine(row, top, full_counts, expected)
        scopes[row['name']] = scope
    assert len(preserved) == union['preserved_proven_point_inputs']
    assert len(report['designs']) == 400
    protected = json.loads((ROOT/'results/novelty_predictor/BOUNDARIES.json').read_text())['protected_local_sha256']
    assert all(n.digest(ROOT/path) == digest for path, digest in protected.items())
    raw_sources = [dict(path=str(path), sha256=n.digest(path), raw_rows=sum(counts[str(path)].values()))
                   for path in (previous, full)]
    report.update(status='adaptive_structural_completion_report', generated_at=datetime.now(timezone.utc).isoformat(),
                  predictor_sha256=n.digest(n.__file__), organizer_score=False,
                  prior_structural_search=old_search,
                  structural_search=dict(mode='mixed_preserved_positive_points_and_full_remaining_candidate_union',
                                         raw_evidence_sources=raw_sources, query_domains=len(domains),
                                         finite_discovery_not_global_reference_exhaustiveness=True),
                  adaptive_completion=dict(plan_sha256=n.digest(HERE/'PLAN.md'),
                                           alignment_checkpoint=str(checkpoint), alignment_checkpoint_sha256=n.digest(checkpoint),
                                           candidate_union_sha256=n.digest(union_path), script_sha256=n.digest(__file__),
                                           rule_compatibility=status['rule_compatibility'], scopes=scopes,
                                           known_points_preserved_unchanged=preserved,
                                           unique_requested_pairs=union['unique_candidate_pairs'],
                                           returned_raw_rows=sum(full_counts.values())), new_gpu_spend_usd=0)
    supplemental = HERE.parent/'antibody_cdr_reference/PREDICTIONS.json'
    if supplemental.exists():
        report['supplemental_antibody_evidence'] = dict(path=str(supplemental), sha256=n.digest(supplemental),
            scope='Conditional reference-loop evidence; organizer structural typer and VH/VHH status unestablished',
            used_to_change_general_scale_levels=False, organizer_scoring_route_established=False)
    report['limitations'] += [
        'Mixed evidence: preserved positive-bound points have no exhaustive absolute best TM; scopes are explicit per input.',
        'Full remaining comparison covers the union of completed domain and whole-parent retrieval, not every reference pair.',
        'Supplied generator coordinates are not independent folds. Domain and antibody-typing policies remain unverified.',
        'Historical PDB snapshot, finite retrieval, missing organizer patent corpus and score-boundary precision remain gaps.',
        'Separate conditional CDR evidence cannot identify the organizer scoring route. No assay or entry qualification is established.']
    n.save_report(report, HERE)
    with (HERE/'EVIDENCE.csv').open('w', newline='') as stream:
        fields = ['name', 'source_arm', 'local_level', 'local_possible_levels', 'scope', 'organizer_level_established']
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        for row, design in zip(report['designs'], designs):
            writer.writerow(dict(name=row['name'], source_arm=design['source_arm'], local_level=row.get('local_predicted_level'),
                                 local_possible_levels='|'.join(map(str, row.get('local_level_interval_with_unfinished_alignments', []))),
                                 scope=scopes[row['name']], organizer_level_established=False))
    saved = json.loads((HERE/'predictions.json').read_text())
    with (HERE/'predictions.csv').open(newline='') as stream:
        csv_rows = list(csv.DictReader(stream))
    assert [r['name'] for r in csv_rows] == [r['name'] for r in saved['designs']]
    assert all((str(r['local_predicted_level']) if r['local_predicted_level'] is not None else '') == c['local_predicted_level']
               for r, c in zip(saved['designs'], csv_rows))
    for old, new in zip(base['designs'], saved['designs']):
        if old['name'] in preserved:
            assert old == new
    by_arm = {}
    for arm in sorted({r['source_arm'] for r in designs}):
        by_arm[arm] = dict(collections.Counter(str(r.get('local_predicted_level')) for r, d in zip(saved['designs'], designs) if d['source_arm'] == arm))
    validation = dict(state='complete', designs=400, by_arm=by_arm,
                      local_level_counts=dict(collections.Counter(str(r.get('local_predicted_level')) for r in saved['designs'])),
                      scope_counts=dict(collections.Counter(scopes.values())), preserved_points=len(preserved),
                      requested_native_pairs=union['unique_candidate_pairs'], exported_native_rows=sum(full_counts.values()),
                      incomplete_domains=[q for q, d in expected.items() if d['unique_union_targets'] == 0 or full_counts[q] != d['unique_union_targets']],
                      no_retrieved_candidate_domains=[q for q, d in expected.items() if d['unique_union_targets'] == 0],
                      organizer_levels_established=0, qualified_designs=0, protected_files_unchanged=True,
                      json_csv_input_order_verified=True, predictions_sha256=n.digest(HERE/'predictions.json'),
                      evidence_csv_sha256=n.digest(HERE/'EVIDENCE.csv'), finalizer_elapsed_seconds=time.monotonic()-start,
                      new_gpu_or_remote_cpu_spend_usd=0)
    n.write_json(HERE/'VALIDATION.json', validation)
    lines = ['# Completed CPU structural evidence', '',
             'Measured the remaining declared candidate union using full native TM-align, exact score export and unchanged general-scale rules.', '',
             '| Arm | Local1 | Local2 | Local3 | Local4 | Unestablished |',
             '| --- | ---: | ---: | ---: | ---: | ---: |']
    for arm, c in by_arm.items():
        lines.append('| '+arm+' | '+' | '.join(str(c.get(k, 0)) for k in ('1', '2', '3', '4', 'None'))+' |')
    lines += ['', 'Requested native pairs: '+str(union['unique_candidate_pairs'])+'; actual exported rows: '+str(sum(full_counts.values()))+'.',
              'Preserved143 positive-bound point records byte-for-value; completed remaining comparisons use all domain/whole-parent retrieved targets, deduplicated before full alignment.',
              'Incomplete domains: '+str(len(validation['incomplete_domains']))+'. Missing output retains conservative completion bounds; no no-hit is recorded as measured TM0.', '',
              'These are provisional general-scale calls on supplied generator coordinates, not organizer scores. Historical/finite reference retrieval, domain policy, antibody typing and missing organizer patent coverage remain unresolved. Coverage-limited local4 is not proof of a new fold. No binding, pH selectivity, assay success or entry qualification measured.', '',
              'See predictions.json, EVIDENCE.csv, VALIDATION.json, immutable ALIGNMENT_CHECKPOINT.json and PLAN.md for evidence and scopes.', '',
              'Recorded own-lane spend: Runpod $56.589988 + Modal estimate $10.718947 = approximately $67.31; new GPU/cloud compute $0.00. Separate generation/account-wide billing excluded.', '']
    (HERE/'RESULTS.md').write_text('\n'.join(lines))
    print(json.dumps(validation))


if __name__ == '__main__':
    main()
