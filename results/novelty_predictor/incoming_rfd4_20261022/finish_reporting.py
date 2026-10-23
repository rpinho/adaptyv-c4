"""Validate and summarize this fixed, unranked CPU novelty batch after completion."""
import collections
import argparse
import copy
import csv
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT / 'scripts'))
import novelty_level_predictor as n


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--workflow-dir',type=Path,default=HERE/'run')
    args=parser.parse_args();run=args.workflow_dir.resolve()
    state = json.loads((run / 'STATUS.json').read_text())
    if state['state'] != 'complete':
        raise RuntimeError('The CPU workflow has not completed: ' + state['state'])
    source = list(csv.DictReader((HERE / 'analysis_inputs.csv').open()))
    path = Path(state['final_report'])
    assert n.digest(path) == state['final_report_sha256']
    report = json.loads(path.read_text())
    assert report['input_csv_sha256'] == n.digest(HERE / 'analysis_inputs.csv')
    assert len(source) == len(report['designs']) == 400
    assert all(s['completed'] for s in state['steps'])
    assert report['organizer_score'] is False
    exported = list(csv.DictReader(path.with_suffix('.csv').open()))
    assert len(exported) == len(source)
    groups = {}
    rows = []
    for original, design, flat in zip(source, report['designs'], exported):
        assert original['name'] == design['name'] == flat['name']
        assert design['length'] == len(original['sequence'])
        assert design['sequence_sha256'] == hashlib.sha256(original['sequence'].encode()).hexdigest()
        assert design['organizer_level_established'] is False
        level = design.get('local_predicted_level')
        assert design.get('predicted_level') == level
        assert flat['local_predicted_level'] == (str(level) if level is not None else '')
        best = design['best_hit'] or {}
        if best:
            assert best['fident'] == best['nident'] / best['alnlen']
        similarity = (design.get('structural_similarity') or {}) if level is not None else {}
        domains = design.get('structure_evidence', {}).get('domains', [])
        if level is not None:
            measured = n.structural_class(design['length'], [dict(
                **n.domain_coordinates(d), tm_score=(d.get('best_hit') or {}).get('qtmscore', 0.),
                normalization='query_domain_length') for d in domains])
            positive=report.get('structural_search',{}).get('positive_evidence_only',False)
            for key, value in measured.items():
                assert similarity['lower_bound_category' if positive and key=='category' else key] == value
            seq_class = design['sequence_evidence_class']
            if positive:
                checked=copy.deepcopy(design)
                checked['structural_similarity']=dict(measured,domain_segmentation_coverage=similarity['domain_segmentation_coverage'])
                n.apply_positive_evidence_bounds(checked)
                assert checked['local_predicted_level']==level
                assert checked['local_level_interval_with_unfinished_alignments']==design['local_level_interval_with_unfinished_alignments']
            else:
                assert level == n.general_level('low' if seq_class == 'no_significant_hit' else seq_class, measured['category'])
        arm = original['source_arm']
        group = groups.setdefault(arm, dict(inputs=0, local_levels=collections.Counter(),
            sequence_classes=collections.Counter(), structural_classes=collections.Counter(),
            best_sequence_hit_databases=collections.Counter()))
        group['inputs'] += 1
        group['local_levels'][str(level) if level is not None else 'unestablished'] += 1
        group['sequence_classes'][design['sequence_evidence_class']] += 1
        group['structural_classes'][similarity.get('category', 'unestablished')] += 1
        group['best_sequence_hit_databases'][best.get('database', 'no_significant_hit')] += 1
        rows.append(dict(name=design['name'], source_arm=arm, length=design['length'],
            provisional_general_level=level, sequence_class=design['sequence_evidence_class'],
            unfinished_alignment_level_interval='|'.join(map(str,design.get('local_level_interval_with_unfinished_alignments',[]))),
            antibody_scoring_mode='not_established',
            antibody_database_in_top_identity_hits=any(h['database'] in ('plabdab', 'therapeutic_antibody') for h in design['top_identity_hits']),
            identity=design['best_detected_identity'], query_coverage=best.get('qcov'),
            sequence_alignment_length=best.get('alnlen'), sequence_matches=best.get('nident'),
            sequence_database=best.get('database'), sequence_target=best.get('target'),
            structural_class=similarity.get('category'), high_coverage=similarity.get('high_coverage'),
            moderate_coverage=similarity.get('moderate_coverage'),
            domain_coverage=similarity.get('domain_segmentation_coverage'),
            domains=len(domains), best_observed_domain_tm=max(((d.get('best_hit') or {}).get('qtmscore', 0.) for d in domains), default=None),
            organizer_level_established=False))
    searches = report['database_searches']
    assert len(searches) == 5 and all(s['search_status'] == 'completed' for s in searches)
    for search in searches:
        assert n.digest(search.get('raw_hits_path',run/'whole_chain'/(search['id']+'.hits.tsv'))) == search['raw_hits_sha256']
    protected = json.loads((ROOT / 'results/novelty_predictor/BOUNDARIES.json').read_text())['protected_local_sha256']
    assert all(n.digest(ROOT / path) == expected for path, expected in protected.items())
    with (HERE / 'EVIDENCE.csv').open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0])); writer.writeheader(); writer.writerows(rows)
    summary = dict(state='complete', validated_at=datetime.now(timezone.utc).isoformat(), inputs=400,
        final_report=str(path), final_report_sha256=n.digest(path), all_order_sequence_hash_and_length_checks_passed=True,
        csv_and_json_levels_agree=True, all_local_rules_recomputed=True,
        all_five_sequence_raw_hashes_verified=True, protected_files_unchanged=True,
        groups=groups, organizer_levels_established=0, binding_or_assay_success_established=False,
        new_gpu_or_remote_cpu_spend_usd=0,
        recorded_spend=dict(runpod=56.5899880444631, own_modal_estimate=10.718947158590089),
        caveat='General scale only; antibody/CDR typing, organizer corpora/domain policy and independently predicted monomer structures are unestablished.')
    (HERE / 'VALIDATION.json').write_text(json.dumps(summary, indent=2) + '\n')
    table = ['| Existing arm | Rows | Local 1 | Local 2 | Local 3 | Local 4 | Unestablished |',
             '| --- | ---: | ---: | ---: | ---: | ---: | ---: |']
    for arm, group in groups.items():
        c = group['local_levels']
        table.append('| ' + ' | '.join([arm, str(group['inputs'])] + [str(c.get(str(k), 0)) for k in range(1, 5)] + [str(c.get('unestablished', 0))]) + ' |')
    measured = '\n'.join(table) + '\n\nActual command and source hashes: SOURCE.json. Final report: ' + str(path.relative_to(HERE)) + '.\n'
    measured += 'Per-design identities, alignment lengths, coverage and domain scores: EVIDENCE.csv.\n'
    measured += 'Input order/sequence/coordinate fingerprints and all five raw sequence-search hashes were verified; JSON/CSV levels agree and the production local rules recompute exactly.\n'
    measured += 'These are general-scale provisional predictions only. Antibody/CDR typing is unimplemented; actual antibody-database match evidence and its consequence are in ANTIBODY_SCOPE.md and ANTIBODY_SCOPE_EXAMPLE.json. No official acceptance/rejection, binding, pH-selectivity or assay result is established.\n'
    if report.get('structural_search',{}).get('positive_evidence_only'):
        measured += 'BOUNDED_POLICY.md applies: this is positive evidence with unfinished alignments, not exhaustive best structural similarity. Unresolved levels are withheld; no local3/4 is certified by missing matches.\n'
    results = HERE / 'RESULTS.md'
    text = results.read_text()
    begin = text.index('## Results\n')
    end = text.index('## Spend\n', begin)
    results.write_text(text[:begin] + '## Results\n\n' + measured + '\n' + text[end:])
    print(json.dumps(dict(inputs=400, groups=groups, organizer_levels_established=0)))


if __name__ == '__main__':
    main()
