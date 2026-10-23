"""Independently check frozen conditional evidence and emit a readable report."""
import collections
import csv
import json
import sqlite3

from score_evidence import identities, exact, sha, HERE, DB, table_hash


def main():
    report = json.loads((HERE/'PREDICTIONS.json').read_text())
    pool = json.loads((HERE/'POOL.json').read_text())
    assert sha(HERE/'POOL.json') == report['pool_manifest_sha256']
    assert sha(report['reference_annotation_checkpoint']) == report['reference_annotation_status_sha256']
    assert sha(HERE/'score_evidence.py') == report['script_sha256']
    annotation = json.loads((HERE/'REFERENCE_ANNOTATION_CHECKPOINT.json').read_text())
    assert sha(HERE/'annotate_references.py') == annotation['script_sha256']
    assert annotation['state'] == 'complete'
    rows = report['designs']
    connection = sqlite3.connect('file:'+str(DB)+'?mode=ro', uri=True)
    assert table_hash(connection, 'SELECT id,sequence_sha256,canonical FROM heavy ORDER BY id') == pool['immutable_heavy_rows_sha256']
    assert table_hash(connection, 'SELECT q,heavy FROM query_heavy ORDER BY q,heavy') == pool['immutable_query_heavy_rows_sha256']
    assert [int(r['query_id'][1:]) for r in rows] == [q for q, in connection.execute('SELECT id FROM query ORDER BY id')]
    compared_best = 0
    for row in rows:
        assert row['organizer_level'] is None and not row['qualifies_for_submission']
        for normalizer, evidence in row['normalizer_evidence'].items():
            best = evidence['best_observed']
            if best is None:
                continue
            sequence, loop, state, kind, flanks, span = connection.execute(
                'SELECT sequence,cdrh3,annotation_state,chain_type,flanks_complete,span_verified FROM heavy WHERE id=?',
                (best['heavy_reference_id'],)).fetchone()
            assert state == 'done' and kind == 'H' and flanks and span
            import hashlib
            assert hashlib.sha256(sequence.encode()).hexdigest() == best['full_reference_sequence_sha256']
            assert loop == best['reference_cdrh3']
            assert exact(identities(row['complete_query_cdrh3'], loop)[normalizer]) == best['identity']
            assert connection.execute('SELECT 1 FROM query_heavy WHERE q=? AND heavy=?',
                                      (int(row['query_id'][1:]), best['heavy_reference_id'])).fetchone()
            compared_best += 1
    connection.close()
    counts = collections.Counter(str(r['conditional_predicted_level']) for r in rows)
    primary_counts = collections.Counter(str(r['normalizer_evidence'][report['primary_normalizer']]['level_over_measured_reference_subset']) for r in rows)
    changed = sum(len({e['level_over_measured_reference_subset'] for e in r['normalizer_evidence'].values()}) > 1 for r in rows)
    with (HERE/'EVIDENCE.csv').open('w', newline='') as stream:
        fields = ['name', 'source_arm', 'returned_heavy_references', 'usable_heavy_references', 'annotation_gaps',
                  'primary_subset_level', 'conditional_level', 'conditional_possible_levels', 'normalizer_changes_subset_level', 'organizer_level_established']
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow(dict(name=row['name'], source_arm=row['source_arm'], returned_heavy_references=row['returned_heavy_reference_candidates'],
                                 usable_heavy_references=row['valid_complete_reference_annotations'], annotation_gaps=sum(row['reference_annotation_gaps'].values()),
                                 primary_subset_level=row['normalizer_evidence'][report['primary_normalizer']]['level_over_measured_reference_subset'],
                                 conditional_level=row['conditional_predicted_level'], conditional_possible_levels='|'.join(map(str, row['conditional_possible_levels'])),
                                 normalizer_changes_subset_level=len({e['level_over_measured_reference_subset'] for e in row['normalizer_evidence'].values()}) > 1,
                                 organizer_level_established=False))
    validation = dict(state='complete', queries=50, annotated_canonical_references=annotation['states']['done'],
                      excluded_noncanonical_references=annotation['states']['excluded_noncanonical'], independently_recomputed_best_identities=compared_best,
                      conditional_level_counts=dict(counts), primary_measured_subset_level_counts=dict(primary_counts),
                      queries_with_normalizer_sensitive_subset_levels=changed, organizer_levels_established=0,
                      source_and_immutable_pool_hashes_verified=True, predictions_sha256=sha(HERE/'PREDICTIONS.json'),
                      evidence_csv_sha256=sha(HERE/'EVIDENCE.csv'), report_builder_sha256=sha(__file__), new_gpu_or_remote_cpu_spend_usd=0)
    (HERE/'VALIDATION.json').write_text(json.dumps(validation, indent=2)+'\n')
    lines = ['# Conditional antibody novelty evidence', '',
             'CPU ANARCII2.0.8 annotation of all returned-hit-associated heavy references completed. This is a supplementary antibody-scale calculation; the organizer structural typer and VH/VHH distinction remain unmeasured.', '',
             '- Inputs:50 complete-flanked H annotations among400 saved designs. The other350 have no CDRH3 result here.',
             '- Restored109,342 unique heavy reference sequences from all raw paired/unpaired source associations;107,440 canonical references annotated,1,902 noncanonical sequences explicitly excluded.105,529 references have usable complete H loops.',
             '- Unit Levenshtein identities were calculated with max-length, query-length and reference-length normalization; normalized indel similarity was tested separately. Boundaries70% and95% are inclusive named constants.22/50 measured-subset levels change with normalization.',
             '- Primary max-length arithmetic over measured reference subsets:11 level2,34 level3,5 level4. These are not predictions under missing annotations.',
             '- Allowing missing-reference completions and all four normalizers:49/50 scalar calls withheld; one conditional local4. That one has only two returned heavy references (one distinct usable loop). Finite full-chain retrieval is a major limitation; this is not a novelty-gate pass.',
             '- All200 exported best-reference identities independently recomputed and source associations verified. Publisher self-loop, unequal-length and boundary controls passed.', '',
             'Reference annotation took1,904.43s; conditional scoring11.30s, on local CPU. Raw per-reference annotations and associations remain in the cache SQLite file named by POOL.json; raw model output hashes and an immutable annotation checkpoint are retained.', '',
             'Not established: organizer antibody class, organizer edit-distance normalization, complete PLAbDab CDR retrieval, organizer reference snapshots/global-identity policy or patent corpus. No binding, pH behavior, assay success, entry ranking or qualification is measured. See PREDICTIONS.json for per-design evidence and gaps, EVIDENCE.csv for compact rows, and VALIDATION.json for checks.', '',
             'Source rule: [Adaptyv novelty](https://www.adaptyvbio.com/blog/novelty). Same public PLAbDab source snapshot as the sequence run; no replacement database was used.', '',
             'Recorded own-lane spend: Runpod$56.589988 + Modal estimate$10.718947 =~$67.31. New GPU/cloud compute$0.00; separate generation/account-wide billing excluded.', '']
    (HERE/'RESULTS.md').write_text('\n'.join(lines))
    print(json.dumps(validation))


if __name__ == '__main__':
    main()
