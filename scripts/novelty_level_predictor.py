#!/usr/bin/env python3
"""CPU MMseqs2 sequence evidence and honest general-novelty level bounds.

This is the sequence half of Proteinbase's general scale. It does not predict
coordinates, segment domains, run Foldseek/TM-align, or detect antibodies.
Missing structural/database evidence never becomes a guessed level4.
"""
import argparse
import csv
from datetime import datetime, timezone
import hashlib
import itertools
import json
import os
from pathlib import Path
import shutil
import subprocess
import time


RULE_SOURCE = 'https://www.adaptyvbio.com/blog/novelty'
REQUIRED_DATABASES = ('swissprot', 'pdb', 'patents', 'therapeutic_antibody', 'plabdab')
AA20 = frozenset('ACDEFGHIKLMNPQRSTVWY')
# Explicit boundary choices. Exact30% stays low; exact70% stays moderate.
SEQUENCE_MODERATE_THRESHOLD = 0.30
SEQUENCE_HIGH_THRESHOLD = 0.70
SEQUENCE_THRESHOLD_EQUALITY_COUNTS_AS_HIGHER = False
# Structural TM cutoffs are inclusive; coverage follows the two published rules.
STRUCTURE_MODERATE_TM = 0.50
STRUCTURE_HIGH_TM = 0.80
STRUCTURE_TM_EQUALITY_COUNTS_AS_MATCH = True
STRUCTURE_MODERATE_COVERAGE = 0.70
STRUCTURE_MODERATE_COVERAGE_INCLUSIVE = False
STRUCTURE_HIGH_COVERAGE = 0.70
STRUCTURE_HIGH_COVERAGE_INCLUSIVE = True
# Search settings are local choices: the organizer did not publish these knobs.
MMSEQS_EVALUE = 0.001
MMSEQS_SENSITIVITY = 7.5
MMSEQS_MAX_SEQS = 10000
MMSEQS_MIN_QUERY_COVERAGE = 0.0
MMSEQS_MEMORY_LIMIT = '4G'
FIELDS = ('query', 'target', 'fident', 'nident', 'alnlen', 'qlen', 'tlen',
          'qstart', 'qend', 'tstart', 'tend', 'evalue', 'bits', 'qcov', 'tcov')


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(1024*1024), b''): h.update(block)
    return h.hexdigest()


def write_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')


def sequence_class(identity):
    if not 0 <= identity <= 1: raise ValueError('Identity must be a fraction in[0,1]')
    if identity > SEQUENCE_HIGH_THRESHOLD: return 'high'
    if identity > SEQUENCE_MODERATE_THRESHOLD: return 'moderate'
    return 'low'


def general_level(seq_class, structure_class):
    """Published general-scale truth table; high structurally includes moderate."""
    if seq_class not in ('low', 'moderate', 'high') or structure_class not in ('low', 'moderate', 'high'):
        raise ValueError('Unknown evidence class')
    if seq_class == 'high' and structure_class in ('moderate', 'high'): return 1
    if seq_class == 'high' or structure_class == 'high' or (seq_class == 'moderate' and structure_class == 'moderate'): return 2
    if seq_class == 'moderate' or structure_class == 'moderate': return 3
    return 4


def structural_class(protein_length, domain_matches):
    """Rule evaluator only, NOT a structural search or domain predictor.

    Inputs are disjoint query-domain intervals (1-based inclusive), with TM-score
    normalized by that query domain's length. Coverage is summed over qualifying
    domains, never over redundant hits and never from an average of TM-scores.
    """
    if not isinstance(protein_length, int) or isinstance(protein_length, bool) or protein_length <= 0:
        raise ValueError('Invalid protein length')
    occupied, high, moderate = set(), set(), set()
    for d in domain_matches:
        start, end, tm = d['start'], d['end'], d['tm_score']
        if not isinstance(start, int) or not isinstance(end, int) or not 1 <= start <= end <= protein_length:
            raise ValueError('Invalid domain coordinates')
        if not 0 <= tm <= 1 or d.get('normalization') != 'query_domain_length':
            raise ValueError('Require query-domain-normalized TM-score in[0,1]')
        span = set(range(start, end+1))
        if span & occupied: raise ValueError('Overlapping domains would double-count coverage')
        occupied |= span
        if tm >= STRUCTURE_MODERATE_TM: moderate |= span
        if tm >= STRUCTURE_HIGH_TM: high |= span
    high_fraction, moderate_fraction = len(high)/protein_length, len(moderate)/protein_length
    category = ('high' if high_fraction >= STRUCTURE_HIGH_COVERAGE else
                'moderate' if moderate_fraction > STRUCTURE_MODERATE_COVERAGE else 'low')
    return dict(category=category, high_coverage=high_fraction, moderate_coverage=moderate_fraction)


def possible_levels(observed_class, missing_databases, structure=None):
    # A missing database may contain a higher-identity hit, never lower an observed hit.
    order = ('low', 'moderate', 'high')
    seq_classes = order[order.index(observed_class):] if missing_databases else (observed_class,)
    structure_classes = (structure,) if structure else order
    return sorted({general_level(s, t) for s, t in itertools.product(seq_classes, structure_classes)})


def load_designs(path):
    with Path(path).open(newline='') as f:
        reader = csv.DictReader(f)
        if not {'name', 'sequence'} <= set(reader.fieldnames or []):
            raise ValueError('Input CSV requires name,sequence headers')
        records = list(reader)
    if not records: raise ValueError('No designs')
    seen = set()
    for r in records:
        if not r['name'] or r['name'] in seen: raise ValueError('Names must be nonempty and unique')
        seen.add(r['name'])
        if not r['sequence'] or set(r['sequence']) - AA20:
            raise ValueError('Sequence must contain only canonical uppercase amino acids; no silent repairs')
        if r.get('molecule_class', 'single_chain') != 'single_chain':
            raise ValueError('Only the general single_chain scale is supported; antibody scale is a documented gap')
    return records


def parse_hits(path, database, queries):
    hits = {q: [] for q in queries}
    with Path(path).open(newline='') as f:
        for fields in csv.reader(f, delimiter='\t'):
            if len(fields) != len(FIELDS): raise ValueError('Unexpected MMseqs output schema')
            r = dict(zip(FIELDS, fields))
            if r['query'] not in hits: raise ValueError('Unknown MMseqs query')
            for k in ('fident', 'evalue', 'bits', 'qcov', 'tcov'): r[k] = float(r[k])
            for k in ('nident', 'alnlen', 'qlen', 'tlen', 'qstart', 'qend', 'tstart', 'tend'): r[k] = int(r[k])
            sequence_class(r['fident'])
            if r['alnlen'] <= 0 or not 0 <= r['nident'] <= r['alnlen']:
                raise ValueError('Invalid identical-residue/alignment counts')
            r['reported_fident'] = r['fident']
            r['fident'] = r['nident'] / r['alnlen']
            if abs(r['reported_fident']-r['fident']) > 0.002:
                raise ValueError('Reported identity disagrees with integer alignment counts')
            if r['qlen'] != len(queries[r['query']]['sequence']): raise ValueError('Query length mismatch')
            if r['evalue'] > MMSEQS_EVALUE * 1.001: raise ValueError('Unqualified alignment')
            r['database'] = database
            hits[r['query']].append(r)
    return hits


def identity_ranked_hits(hits):
    return sorted(hits, key=lambda h: (-h['fident'], -h['bits'], h['database'], h['target']))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('designs', type=Path)
    ap.add_argument('--databases', type=Path, required=True, help='Explicit provenance/scope manifest; no substitutions')
    ap.add_argument('--out', type=Path, required=True, help='Fresh output directory')
    ap.add_argument('--mmseqs', default='mmseqs')
    ap.add_argument('--threads', type=int, default=4)
    ap.add_argument('--timeout-per-database', type=int, default=600)
    args = ap.parse_args(argv)
    if not 1 <= args.threads <= 4: ap.error('CPU-only task budget:1–4 threads')
    if args.timeout_per_database <= 0: ap.error('Timeout must be positive')
    records = load_designs(args.designs)
    manifest = json.loads(args.databases.read_text())
    dbs = manifest['databases']
    if {d['id'] for d in dbs} != set(REQUIRED_DATABASES) or len(dbs) != len(REQUIRED_DATABASES):
        raise ValueError('Manifest must explicitly account for every named database exactly once')
    if args.out.exists(): raise FileExistsError('Fresh output directory required')
    binary = shutil.which(args.mmseqs)
    if not binary: raise FileNotFoundError('MMseqs2 not installed; no fallback aligner is substituted')
    args.out.mkdir(parents=True)
    queries = {f'q{i:06d}': r for i, r in enumerate(records, 1)}
    query_fasta = args.out / 'queries.fasta'
    query_fasta.write_text(''.join(f'>{q}\n{r["sequence"]}\n' for q, r in queries.items()))
    version = subprocess.check_output([binary, 'version'], text=True).strip()
    all_hits = {q: [] for q in queries}
    search_records, unavailable = [], []
    for d in dbs:
        if d['status'] != 'available':
            unavailable.append(d['id']); search_records.append(dict(**d, search_status='not_run')); continue
        rec = dict(d)
        fasta = Path(d['fasta_path']).expanduser()
        if not fasta.is_absolute(): fasta = args.databases.resolve().parent / fasta
        try:
            if digest(fasta) != d['fasta_sha256']: raise ValueError('Reference FASTA fingerprint differs')
            tsv = args.out / (d['id'] + '.hits.tsv')
            cmd = [binary, 'easy-search', str(query_fasta), str(fasta), str(tsv),
                   str(args.out / ('tmp_' + d['id'])), '--threads', str(args.threads),
                   '--gpu', '0', '--alignment-mode', '3', '-a', '-s', str(MMSEQS_SENSITIVITY),
                   '-e', str(MMSEQS_EVALUE), '--max-seqs', str(MMSEQS_MAX_SEQS),
                   '--min-seq-id', '0', '-c', str(MMSEQS_MIN_QUERY_COVERAGE),
                   '--split-memory-limit', MMSEQS_MEMORY_LIMIT, '--remove-tmp-files', '1',
                   '--format-output', ','.join(FIELDS)]
            rec['command'] = cmd; start = time.monotonic()
            with (args.out / (d['id'] + '.search.log')).open('w') as log:
                subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT, check=True,
                    timeout=args.timeout_per_database, env={**os.environ, 'CUDA_VISIBLE_DEVICES': ''})
            matched = parse_hits(tsv, d['id'], queries)
            for q in queries: all_hits[q].extend(matched[q])
            rec.update(search_status='completed', elapsed_seconds=time.monotonic()-start,
                       hit_rows=sum(map(len, matched.values())), raw_hits_sha256=digest(tsv))
            if not d.get('scope_complete', False): unavailable.append(d['id'])
        except (OSError, ValueError, subprocess.SubprocessError) as exc:
            rec.update(search_status='failed', error=repr(exc)); unavailable.append(d['id'])
        search_records.append(rec)
        write_json(args.out / 'search_status.json', search_records)
    results = []
    for q, design in queries.items():
        hs = identity_ranked_hits(all_hits[q])
        best = hs[0] if hs else None
        observed = sequence_class(best['fident']) if best else 'low'
        saturated = any(sum(h['database'] == d for h in hs) >= MMSEQS_MAX_SEQS for d in REQUIRED_DATABASES)
        gaps = sorted(set(unavailable))
        if saturated: gaps.append('MMseqs_max_seqs_cap')
        levels = possible_levels(observed, gaps)
        results.append(dict(name=design['name'], length=len(design['sequence']),
            sequence_sha256=hashlib.sha256(design['sequence'].encode()).hexdigest(),
            best_detected_identity=best['fident'] if best else None,
            sequence_evidence_class=observed if best else 'no_significant_hit',
            best_hit=best, hit_count=len(hs), top_identity_hits=hs[:10],
            sequence_database_gaps=gaps, structure_status='not_measured',
            predicted_level=levels[0] if len(levels) == 1 else None,
            possible_levels=levels,
            levels_if_complete_sequence_search_and_structure_class={s: general_level(observed, s) for s in ('low', 'moderate', 'high')},
            local_gate_assessment='sequence_high_at_risk_of_level_below3' if observed == 'high' else 'unestablished',
            explanation='Level withheld: structural prediction, consensus domain segmentation, Foldseek and TM-align not run. Missing databases can raise detected sequence similarity. Local MMseqs parameters may differ from ProteinTyper.'))
    report = dict(status='sequence_half_complete' if not unavailable else 'sequence_half_partial_databases',
        rule_source=RULE_SOURCE, generated_at=datetime.now(timezone.utc).isoformat(),
        organizer_score=False, general_scale_only=True, antibody_detection='not_implemented',
        input_csv_sha256=digest(args.designs), database_manifest_sha256=digest(args.databases),
        predictor_sha256=digest(__file__), mmseqs_version=version,
        constants={k: v for k, v in globals().items() if k.startswith(('SEQUENCE_', 'STRUCTURE_', 'MMSEQS_')) and isinstance(v, (str,int,float,bool))},
        database_searches=search_records, designs=results, new_gpu_spend_usd=0,
        limitations=['No structure prediction from sequences on CPU within this task budget',
                     'No consensus domain segmentation or Foldseek -> TM-align search',
                     'No organizer database snapshots or alignment parameters published',
                     'No antibody typing or antibody-specific CDR novelty scoring',
                     'No-hit search is not measured zero identity or proof of novelty'])
    write_json(args.out / 'predictions.json', report)
    with (args.out / 'predictions.csv').open('w', newline='') as f:
        columns = ['name','length','best_detected_identity','sequence_evidence_class',
                   'predicted_level','possible_levels','local_gate_assessment','structure_status',
                   'sequence_database_gaps','best_hit_database','best_hit_target']
        w = csv.DictWriter(f, fieldnames=columns); w.writeheader()
        for r in results:
            row = {k: r[k] for k in columns if k in r}
            row['possible_levels'] = '|'.join(map(str, r['possible_levels']))
            row['sequence_database_gaps'] = '|'.join(r['sequence_database_gaps'])
            row['best_hit_database'] = (r['best_hit'] or {}).get('database', '')
            row['best_hit_target'] = (r['best_hit'] or {}).get('target', '')
            w.writerow(row)
    print(json.dumps(dict(status=report['status'], designs=len(results),
        databases_searched=sum(r['search_status']=='completed' for r in search_records),
        exact_levels_established=sum(r['predicted_level'] is not None for r in results),
        out=str(args.out), new_gpu_spend_usd=0)))


if __name__ == '__main__': main()
