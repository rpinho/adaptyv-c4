#!/usr/bin/env python3
"""CPU MMseqs2 + Foldseek -> TM-align provisional general-novelty predictor.

Requires existing sequence-verified coordinates for structural searches. Reports
local levels with explicit confidence/scope, never an established organizer score.
No structure inference, consensus domain prediction, antibody typing or GPU use.
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
FOLDSEEK_FIELDS = ('query', 'target', 'qlen', 'tlen', 'qstart', 'qend',
                   'tstart', 'tend', 'alnlen', 'qtmscore', 'ttmscore', 'qaln', 'taln')
FOLDSEEK_MAX_SEQS = 50000
FOLDSEEK_EVALUE = 10.0  # permissive candidate discovery; final rule uses TM only
FOLDSEEK_SENSITIVITY = 9.5
FOLDSEEK_EXACT_TM_SCORE = 1  # convertalis otherwise recomputes approximate TM, even after tmalign


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


def domain_intervals(domain):
    """One domain may contain multiple ordered, disjoint sequence fragments."""
    if 'segments' in domain:
        if 'start' in domain or 'end' in domain:
            raise ValueError('Use segments OR start/end, not both')
        intervals = domain['segments']
        if not isinstance(intervals, (list, tuple)) or not intervals:
            raise ValueError('Domain segments must be nonempty')
    else:
        if 'start' not in domain or 'end' not in domain:
            raise ValueError('Domain requires start/end or segments')
        intervals = [(domain['start'], domain['end'])]
    previous_end = 0
    checked = []
    for interval in intervals:
        if not isinstance(interval, (list, tuple)) or len(interval) != 2:
            raise ValueError('Invalid domain coordinates')
        start, end = interval
        if (any(not isinstance(x, int) or isinstance(x, bool) for x in (start, end))
                or not 1 <= start <= end or start <= previous_end):
            raise ValueError('Invalid, overlapping or unordered domain coordinates')
        checked.append((start, end))
        previous_end = end
    return checked


def domain_coordinates(domain):
    intervals = domain_intervals(domain)
    if 'segments' in domain:
        return dict(segments=[list(pair) for pair in intervals])
    return dict(start=intervals[0][0], end=intervals[0][1])


def domain_positions(domain, protein_length):
    positions = [i for start, end in domain_intervals(domain) for i in range(start, end+1)]
    if positions[-1] > protein_length: raise ValueError('Invalid domain coordinates')
    return positions


def domain_query_pdb(residues, domain):
    """Copy existing CAs, retaining full-query residue IDs across domain gaps."""
    positions = domain_positions(domain, len(residues))
    lines = []
    for serial, position in enumerate(positions, 1):
        _, _, (x,y,z), comp = residues[position-1]
        lines.append(f'ATOM  {serial:5d}  CA  {comp:3s} A{position:4d}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           C\n')
    return ''.join(lines)+'TER\nEND\n'


def structural_class(protein_length, domain_matches):
    """Rule evaluator only, NOT a structural search or domain predictor.

    Inputs are disjoint query domains (1-based inclusive intervals or segments), with TM-score
    normalized by that query domain's length. Coverage is summed over qualifying
    domains, never over redundant hits and never from an average of TM-scores.
    """
    if not isinstance(protein_length, int) or isinstance(protein_length, bool) or protein_length <= 0:
        raise ValueError('Invalid protein length')
    occupied, high, moderate = set(), set(), set()
    for d in domain_matches:
        span = set(domain_positions(d, protein_length))
        tm = d['tm_score']
        if not 0 <= tm <= 1 or d.get('normalization') != 'query_domain_length':
            raise ValueError('Require query-domain-normalized TM-score in[0,1]')
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


def ca_residues(path, chain):
    """Read first-model CAs, including minimal BindCraft CIFs lacking occupancy.

    No sequence repair or target-chain inference. Author numbering is preserved
    for extraction; domain coordinates subsequently index the verified sequence.
    """
    from Bio.Data.PDBData import protein_letters_3to1 as aa
    import math
    if path.suffix.lower() == '.cif':
        from Bio.PDB.MMCIF2Dict import MMCIF2Dict
        data = MMCIF2Dict(str(path))
        def column(auth, label):
            return data.get('_atom_site.'+auth, data.get('_atom_site.'+label))
        names = column('auth_atom_id', 'label_atom_id')
        chains = column('auth_asym_id', 'label_asym_id')
        comps = column('auth_comp_id', 'label_comp_id')
        nums = column('auth_seq_id', 'label_seq_id')
        n = len(names)
        models = data.get('_atom_site.pdbx_PDB_model_num', ['1']*n)
        alt = data.get('_atom_site.label_alt_id', ['.']*n)
        ins = data.get('_atom_site.pdbx_PDB_ins_code', ['?']*n)
        xyz = zip(*(data['_atom_site.Cartn_'+d] for d in 'xyz'))
        records = [(int(number), comp, tuple(map(float, coords)), insertion)
                   for name, cid, comp, number, coords, model, al, insertion in
                   zip(names, chains, comps, nums, xyz, models, alt, ins)
                   if name == 'CA' and cid == chain and model == models[0] and al in ('.','?','A')]
    else:
        from Bio.PDB import PDBParser
        model = PDBParser(QUIET=True).get_structure('query', str(path))[0]
        records = [(r.id[1], r.resname, tuple(map(float, r['CA'].coord)), r.id[2])
                   for r in model[chain] if r.id[0] == ' ' and 'CA' in r]
    if not records or len({r[0] for r in records}) != len(records):
        raise ValueError('Empty/duplicate CA residue numbering')
    if any(r[1] not in aa or r[3] not in (' ', '.', '?') or
           not all(math.isfinite(x) and -999.999 <= x <= 9999.999 for x in r[2]) for r in records):
        raise ValueError('Noncanonical residue, insertion code or invalid CA coordinates')
    return [(n, aa[comp], xyz, comp) for n, comp, xyz, _ in records]


def parse_structure_hits(path, domains):
    hits = {q: [] for q in domains}
    with Path(path).open(newline='') as f:
        for fields in csv.reader(f, delimiter='\t'):
            if len(fields) != len(FOLDSEEK_FIELDS):
                raise ValueError('Unexpected Foldseek TM-align output schema')
            r = dict(zip(FOLDSEEK_FIELDS, fields))
            if r['query'] not in domains: raise ValueError('Unknown structural query')
            for k in ('qlen','tlen','qstart','qend','tstart','tend','alnlen'): r[k] = int(r[k])
            for k in ('qtmscore','ttmscore'): r[k] = float(r[k])
            d = domains[r['query']]
            expected_length = sum(end-start+1 for start,end in domain_intervals(d))
            if r['qlen'] != expected_length or not 0 <= r['qtmscore'] <= 1:
                raise ValueError('Wrong query domain length or TM-score normalization')
            if not 1 <= r['qstart'] <= r['qend'] <= r['qlen']:
                raise ValueError('Invalid structural alignment coordinates')
            hits[r['query']].append(r)
    # Two discovery modes may return the same target. A repeated hit is not
    # independent evidence or additional domain coverage; retain its best qTM.
    for q, hs in hits.items():
        unique = {}
        for h in hs:
            if h['target'] not in unique or h['qtmscore'] > unique[h['target']]['qtmscore']:
                unique[h['target']] = h
        hits[q] = sorted(unique.values(), key=lambda h: (-h['qtmscore'], h['target']))
    return hits


def run_structural_half(records, report, args):
    """Two real CPU stages: permissive Foldseek discovery, full TM-align rerun.

    Database files and every coordinate file must match provenance fingerprints.
    Explicit disjoint query domains are counted once each, using their best qTM.
    Whole-chain domains are allowed only as a stated local assumption.
    """
    db = json.loads(args.structure_database.read_text())
    structures = json.loads(args.structures.read_text())
    binary = shutil.which(args.foldseek)
    if not binary: raise FileNotFoundError('CPU Foldseek not installed')
    base = Path(db['prefix']).expanduser()
    if not base.is_absolute(): base = args.structure_database.resolve().parent/base
    if not db.get('files'): raise ValueError('Structure DB requires file fingerprints')
    for suffix, fingerprint in db['files'].items():
        if digest(str(base)+suffix) != fingerprint: raise ValueError('Structure DB fingerprint differs: '+suffix)
    out = args.out/'structure'; out.mkdir()
    folder = out/'queries'; folder.mkdir()
    specs = {r['name']: r for r in structures['structures']}
    if len(specs) != len(structures['structures']): raise ValueError('Duplicate structure name')
    if set(specs)-{r['name'] for r in records}: raise ValueError('Unknown structure name')
    domains, evidence = {}, {}
    for i, design in enumerate(records, 1):
        spec = specs.get(design['name'])
        if not spec or spec.get('status') != 'available':
            evidence[design['name']] = dict(status='coordinates_unavailable', reason=(spec or {}).get('reason','No coordinate manifest entry'))
            continue
        path = Path(spec['path']).expanduser()
        if not path.is_absolute(): path = args.structures.resolve().parent/path
        if digest(path) != spec['sha256']: raise ValueError('Query coordinate fingerprint differs')
        rr = ca_residues(path, spec['chain'])
        if 'author_residue_range' in spec:
            lo, hi = spec['author_residue_range']; rr = [r for r in rr if lo <= r[0] <= hi]
        if ''.join(r[1] for r in rr) != design['sequence']:
            raise ValueError('Coordinate/sequence mismatch: '+design['name'])
        intervals = spec['domains']
        if not intervals: raise ValueError('At least one declared domain is required')
        # Reuse interval validation, rejecting overlaps before any external search.
        structural_class(len(rr), [dict(**domain_coordinates(d), tm_score=0., normalization='query_domain_length') for d in intervals])
        for j, d in enumerate(intervals, 1):
            if not d.get('source'): raise ValueError('Domain provenance is required')
            q = f'q{i:06d}d{j:03d}'
            domains[q] = dict(**d, name=design['name'],
                query_residue_positions=domain_positions(d, len(rr)))
            (folder/(q+'.pdb')).write_text(domain_query_pdb(rr, d))
        evidence[design['name']] = dict(status='coordinates_verified', input=spec, domains=[])
    write_json(out/'query_domains.json', domains)
    if not domains: raise ValueError('No verified coordinate domains to search')
    querydb, found, realigned = out/'querydb', out/'discovered', out/'tmaligned'
    # No GPU inference or GPU search; explicit 4G search memory and <=4 threads.
    modes = [0,1] if args.structure_prefilter_mode == 'union' else [int(args.structure_prefilter_mode)]
    commands = [[binary,'createdb',str(folder),str(querydb),'--threads',str(args.threads)]]
    for mode in modes:
        discovery = out/f'discovered_mode{mode}'
        commands.append([binary,'search',str(querydb),str(base),str(discovery),str(out/f'tmp_mode{mode}'),
         '--threads',str(args.threads),'--gpu','0','--split-memory-limit','4G',
         '-s',str(FOLDSEEK_SENSITIVITY),'-e',str(FOLDSEEK_EVALUE),
         '--max-seqs',str(args.structure_max_seqs),
         '--prefilter-mode',str(mode),'-a','1'])
    commands.append([binary,'mergedbs',str(querydb),str(found),
                     *[str(out/f'discovered_mode{mode}') for mode in modes]])
    commands.extend([
        [binary,'tmalign',str(querydb),str(base),str(found),str(realigned),
         '--threads',str(args.threads),'--tmalign-fast','0','--tmscore-threshold','0',
         '--tmscore-threshold-mode','1','--tmalign-hit-order','1','-a','1'],
        [binary,'convertalis',str(querydb),str(base),str(realigned),str(out/'hits.tsv'),
         '--threads',str(args.threads),'--exact-tmscore',str(FOLDSEEK_EXACT_TM_SCORE),
         '--format-output',','.join(FOLDSEEK_FIELDS)]])
    start = time.monotonic()
    for i, cmd in enumerate(commands):
        with (out/f'stage_{i}.log').open('w') as log:
            subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT, check=True,
                timeout=args.timeout_structure, env={**os.environ,'CUDA_VISIBLE_DEVICES':'','OMP_NUM_THREADS':str(args.threads)})
    hits = parse_structure_hits(out/'hits.tsv', domains)
    for q, d in domains.items():
        hs = hits[q]
        evidence[d['name']]['domains'].append(dict(**d, query_id=q,
            best_hit=hs[0] if hs else None, hit_count=len(hs), top_hits=hs[:10]))
    for design, r in zip(records, report['designs']):
        if r['name'] != design['name']: raise ValueError('Sequence-report order mismatch')
        e = evidence[r['name']]; r['structure_evidence'] = e
        if e['status'] == 'coordinates_unavailable':
            r.update(structure_status=e['status'], local_predicted_level=None,
                     prediction_confidence='unestablished: coordinates unavailable', organizer_level_established=False)
            continue
        matches = [dict(**domain_coordinates(d),tm_score=(d['best_hit'] or {}).get('qtmscore',0.),
                        normalization='query_domain_length') for d in e['domains']]
        measured = structural_class(len(design['sequence']), matches)
        measured['evidence_status'] = ('matches_measured' if any(d['best_hit'] for d in e['domains'])
                                       else 'no_candidate_match_detected_not_zero_TM')
        measured['domain_segmentation_coverage'] = sum(len(domain_positions(d, len(design['sequence']))) for d in e['domains'])/len(design['sequence'])
        measured['tm_boundary_precision_caveat'] = any(
            (d['best_hit'] or {}).get('qtmscore') in (.5,.8) for d in e['domains'])
        observed = r['sequence_evidence_class']
        if observed == 'no_significant_hit': observed = 'low'
        local = general_level(observed, measured['category'])
        # Missing sequence collections and domain/structure-policy differences
        # prevent an exact organizer claim even when the local class is singular.
        r.update(structure_status='foldseek_then_tmalign_completed',
            structural_similarity=measured, local_predicted_level=local,
            predicted_level=local, organizer_level_established=False,
            level_status='provisional_local_prediction',
            possible_levels=possible_levels(observed,r['sequence_database_gaps']),
            possible_levels_if_structure_policy_accepted=possible_levels(observed,r['sequence_database_gaps'],measured['category']),
            prediction_confidence='low: public historical PDB snapshot; provided coordinates; declared domains; organizer consensus/versions and patent corpus unverified; no-hit sequence class is provisional',
            local_gate_assessment='provisional_local_level_'+str(local),
            explanation='Level applies only to searched databases, supplied coordinates and declared domain segmentation. Not an organizer score; see confidence and possible_levels.')
    report.update(status='cpu_sequence_and_structure_screen_complete_with_scope_caveats',
        structure_database=db, structure_manifest_sha256=digest(args.structures),
        domain_segmentation_provenance=structures.get('segmentation_provenance'),
        foldseek_version=subprocess.check_output([binary,'version'],text=True).strip(),
        structural_search=dict(commands=commands, elapsed_seconds=time.monotonic()-start,
            query_domains=len(domains), hit_rows=sum(map(len,hits.values())),
            raw_alignment_rows=sum(1 for _ in (out/'hits.tsv').open()),
            repeated_targets='Deduplicated per query, retaining maximum query-normalized TM',
            exact_tmscore_export=FOLDSEEK_EXACT_TM_SCORE,
            raw_hits_sha256=digest(out/'hits.tsv'),
            alignment='Foldseek 3Di+AA search -> Foldseek built-in TM-align, fast mode OFF; scores normalized by complete query domain length',
            hit_cap_per_discovery_mode=args.structure_max_seqs, prefilter_modes=modes,
            evalue=FOLDSEEK_EVALUE),
        limitations=['Public historical PDB snapshot, not the unpublished organizer structural corpus/version',
            'Existing verified coordinates used; organizer ESMFold/ESMFold2-Fast structures not reproduced',
            'Domain provenance is supplied; organizer domain-consensus implementation/version not verified',
            'Foldseek finite candidate search may miss matches; no-hit is not proof of new fold',
            'Organizer patent collection/snapshot and unpublished search settings remain unverified',
            'Antibody detection and separate CDR novelty scale not implemented',
            'All numeric levels are local provisional predictions; organizer levels remain unestablished'])


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('designs', type=Path)
    ap.add_argument('--databases', type=Path, help='Explicit sequence provenance/scope manifest')
    ap.add_argument('--sequence-report', type=Path, help='Reuse fingerprint-verified existing sequence evidence')
    ap.add_argument('--structures', type=Path, help='Sequence-matched coordinate + disjoint domain manifest')
    ap.add_argument('--structure-database', type=Path, help='Foldseek DB prefix, date, size and fingerprints')
    ap.add_argument('--foldseek', default='foldseek')
    ap.add_argument('--structure-prefilter-mode', choices=('0','1','union'), default='union',
                    help='Default union: both kmer0 and all-reference ungapped1 discovery; CPU only')
    ap.add_argument('--structure-max-seqs', type=int, default=FOLDSEEK_MAX_SEQS,
                    help='Structural candidate cap per query, maximum50000')
    ap.add_argument('--timeout-structure', type=int, default=3600, help='Seconds per CPU structural stage')
    ap.add_argument('--out', type=Path, required=True, help='Fresh output directory')
    ap.add_argument('--mmseqs', default='mmseqs')
    ap.add_argument('--threads', type=int, default=4)
    ap.add_argument('--timeout-per-database', type=int, default=600)
    args = ap.parse_args(argv)
    if not 1 <= args.threads <= 4: ap.error('CPU-only task budget:1–4 threads')
    if args.timeout_per_database <= 0: ap.error('Timeout must be positive')
    if not 1 <= args.structure_max_seqs <= 50000: ap.error('Structural candidate cap must be1–50000')
    if args.timeout_structure <= 0: ap.error('Structural stage timeout must be positive')
    if bool(args.structures) != bool(args.structure_database): ap.error('Provide both structural manifests')
    records = load_designs(args.designs)
    if args.sequence_report:
        report = json.loads(args.sequence_report.read_text())
        if report['input_csv_sha256'] != digest(args.designs): raise ValueError('Sequence report input differs')
        if [(r['name'],r['sequence_sha256']) for r in report['designs']] != [(r['name'],hashlib.sha256(r['sequence'].encode()).hexdigest()) for r in records]:
            raise ValueError('Sequence evidence does not match input')
        if args.out.exists(): raise FileExistsError('Fresh output directory required')
        if not args.structures: ap.error('Sequence report reuse requires structural manifests')
        args.out.mkdir(parents=True)
        report.update(sequence_report_sha256=digest(args.sequence_report), predictor_sha256=digest(__file__),
                      generated_at=datetime.now(timezone.utc).isoformat())
        run_structural_half(records,report,args)
        save_report(report,args.out)
        return
    if not args.databases: ap.error('Provide --databases or --sequence-report')
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
    if args.structures: run_structural_half(records,report,args)
    save_report(report,args.out)


def save_report(report, out):
    results = report['designs']
    write_json(out / 'predictions.json', report)
    with (out / 'predictions.csv').open('w', newline='') as f:
        columns = ['name','length','best_detected_identity','sequence_evidence_class',
                   'predicted_level','possible_levels','local_gate_assessment','structure_status',
                   'sequence_database_gaps','best_hit_database','best_hit_target',
                   'local_predicted_level','prediction_confidence','organizer_level_established']
        w = csv.DictWriter(f, fieldnames=columns); w.writeheader()
        for r in results:
            row = {k: r[k] for k in columns if k in r}
            row['possible_levels'] = '|'.join(map(str, r['possible_levels']))
            row['sequence_database_gaps'] = '|'.join(r['sequence_database_gaps'])
            row['best_hit_database'] = (r['best_hit'] or {}).get('database', '')
            row['best_hit_target'] = (r['best_hit'] or {}).get('target', '')
            w.writerow(row)
    print(json.dumps(dict(status=report['status'], designs=len(results),
        databases_searched=sum(r['search_status']=='completed' for r in report['database_searches']),
        exact_levels_established=sum(r.get('organizer_level_established',False) for r in results),
        local_levels_predicted=sum(r.get('local_predicted_level') is not None for r in results),
        out=str(out), new_gpu_spend_usd=0)))


if __name__ == '__main__': main()
