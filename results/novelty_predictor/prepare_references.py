"""Normalize the downloaded named databases; preserve scope and fingerprints."""
import csv
import gzip
import hashlib
import json
from pathlib import Path
import re

OUT = Path(__file__).resolve().parent
CACHE = Path.home() / '.cache/adaptyv-c4/novelty'


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(1024*1024), b''): h.update(block)
    return h.hexdigest()


def fasta_rows(path):
    with gzip.open(path, 'rt') as f:
        name, seq = None, []
        for line in f:
            if line.startswith('>'):
                if name is not None: yield name, ''.join(seq)
                name, seq = line[1:].strip(), []
            else: seq.append(line.strip())
        if name is not None: yield name, ''.join(seq)


def csv_rows(path, family):
    op = gzip.open if str(path).endswith('.gz') else open
    with op(path, 'rt', encoding='utf-8-sig', newline='') as f:
        reader = csv.DictReader(f)
        if family == 'therapeutic_antibody':
            fields = ['HeavySequence', 'LightSequence', 'HeavySequence(ifbispec)', 'LightSequence(ifbispec)']
        elif 'GBSeq_sequence' in reader.fieldnames:
            fields = ['GBSeq_sequence']
        else:
            fields = [c for c in reader.fieldnames if c.lower() in ('sequence', 'heavy_sequence', 'light_sequence')]
        if not fields or not set(fields) <= set(reader.fieldnames): raise ValueError('Unexpected sequence column schema')
        for i, row in enumerate(reader, 1):
            accession = row.get('Therapeutic', row.get('ID', row.get('GBSeq_accession-version', str(i))))
            accession = re.sub(r'[^A-Za-z0-9_.-]', '_', accession)
            for field in fields:
                yield f'{family}:{i}:{accession}:{field}', row[field]


def normalize(database, download_records):
    dest = CACHE / (database + '.fasta')
    stats = dict(emitted=0, duplicates_removed=0, placeholder_fields=0, nonprotein_records=0, invalid_sequences=0,
                 uppercased=0)
    seen = set()
    with dest.open('w') as f:
        for download in download_records:
            path = Path(download['path'])
            assert sha(path) == download['sha256']
            entries = fasta_rows(path) if database in ('swissprot', 'pdb') else csv_rows(path, database)
            for identifier, sequence in entries:
                if database == 'pdb' and 'mol:protein' not in identifier:
                    stats['nonprotein_records'] += 1; continue
                if sequence is None or sequence.strip().lower() in ('', 'na', 'n/a', 'nan', 'none'):
                    stats['placeholder_fields'] += 1; continue
                # Reference ambiguity codes are retained. Queries still require canonical AA20.
                if sequence != sequence.upper(): stats['uppercased'] += 1
                sequence = sequence.upper()
                if not re.fullmatch('[A-Z]+', sequence):
                    stats['invalid_sequences'] += 1; continue
                if sequence in seen:
                    stats['duplicates_removed'] += 1; continue
                seen.add(sequence)
                f.write('>' + identifier.split()[0] + '\n' + sequence + '\n')
                stats['emitted'] += 1
    return dict(id=database, status='available', fasta_path=str(dest), fasta_sha256=sha(dest),
        scope_complete=stats['invalid_sequences'] == 0,
        scope='Public downloaded snapshot only; organizer snapshot/version unknown',
        normalization=stats, sources=download_records)


if __name__ == '__main__':
    sources = {r['database']: r for r in json.loads((OUT / 'DOWNLOADS.json').read_text())['sources']}
    dbs = []
    for database, keys in [('swissprot', ['swissprot']), ('pdb', ['pdb']),
                           ('therapeutic_antibody', ['therapeutic_antibody']),
                           ('plabdab', ['plabdab_paired', 'plabdab_unpaired'])]:
        available = [sources[k] for k in keys if sources[k]['status'] == 'downloaded']
        if available:
            d = normalize(database, available)
            if len(available) != len(keys): d['scope_complete'] = False
            dbs.append(d)
        else: dbs.append(dict(id=database, status='unavailable', reason='Named public source download failed'))
    dbs.append(dict(id='patents', status='unavailable',
        reason='Organizer names patent sequences without identifying the database/provider/release. No exact corpus obtained; PLAbDab is not substituted.'))
    result = dict(schema_version=1, databases=dbs,
                  omitted_scope='No arbitrary local binder reference set is substituted for any named database')
    (OUT / 'DATABASES.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps([dict(database=d['id'], status=d['status'], counts=d.get('normalization')) for d in dbs]))
