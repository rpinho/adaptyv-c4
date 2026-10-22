"""Obtain, verify and inventory explicitly named NCBI pataa public evidence."""
import csv
from collections import Counter
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import tarfile
import time
import urllib.request

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
CACHE = Path.home()/'.cache/adaptyv-c4/novelty/patents_ncbi_pataa_20261021'
DB_URL = 'https://ftp.ncbi.nlm.nih.gov/blast/db/pataa.tar.gz'
TOOL_URL = ('https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.17.0/'
            'ncbi-blast-2.17.0+-aarch64-macosx.tar.gz')
AA20 = frozenset('ACDEFGHIKLMNPQRSTVWY')


def sha(path, algorithm='sha256'):
    h = hashlib.new(algorithm)
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(4*1024*1024), b''):
            h.update(block)
    return h.hexdigest()


def save(name, value):
    (HERE/name).write_text(json.dumps(value, indent=2, allow_nan=False)+'\n')


def fetch_small(url):
    with urllib.request.urlopen(url, timeout=40) as response:
        return response.read(1000000), dict(response.headers)


def download(url, path, expected_md5, expected_bytes=None):
    if path.exists():
        if sha(path, 'md5') != expected_md5:
            raise ValueError('Existing archive does not match publisher: '+str(path))
        return dict(url=url, path=str(path), reused=True,
                    bytes=path.stat().st_size, md5=expected_md5, sha256=sha(path))
    partial = path.with_suffix(path.suffix+'.partial')
    if partial.exists():
        raise FileExistsError('Inspect interrupted partial before retry: '+str(partial))
    start = time.monotonic()
    with urllib.request.urlopen(url, timeout=40) as response, partial.open('xb') as output:
        headers = dict(response.headers)
        length = int(response.headers.get('Content-Length', 0))
        if expected_bytes is not None and length != expected_bytes:
            raise ValueError('Archive size changed from pinned metadata')
        received = 0
        last_checkpoint = 0
        while True:
            block = response.read(4*1024*1024)
            if not block:
                break
            output.write(block)
            received += len(block)
            if received-last_checkpoint >= 64*1024*1024:
                save('DOWNLOAD_PROGRESS.json', dict(url=url, received_bytes=received,
                    expected_bytes=length, elapsed_seconds=time.monotonic()-start))
                last_checkpoint = received
    if length and partial.stat().st_size != length:
        raise ValueError('Truncated archive')
    if sha(partial, 'md5') != expected_md5:
        raise ValueError('Publisher MD5 verification failed')
    partial.rename(path)
    return dict(url=url, path=str(path), headers=headers, bytes=path.stat().st_size,
        md5=expected_md5, sha256=sha(path), elapsed_seconds=time.monotonic()-start)


def main():
    if (HERE/'ACQUISITION.json').exists():
        raise FileExistsError('Do not replace frozen acquisition')
    CACHE.mkdir(parents=True, exist_ok=True)
    if shutil.disk_usage(CACHE).free < 15*1024**3:
        raise OSError('Need at least 15GiB free for local bounded acquisition/search')
    metadata_raw, headers = fetch_small('https://ftp.ncbi.nlm.nih.gov/blast/db/pataa-prot-metadata.json')
    metadata = json.loads(metadata_raw)
    if (metadata['last-updated'], metadata['number-of-sequences'], metadata['number-of-letters']) != (
            '2026-10-21T00:00:00', 3982722, 787058342):
        raise ValueError('Publisher snapshot differs from registered source')
    (HERE/'pataa-prot-metadata.json').write_bytes(metadata_raw)
    receipts = {}
    for kind, url, target in [('exporter', TOOL_URL, CACHE/'blast_tools.tar.gz'),
                              ('patents', DB_URL, CACHE/'pataa.tar.gz')]:
        checksum, checksum_headers = fetch_small(url+'.md5')
        expected = checksum.decode().split()[0]
        (HERE/(kind+'.publisher.md5')).write_bytes(checksum)
        receipts[kind] = download(url, target, expected,
            metadata['bytes-total-compressed'] if kind == 'patents' else None)
        receipts[kind]['checksum_http_headers'] = checksum_headers
        save('DOWNLOADS.json', receipts)
    binary = CACHE/'blastdbcmd'
    with tarfile.open(CACHE/'blast_tools.tar.gz', 'r:gz') as archive:
        selected = [m for m in archive if m.name.endswith('/bin/blastdbcmd')]
        if len(selected) != 1 or not selected[0].isfile():
            raise ValueError('Expected one regular exporter binary')
        with archive.extractfile(selected[0]) as src:
            contents = src.read()
        if binary.exists():
            if binary.read_bytes() != contents:
                raise ValueError('Existing exporter differs from verified archive')
        else:
            binary.write_bytes(contents)
    binary.chmod(0o755)
    version = subprocess.check_output([str(binary), '-version'], text=True).strip()
    dbdir = CACHE/'db'
    dbdir.mkdir(exist_ok=True)
    extracted = []
    with tarfile.open(CACHE/'pataa.tar.gz', 'r:gz') as archive:
        for member in archive:
            if (not member.isfile() or Path(member.name).name != member.name
                    or not (member.name.startswith('pataa.') or member.name in (
                        'taxdb.btd', 'taxdb.bti', 'taxonomy4blast.sqlite3'))):
                raise ValueError('Unexpected archive entry: '+member.name)
            target = dbdir/member.name
            h = hashlib.sha256()
            if target.exists():
                with archive.extractfile(member) as src:
                    for block in iter(lambda: src.read(4*1024*1024), b''):
                        h.update(block)
                if target.stat().st_size != member.size or sha(target) != h.hexdigest():
                    raise ValueError('Existing extracted file differs: '+member.name)
            else:
                with archive.extractfile(member) as src, target.open('xb') as dst:
                    for block in iter(lambda: src.read(4*1024*1024), b''):
                        dst.write(block)
                        h.update(block)
            extracted.append(dict(name=member.name, bytes=member.size, sha256=h.hexdigest()))
    info = subprocess.check_output([str(binary), '-db', str(dbdir/'pataa'), '-info'], text=True)
    (HERE/'blastdbcmd_info.txt').write_text(info)
    fasta = CACHE/'pataa.fasta'
    command = [str(binary), '-db', str(dbdir/'pataa'), '-entry', 'all', '-outfmt', '%f', '-out', str(fasta)]
    with (HERE/'export.log').open('w') as log:
        subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, check=True, timeout=600)
    count = letters = canonical = 0
    alphabet = Counter()
    positive = None
    name = None
    parts = []

    def account(header, sequence):
        nonlocal count, letters, canonical, positive
        count += 1
        letters += len(sequence)
        alphabet.update(sequence)
        is_canonical = not set(sequence)-AA20
        canonical += is_canonical
        if positive is None and is_canonical and 60 <= len(sequence) <= 250:
            positive = dict(source_header=header, sequence=sequence, source_ordinal=count)

    with fasta.open() as stream:
        for line in stream:
            if line.startswith('>'):
                if name is not None:
                    account(name, ''.join(parts))
                name, parts = line[1:].strip(), []
            else:
                parts.append(line.strip())
    if name is not None:
        account(name, ''.join(parts))
    if count != metadata['number-of-sequences'] or letters != metadata['number-of-letters']:
        raise ValueError('Exported record/residue counts disagree with publisher metadata')
    if positive is None:
        raise ValueError('No preregistered positive control found')
    with (HERE/'PATENT_POSITIVE_CONTROL.csv').open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=['name','sequence'])
        writer.writeheader()
        writer.writerow(dict(name='KNOWN_PATENT_SOURCE_POSITIVE_CONTROL', sequence=positive['sequence']))
    receipt = dict(generated_at=datetime.now(timezone.utc).isoformat(),
        source_metadata=metadata, metadata_http_headers=headers, downloads=receipts, extracted_files=extracted,
        fasta_path=str(fasta), fasta_sha256=sha(fasta), fasta_bytes=fasta.stat().st_size,
        exported_records=count, exported_residues=letters, canonical_records=canonical,
        alphabet=dict(sorted(alphabet.items())), duplicates='Preserved; no unique-sequence count claimed',
        sequence_transformations='None; official blastdbcmd FASTA export passed unchanged to MMseqs2',
        exporter=dict(path=str(binary), sha256=sha(binary), version=version, command=command),
        positive_control=positive, scope_complete=False,
        scope='NCBI pataa Patent-division protein records only; not verified as organizer patent corpus')
    save('ACQUISITION.json', receipt)
    print(json.dumps(dict(status='acquired_and_export_verified', records=count,
        residues=letters, canonical_records=canonical, new_gpu_spend_usd=0)))


if __name__ == '__main__':
    main()
