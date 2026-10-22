"""Download named public sequence snapshots only; no database substitutions."""
import concurrent.futures
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import urllib.request

OUT = Path(__file__).resolve().parent
CACHE = Path.home() / '.cache/adaptyv-c4/novelty'
CACHE.mkdir(parents=True, exist_ok=True)
SOURCES = {
 'swissprot': 'https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz',
 'pdb': 'https://files.rcsb.org/pub/pdb/derived_data/pdb_seqres.txt.gz',
 'therapeutic_antibody': 'https://opig.stats.ox.ac.uk/webapps/sabdab-sabpred/static/downloads/TheraSAbDab_SeqStruc_OnlineDownload.csv',
 'plabdab_paired': 'https://opig.stats.ox.ac.uk/webapps/plabdab/static/downloads/paired_sequences.csv.gz',
 'plabdab_unpaired': 'https://opig.stats.ox.ac.uk/webapps/plabdab/static/downloads/unpaired_sequences.csv.gz',
}

def fetch(item):
    key, url = item
    dest = CACHE / (key + ('.csv.gz' if 'sequences.csv.gz' in url else '.fasta.gz' if url.endswith('.gz') else '.csv'))
    rec = dict(database=key, url=url, path=str(dest))
    try:
        with urllib.request.urlopen(url, timeout=120) as response, dest.with_suffix(dest.suffix + '.partial').open('wb') as fh:
            rec['http_headers'] = dict(response.headers)
            size = 0
            while True:
                chunk = response.read(1024*1024)
                if not chunk: break
                size += len(chunk)
                if size > 300*1024*1024: raise ValueError('CPU download budget300MiB per source exceeded')
                fh.write(chunk)
        dest.with_suffix(dest.suffix + '.partial').replace(dest)
        rec.update(status='downloaded', bytes=size, sha256=hashlib.sha256(dest.read_bytes()).hexdigest())
    except Exception as e:
        rec.update(status='unavailable', error=repr(e))
    (OUT / ('download_' + key + '.json')).write_text(json.dumps(rec, indent=2)+'\n')
    return rec

if __name__ == '__main__':
    with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:
        records = list(pool.map(fetch, SOURCES.items()))
    (OUT / 'DOWNLOADS.json').write_text(json.dumps(dict(checked_at=datetime.now(timezone.utc).isoformat(), sources=records), indent=2)+'\n')
    print(json.dumps([dict(database=r['database'],status=r['status'],bytes=r.get('bytes'),error=r.get('error')) for r in records]))
