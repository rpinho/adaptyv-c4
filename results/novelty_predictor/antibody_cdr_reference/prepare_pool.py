"""Restore all heavy-reference associations for actual returned PLAbDab hits."""
import collections
import csv
from datetime import datetime, timezone
import gzip
import hashlib
import json
from pathlib import Path
import re
import sqlite3
import sys

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
CACHE=Path.home()/'.cache/adaptyv-c4/novelty_antibody_scope'
DB=CACHE/'cdr_reference.sqlite3'
INPUT=ROOT/'results/novelty_predictor/incoming_rfd4_20261022'
AUDIT=ROOT/'results/novelty_predictor/antibody_scope_audit'


def sha(path):
    result=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(1024*1024),b''):result.update(block)
    return result.hexdigest()


def table_hash(connection,sql):
    result=hashlib.sha256()
    for row in connection.execute(sql):
        result.update((json.dumps(row,separators=(',',':'))+'\n').encode())
    return result.hexdigest()


def fasta(path):
    name=None;parts=[]
    with path.open() as stream:
        for line in stream:
            if line.startswith('>'):
                if name is not None:yield name,''.join(parts)
                name=line[1:].split()[0];parts=[]
            else:parts.append(line.strip())
        if name is not None:yield name,''.join(parts)


def main():
    if DB.exists() or (HERE/'POOL.json').exists():raise FileExistsError('Fresh reference pool required')
    state=dict(state='preparing_pool',started_at=datetime.now(timezone.utc).isoformat(),new_gpu_or_remote_cpu_spend_usd=0)
    def save():(HERE/'PREPARATION_STATUS.json').write_text(json.dumps(state,indent=2)+'\n')
    save()
    connection=None
    try:
        manifest_path=ROOT/'results/novelty_predictor/patent_extension/DATABASES.json'
        source=next(d for d in json.loads(manifest_path.read_text())['databases'] if d['id']=='plabdab')
        reference_fasta=Path(source['fasta_path']);assert sha(reference_fasta)==source['fasta_sha256']
        sources={d['database']:d for d in source['sources']}
        for item in sources.values():assert sha(item['path'])==item['sha256']
        rows=list(csv.DictReader((INPUT/'analysis_inputs.csv').open()))
        raw_path=AUDIT/'batch_raw.json';raw=json.loads(raw_path.read_text())
        validation=json.loads((AUDIT/'VALIDATION.json').read_text())
        assert sha(raw_path)==validation['raw_sha256']
        assert json.loads((AUDIT/'ANNOTATIONS.json').read_text())['input_sha256']==sha(INPUT/'analysis_inputs.csv')
        global_report=json.loads((INPUT/'run_bounded/whole_chain/predictions.json').read_text())
        assert global_report['input_csv_sha256']==sha(INPUT/'analysis_inputs.csv')
        globals_by_name={r['name']:r for r in global_report['designs']}
        query_rows={}
        for i,row in enumerate(rows,1):
            annotation=raw[row['name']]
            if annotation['chain_type']!='H' or annotation.get('error'):continue
            mapped={position[0]:aa for position,aa in annotation['numbering'] if not position[1].strip()}
            if mapped.get(104,'-')=='-' or mapped.get(118,'-')=='-':continue
            loop=''.join(aa for (position,insertion),aa in annotation['numbering'] if 105<=position<=117 and aa!='-')
            assert loop and re.fullmatch('[ACDEFGHIKLMNPQRSTVWY]+',loop)
            query_rows[i]=dict(row,cdrh3=loop,global_evidence=globals_by_name[row['name']])
        assert len(query_rows)==50
        raw_hits=INPUT/'run/whole_chain/plabdab.hits.tsv'
        recorded=json.loads((INPUT/'run/whole_chain/search_status.json').read_text())
        assert sha(raw_hits)==next(r['raw_hits_sha256'] for r in recorded if r['id']=='plabdab')
        candidates={q:set() for q in query_rows};scanned=0
        with raw_hits.open() as stream:
            for fields in csv.reader(stream,delimiter='\t'):
                assert len(fields)==15;scanned+=1
                number=int(fields[0][1:])
                if number in candidates:candidates[number].add(fields[1])
        wanted=set().union(*candidates.values())
        targets={name:sequence for name,sequence in fasta(reference_fasta) if name in wanted}
        assert set(targets)==wanted,'Returned reference ID absent from fingerprinted FASTA'
        by_sequence=collections.defaultdict(list)
        for name,sequence in targets.items():by_sequence[sequence].append(name)
        connection=sqlite3.connect(DB)
        connection.executescript('''
            PRAGMA journal_mode=WAL;
            CREATE TABLE query(id INTEGER PRIMARY KEY,name TEXT UNIQUE,arm TEXT,cdrh3 TEXT,global_json TEXT);
            CREATE TABLE target(id TEXT PRIMARY KEY,sequence TEXT);
            CREATE TABLE candidate(q INTEGER,target TEXT,PRIMARY KEY(q,target));
            CREATE TABLE heavy(id INTEGER PRIMARY KEY,sequence TEXT UNIQUE,sequence_sha256 TEXT,canonical INTEGER,
                annotation_state TEXT DEFAULT 'pending',chain_type TEXT,score REAL,error TEXT,raw BLOB,raw_sha256 TEXT,
                cdrh3 TEXT,flanks_complete INTEGER,span_verified INTEGER);
            CREATE TABLE association(target TEXT,heavy INTEGER,source TEXT,PRIMARY KEY(target,heavy,source));
            CREATE TABLE query_heavy(q INTEGER,heavy INTEGER,PRIMARY KEY(q,heavy));
        ''')
        connection.executemany('INSERT INTO query VALUES(?,?,?,?,?)',
            [(q,r['name'],r['source_arm'],r['cdrh3'],json.dumps(r['global_evidence'])) for q,r in query_rows.items()])
        connection.executemany('INSERT INTO target VALUES(?,?)',targets.items())
        for q,refs in candidates.items():
            connection.executemany('INSERT INTO candidate VALUES(?,?)',((q,t) for t in sorted(refs)))
        heavy_ids={}
        def associate(target_sequence,heavy_sequence,kind):
            if target_sequence not in by_sequence:return
            if heavy_sequence not in heavy_ids:
                canonical=bool(re.fullmatch('[ACDEFGHIKLMNPQRSTVWY]+',heavy_sequence))
                cursor=connection.execute('INSERT INTO heavy(sequence,sequence_sha256,canonical) VALUES(?,?,?)',
                    (heavy_sequence,hashlib.sha256(heavy_sequence.encode()).hexdigest(),int(canonical)))
                heavy_ids[heavy_sequence]=cursor.lastrowid
            key=heavy_ids[heavy_sequence]
            connection.executemany('INSERT OR IGNORE INTO association VALUES(?,?,?)',
                ((target,key,kind) for target in by_sequence[target_sequence]))
        source_counts=collections.Counter()
        with gzip.open(sources['plabdab_paired']['path'],'rt') as stream:
            for row in csv.DictReader(stream):
                heavy=row['heavy_sequence'].upper();light=row['light_sequence'].upper()
                associate(heavy,heavy,'paired_heavy');associate(light,heavy,'paired_light')
                source_counts['paired_rows']+=1
        with gzip.open(sources['plabdab_unpaired']['path'],'rt') as stream:
            for row in csv.DictReader(stream):
                source_counts['unpaired_rows']+=1
                chain=row['chain'];source_counts['unpaired_chain_'+chain]+=1
                if chain=='H':
                    sequence=row['GBSeq_sequence'].upper();associate(sequence,sequence,'unpaired_H')
        connection.execute('INSERT OR IGNORE INTO query_heavy SELECT c.q,a.heavy FROM candidate c JOIN association a ON c.target=a.target')
        connection.commit()
        immutable_heavy=table_hash(connection,'SELECT id,sequence_sha256,canonical FROM heavy ORDER BY id')
        immutable_pairs=table_hash(connection,'SELECT q,heavy FROM query_heavy ORDER BY q,heavy')
        counts={name:connection.execute('SELECT count(*) FROM '+name).fetchone()[0]
                for name in ('query','target','candidate','heavy','association','query_heavy')}
        associations=dict(connection.execute('SELECT source,count(*) FROM association GROUP BY source'))
        connection.execute('PRAGMA wal_checkpoint(TRUNCATE)');connection.close();connection=None
        pool=dict(state='pool_prepared',database=str(DB),counts=counts,source_counts=source_counts,
            association_routes=associations,raw_plab_alignment_rows_scanned=scanned,
            input_csv_sha256=sha(INPUT/'analysis_inputs.csv'),query_annotation_raw_sha256=sha(raw_path),
            sequence_hit_sha256=sha(raw_hits),reference_manifest_sha256=sha(manifest_path),
            source_fasta_sha256=source['fasta_sha256'],source_downloads=sources,
            immutable_heavy_rows_sha256=immutable_heavy,immutable_query_heavy_rows_sha256=immutable_pairs,
            prepared_sqlite_sha256=sha(DB),script_sha256=sha(__file__),
            selection='All returned PLAbDab hits for all50 complete H annotations; all associated raw-source heavy sequences restored',
            organizer_typing_or_cdr_level_established=False,new_gpu_or_remote_cpu_spend_usd=0)
        (HERE/'POOL.json').write_text(json.dumps(pool,indent=2)+'\n')
        state.update(state='complete',counts=counts,finished_at=datetime.now(timezone.utc).isoformat());save()
        print(json.dumps(dict(counts=counts,source_counts=source_counts,association_routes=associations)))
    except Exception as error:
        if connection is not None:connection.close()
        state.update(state='failed',error=repr(error));save();raise


if __name__=='__main__':main()
