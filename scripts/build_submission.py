"""Strict Proteinbase CSV export for this project's single-chain protein lane."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import re

AA=set('ACDEFGHIKLMNPQRSTVWY')

def validate(entries,*,dry_run=False):
    if not 1<=len(entries)<=20:raise ValueError('Track3 requires1-20 selected designs')
    names=set();sequences=set();rows=[]
    for entry in entries:
        name=entry['name'];seq=entry['sequence'];kind=entry['molecule_class']
        if not isinstance(name,str) or not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]{0,127}',name):raise ValueError('Invalid or unsafe design name')
        if name in names or seq in sequences:raise ValueError('Duplicate design name or sequence')
        if not isinstance(seq,str) or not 10<=len(seq)<=250 or set(seq)-AA:raise ValueError('Protein must contain10-250 canonical uppercase amino acids')
        if kind!='single_chain':raise ValueError('Live Proteinbase template requires single_chain here; antibody classes need separately validated annotations')
        if not dry_run and entry.get('recommendation_status')!='recommended':raise ValueError('Unreviewed design: use an explicitly labeled dry run')
        names.add(name);sequences.add(seq);rows.append({'name':name,'sequence':seq,'molecule_class':kind})
    return rows

def build(selection,output,*,dry_run=False):
    selection=Path(selection);output=Path(output)
    if output.exists() or output.with_suffix('.audit.json').exists():raise FileExistsError('Use fresh reviewable output paths')
    document=json.loads(selection.read_text());rows=validate(document['designs'],dry_run=dry_run)
    output.parent.mkdir(parents=True,exist_ok=True)
    with output.open('w',newline='') as fh:
        writer=csv.DictWriter(fh,fieldnames=['name','sequence','molecule_class']);writer.writeheader();writer.writerows(rows)
    with output.open(newline='') as fh:
        reader=csv.DictReader(fh)
        if reader.fieldnames!=['name','sequence','molecule_class'] or list(reader)!=rows:raise RuntimeError('CSV roundtrip differs')
    audit=dict(status='format_dry_run_not_for_upload' if dry_run else 'reviewed_export_for_Ricardo',
        row_count=len(rows),rank_order='Input order preserved; best first',
        selection_sha256=hashlib.sha256(selection.read_bytes()).hexdigest(),
        csv_sha256=hashlib.sha256(output.read_bytes()).hexdigest(),
        eligibility_scope='Format/count/sequence checks only. Novelty, diversity, evidence and organizer eligibility require separate review.',
        submitted=False,account_created=False)
    output.with_suffix('.audit.json').write_text(json.dumps(audit,indent=2)+'\n')
    return audit

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('selection');ap.add_argument('output');ap.add_argument('--dry-run',action='store_true');a=ap.parse_args()
    print(json.dumps(build(a.selection,a.output,dry_run=a.dry_run)))
