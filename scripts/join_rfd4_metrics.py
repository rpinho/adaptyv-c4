"""Join existing RFD4 sequence/metric exports for analysis, without ranking.

Not a submission emitter. Missing metric rows are retained and explicitly marked.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path

AA20=set('ACDEFGHIKLMNPQRSTVWY')


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def indexed(path,key,required):
    with Path(path).open(newline='') as stream:
        reader=csv.DictReader(stream)
        fields=reader.fieldnames or []
        if not required<=set(fields):raise ValueError('Required columns missing: '+str(required-set(fields)))
        result={}
        for row in reader:
            name=row[key]
            if not name or name in result:raise ValueError('Empty/duplicate design ID in '+str(path))
            result[name]=row
    return fields,result


def join(sequences,metrics,arm,out):
    out=Path(out)
    receipt=out.with_suffix('.provenance.json')
    if out.exists() or receipt.exists():raise FileExistsError('Preserve the existing join and receipt')
    sequence_fields,seq=indexed(sequences,'name',{'name','sequence'})
    metric_fields,measured=indexed(metrics,'example_id',{'example_id'})
    if not seq:raise ValueError('No sequences')
    if set(measured)-set(seq):raise ValueError('Metric IDs have no matching sequence record')
    metadata=['source_arm','source_design_id','sequence_sha256','metric_record_status']
    if set(metric_fields)&(set(sequence_fields)|set(metadata)):
        raise ValueError('Metric columns collide with sequence/analysis columns')
    rows=[]
    for name,sequence_row in seq.items():
        sequence=sequence_row['sequence']
        if not sequence or set(sequence)-AA20:raise ValueError('Noncanonical or empty sequence; no silent repair')
        row=dict(sequence_row)
        row.update(name=arm+'__'+name,source_arm=arm,source_design_id=name,
                   sequence_sha256=hashlib.sha256(sequence.encode()).hexdigest(),
                   metric_record_status='present' if name in measured else 'missing')
        row.update(measured.get(name,{}));rows.append(row)
    out.parent.mkdir(parents=True,exist_ok=True)
    with out.open('w',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=sequence_fields+metadata+metric_fields,lineterminator='\n')
        writer.writeheader();writer.writerows(rows)
    result=dict(state='complete',arm=arm,inputs=len(rows),metric_records=len(measured),
                missing_metric_ids=[name for name in seq if name not in measured],
                input_sequence_sha256=sha(sequences),input_metrics_sha256=sha(metrics),output_sha256=sha(out),
                script_sha256=sha(__file__),sequence_order_preserved=True,ranking_performed=False,
                organizer_acceptance_measured=False,submission_format_emitted=False)
    receipt.write_text(json.dumps(result,indent=2)+'\n')
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--sequences',type=Path,required=True)
    parser.add_argument('--metrics',type=Path,required=True)
    parser.add_argument('--arm',required=True)
    parser.add_argument('--out',type=Path,required=True)
    args=parser.parse_args()
    print(json.dumps(join(args.sequences,args.metrics,args.arm,args.out)))


if __name__=='__main__':main()
