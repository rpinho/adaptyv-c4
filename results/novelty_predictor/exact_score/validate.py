"""Validate exact score export, fresh execution and frozen-alignment comparison."""
from collections import Counter
import csv
import importlib.util
import json
from pathlib import Path

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
spec=importlib.util.spec_from_file_location('novelty',ROOT/'scripts/novelty_level_predictor.py')
novelty=importlib.util.module_from_spec(spec);spec.loader.exec_module(novelty)


def best(row):
    hs=[d['best_hit'] for d in row.get('structure_evidence',{}).get('domains',[]) if d.get('best_hit')]
    return max(hs,key=lambda h:h['qtmscore']) if hs else None


def main():
    oldbase=HERE.parent/'sensitivity/full_pool_union'
    old=json.loads((oldbase/'predictions.json').read_text())
    new=json.loads((HERE/'run/predictions.json').read_text())
    domains=json.loads((oldbase/'structure/query_domains.json').read_text())
    cached=Path(json.loads((HERE/'CACHED_EXPORT_COMMAND.json').read_text())['output'])
    assert new['predictor_sha256']==novelty.digest(ROOT/'scripts/novelty_level_predictor.py')
    assert new['input_csv_sha256']==old['input_csv_sha256']
    assert new['structure_manifest_sha256']==old['structure_manifest_sha256']
    assert new['structure_database']['files']==old['structure_database']['files']
    assert new['structural_search']['exact_tmscore_export']==1
    export_cmd=new['structural_search']['commands'][-1]
    assert export_cmd[export_cmd.index('--exact-tmscore')+1]=='1'
    def without_scores(path):
        with path.open(newline='') as f:
            for row in csv.reader(f,delimiter='\t'):
                yield tuple(value for key,value in zip(novelty.FOLDSEEK_FIELDS,row)
                            if key not in ('qtmscore','ttmscore'))
    # The cache-only A/B changes no correspondence, sequence or candidate row.
    assert list(without_scores(oldbase/'structure/hits.tsv'))==list(without_scores(cached))
    a=novelty.parse_structure_hits(cached,domains)
    b=novelty.parse_structure_hits(HERE/'run/structure/hits.tsv',domains)
    flat=lambda hits:{(q,h['target']):h['qtmscore'] for q,hs in hits.items() for h in hs}
    assert flat(a)==flat(b),'Fresh exact pipeline differs from isolated exact rescore'
    previous={r['name']:r for r in old['designs']}
    compared=[]
    for r in new['designs']:
        assert not r.get('organizer_level_established',False)
        before=previous[r['name']];oh,nh=best(before),best(r)
        compared.append(dict(name=r['name'],approximate_local_level=before.get('local_predicted_level'),
            exact_local_level=r.get('local_predicted_level'),approximate_best_TM=oh['qtmscore'] if oh else None,
            exact_best_TM=nh['qtmscore'] if nh else None,exact_best_target=nh['target'] if nh else None))
    retired=[r for r in compared if r['name']!='KNOWN_TFR1_POSITIVE_CONTROL']
    control=next(r for r in compared if r['name']=='KNOWN_TFR1_POSITIVE_CONTROL')
    assert control['exact_local_level']==1 and control['exact_best_TM']>=.98
    agreement=json.loads((HERE/'CACHED_EXACT_VALIDATION.json').read_text())
    saved=json.loads((HERE.parent/'independent_tmalign/SAVED_ALIGNMENT_RESULTS.json').read_text())
    summary=dict(status='validated_exact_export',retired_sequences=len(retired),
        approximate_counts=dict(Counter(str(r['approximate_local_level']) for r in retired)),
        exact_counts=dict(Counter(str(r['exact_local_level']) for r in retired)),
        changed_local_calls=sum(r['approximate_local_level']!=r['exact_local_level'] for r in retired),
        score_only_AB_alignment_records_identical=True,fresh_pipeline_equals_cached_exact_scores=True,
        unique_query_target_pairs=len(flat(b)),raw_alignment_rows=new['structural_search']['raw_alignment_rows'],
        full_pipeline_seconds=new['structural_search']['elapsed_seconds'],positive_control=control,
        standalone_fixed_before_max_abs_delta=saved['fixed_maximum_absolute_difference'],
        standalone_fixed_exact_max_abs_delta=agreement['agreement_max_absolute_delta'],
        standalone_fixed_agreement_tolerance=agreement['tolerance'],
        standalone_fixed_agreement_pass=agreement['agreement_pass'],comparisons=compared,
        new_gpu_spend_usd=0)
    (HERE/'VALIDATION.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps({k:v for k,v in summary.items() if k not in ('comparisons','positive_control')}))


if __name__=='__main__':main()
