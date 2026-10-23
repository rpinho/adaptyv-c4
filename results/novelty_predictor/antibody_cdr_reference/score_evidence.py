"""Conditional antibody-scale arithmetic on measured reference-loop identities."""
import collections
from datetime import datetime,timezone
from fractions import Fraction
from functools import lru_cache
import hashlib
import json
from pathlib import Path
import sqlite3
import time

from prepare_pool import HERE,ROOT,DB,INPUT,AUDIT,sha,table_hash

CDR_FAMILIAR_MIN_INCLUSIVE=Fraction(7,10)
CDR_KNOWN_MIN_INCLUSIVE=Fraction(19,20)
GLOBAL_FAMILIAR_MIN_INCLUSIVE=Fraction(7,10)
GLOBAL_KNOWN_MIN_INCLUSIVE=Fraction(19,20)
PRIMARY_IDENTITY='levenshtein_max_length'
NORMALIZERS=(PRIMARY_IDENTITY,'levenshtein_query_length','levenshtein_reference_length','normalized_indel')


@lru_cache(maxsize=200000)
def edit_distance(a,b,substitution_cost=1):
    previous=list(range(len(b)+1))
    for i,letter in enumerate(a,1):
        current=[i]
        for j,other in enumerate(b,1):
            current.append(min(current[-1]+1,previous[j]+1,previous[j-1]+(0 if letter==other else substitution_cost)))
        previous=current
    return previous[-1]


def identities(a,b):
    if not a or not b:raise ValueError('A missing/empty loop is not measured zero identity')
    distance=edit_distance(a,b)
    indel=edit_distance(a,b,2)
    return dict(zip(NORMALIZERS,(
        Fraction(max(len(a),len(b))-distance,max(len(a),len(b))),
        max(Fraction(0),Fraction(len(a)-distance,len(a))),
        max(Fraction(0),Fraction(len(b)-distance,len(b))),
        Fraction(len(a)+len(b)-indel,len(a)+len(b)))))


def conditional_level(cdr_identity,global_identity):
    if cdr_identity>=CDR_KNOWN_MIN_INCLUSIVE:return 1
    if cdr_identity>=CDR_FAMILIAR_MIN_INCLUSIVE:
        return 1 if global_identity>=GLOBAL_KNOWN_MIN_INCLUSIVE else 2
    return 3 if global_identity>=GLOBAL_FAMILIAR_MIN_INCLUSIVE else 4


def completion_levels(observed,global_identity,missing_reference_cdrs):
    if not missing_reference_cdrs and observed is not None:
        return [conditional_level(observed,global_identity)]
    lower=observed if observed is not None else Fraction(0)
    candidates={lower,Fraction(1)}|{x for x in (CDR_FAMILIAR_MIN_INCLUSIVE,CDR_KNOWN_MIN_INCLUSIVE) if x>=lower}
    return sorted({conditional_level(c,global_identity) for c in candidates})


def exact(value):
    return None if value is None else dict(numerator=value.numerator,denominator=value.denominator,value=float(value))


def controls():
    raw=json.loads((AUDIT/'control_raw.json').read_text())['CONTROL_12E8_H']
    loop=''.join(aa for (position,insertion),aa in raw['numbering'] if 105<=position<=117 and aa!='-')
    self_scores=identities(loop,loop)
    assert all(v==1 and conditional_level(v,Fraction(1))==1 for v in self_scores.values())
    unequal=identities('AAAAAAAAAA','AAAAAAAAAAAAAAA')
    assert unequal[PRIMARY_IDENTITY]==Fraction(2,3) and unequal['normalized_indel']==Fraction(4,5)
    assert conditional_level(Fraction(7,10),Fraction(7,10))==2
    assert conditional_level(Fraction(19,20),Fraction(0))==1
    assert conditional_level(Fraction(69,100),Fraction(7,10))==3
    assert conditional_level(Fraction(69,100),Fraction(69,100))==4
    assert conditional_level(Fraction(7,10),Fraction(19,20))==1
    assert completion_levels(Fraction(3,5),Fraction(4,5),True)==[1,2,3]
    return dict(publisher_loop_self_identity={k:exact(v) for k,v in self_scores.items()},
        self_conditional_level=1,unequal_length_normalization_control={k:exact(v) for k,v in unequal.items()},
        control_annotation_sha256=sha(AUDIT/'control_raw.json'),boundary_controls_passed=True)


def main():
    start=time.monotonic()
    annotation_path=HERE/'REFERENCE_ANNOTATION_CHECKPOINT.json'
    if not annotation_path.exists():
        original_status=(HERE/'ANNOTATION_STATUS.json').read_bytes()
        if json.loads(original_status)['state']!='complete':raise RuntimeError('Reference annotation is not complete')
        annotation_path.write_bytes(original_status)
    annotation=json.loads(annotation_path.read_text())
    if annotation['state']!='complete':raise RuntimeError('Reference annotation is not complete')
    if (HERE/'PREDICTIONS.json').exists():raise FileExistsError('Conditional evidence already frozen')
    pool=json.loads((HERE/'POOL.json').read_text())
    assert sha(INPUT/'analysis_inputs.csv')==pool['input_csv_sha256']
    assert sha(AUDIT/'batch_raw.json')==pool['query_annotation_raw_sha256']
    connection=sqlite3.connect(DB)
    assert table_hash(connection,'SELECT id,sequence_sha256,canonical FROM heavy ORDER BY id')==pool['immutable_heavy_rows_sha256']
    assert table_hash(connection,'SELECT q,heavy FROM query_heavy ORDER BY q,heavy')==pool['immutable_query_heavy_rows_sha256']
    software_controls=controls()
    results=[]
    for number,name,arm,query_loop,global_json in connection.execute('SELECT id,name,arm,cdrh3,global_json FROM query ORDER BY id'):
        original=json.loads(global_json);hit=original.get('best_hit')
        global_identity=Fraction(hit['nident'],hit['alnlen']) if hit else None
        references=connection.execute('''SELECT h.id,h.sequence_sha256,h.cdrh3,h.annotation_state,
            h.chain_type,h.flanks_complete,h.span_verified,h.error FROM query_heavy q JOIN heavy h ON h.id=q.heavy
            WHERE q.q=? ORDER BY h.id''',(number,)).fetchall()
        usable={};gaps=collections.Counter();valid_reference_count=0
        for key,sequence_sha,loop,state,kind,flanks,span,error in references:
            if state!='done':gaps[state]+=1
            elif kind!='H':gaps['not_an_error_free_H_annotation']+=1
            elif not span:gaps['source_span_not_verified']+=1
            elif not flanks:gaps['missing_cdr3_flank']+=1
            elif not loop:gaps['empty_cdr3_annotation']+=1
            else:
                valid_reference_count+=1
                usable.setdefault(loop,dict(heavy_reference_id=key,full_reference_sequence_sha256=sequence_sha))
        best={key:None for key in NORMALIZERS};compared=0
        for loop,ref in usable.items():
            scores=identities(query_loop,loop);compared+=1
            for convention,value in scores.items():
                if best[convention] is None or value>best[convention]['fraction']:
                    best[convention]=dict(**ref,fraction=value,reference_cdrh3=loop,
                        reference_cdrh3_sha256=hashlib.sha256(loop.encode()).hexdigest(),
                        query_cdr_length=len(query_loop),reference_cdr_length=len(loop),
                        unit_edit_distance=edit_distance(query_loop,loop),indel_distance=edit_distance(query_loop,loop,2))
            if all(b['fraction']==1 for b in best.values()):break
        evidence={};all_completion_levels=set()
        missing=bool(gaps) or not usable
        for convention,match in best.items():
            fraction=match['fraction'] if match else None
            observed_level=conditional_level(fraction,global_identity) if fraction is not None and global_identity is not None else None
            possible=completion_levels(fraction,global_identity,missing) if global_identity is not None else [1,2,3,4]
            all_completion_levels.update(possible)
            exported={k:v for k,v in match.items() if k!='fraction'} if match else None
            if exported is not None:
                exported['identity']=exact(fraction)
                exported['actual_query_hit_associations']=[dict(returned_target=t,source_route=s) for t,s in connection.execute(
                    'SELECT a.target,a.source FROM association a JOIN candidate c ON a.target=c.target WHERE c.q=? AND a.heavy=? ORDER BY a.target,a.source LIMIT 5',
                    (number,exported['heavy_reference_id']))]
            evidence[convention]=dict(best_observed=exported,
                level_over_measured_reference_subset=observed_level,
                conditional_levels_allowing_unmeasured_reference_cdrs=possible)
        point=next(iter(all_completion_levels)) if len(all_completion_levels)==1 else None
        results.append(dict(name=name,source_arm=arm,query_id=f'q{number:06d}',
            complete_query_cdrh3=query_loop,query_cdrh3_sha256=hashlib.sha256(query_loop.encode()).hexdigest(),
            observed_global_identity=exact(global_identity),best_global_hit=hit,
            returned_heavy_reference_candidates=len(references),valid_complete_reference_annotations=valid_reference_count,
            distinct_usable_reference_cdrs=len(usable),reference_cdrs_actually_compared=compared,
            maximum_proven_by_identity_one=bool(usable) and all(b and b['fraction']==1 for b in best.values()),
            reference_annotation_gaps=dict(gaps),normalizer_evidence=evidence,
            conditional_predicted_level=point,conditional_possible_levels=sorted(all_completion_levels),
            organizer_level=None,organizer_antibody_type_established=False,qualifies_for_submission=False))
        (HERE/'SCORING_STATUS.json').write_text(json.dumps(dict(state='running',completed_queries=len(results),total_queries=50,elapsed_seconds=time.monotonic()-start),indent=2)+'\n')
    connection.close()
    assert len(results)==50
    report=dict(state='complete_conditional_evidence',generated_at=datetime.now(timezone.utc).isoformat(),
        query_count=50,total_design_inputs=400,primary_normalizer=PRIMARY_IDENTITY,
        thresholds={k:exact(v) for k,v in dict(CDR_FAMILIAR_MIN_INCLUSIVE=CDR_FAMILIAR_MIN_INCLUSIVE,
            CDR_KNOWN_MIN_INCLUSIVE=CDR_KNOWN_MIN_INCLUSIVE,GLOBAL_FAMILIAR_MIN_INCLUSIVE=GLOBAL_FAMILIAR_MIN_INCLUSIVE,
            GLOBAL_KNOWN_MIN_INCLUSIVE=GLOBAL_KNOWN_MIN_INCLUSIVE).items()},
        controls=software_controls,designs=results,pool_manifest_sha256=sha(HERE/'POOL.json'),
        reference_annotation_checkpoint=str(annotation_path),reference_annotation_status_sha256=sha(annotation_path),
        script_sha256=sha(__file__),elapsed_seconds=time.monotonic()-start,
        organizer_levels_established=0,new_gpu_or_remote_cpu_spend_usd=0,
        limitations=['Conditional on antibody scoring; not organizer structural typing or VH/VHH determination',
            'Organizer edit-distance convention, search settings and PLAbDab snapshot unknown',
            'All associated heavy sequences restored for actually returned hits only; finite MMseqs retrieval/caps remain',
            'One model-selected H variable domain per reference, with source span and flanks checked; additional domains not exhaustively numbered',
            'Best global identity is the existing observed alignment identity; organism/coverage policy and missing patent corpus remain gaps',
            'A level over usable reference subsets is arithmetic on measured data, not a prediction under missing CDR annotations',
            'Numeric conditional prediction requires invariance across declared normalizers and possible missing-reference completions',
            'No binding, pH selectivity, assay success, ranked submission or organizer acceptance established'])
    (HERE/'PREDICTIONS.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    (HERE/'SCORING_STATUS.json').write_text(json.dumps(dict(state='complete',completed_queries=50,elapsed_seconds=time.monotonic()-start),indent=2)+'\n')
    print(json.dumps(dict(queries=50,conditional_point_levels=dict(collections.Counter(str(d['conditional_predicted_level']) for d in results)),organizer_levels_established=0)))


if __name__=='__main__':main()
