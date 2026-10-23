"""Behavioral check: duplicated original rows versus unique zero-seed rows."""
import json
from native import HERE,DOMAIN,BINARY,CACHE,n,CandidateDB,lookup,run,source_audit


def main():
    out=HERE/'seed_check'
    if out.exists():raise FileExistsError('Fresh behavioral check required')
    out.mkdir();source=source_audit();steps=[]
    query=DOMAIN/'querydb';db=json.loads((n.Path('results/novelty_predictor/structural/DATABASE.json')).read_text())
    target=n.Path(db['prefix']);domains=json.loads((DOMAIN/'query_domains.json').read_text())
    keys=lookup(query);candidates=CandidateDB(DOMAIN/'discovered')
    chosen=[]
    for q in sorted(domains):
        rows=candidates.records(keys[q])[:12]
        if rows:chosen.append((q,rows))
        if len(chosen)==8:break
    assert len(chosen)==8
    original=out/'duplicated_original.tsv';unique=out/'unique_zero_seed.tsv'
    with original.open('wb') as a,unique.open('w') as b:
        for q,rows in sorted(chosen,key=lambda item:keys[item[0]]):
            for row in rows+rows:a.write(str(keys[q]).encode()+b'\t'+row+b'\n')
            for ref in sorted({int(row.split(b'\t',1)[0]) for row in rows}):b.write(f'{keys[q]}\t{ref}\t0\t0\n')
    candidates.close()
    observed={}
    for label,seeds,kind in [('original',original,'5'),('deduplicated',unique,'7')]:
        prefix=out/(label+'_prefilter');aligned=out/(label+'_aligned');tsv=out/(label+'_hits.tsv')
        commands=[
            [BINARY,'tsv2db',seeds,prefix,'--output-dbtype',kind],
            [BINARY,'tmalign',query,target,prefix,aligned,'--threads','4','--tmalign-fast','0',
             '--tmscore-threshold','0','--tmscore-threshold-mode','1','--tmalign-hit-order','1','-a','1'],
            [BINARY,'convertalis',query,target,aligned,tsv,'--threads','4','--exact-tmscore','1',
             '--format-output',','.join(n.FOLDSEEK_FIELDS)]]
        for i,command in enumerate(commands):run(command,out/(label+'_'+str(i)+'.log'),180,steps)
        hits=n.parse_structure_hits(tsv,domains)
        observed[label]={(q,h['target']):(h['qtmscore'],h['ttmscore']) for q,records in hits.items() for h in records}
    assert observed['original'] and observed['original']==observed['deduplicated'],'Seed/deduplication changed actual exact scores'
    result=dict(state='passed',source=source,queries=[q for q,rows in chosen],unique_pairs=len(observed['original']),
        exact_query_and_target_tm_scores_equal=True,tolerance=0,steps=steps,
        raw_hits_sha256={label:n.digest(out/(label+'_hits.tsv')) for label in observed},
        scope='Deterministic first8 nonempty domain queries, first12 candidate rows, deliberately duplicated; no global numerical-equivalence claim beyond this fixture and source audit',
        new_gpu_or_remote_cpu_spend_usd=0)
    n.write_json(HERE/'SEED_CHECK.json',result)
    print(json.dumps({k:result[k] for k in ('state','unique_pairs','exact_query_and_target_tm_scores_equal','tolerance')}))


if __name__=='__main__':main()
