"""Materialize all checked high/medium domains, preserving discontinuities."""
import csv
import argparse
import hashlib
import json
from pathlib import Path
import sys

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
CACHE=Path.home()/'.cache/adaptyv-c4/novelty_domain_source_audit'
sys.path.insert(0,str(ROOT/'scripts'))
import novelty_level_predictor as n


def domains(chopping):
    if chopping=='na': return []
    return [dict(segments=[list(map(int,fragment.split('-'))) for fragment in domain.split('_')])
            for domain in chopping.split(',')]


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--input-manifest',type=Path,default=HERE/'INPUTS.json')
    ap.add_argument('--model-run',type=Path,default=HERE/'all')
    ap.add_argument('--designs',type=Path)
    ap.add_argument('--parent-structures',type=Path)
    ap.add_argument('--out',type=Path,default=HERE/'STRUCTURES.json')
    args=ap.parse_args();benchmark=args.designs is None
    if args.out.exists():raise FileExistsError('Domain manifest already frozen')
    state=json.loads((args.model_run/'STATUS.json').read_text());assert state['state']=='complete'
    assert n.digest(args.model_run/'consensus.tsv')==state['consensus_sha256']
    indexing=json.loads((HERE/'INDEXING_CHECKS.json').read_text());assert indexing['all_passed']
    inventory=json.loads(args.input_manifest.read_text())
    assert n.digest(HERE/'PROTOCOL.md')==inventory['protocol_sha256']
    sources={r['query']:r for r in inventory['inputs'] if r['status']=='full_backbone_verified'}
    rows=[line.split('\t') for line in (args.model_run/'consensus.tsv').read_text().splitlines()]
    assert len(rows)==len(sources) and {r[0] for r in rows}==set(sources)
    input_csv=args.designs or ROOT/'results/novelty_predictor/patent_extension/analysis_inputs.csv'
    sequences={r['name']:r['sequence'] for r in csv.DictReader(input_csv.open())}
    specs=json.loads((args.parent_structures or ROOT/'results/novelty_predictor/patent_extension/STRUCTURES.json').read_text())
    by_name={r['name']:r for r in specs['structures']}
    census=[];total=discontinuous=0
    for row in sorted(rows):
        query,md5,nres,nhigh,nmedium,nlow,high,medium,low=row
        source=sources[query]; sequence=sequences[source['name']]
        assert md5==hashlib.md5(sequence.encode()).hexdigest() and int(nres)==len(sequence)
        assert n.digest(source['path'])==source['sha256']
        spec=by_name[source['name']]
        assert spec['status']=='available' and spec['sha256']==source['source']['sha256']
        accepted=[]
        for quality,text,expected_count in [('high',high,int(nhigh)),('medium',medium,int(nmedium))]:
            parsed=domains(text);assert len(parsed)==expected_count
            for d in parsed:
                if len(d['segments'])==1:
                    start,end=d.pop('segments')[0];d.update(start=start,end=end)
                d.update(source='Local three-method TED-style consensus; published indexing bug explicitly corrected; organizer versions unverified',
                    consensus_quality=quality,consensus_query=query,raw_consensus=text)
                accepted.append(d)
        n.structural_class(len(sequence),[dict(**n.domain_coordinates(d),tm_score=0.,normalization='query_domain_length') for d in accepted])
        covered=sum(len(n.domain_positions(d,len(sequence))) for d in accepted)
        discont=any('segments' in d for d in accepted)
        total+=len(accepted);discontinuous+=discont
        spec['domains']=accepted
        if not accepted:
            spec.update(status='unavailable',reason='No high/medium consensus domain assigned; structural level unestablished')
        census.append(dict(name=source['name'],query=query,length=len(sequence),accepted_domains=len(accepted),
            accepted_residues=covered,domain_coverage=covered/len(sequence),discontinuous_domain=discont,
            high=high,medium=medium,low_excluded=low))
    if benchmark:assert len(sources)==68 and total==70 and discontinuous==4
    weights=json.loads((HERE/'WEIGHTS.json').read_text())
    specs['segmentation_provenance']=dict(method='Local three-method TED-style consensus',
        ted_source_commit='f07e24ef4a19a722f2c3c3001cbd03d477ba3f31',
        merizo_weights_commit=weights['merizo_weights_commit'],
        chainsaw_weights_sha256=n.digest(CACHE/'ted-pinned/programs/chainsaw/saved_models/secondary_structure_epoch17/version_2/weights.pt'),
        unidoc_source_url='https://yanglab.qd.sdu.edu.cn/UniDoc/download/UniDoc_20250514.tgz',
        native_build_receipt_sha256=n.digest(HERE/'NATIVE_BUILD.json'),
        indexing_check_receipt_sha256=n.digest(HERE/'INDEXING_CHECKS.json'),
        consensus_sha256=n.digest(args.model_run/'consensus.tsv'),
        predictor_tables_sha256={method:n.digest(args.model_run/('chopping_'+method+'.txt')) for method in ['chainsaw','merizo','unidoc']},
        accepted_consensus_levels=['high_3_methods','medium_2_methods'],
        excluded='Low one-method domains',model_inputs='Provided sequence-verified coordinates; no structure inference',
        merizo_unidoc_dependency='UniDoc receives Merizo non-domain-residue filtering',
        input_models=len(sources),query_domains=total,inputs_with_discontinuous_domains=discontinuous,
        organizer_implementation_equivalence_established=False)
    n.write_json(args.out,specs)
    n.write_json(args.out.parent/'DOMAIN_CENSUS.json',dict(inputs=census,provenance=specs['segmentation_provenance']))
    with (args.out.parent/'DOMAIN_CENSUS.csv').open('w',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(census[0]));writer.writeheader();writer.writerows(census)
    print(json.dumps(dict(input_models=len(sources),query_domains=total,discontinuous_inputs=discontinuous,
        unassigned_inputs=sum(r['accepted_domains']==0 for r in census))))


if __name__=='__main__':main()
