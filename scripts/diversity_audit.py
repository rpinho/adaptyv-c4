"""Within-pool sequence diversity; neither external nor structural novelty proof."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
from Bio.Align import PairwiseAligner

def audit(rows,threshold=.8):
    names=[r['name'] for r in rows];seqs=[r['sequence'] for r in rows]
    if len(set(names))!=len(names) or not names:raise ValueError('Need distinct nonempty names')
    aligner=PairwiseAligner(mode='global',match_score=2,mismatch_score=-1,open_gap_score=-3,extend_gap_score=-.5)
    parent=list(range(len(rows)));pairs=[]
    def root(i):
        while parent[i]!=i:parent[i]=parent[parent[i]];i=parent[i]
        return i
    for i in range(len(rows)):
        for j in range(i):
            alignment=aligner.align(seqs[i],seqs[j])[0];idx=alignment.indices
            identical=sum(a>=0 and b>=0 and seqs[i][a]==seqs[j][b] for a,b in zip(*idx))
            identity=identical/idx.shape[1]
            pairs.append(dict(name_a=names[j],name_b=names[i],global_identity=identity,alignment_columns=idx.shape[1]))
            if identity>=threshold:parent[root(i)]=root(j)
    clusters={}
    for i,name in enumerate(names):clusters.setdefault(root(i),[]).append(name)
    families=sorted(clusters.values(),key=lambda x:(-len(x),x[0]))
    return dict(designs=len(rows),unique_sequences=len(set(seqs)),sequence_clusters=families,
        cluster_count=len(families),largest_cluster=len(families[0]),pairs=pairs,
        operational_identity_threshold=threshold,organizer_threshold=False,
        alignment='Global match2/mismatch-1/gapopen-3/extend-.5; first optimal alignment; identities/all columns',
        external_known_protein_novelty_assessed=False,structural_novelty_assessed=False,
        interpretation='Sequence clusters identify related project families; no eligibility or binding guarantee. Structural and external database comparisons remain required for a real shortlist.')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('input');ap.add_argument('output');a=ap.parse_args()
    p=Path(a.input)
    with p.open() as fh:rows=list(csv.DictReader(fh))
    result=audit(rows);result['input_sha256']=hashlib.sha256(p.read_bytes()).hexdigest()
    Path(a.output).write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:result[k] for k in ['designs','unique_sequences','cluster_count','largest_cluster']}))
