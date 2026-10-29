"""The bounded-memory reader must preserve exact top scores and duplicate ties."""
import csv
from pathlib import Path
import random
import sys
import tempfile
import unittest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'results/novelty_predictor/adaptive_completion'))
from native import n,stream_best


class Streaming(unittest.TestCase):
    def test_matches_materialized_parser_with_reappearing_targets_and_ties(self):
        domains={'q1':dict(start=1,end=10),'q2':dict(segments=[[1,4],[7,12]])}
        rng=random.Random(0)
        with tempfile.TemporaryDirectory() as folder:
            path=Path(folder)/'hits.tsv'
            with path.open('w') as stream:
                writer=csv.writer(stream,delimiter='\t',lineterminator='\n')
                for i in range(300):
                    query=rng.choice(list(domains));target='target'+str(rng.randrange(25))
                    score=rng.choice([.1,.2,.4,.5,.6,.8,.9,1.])
                    row=dict(query=query,target=target,qlen=10,tlen=10,qstart=1,qend=10,tstart=1,tend=10,
                        alnlen=10,qtmscore=score,ttmscore=score,qaln='A'*10,taln='A'*10)
                    writer.writerow([row[key] for key in n.FOLDSEEK_FIELDS])
            expected=n.parse_structure_hits(path,domains)
            actual,counts=stream_best([path],domains,keep=3)
            self.assertEqual(actual,{q:records[:3] for q,records in expected.items()})
            self.assertEqual(sum(counts[str(path)].values()),300)


if __name__=='__main__':unittest.main()
