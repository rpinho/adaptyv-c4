import csv
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
import rank_candidates as ranker


class GeometryRanking(unittest.TestCase):
    def test_hotspot_order_does_not_synthesize_confidence_or_weighted_scores(self):
        fields=[key for key,_ in ranker.GEOMETRY_ORDER]
        records=[dict(source_file='same',source_id='a',metrics=dict(zip(fields,['1','.9','.8']))),
                 dict(source_file='same',source_id='b',metrics=dict(zip(fields,['.5','0','0']))),
                 dict(source_file='same',source_id='c',metrics=dict(zip(fields,['NaN','0','0'])))]
        ranker.rank_geometry(records);records.sort(key=ranker.geometry_sort_key)
        self.assertEqual([r['source_id'] for r in records],['a','b','c'])
        self.assertTrue(all(r['rank_score'] is None for r in records))
        self.assertEqual(records[-1]['ranking_missing_fields'],[fields[0]])

    def test_missing_geometry_is_retained_beyond_fifty_deferrals(self):
        with tempfile.TemporaryDirectory() as folder:
            folder=Path(folder);source=folder/'input.csv';out=folder/'shortlist.csv'
            fields=['name','sequence']+[key for key,_ in ranker.GEOMETRY_ORDER]
            with source.open('w',newline='') as stream:
                writer=csv.DictWriter(stream,fieldnames=fields);writer.writeheader()
                for i in range(61):writer.writerow(dict(name='missing'+str(i),sequence='A'*(i+1)))
                writer.writerow(dict(zip(fields,['measured','W'*30,'1','0','0'])))
            result=subprocess.run([sys.executable,str(ROOT/'scripts/rank_candidates.py'),str(source),'-o',str(out),'--ranking-policy','rfd4-geom'],capture_output=True,text=True)
            self.assertEqual(result.returncode,0,result.stderr)
            report=json.loads(out.with_suffix('.provenance.json').read_text())
            self.assertEqual(report['counts']['ranking_incomplete'],61)
            self.assertEqual(len(report['deferred']),61)
            self.assertEqual(report['counts']['shortlisted'],1)
            self.assertEqual(report['ranking_weights'],{})
            self.assertFalse(report['ranking_is_official_acceptance'])
            self.assertEqual(report['shortlist'][0]['rank_score'],None)
            self.assertNotIn('diversity cap held back',result.stdout)

    def test_legacy_finite_confidence_scores_still_use_the_existing_weights(self):
        records=[dict(metrics={key:0 for key in ranker.RANK_WEIGHTS}),
                 dict(metrics={key:1 for key in ranker.RANK_WEIGHTS})]
        ranker.rank_scores(records)
        self.assertAlmostEqual(records[0]['rank_score'],.2)
        self.assertAlmostEqual(records[1]['rank_score'],.8)
        self.assertIsNone(ranker.fnum({'value':'inf'},'value'))


if __name__=='__main__':unittest.main()
