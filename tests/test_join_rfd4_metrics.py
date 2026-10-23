import csv
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from join_rfd4_metrics import join


class JoinMetrics(unittest.TestCase):
    def test_exact_id_join_preserves_order_missing_rows_and_sequence_values(self):
        with tempfile.TemporaryDirectory() as folder:
            folder=Path(folder);seq=folder/'seq.csv';metric=folder/'metrics.csv';out=folder/'joined.csv'
            seq.write_text('name,sequence\nb,WWWW\na,AAAA\nc,CCCC\n')
            metric.write_text('example_id,measurement\na,1\nb,0\n')
            result=join(seq,metric,'arm',out)
            with out.open() as stream:rows=list(csv.DictReader(stream))
            self.assertEqual([r['source_design_id'] for r in rows],['b','a','c'])
            self.assertEqual([r['sequence'] for r in rows],['WWWW','AAAA','CCCC'])
            self.assertEqual(rows[-1]['measurement'],'')
            self.assertEqual(rows[-1]['metric_record_status'],'missing')
            self.assertEqual(result['missing_metric_ids'],['c'])
            self.assertFalse(result['ranking_performed'])

    def test_duplicate_or_unmatched_metric_ids_never_select_a_favorable_row(self):
        with tempfile.TemporaryDirectory() as folder:
            folder=Path(folder);seq=folder/'seq.csv';metric=folder/'metric.csv'
            seq.write_text('name,sequence\na,AAAA\n')
            for text in ('example_id,value\na,1\na,0\n','example_id,value\nunknown,1\n'):
                metric.write_text(text)
                with self.assertRaises(ValueError):join(seq,metric,'arm',folder/'out.csv')
                self.assertFalse((folder/'out.csv').exists())


if __name__=='__main__':unittest.main()
