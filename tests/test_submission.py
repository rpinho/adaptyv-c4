import csv
import json
from pathlib import Path
import sys
import tempfile
import unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from build_submission import validate,build

class SubmissionTests(unittest.TestCase):
    def entry(self,**kw):return dict(name='original_1',sequence='ACDEFGHIKL',molecule_class='single_chain',**kw)
    def test_real_rejected_sequence_needs_dry_run(self):
        with self.assertRaises(ValueError):validate([self.entry()])
        self.assertEqual(len(validate([self.entry()],dry_run=True)),1)
    def test_format_and_rank_survive_roundtrip(self):
        with tempfile.TemporaryDirectory() as folder:
            p=Path(folder);selection=p/'selection.json';output=p/'dry_run.csv'
            first=self.entry();second=dict(name='original_2',sequence='LMNPQRSTVW',molecule_class='single_chain')
            selection.write_text(json.dumps({'designs':[second,first]}));audit=build(selection,output,dry_run=True)
            with output.open() as fh:rows=list(csv.DictReader(fh))
            self.assertEqual([r['name'] for r in rows],['original_2','original_1']);self.assertFalse(audit['submitted'])
    def test_duplicates_and_noncanonical_input_rejected(self):
        for entries in [[self.entry(),self.entry()],[dict(self.entry(),sequence='ACDEFGHIKX')],[],[dict(self.entry(),molecule_class='nanobody')]]:
            with self.subTest(entries=entries),self.assertRaises(ValueError):validate(entries,dry_run=True)
    def test_reviewed_export_and_count_limit(self):
        self.assertEqual(len(validate([self.entry(recommendation_status='recommended')])),1)
        with self.assertRaises(ValueError):validate([self.entry()]*21,dry_run=True)

if __name__=='__main__':unittest.main()
