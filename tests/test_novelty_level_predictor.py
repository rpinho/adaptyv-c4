"""Boundary and missing-evidence tests, independent of GPU/generator work."""
import csv
import importlib.util
from pathlib import Path
import unittest

SPEC = importlib.util.spec_from_file_location('novelty_predictor', Path(__file__).resolve().parents[1] / 'scripts/novelty_level_predictor.py')
novelty = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(novelty)


def domain(end, tm):
    return dict(start=1, end=end, tm_score=tm, normalization='query_domain_length')


class NoveltyTests(unittest.TestCase):
    def test_sequence_boundaries(self):
        for identity, expected in [(0.30,'low'), (0.300001,'moderate'), (0.70,'moderate'), (0.700001,'high')]:
            with self.subTest(identity=identity): self.assertEqual(novelty.sequence_class(identity), expected)

    def test_published_general_scale(self):
        for sequence, expected in [('low',[4,3,2]), ('moderate',[3,2,2]), ('high',[2,1,1])]:
            with self.subTest(sequence=sequence):
                self.assertEqual([novelty.general_level(sequence,s) for s in ('low','moderate','high')], expected)

    def test_tm_and_coverage_boundaries(self):
        for end,tm,expected in [(70,.5,'low'), (71,.5,'moderate'), (71,.499999,'low'), (70,.8,'high'), (70,.799999,'low'), (69,1.,'low')]:
            with self.subTest(end=end,tm=tm): self.assertEqual(novelty.structural_class(100,[domain(end,tm)])['category'], expected)

    def test_no_double_counting_redundant_domain_hits(self):
        with self.assertRaisesRegex(ValueError,'Overlapping'):
            novelty.structural_class(100,[domain(40,.9),domain(40,.9)])

    def test_tm_normalization_must_be_query_domain(self):
        with self.assertRaisesRegex(ValueError,'normalized'):
            novelty.structural_class(100,[dict(start=1,end=80,tm_score=.9,normalization='reference_length')])

    def test_missing_database_can_only_increase_similarity(self):
        self.assertEqual(novelty.possible_levels('low',['patents']),[1,2,3,4])
        self.assertEqual(novelty.possible_levels('moderate',['patents']),[1,2,3])
        self.assertEqual(novelty.possible_levels('high',['patents']),[1,2])
        self.assertEqual(novelty.possible_levels('low',[]),[2,3,4])
        self.assertEqual(novelty.possible_levels('low',[],structure='low'),[4])

    def test_reference_hit_identity_is_a_fraction_not_percent(self):
        import tempfile
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/'hits.tsv'
            path.write_text('q000001\tref\t80\t40\t50\t50\t50\t1\t50\t1\t50\t0\t100\t1\t1\n')
            with self.assertRaisesRegex(ValueError,'fraction'):
                novelty.parse_hits(path,'swissprot',{'q000001':dict(sequence='A'*50)})

    def test_identity_uses_integer_counts_before_boundary_comparison(self):
        import tempfile
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/'hits.tsv'
            path.write_text('q000001\tref\t0.7\t701\t1000\t1000\t1000\t1\t1000\t1\t1000\t0\t100\t1\t1\n')
            hit=novelty.parse_hits(path,'swissprot',{'q000001':dict(sequence='A'*1000)})['q000001'][0]
            self.assertEqual(hit['fident'],.701)
            self.assertEqual(novelty.sequence_class(hit['fident']),'high')

    def test_highest_identity_is_not_highest_score(self):
        hits=[dict(fident=.4,bits=200,database='pdb',target='long'),dict(fident=.8,bits=100,database='pdb',target='closer')]
        self.assertEqual(novelty.identity_ranked_hits(hits)[0]['target'],'closer')

    def test_input_does_not_silently_strip_colons_or_fix_case(self):
        import tempfile
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/'designs.csv'
            with path.open('w',newline='') as f:
                w=csv.writer(f);w.writerow(['name','sequence']);w.writerow(['bad','ACD:efg'])
            with self.assertRaisesRegex(ValueError,'silent repairs'): novelty.load_designs(path)


if __name__ == '__main__': unittest.main()
