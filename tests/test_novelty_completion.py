"""Completed and missing native outputs must support different conclusions."""
import copy
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'results/novelty_predictor/adaptive_completion'))
from finalize import refine


def example(end=100):
    return dict(length=100, sequence_evidence_class='low', sequence_database_gaps=['patents'],
                structure_evidence=dict(domains=[dict(start=1, end=end, query_id='q1')]))


class Completion(unittest.TestCase):
    def test_same_observed_match_is_unknown_until_all_pairs_return(self):
        for returned, expected_level in ((9, None), (10, 3)):
            row = example()
            refine(row, {'q1': [dict(qtmscore=.6)]}, {'q1': returned},
                   {'q1': dict(unique_union_targets=10)})
            self.assertEqual(row['predicted_level'], expected_level)
            self.assertEqual(row['local_level_interval_with_unfinished_alignments'], [2, 3] if returned == 9 else [3])

    def test_missing_evidence_and_low_coverage_do_not_manufacture_four(self):
        row = example(55)
        refine(row, {'q1': []}, {'q1': 0}, {'q1': dict(unique_union_targets=10)})
        self.assertIsNone(row['predicted_level'])
        self.assertIsNone(row['structure_evidence']['domains'][0]['best_hit'])
        self.assertFalse(row['adaptive_completion']['all_requested_pairs_exported'])

    def test_completed_coverage_limited_four_is_explicit_despite_known_domain(self):
        row = example(55)
        refine(row, {'q1': [dict(qtmscore=.95)]}, {'q1': 10}, {'q1': dict(unique_union_targets=10)})
        self.assertEqual(row['predicted_level'], 4)
        self.assertTrue(row['structural_similarity']['coverage_limited_low_category'])

    def test_exact_seventy_percent_high_coverage_keeps_published_asymmetry(self):
        for score, level in ((.5, 4), (.8, 2)):
            row = example(70)
            refine(row, {'q1': [dict(qtmscore=score)]}, {'q1': 10}, {'q1': dict(unique_union_targets=10)})
            self.assertEqual(row['predicted_level'], level)
            self.assertTrue(row['structural_similarity']['tm_boundary_precision_caveat'])

    def test_complete_nonempty_set_requires_real_evidence_and_no_extra_rows(self):
        for hits, count in (([], 10), ([dict(qtmscore=.6)], 11)):
            with self.assertRaises(ValueError):
                refine(example(), {'q1': hits}, {'q1': count}, {'q1': dict(unique_union_targets=10)})


if __name__ == '__main__':
    unittest.main()
