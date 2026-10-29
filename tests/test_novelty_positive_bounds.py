"""Partial structural evidence must never invent an eligible novelty score."""
import copy
import itertools
from pathlib import Path
import sys
import unittest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
import novelty_level_predictor as n


def partial(sequence_class,domains,scores):
    matches=[dict(**d,tm_score=s,normalization='query_domain_length') for d,s in zip(domains,scores)]
    measured=n.structural_class(100,matches)
    measured['domain_segmentation_coverage']=sum(len(n.domain_positions(d,100)) for d in domains)/100
    return dict(sequence_evidence_class=sequence_class,sequence_database_gaps=['patents'],
        structure_evidence=dict(domains=domains),structural_similarity=measured,
        local_predicted_level=4,predicted_level=4)


class PositiveBounds(unittest.TestCase):
    def test_points_survive_every_critical_value_completion(self):
        # Similarity categories change only at .5/.8. Enumerate those critical
        # completions, including exactly70% coverage and discontinuous domains.
        partitions=[[dict(start=1,end=100)],[dict(start=1,end=55)],
            [dict(start=1,end=70)],
            [dict(start=1,end=80),dict(start=81,end=100)],
            [dict(start=1,end=35),dict(start=36,end=70)],
            [dict(segments=[[1,30],[71,100]]),dict(start=31,end=70)]]
        for domains in partitions:
            for observed in itertools.product((0.,.49,.5,.79,.8,.95),repeat=len(domains)):
                for sequence in ('low','moderate','high'):
                    r=partial(sequence,domains,observed);n.apply_positive_evidence_bounds(r)
                    if r['local_predicted_level'] is None:continue
                    self.assertNotIn(r['local_predicted_level'],(3,4))
                    choices=[sorted({score,*[t for t in (.5,.8,1.) if t>=score]}) for score in observed]
                    for completion in itertools.product(*choices):
                        actual=partial(sequence,domains,completion)['structural_similarity']['category']
                        self.assertEqual(r['local_predicted_level'],n.general_level(sequence,actual))

    def test_absence_and_coverage_only_do_not_emit_novelty_four(self):
        for domains in ([dict(start=1,end=100)],[dict(start=1,end=55)]):
            r=partial('low',domains,[0.]);n.apply_positive_evidence_bounds(r)
            self.assertIsNone(r['predicted_level'])
            self.assertIsNone(r['local_predicted_level'])

    def test_moderate_positive_for_low_sequence_is_still_unresolved(self):
        r=partial('low',[dict(start=1,end=100)],[.6]);n.apply_positive_evidence_bounds(r)
        self.assertEqual(r['local_level_interval_with_unfinished_alignments'],[2,3])
        self.assertIsNone(r['predicted_level'])
        self.assertEqual(r['structural_similarity']['category'],'unestablished')


if __name__=='__main__':unittest.main()
