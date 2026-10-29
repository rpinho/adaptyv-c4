"""Regression checks for chain mapping, native histidines and ranking provenance."""
from pathlib import Path
import sys
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
from install_switch import make_plan, map_sites
from run_switch_campaign import substrates
from switch_sites import score_design

BEST = ROOT / 'results/v3_all_pdb/c4_tfr1_groove_human_forced_l58_9f94b53776050aad_candidate3.pdb'


class RankingTests(unittest.TestCase):
    def test_real_c4_keeps_native_histidine_and_scorer_order(self):
        plan = make_plan(BEST, 'C', 'AB', max_pairs=4)
        ranked = score_design(str(BEST), 'C', 'AB')
        eligible = [r for r in ranked if not r['note'] and r['est_contribution'] < 0
                    and r['d_burial'] > 0]
        self.assertEqual(plan['shortlist'], eligible[:8])
        self.assertEqual(plan['pairs'][0]['sites'][0]['resnum'], 24)
        self.assertEqual(plan['pairs'][0]['sites'][0]['wt'], 'H')
        sums = [p['estimated_site_sum'] for p in plan['pairs']]
        self.assertEqual(sums, sorted(sums))
        self.assertTrue(all(len(p['sites']) == 2 for p in plan['pairs']))
        self.assertTrue(all(not r['note'] for p in plan['pairs'] for r in p['sites']))

    def test_missing_fragment_assignment_fails(self):
        with self.assertRaisesRegex(ValueError, 'assign all'):
            make_plan(BEST, 'C', 'A')
        with self.assertRaises(ValueError):
            make_plan(BEST, 'B', 'AB')

    def test_no_pair_is_not_a_successful_mutation(self):
        with patch('install_switch.score_design', return_value=[]):
            self.assertEqual(make_plan(BEST, 'C', 'AB')['pairs'], [])

    def test_maps_pdb_number_instead_of_pose_index(self):
        class Info:
            def pdb2pose(self, chain, number):
                return {('C', 24): 97}.get((chain, number), 0)
            def chain(self, i):
                return 'C'
            def number(self, i):
                return 24
        class Residue:
            def name1(self):
                return 'H'
        class Pose:
            def pdb_info(self):
                return Info()
            def residue(self, i):
                return Residue()
        self.assertEqual(map_sites(Pose(), [dict(chain='C', resnum=24, wt='H')]), [97])
        with self.assertRaisesRegex(ValueError, 'identity mismatch'):
            map_sites(Pose(), [dict(chain='C', resnum=24, wt='Q')])
        with self.assertRaisesRegex(ValueError, 'not present'):
            map_sites(Pose(), [dict(chain='C', resnum=25, wt='H')])

    def test_campaign_covers_25_ranked_backbones_best_L_first(self):
        selected = substrates(ROOT / 'results/pods/c4_E1_human/3_Ranked/!_Ranked.csv',
                              ROOT / 'results/v3_all_linkage.csv', ROOT / 'results/v3_all_pdb')
        self.assertEqual(len(selected), 25)
        self.assertEqual(len({s['backbone'] for s in selected}), 25)
        values = [s['historical_L'] for s in selected]
        self.assertEqual(values, sorted(values))
        self.assertEqual(values[0], -1.175)


class LinkageFailureTests(unittest.TestCase):
    def test_suppressed_histidine_never_becomes_a_complete_gate_score(self):
        from run_switch_campaign import score_model
        result = dict(binder_chain='C', target_chain='AB',
                      pair_references=[dict(pair=1, pdb='ref.pdb')],
                      variants=[dict(pair=1, pdb='HH.pdb', reference='ref.pdb',
                                     kept=True, kind='switch',
                                     sites=[dict(chain='C', resnum=1), dict(chain='C', resnum=7)])])
        partial = dict(L=-2.3, switch_ok=True,
                       rows=[dict(chain='C', resnum=7, resname='HIS')])
        with patch('run_switch_campaign.score_one', side_effect=[dict(L=-.5), partial]):
            checked = score_model(result, Path('/tmp'))['variants'][0]
        self.assertTrue(checked['kept'])  # energy retention is a separate fact
        self.assertNotIn('linkage', checked)
        self.assertNotIn('delta_L_vs_reference', checked)
        self.assertEqual(checked['missing_histidines'], [dict(chain='C', resnum=1)])
        self.assertEqual(checked['linkage_status'], 'incomplete_selected_histidines')
        self.assertEqual(checked['partial_linkage']['L'], -2.3)

    def test_free_state_failure_cannot_produce_partial_score(self):
        from ph_linkage import score_one
        from tempfile import TemporaryDirectory
        with TemporaryDirectory() as td:
            pdb = Path(td) / 'input.pdb'
            pdb.write_text(BEST.read_text())
            cx = {('C', 24, 'HIS'): 5.67}
            for missing_binder in (True, False):
                free = [None, cx] if missing_binder else [cx, None]
                with patch('ph_linkage.run_propka', side_effect=[cx] + free):
                    self.assertIsNone(score_one(str(pdb), 'C', 'AB', 6.0, 7.4, 0.1, False))


if __name__ == '__main__':
    unittest.main()
