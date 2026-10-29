"""Regression tests for the pH direction of the linkage scorers.

C4 is acid-OFF: binding at pH 7.4 and none at pH 6.0. With
L = log10 Ka(pH_low) - log10 Ka(pH_high), that means L <= -2 (100-fold weaker at pH 6.0).

Commit 75c834b (2026-10-23) flipped the pass flag, banner and sort order in ph_linkage.py to
L >= +2 and negated the value in complex_propka_switch.py, which is the opposite direction.
These tests fail on that state and pass on the corrected one. They use synthetic pKa sets, so
they do not depend on PROPKA being installed.
"""
from pathlib import Path
import math
import sys
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))

BEST = ROOT / 'results/v3_all_pdb/c4_tfr1_groove_human_forced_l58_9f94b53776050aad_candidate3.pdb'


def _score(bound, free_binder, free_target):
    from ph_linkage import score_one
    from tempfile import TemporaryDirectory
    with TemporaryDirectory() as td:
        pdb = Path(td) / 'input.pdb'
        pdb.write_text(BEST.read_text())
        with patch('ph_linkage.run_propka', side_effect=[bound, free_binder, free_target]):
            return score_one(str(pdb), 'C', 'AB', 6.0, 7.4, 0.1, False)


class PhLinkageDirection(unittest.TestCase):
    TARGET = {('A', 1, 'LYS'): 10.0}

    def test_depressed_pka_on_binding_is_acid_off_and_passes(self):
        # three His, pKa 7.0 free -> 5.0 bound: weaker binding at pH 6.0, the C4 direction
        keys = [('C', 24, 'HIS'), ('C', 30, 'HIS'), ('C', 40, 'HIS')]
        bound = {k: 5.0 for k in keys}
        free = {k: 7.0 for k in keys}
        r = _score({**bound, **self.TARGET}, free, self.TARGET)
        self.assertLess(r['L'], -2.0)
        self.assertTrue(r['switch_ok'], 'an acid-OFF design must pass the C4 gate')

    def test_elevated_pka_on_binding_is_acid_on_and_fails(self):
        # three His, pKa 6.0 free -> 8.0 bound: tighter binding at pH 6.0, the wrong direction
        keys = [('C', 24, 'HIS'), ('C', 30, 'HIS'), ('C', 40, 'HIS')]
        bound = {k: 8.0 for k in keys}
        free = {k: 6.0 for k in keys}
        r = _score({**bound, **self.TARGET}, free, self.TARGET)
        self.assertGreater(r['L'], 2.0)
        self.assertFalse(r['switch_ok'], 'an acid-ON design must NOT pass the C4 gate')


class ComplexProPkaSwitchDirection(unittest.TestCase):
    def test_wyman_ddg_positive_means_pka_depressed_on_binding(self):
        from complex_propka_switch import wyman_linkage_ddg as f
        self.assertGreater(f(5.0) - f(7.0), 0.0, 'depressed pKa must give a positive switch value')
        self.assertLess(f(8.0) - f(6.5), 0.0, 'elevated pKa must give a negative switch value')

    def test_source_does_not_negate_the_interface_total(self):
        src = (ROOT / 'scripts/complex_propka_switch.py').read_text()
        self.assertNotIn('round(-total_ddg', src,
                         'negating total_ddg ranks the acid-ON direction first')


if __name__ == '__main__':
    unittest.main()
