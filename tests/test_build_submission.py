"""Format-contract fixtures for the submission emitter; none is a binding candidate.

Every sequence here is synthetic filler chosen to exercise one validation rule. A
validator that has only ever seen well-formed input proves nothing, so the negative
cases are the point of this file: each asserts that build_submission.py REFUSES to
write, and names the reason.
"""
from pathlib import Path
import csv
import json
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
BUILDER = ROOT / 'scripts' / 'build_submission.py'

OK = 'MKEELLRRLEELARRHEELAKGDLSEEEYRRIVEELGRQLREMAEWLAERLG'  # 52 aa, canonical 20


def shortlist(rows, columns=('name', 'sequence', 'molecule_class')):
    """Write a ranked shortlist CSV and return its path."""
    fh = tempfile.NamedTemporaryFile('w', suffix='.csv', delete=False, newline='')
    w = csv.DictWriter(fh, fieldnames=list(columns), lineterminator='\n')
    w.writeheader()
    for r in rows:
        w.writerow({k: r.get(k, '') for k in columns})
    fh.close()
    return fh.name


def build(path, *extra):
    """Run the emitter; return (returncode, combined output, out_path)."""
    out = Path(tempfile.mkdtemp()) / 'submission.csv'
    p = subprocess.run([sys.executable, str(BUILDER), path, '-o', str(out)] + list(extra),
                       capture_output=True, text=True)
    return p.returncode, p.stdout + p.stderr, out


def row(name, seq, cls='protein'):
    return {'name': name, 'sequence': seq, 'molecule_class': cls}


class Conformant(unittest.TestCase):
    def test_valid_file_is_written_in_rank_order(self):
        # distinct sequences: a truncation alone would collide once it hits len(OK)
        rows = [row('c4_tfr1_%02d' % i, OK[:-1] + aa) for i, aa in enumerate('GAV', 1)]
        code, out, path = build(shortlist(rows))
        self.assertEqual(code, 0, out)
        got = list(csv.DictReader(open(path)))
        self.assertEqual([r['name'] for r in got],
                         ['c4_tfr1_01', 'c4_tfr1_02', 'c4_tfr1_03'])
        self.assertEqual(list(got[0].keys()), ['name', 'sequence', 'molecule_class'])

    def test_manifest_records_the_unverified_tag_state(self):
        code, out, path = build(shortlist([row('d1', OK)]))
        self.assertEqual(code, 0, out)
        man = json.load(open(str(path).replace('.csv', '.manifest.json')))
        self.assertEqual(man['twin_strep_tag'], 'excluded')
        self.assertIn('UNVERIFIED', man['twin_strep_status'])
        self.assertTrue(man['row_order_is_rank'])
        self.assertFalse(man['is_format_test'])

    def test_format_test_flag_is_recorded_and_announced(self):
        code, out, path = build(shortlist([row('d1', OK)]), '--format-test')
        self.assertEqual(code, 0, out)
        self.assertIn('FORMAT TEST', out)
        self.assertTrue(json.load(open(str(path).replace('.csv', '.manifest.json')))
                        ['is_format_test'])

    def test_extra_columns_pass_through(self):
        rows = [dict(row('d1', OK), i_pTM='0.81')]
        code, out, path = build(shortlist(rows, ('name', 'sequence', 'molecule_class', 'i_pTM')),
                                '--extra-columns', 'i_pTM')
        self.assertEqual(code, 0, out)
        got = list(csv.DictReader(open(path)))
        self.assertEqual(got[0]['i_pTM'], '0.81')

    def test_fab_may_use_the_chain_separator(self):
        code, out, _ = build(shortlist([row('fab1', OK + ':' + OK, 'fab_kappa')]))
        self.assertEqual(code, 0, out)


class Rejected(unittest.TestCase):
    def assert_refused(self, rows, needle, *extra, **kw):
        code, out, path = build(shortlist(rows, **kw), *extra)
        self.assertNotEqual(code, 0, 'expected refusal, got success:\n' + out)
        self.assertIn(needle, out)
        self.assertFalse(Path(path).exists(), 'a rejected run must write no file')

    def test_over_the_track_row_cap(self):
        self.assert_refused([row('d%02d' % i, OK) for i in range(21)], 'exceeds')

    def test_sequence_too_short(self):
        self.assert_refused([row('d1', 'MKEELLRRL')], 'outside the documented')

    def test_sequence_too_long(self):
        self.assert_refused([row('d1', 'A' * 251)], 'outside the documented')

    def test_noncanonical_residue(self):
        self.assert_refused([row('d1', OK[:-1] + 'X')], 'non-standard residues X')

    def test_lowercase_is_not_silently_upcased(self):
        self.assert_refused([row('d1', OK.lower())], 'lowercase not accepted')

    def test_colon_in_a_protein_row(self):
        self.assert_refused([row('d1', OK + ':' + OK)], 'malformed row')

    def test_duplicate_name(self):
        self.assert_refused([row('same', OK), row('same', OK[:-1])], 'duplicate name')

    def test_duplicate_sequence(self):
        self.assert_refused([row('d1', OK), row('d2', OK)], 'duplicate sequence')

    def test_empty_sequence(self):
        self.assert_refused([row('d1', '')], 'empty sequence')

    def test_bad_molecule_class_is_fatal_under_strict(self):
        self.assert_refused([row('d1', OK, 'minibinder')], 'not in', '--strict-class')

    def test_name_violating_our_convention(self):
        self.assert_refused([row('bad name!', OK)], 'convention')

    def test_missing_required_column(self):
        code, out, path = build(shortlist([{'name': 'd1', 'sequence': OK}],
                                          ('name', 'sequence')))
        self.assertNotEqual(code, 0, out)
        self.assertIn('missing required column', out)
        self.assertFalse(Path(path).exists())

    def test_unknown_molecule_class_only_warns_by_default(self):
        """The enum spelling is paraphrased from a capped fetch, so it must not be fatal
        unless asked -- encoding a guess as a hard rule is the failure mode here."""
        code, out, path = build(shortlist([row('d1', OK, 'minibinder')]))
        self.assertEqual(code, 0, out)
        self.assertIn('warning', out)
        self.assertTrue(Path(path).exists())


if __name__ == '__main__':
    unittest.main()
