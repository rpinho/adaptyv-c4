"""Scientific failures and sequence siblings must not become scaffold successes."""
import csv
import hashlib
import json
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from campaign_admission import read_outcome,trajectory_id,unique_parents

HASH='0123456789abcdef'
SECOND='fedcba9876543210'


class AdmissionTests(unittest.TestCase):
    def fixture(self,root,rows,status):
        campaign=root/'campaign';ranked=campaign/'3_Ranked';ranked.mkdir(parents=True)
        with (ranked/'!_Ranked.csv').open('w') as fh:
            writer=csv.DictWriter(fh,fieldnames=['design','hash']);writer.writeheader();writer.writerows(rows)
        followup=root/'followup';followup.mkdir()
        (followup/'status.json').write_text(json.dumps(status))
        return campaign,followup

    def valid(self,name):
        return dict(design=name,status='no_eligible_pair',audit=dict(target_sequence_verified=True,
                    canonical_sha256=hashlib.sha256(b'canonical test structure').hexdigest(),
                    minimum_binder_target_backbone_distance=2.2,target_subunit_backbone_min_distance=2.0,
                    observed_NAG_min_distance=2.5))

    def evaluate(self,rows,records,**status_changes):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);status=dict(status='completed',errors=[],designs=records);status.update(status_changes)
            campaign,followup=self.fixture(root,rows,status)
            for row in rows:(followup/(row['design']+'_canonical.pdb')).write_bytes(b'canonical test structure')
            return read_outcome(campaign,followup)

    def test_sequence_variants_are_one_parent_before_cpu_selection(self):
        rows=[dict(design=f'run_{HASH}_seq0'),dict(design=f'run_{HASH}_seq1'),dict(design=f'run_{SECOND}_seq0')]
        self.assertEqual([trajectory_id(r) for r in unique_parents(rows)],[HASH,SECOND])

    def test_missing_or_conflicting_identity_is_not_guessed(self):
        for row in [dict(design='anonymous'),dict(design=f'run_{HASH}',hash=SECOND),dict(design=f'../run_{HASH}')]:
            with self.assertRaises(ValueError):trajectory_id(row)

    def test_multiple_successful_sequence_rows_count_one_trajectory(self):
        rows=[dict(design=f'run_{HASH}_seq{i}',hash=HASH) for i in range(2)]
        out=self.evaluate(rows,[self.valid(r['design']) for r in rows])
        self.assertEqual(out['generator_accepted_rows'],2)
        self.assertEqual(out['native_context_admitted_trajectory_ids'],[HASH])
        self.assertEqual(out['engineering_advance_trajectory_ids'],[])
        self.assertTrue(out['followup_completed_without_errors'])

    def test_caught_scientific_errors_block_despite_completed_process(self):
        name=f'run_{HASH}_seq0';row=dict(design=name,hash=HASH)
        out=self.evaluate([row],[dict(design=name,status='failed',error='native clash')],errors=[{'error':'native clash'}])
        self.assertFalse(out['followup_completed_without_errors'])
        self.assertEqual(out['native_context_admitted_trajectory_ids'],[])

    def test_partial_followup_cannot_advance(self):
        name=f'run_{HASH}_seq0';row=dict(design=name,hash=HASH)
        out=self.evaluate([row],[self.valid(name)],status='time_cap_partial')
        self.assertFalse(out['followup_completed_without_errors'])

    def test_native_sequence_clearance_and_file_hash_are_required(self):
        name=f'run_{HASH}_seq0';rows=[dict(design=name,hash=HASH)]
        mutations=[('target_sequence_verified',False),('minimum_binder_target_backbone_distance',1.9),
                   ('target_subunit_backbone_min_distance',1.7),('observed_NAG_min_distance',1.9),
                   ('observed_NAG_min_distance',float('nan')),('canonical_sha256','incorrect')]
        for key,value in mutations:
            record=self.valid(name);record['audit'][key]=value
            with self.subTest(key=key,value=value):
                out=self.evaluate(rows,[record])
                self.assertFalse(out['followup_completed_without_errors'])
                self.assertEqual(out['native_context_admitted_trajectory_ids'],[])

    def test_missing_followup_does_not_inherit_raw_generator_success(self):
        with tempfile.TemporaryDirectory() as tmp:
            out=read_outcome(Path(tmp)/'missing_campaign',Path(tmp)/'missing_followup')
            self.assertFalse(out['followup_completed_without_errors'])

    def test_empty_batch_requires_explicit_cpu_no_acceptance_status(self):
        with tempfile.TemporaryDirectory() as tmp:
            campaign,followup=self.fixture(Path(tmp),[],dict(status='no_accepted_scaffolds',errors=[]))
            self.assertTrue(read_outcome(campaign,followup)['followup_completed_without_errors'])
            (followup/'status.json').write_text(json.dumps(dict(status='preflight',errors=[])))
            self.assertFalse(read_outcome(campaign,followup)['followup_completed_without_errors'])


if __name__=='__main__':unittest.main()
