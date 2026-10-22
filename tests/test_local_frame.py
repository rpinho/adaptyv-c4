import json
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
import prepare_local_frame as frame
import modal_search_policy as policy

ROOT=Path(__file__).resolve().parents[1]


class FrameAndPolicyTests(unittest.TestCase):
    def test_full_receptor_geometry_and_atom_identity_survive_placement(self):
        with tempfile.TemporaryDirectory() as tmp:
            report=frame.prepare(ROOT/'results/gpu_campaign_v2/human_11qc_full_monomer.pdb',
                                 ROOT/'results/gpu_campaign/human_11qc_full_dimer.pdb',Path(tmp)/'local.pdb')
            self.assertEqual(report['residues'],639)
            self.assertTrue(report['non_coordinate_columns_identical'])
            self.assertLess(report['max_CA_pair_distance_change_A'],0.002)
            self.assertGreater(report['new_origin_full_context_atom_clearance_A'],6)
            self.assertAlmostEqual(report['old_origin_to_hotspot_A'],188.816,places=2)

    def test_initial_point_inside_full_receptor_is_rejected(self):
        with tempfile.TemporaryDirectory() as tmp,self.assertRaisesRegex(ValueError,'inside the restored context'):
            frame.prepare(ROOT/'results/gpu_campaign_v2/human_11qc_full_monomer.pdb',
                          ROOT/'results/gpu_campaign/human_11qc_full_dimer.pdb',Path(tmp)/'bad.pdb',offset=0)

    def test_declared_changes_preserve_all_gates_and_target_hotspots(self):
        base=json.loads((ROOT/'results/gpu_campaign_v3/campaign.json').read_text())
        for profile in policy.PROFILES:
            with self.subTest(profile=profile):
                changed=policy.make_config(base,profile,{'campaign_seed':123},'/results/test')
                self.assertEqual(changed['filters'],base['filters'])
                self.assertEqual(changed['targets'][0]['hotspots'],base['targets'][0]['hotspots'])
                self.assertEqual(changed['targets'][0]['chains'],'A')
                self.assertFalse(changed['forced_targeting'])
                self.assertEqual(changed['max_trajectories'],2)
                self.assertEqual(base['targets'][0]['target_path'],'/workspace/runs/human_11qc_full_monomer.pdb')

    def test_screen_and_refine_labels_do_not_prove_a_fold_confidence_failure(self):
        for profile in policy.PROFILES:
            selected,reason=policy.next_profile(profile,0,['screen','refine'])
            self.assertIsNone(selected)
            self.assertEqual(reason,'physical_failure_review_required')

    def test_repaired_test_changes_only_initialization_and_trajectory_limit(self):
        base=json.loads((ROOT/'results/gpu_campaign_v3/campaign.json').read_text())
        changed=policy.make_config(base,policy.REPAIRED,{'campaign_seed':20261021},'/results/test')
        expected=base|dict(large_complex_initialization=True,max_trajectories=1,
                           campaign_name=policy.REPAIRED,project_folder='/results/test')
        expected['targets']=[dict(base['targets'][0],target_path='/opt/c4/local-frame-target.pdb')]
        self.assertEqual(changed,expected)

    def test_interface8_changes_only_the_contact_objective_from_repaired_control(self):
        base=json.loads((ROOT/'results/gpu_campaign_v3/campaign.json').read_text())
        reference=policy.make_config(base,policy.REPAIRED,{'campaign_seed':20261021},'/results/test')
        changed=policy.make_config(base,policy.INTERFACE8,{'campaign_seed':20261021},'/results/test')
        reference['campaign_name']=policy.INTERFACE8
        reference['losses']['interface_contacts']['params']['cutoff']=8.0
        self.assertEqual(changed,reference)
        self.assertEqual(base['losses']['interface_contacts']['params']['cutoff'],20.0)

    def test_forced_search_changes_only_targeting_mode_from_native_control(self):
        base=json.loads((ROOT/'results/gpu_campaign_v3/campaign.json').read_text())
        reference=policy.make_config(base,policy.REPAIRED,{'campaign_seed':20261021},'/results/test')
        changed=policy.make_config(base,policy.FORCED,{'campaign_seed':20261021},'/results/test')
        reference.update(campaign_name=policy.FORCED,forced_targeting=True)
        self.assertEqual(changed,reference)
        self.assertEqual(changed['losses']['interface_contacts']['params']['cutoff'],20.0)
        self.assertFalse(base['forced_targeting'])

    def test_structural_admission_alone_never_triggers_more_gpu_sampling(self):
        selected,reason=policy.next_profile(policy.PROFILES[0],1,['screen'])
        self.assertIsNone(selected)
        self.assertEqual(reason,'native_context_outputs_require_acquisition_review')

    def test_confidence_profile_changes_only_effective_weight_from_forced_control(self):
        base=json.loads((ROOT/'results/gpu_campaign_v3/campaign.json').read_text())
        reference=policy.make_config(base,policy.FORCED,{'campaign_seed':20261021},'/results/test')
        changed=policy.make_config(base,policy.CONFIDENCE,{'campaign_seed':20261021},'/results/test')
        reference.update(campaign_name=policy.CONFIDENCE,weights_plddt_loss=.5)
        self.assertEqual(changed,reference)
        self.assertEqual(base['weights_plddt_loss'],.1)

    def test_later_stage_failure_requires_different_diagnosis(self):
        selected,reason=policy.next_profile(policy.PROFILES[0],0,['anneal'])
        self.assertIsNone(selected)
        self.assertEqual(reason,'needs_interface_or_later_stage_diagnosis')

    def test_profile_names_and_gate_override_are_not_arbitrary(self):
        with self.assertRaises(ValueError):policy.profile_for_name('unbounded_run')
        base=json.loads((ROOT/'results/gpu_campaign_v3/campaign.json').read_text())
        with self.assertRaisesRegex(ValueError,'campaign seed'):
            policy.make_config(base,policy.PROFILES[0],{'filters':{}},'/results/test')


if __name__=='__main__':unittest.main()
