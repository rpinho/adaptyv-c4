import importlib.util
import json
from pathlib import Path
import sys
import tempfile
import types
import unittest
from unittest.mock import Mock,patch

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))


def atomic(path,data):path.write_text(json.dumps(data))


class AdaptiveSafetyTests(unittest.TestCase):
    def load(self):
        driver=types.SimpleNamespace(HOURLY_MAX=4.649328,atomic=atomic,main=Mock(),modal=Mock())
        spec=importlib.util.spec_from_file_location('adaptive_test',Path(__file__).resolve().parents[1]/'scripts/drive_modal_adaptive.py')
        module=importlib.util.module_from_spec(spec)
        with patch.dict(sys.modules,{'drive_modal_campaign':driver}):spec.loader.exec_module(module)
        return module,driver

    def test_budget_reserves_complete_call_before_any_gpu_allocation(self):
        module,driver=self.load()
        with tempfile.TemporaryDirectory() as tmp:
            state=Path(tmp)
            (state/'search.json').write_text(json.dumps(dict(status='active',profile='c4_local_frame_v5',completed=[],accepted=0,runtime_estimate_upper_usd=25)))
            with patch.object(module,'STATE',state):module.main()
            driver.main.assert_not_called()
            self.assertEqual(json.loads((state/'search.json').read_text())['status'],'runtime_budget_reserve_reached')

    def test_finished_search_never_restarts_itself(self):
        module,driver=self.load()
        with tempfile.TemporaryDirectory() as tmp:
            state=Path(tmp)
            (state/'search.json').write_text(json.dumps(dict(status='declared_fold_hypotheses_exhausted')))
            with patch.object(module,'STATE',state):module.main()
            driver.main.assert_not_called()

    def test_driver_exception_cancels_exact_call_and_verifies_actual_release(self):
        module,driver=self.load()
        with tempfile.TemporaryDirectory() as tmp:
            state=Path(tmp)
            batch=state/'batch_1';batch.mkdir()
            (batch/'supervisor.json').write_text(json.dumps(dict(call_id='fc-exact-test')))
            driver.main.side_effect=RuntimeError('mirror failure')
            driver.modal.Function.from_name.return_value.get_current_stats.return_value=types.SimpleNamespace(
                num_total_runners=0,num_running_inputs=0,backlog=0)
            with patch.object(module,'STATE',state):module.main()
            driver.modal.FunctionCall.from_id.assert_called_once_with('fc-exact-test')
            driver.modal.FunctionCall.from_id.return_value.cancel.assert_called_once_with(terminate_containers=True)
            self.assertTrue(json.loads((state/'search.json').read_text())['emergency_release']['gpu_released_verified'])

    def test_process_exit_zero_does_not_override_scientific_errors(self):
        module,driver=self.load()
        driver.main.return_value=dict(gpu_released_verified=True,status='returned',
            result={'status':'completed'},cpu_followup_returncode=0,submitted=1,released=2)
        evidence=dict(followup_completed_without_errors=False,error='native-context failure',
                      generator_accepted_rows=3,native_context_admitted_trajectory_ids=[])
        with tempfile.TemporaryDirectory() as tmp:
            state=Path(tmp)
            with patch.object(module,'STATE',state),patch.object(module,'outcome',return_value=(evidence,['relax'])):
                module.main()
            result=json.loads((state/'search.json').read_text())
            self.assertEqual(result['status'],'stopped_on_scientific_followup_failure')
            self.assertEqual(result['accepted'],0)
            self.assertEqual(driver.main.call_count,1)

    def test_global_identity_deduplicates_previous_and_current_batch(self):
        module,driver=self.load();identity='0123456789abcdef'
        driver.main.return_value=dict(gpu_released_verified=True,status='returned',
            result={'status':'completed'},cpu_followup_returncode=0,submitted=1,released=2)
        evidence=dict(followup_completed_without_errors=True,generator_accepted_rows=2,
                      native_context_admitted_trajectory_ids=[identity])
        with tempfile.TemporaryDirectory() as tmp:
            state=Path(tmp)
            (state/'search.json').write_text(json.dumps(dict(status='active',profile='c4_local_frame_v5',
                completed=[],accepted=99,native_context_trajectory_ids=[identity],runtime_estimate_upper_usd=1)))
            with patch.object(module,'STATE',state),patch.object(module,'outcome',return_value=(evidence,['relax'])):
                module.main()
            result=json.loads((state/'search.json').read_text())
            self.assertEqual(result['accepted'],1)
            self.assertEqual(result['status'],'native_context_outputs_require_acquisition_review')
            self.assertEqual(driver.main.call_count,1)

    def test_legacy_raw_row_counts_cannot_resume_into_new_allocations(self):
        module,driver=self.load()
        with tempfile.TemporaryDirectory() as tmp:
            state=Path(tmp)
            (state/'search.json').write_text(json.dumps(dict(status='active',completed=[{'accepted':3}])))
            with patch.object(module,'STATE',state):module.main()
            driver.main.assert_not_called()
            self.assertEqual(json.loads((state/'search.json').read_text())['status'],'legacy_acceptance_accounting_requires_review')


if __name__=='__main__':unittest.main()
