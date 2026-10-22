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


if __name__=='__main__':unittest.main()
