"""A seed study must stop on evidence/release failures and never resubmit finished calls."""
import importlib.util
import json
from pathlib import Path
import sys
import tempfile
import types
import unittest
from unittest.mock import Mock,patch


class NativeSeedQueueTests(unittest.TestCase):
    def exercise(self,good=True,admitted=False,error=False,prior=None,idle=True):
        def atomic(path,value):path.write_text(json.dumps(value))
        result=dict(status='returned',gpu_released_verified=good,result=dict(status='completed'),
                    cpu_followup_returncode=0,submitted=1,released=451)
        fake=types.SimpleNamespace(HOURLY_MAX=4.649328,atomic=atomic,main=Mock(return_value=result))
        stats=types.SimpleNamespace(num_total_runners=0 if idle else 1,num_running_inputs=0,backlog=0)
        fake.modal=types.SimpleNamespace(Function=types.SimpleNamespace(from_name=Mock(return_value=types.SimpleNamespace(get_current_stats=Mock(return_value=stats)))))
        spec=importlib.util.spec_from_file_location('native_seed_test',Path(__file__).resolve().parents[1]/'scripts/run_native_replicates.py')
        module=importlib.util.module_from_spec(spec)
        with patch.dict(sys.modules,{'drive_modal_campaign':fake,'drive_modal_adaptive':types.SimpleNamespace(abort_after_driver_error=Mock())}):
            spec.loader.exec_module(module)
        outcome=dict(followup_completed_without_errors=not error,native_context_admitted_trajectory_ids=['0123456789abcdef'] if admitted else [])
        with tempfile.TemporaryDirectory() as tmp:
            state=Path(tmp)
            if prior:(state/'series.json').write_text(json.dumps(prior))
            with patch.object(module,'STATE',state),patch.object(module,'read_outcome',return_value=outcome):
                final=module.main()
        return final,fake.main.call_args_list

    def test_five_unique_native_seeds_are_finite(self):
        result,calls=self.exercise()
        self.assertEqual(result['status'],'completed');self.assertEqual(len(calls),5)
        self.assertEqual([call.args[1]['campaign_seed'] for call in calls],list(range(20261022,20261027)))
        self.assertTrue(all(call.kwargs['max_runtime_seconds']==2100 for call in calls))

    def test_actual_release_and_scientific_errors_stop_successors(self):
        for kwargs in [dict(good=False),dict(error=True)]:
            with self.subTest(kwargs=kwargs):
                result,calls=self.exercise(**kwargs)
                self.assertEqual(result['status'],'stopped_on_failure');self.assertEqual(len(calls),1)

    def test_admitted_output_stops_for_review(self):
        result,calls=self.exercise(admitted=True)
        self.assertEqual(result['status'],'stopped_for_review');self.assertEqual(len(calls),1)

    def test_completed_series_never_restarts(self):
        result,calls=self.exercise(prior=dict(status='completed'))
        self.assertEqual(calls,[])

    def test_fresh_cloud_state_blocks_spawn_even_without_prior_local_job(self):
        result,calls=self.exercise(idle=False)
        self.assertEqual(result['status'],'stopped_on_failure');self.assertEqual(calls,[])


if __name__=='__main__':unittest.main()
