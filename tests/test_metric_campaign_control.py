"""Cloud spending controls verified without loading Modal or creating an app."""
import contextlib
import hashlib
import json
from pathlib import Path
from types import SimpleNamespace
import sys
import tempfile
import unittest
from unittest.mock import patch,Mock

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'results/generation_pivot'))
import run_metric_arms as runner


class MetricCampaignControl(unittest.TestCase):
    def test_release_requires_app_stop_as_well_as_zero_tasks(self):
        self.assertTrue(runner.apps_released([dict(app_id='owned',state='stopped',tasks='0')],['owned']))
        self.assertFalse(runner.apps_released([dict(app_id='owned',state='ephemeral',tasks='0')],['owned']))
        self.assertFalse(runner.apps_released([dict(app_id='owned',state='stopped',tasks='1')],['owned']))
        self.assertTrue(runner.apps_released([dict(app_id='other',state='ephemeral',tasks='9')],['owned']))
        self.assertFalse(runner.apps_released([],[]))

    def test_remote_code_is_not_loaded_before_completed_report_is_committed(self):
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder);source=root/'source.py';source.write_text('source only')
            with patch.object(runner,'ROOT',root),patch.object(runner,'SOURCE',source),patch.object(runner,'EXPECTED_SOURCE',runner.sha(source)),patch.object(sys,'argv',['runner','--execute']),patch.object(runner.importlib.util,'spec_from_file_location') as load:
                with self.assertRaises(FileNotFoundError):runner.main()
                load.assert_not_called()

    def test_function_failure_explicitly_stops_only_its_own_app(self):
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder);source=root/'source.py';source.write_text('source only')
            report=root/'results/novelty_predictor/adaptive_completion/predictions.json';report.parent.mkdir(parents=True);report.write_text('{}')
            (report.parent/'VALIDATION.json').write_text(json.dumps(dict(state='complete',predictions_sha256=runner.sha(report))))
            app=SimpleNamespace(app_id='ap-owned-test',run=lambda **kwargs:contextlib.nullcontext())
            call=Mock();call.object_id='fc-test';call.get.side_effect=RuntimeError('injected remote failure')
            function=Mock();function.with_options.return_value=function;function.spawn.return_value=call
            module=SimpleNamespace(app=app,generate=function,modal=SimpleNamespace(enable_output=lambda:contextlib.nullcontext()))
            spec=SimpleNamespace(loader=SimpleNamespace(exec_module=lambda target:None))
            def command(argv,**kwargs):
                return SimpleNamespace(returncode=0,stdout='[]' if 'list' in argv else 'Stopped owned app',stderr='')
            with patch.object(runner,'ROOT',root),patch.object(runner,'SOURCE',source),patch.object(runner,'EXPECTED_SOURCE',runner.sha(source)),patch.object(sys,'argv',['runner','--execute']),patch.object(runner.importlib.util,'spec_from_file_location',return_value=spec),patch.object(runner.importlib.util,'module_from_spec',return_value=module),patch.object(runner.subprocess,'check_output',return_value=report.read_bytes()),patch.object(runner.subprocess,'run',side_effect=command) as remote:
                with self.assertRaisesRegex(RuntimeError,'injected remote failure'):runner.main()
                self.assertEqual(function.spawn.call_count,2)
                self.assertEqual(remote.call_args_list[0].args[0][-2:],['stop','ap-owned-test'])
                state=json.loads((root/'results/rfd4_geom_20261023/STATUS.json').read_text())
                self.assertEqual(state['state'],'failed')
                self.assertTrue(state['app_stop_confirmed'])
                self.assertTrue(state['app_absent_from_active_listing'])


if __name__=='__main__':unittest.main()
