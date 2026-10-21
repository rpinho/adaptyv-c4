"""The continuation may start before its monitor has written the first status."""
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
import continue_after_gpu as queue


class ContinuationTests(unittest.TestCase):
    def test_missing_initial_monitor_file_waits_then_runs_after_verified_stop(self):
        with tempfile.TemporaryDirectory() as tmp:
            out=Path(tmp)
            def publish_monitor(_seconds):
                (out/'supervisor.json').write_text(json.dumps({'status':'pod_stopped','stop_verified':True}))
            with patch.object(queue.time,'sleep',side_effect=publish_monitor) as sleep, \
                 patch.object(queue.subprocess,'run',return_value=subprocess.CompletedProcess([],0)) as run:
                queue.main(out,'results/test_followup')
            self.assertEqual(sleep.call_count,1)
            self.assertEqual(run.call_count,1)
            self.assertEqual(json.loads((out/'continuation.json').read_text())['status'],'cpu_followup_finished')

    def test_unverified_stop_never_launches_scoring(self):
        with tempfile.TemporaryDirectory() as tmp:
            out=Path(tmp);(out/'supervisor.json').write_text(json.dumps(
                {'status':'STOP_FAILED_REQUIRES_ATTENTION','stop_verified':False}))
            with patch.object(queue.subprocess,'run') as run:queue.main(out)
            run.assert_not_called()
            self.assertEqual(json.loads((out/'continuation.json').read_text())['status'],'blocked_on_verified_pod_stop')


if __name__=='__main__':unittest.main()
