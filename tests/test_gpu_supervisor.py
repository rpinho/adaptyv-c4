"""Critical shutdown behavior, with all network/process actions mocked."""
import datetime
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import time
import unittest
from unittest.mock import patch

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
import supervise_gpu_campaign as supervisor


class ShutdownTests(unittest.TestCase):
    def run_case(self, inspections, worker_output='0 %, 2500 MiB\n', mirror_error=False):
        with tempfile.TemporaryDirectory() as tmp:
            with patch.object(supervisor,'OUT',Path(tmp)), \
                 patch.object(supervisor,'inspect',side_effect=inspections), \
                 patch.object(supervisor,'request',return_value={'content':[]}) as req, \
                 patch.object(supervisor.subprocess,'run',return_value=subprocess.CompletedProcess([],0,worker_output,'')), \
                 patch.object(supervisor,'mirror',side_effect=RuntimeError('mirror failure') if mirror_error else None,return_value={}), \
                 patch.object(supervisor.time,'sleep'):
                supervisor.main()
                result=json.loads((Path(tmp)/'supervisor.json').read_text())
                calls=req.call_args_list
            return result,calls

    def pod(self,status='RUNNING',age=60):
        return dict(status=status,cost=1.79,startedAt=datetime.datetime.fromtimestamp(
            time.time()-age,datetime.timezone.utc).isoformat())

    def test_worker_exit_stops_pod_and_reads_back(self):
        result,calls=self.run_case([self.pod(),self.pod(),self.pod('EXITED')])
        self.assertTrue(result['stop_verified'])
        self.assertEqual(result['stop_reason'],'worker_exited')
        self.assertEqual(calls[0].args[1]['arguments']['body'],{'action':'stop'})

    def test_budget_cap_stops_even_with_worker_alive(self):
        result,calls=self.run_case([self.pod(age=7300),self.pod(),self.pod('EXITED')],worker_output='Rl\n90 %, 30000 MiB')
        self.assertTrue(result['stop_verified'])
        self.assertEqual(result['stop_reason'],'compute_cap')

    def test_mirror_failure_does_not_prevent_stop(self):
        result,calls=self.run_case([self.pod(),self.pod(),self.pod('EXITED')],mirror_error=True)
        self.assertTrue(result['stop_verified'])
        self.assertTrue(result['errors'])

    def test_initial_read_failure_still_attempts_stop(self):
        result,calls=self.run_case([RuntimeError('initial read failed'),
                                    RuntimeError('pre-stop read failed'),self.pod('EXITED')])
        self.assertTrue(result['stop_verified'])
        self.assertEqual(result['stop_reason'],'supervisor_exception')
        self.assertEqual(len(calls),1)


if __name__=='__main__':unittest.main()
