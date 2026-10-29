"""Paid launch guards; no real API, SSH or GPU actions in these tests."""
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
import wait_for_gpu_capacity as queue


class CapacityQueueTests(unittest.TestCase):
    def test_capacity_failure_never_launches_workers(self):
        def call(name,args):
            if name=='get-pod':return {'status':'EXITED'}
            if name=='get-gpu-type':return {'price':{'secure':1.79}}
            raise queue.CapacityUnavailable('occupied')
        with tempfile.TemporaryDirectory() as tmp, patch.object(queue,'OUT',Path(tmp)), \
             patch.object(queue,'call',side_effect=call),patch.object(queue.subprocess,'run') as run, \
             patch.object(queue,'stop_verified') as stop:
            with self.assertRaises(queue.CapacityUnavailable):queue.launch({})
            run.assert_not_called();stop.assert_not_called()

    def test_existing_running_pod_is_left_alone(self):
        with patch.object(queue,'call',return_value={'status':'RUNNING'}) as call:
            self.assertFalse(queue.launch({}))
            self.assertEqual(call.call_count,1)

    def test_price_above_quote_never_starts_pod(self):
        with patch.object(queue,'call',side_effect=[{'status':'EXITED'},{'price':{'secure':2.5}}]) as call:
            with self.assertRaisesRegex(RuntimeError,'prequoted'):queue.launch({})
            self.assertEqual(call.call_count,2)

    def test_missing_gpu_after_start_stops_the_pod(self):
        def call(name,args):
            if name=='get-pod':
                if call.started:return {'status':'RUNNING','gpu':{'count':0},'ssh':{'direct':{'host':'test','port':22}}}
                return {'status':'EXITED'}
            if name=='get-gpu-type':return {'price':{'secure':1.79}}
            call.started=True;return {}
        call.started=False
        with tempfile.TemporaryDirectory() as tmp,patch.object(queue,'OUT',Path(tmp)), \
             patch.object(queue,'call',side_effect=call),patch.object(queue,'stop_verified',return_value=True) as stop:
            state={}
            with self.assertRaisesRegex(RuntimeError,'allocation missing'):queue.launch(state)
            stop.assert_called_once();self.assertTrue(state['stop_verified'])


if __name__=='__main__':unittest.main()
