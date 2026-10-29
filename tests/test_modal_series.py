"""The rolling seed series has finite spend and never advances past failed release."""
import importlib.util
import json
from pathlib import Path
import sys
import tempfile
import types
import unittest
from unittest.mock import Mock, patch


def atomic(path, value):
    path.write_text(json.dumps(value))


class SeriesTests(unittest.TestCase):
    def exercise(self, result, prior=None):
        driver = types.SimpleNamespace(HOURLY_MAX=4.649328, RUN="test_run", atomic=atomic,
                                       main=Mock(return_value=result))
        spec = importlib.util.spec_from_file_location("modal_series_test", Path(__file__).resolve().parents[1] / "scripts/drive_modal_series.py")
        series = importlib.util.module_from_spec(spec)
        with patch.dict(sys.modules, {"drive_modal_campaign": driver}):
            spec.loader.exec_module(series)
        with tempfile.TemporaryDirectory() as tmp:
            state = Path(tmp)
            if prior:
                (state / "series.json").write_text(json.dumps(prior))
            with patch.object(series, "STATE", state):
                series.main()
            return json.loads((state / "series.json").read_text()), driver.main.call_args_list

    def good(self):
        return {"gpu_released_verified": True, "status": "returned", "result": {"status": "completed"},
                "cpu_followup_returncode": 0, "compute_estimate_upper_usd": 1}

    def test_exactly_three_unique_seed_batches_then_stop(self):
        state, calls = self.exercise(self.good())
        self.assertEqual(state["status"], "series_completed")
        self.assertEqual(len(calls), 3)
        self.assertEqual([c.args[1]["campaign_seed"] for c in calls], [20261021, 20261022, 20261023])
        self.assertEqual(len({c.args[0] for c in calls}), 3)
        self.assertEqual(state["runtime_estimate_upper_usd"], 3)

    def test_no_next_allocation_without_verified_release(self):
        result = self.good()
        result["gpu_released_verified"] = False
        state, calls = self.exercise(result)
        self.assertEqual(len(calls), 1)
        self.assertEqual(state["status"], "stopped_on_execution_failure")

    def test_resuming_never_repeats_completed_batch(self):
        state, calls = self.exercise(self.good(), {"finished_batches": [0]})
        self.assertEqual(len(calls), 2)
        self.assertTrue(all(c.args[0] != "test_run" for c in calls))
        self.assertEqual(state["finished_batches"], [0, 1, 2])


if __name__ == "__main__":
    unittest.main()
