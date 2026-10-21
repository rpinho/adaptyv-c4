"""Paid containers must disappear before CPU follow-up; uncertain spawns never repeat."""
import dataclasses
import importlib.util
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import types
import unittest
from unittest.mock import Mock, patch


@dataclasses.dataclass
class Stats:
    backlog: int = 0
    num_total_runners: int = 0
    num_running_inputs: int = 0
    input_headroom: int = 0


class ModalDriverTests(unittest.TestCase):
    def driver(self):
        fake = types.SimpleNamespace(Function=Mock(), FunctionCall=Mock(), Volume=Mock())
        spec = importlib.util.spec_from_file_location("modal_driver_test", Path(__file__).resolve().parents[1] / "scripts/drive_modal_campaign.py")
        module = importlib.util.module_from_spec(spec)
        with patch.dict(sys.modules, {"modal": fake}):
            spec.loader.exec_module(module)
        return module, fake

    def exercise(self, stats, cancel=False):
        driver, fake = self.driver()
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            state = root / "results/modal_campaign_v3"
            state.mkdir(parents=True)
            if cancel:
                (state / "CANCEL").touch()
            fn = fake.Function.from_name.return_value
            call = fn.spawn.return_value
            call.object_id = "fc-test"
            call.get.return_value = {"status": "completed"}
            fn.get_current_stats.side_effect = stats
            with patch.object(driver, "ROOT", root), patch.object(driver, "STATE", state), \
                    patch.object(driver, "mirror"), patch.object(driver.time, "sleep"), \
                    patch.object(driver.subprocess, "run", return_value=subprocess.CompletedProcess([], 0)) as followup:
                driver.main()
            return json.loads((state / "supervisor.json").read_text()), call, followup

    def test_result_waits_for_container_shutdown_before_cpu(self):
        status, call, followup = self.exercise([Stats(num_total_runners=1), Stats(num_total_runners=1), Stats()])
        self.assertTrue(status["gpu_released_verified"])
        followup.assert_called_once()
        call.cancel.assert_not_called()

    def test_cancel_terminates_container_and_failed_release_blocks_cpu(self):
        status, call, followup = self.exercise([Stats(num_total_runners=1)] * 31, cancel=True)
        call.cancel.assert_called_once_with(terminate_containers=True)
        self.assertFalse(status["gpu_released_verified"])
        followup.assert_not_called()

    def test_uncertain_previous_spawn_never_spawns_again(self):
        driver, fake = self.driver()
        with tempfile.TemporaryDirectory() as tmp:
            state = Path(tmp)
            (state / "supervisor.json").write_text(json.dumps({"spawn_intent": 123}))
            with patch.object(driver, "STATE", state), self.assertRaisesRegex(RuntimeError, "uncertain outcome"):
                driver.main()
            fake.Function.from_name.return_value.spawn.assert_not_called()


if __name__ == "__main__":
    unittest.main()
