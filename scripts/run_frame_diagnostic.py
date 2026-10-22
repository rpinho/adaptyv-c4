"""Persistent monitor for exactly one repaired forward/gradient diagnostic."""
from pathlib import Path
import drive_modal_campaign as driver

if __name__=='__main__':
    driver.STATE=Path(__file__).resolve().parents[1]/'results/modal_frame_diagnostic_v8'
    driver.CANCEL_ROOT=driver.STATE
    driver.main('c4_frame_guard_diagnostic_v8',{},function_name='frame_diagnostic',perform_cpu_followup=False)
