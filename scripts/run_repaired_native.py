"""Exactly one matched repaired trajectory, no automatic successor."""
from pathlib import Path
import drive_modal_campaign as driver

ROOT=Path(__file__).resolve().parents[1]


if __name__=='__main__':
    driver.STATE=ROOT/'results/modal_repaired_native_v9'
    driver.CANCEL_ROOT=driver.STATE
    driver.main('c4_repaired_native_v9',{'campaign_seed':20261021},max_runtime_seconds=2100)
