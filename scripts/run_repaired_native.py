"""Exactly one matched repaired trajectory, no automatic successor."""
from pathlib import Path
import sys
import drive_modal_campaign as driver
import modal_search_policy as policy

ROOT=Path(__file__).resolve().parents[1]


if __name__=='__main__':
    run=sys.argv[1] if len(sys.argv)>1 else policy.REPAIRED
    folders={policy.REPAIRED:'modal_repaired_native_v9',policy.INTERFACE8:'modal_interface8_v10',policy.FORCED:'modal_forced_native_restore_v11',policy.CONFIDENCE:'modal_forced_confidence_v12'}
    if run not in folders:raise ValueError('Only the declared single-trajectory tests are allowed')
    driver.STATE=ROOT/'results'/folders[run]
    driver.CANCEL_ROOT=driver.STATE
    driver.main(run,{'campaign_seed':20261021},max_runtime_seconds=2100)
