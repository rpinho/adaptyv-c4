"""Finite, declared failure-dependent search; no acceptance threshold relaxation."""
import copy
import re

REFERENCE = 'c4_native_apical_ablation_v3'
UNCENTERED = 'c4_native_largeinit_v4'
REPAIRED = 'c4_repaired_native_v9'
PROFILES = (
    'c4_local_frame_v5',
    'c4_local_frame_longer_v6',
    'c4_local_frame_foldweight_v7',
)
CHANGES = {
    REFERENCE: {},
    UNCENTERED: {'large_complex_initialization': True},
    REPAIRED: {'large_complex_initialization': True,'max_trajectories':1},
    PROFILES[0]: {'large_complex_initialization': True},
    PROFILES[1]: {'large_complex_initialization': True, 'binder_lengths': [90,110]},
    PROFILES[2]: {'large_complex_initialization': True, 'binder_lengths': [90,110], 'weights_plddt_loss': 0.5},
}


def profile_for_name(run_name):
    profile = re.sub(r'_replicate_[2-6]$', '', run_name)
    if profile not in CHANGES:
        raise ValueError('Unknown bounded campaign name')
    return profile


def make_config(base, run_name, overrides, project_folder):
    if set(overrides)-{'campaign_seed'}:
        raise ValueError('Only a campaign seed may be overridden outside the declared profiles')
    profile = profile_for_name(run_name)
    config = copy.deepcopy(base)
    config.update(copy.deepcopy(CHANGES[profile]))
    config.update(overrides)
    config.update(campaign_name=run_name, project_folder=project_folder)
    if profile in PROFILES or profile==REPAIRED:
        config['targets'][0]['target_path']='/opt/c4/local-frame-target.pdb'
    else:
        config['targets'][0]['target_path']='/opt/c4/target.pdb'
    if config['filters']!=base['filters']:
        raise ValueError('A search profile changed an acceptance gate')
    return config


def next_profile(current, accepted, terminated_stages):
    if accepted:
        return None, 'native_context_outputs_require_acquisition_review'
    if not terminated_stages or any(stage not in {'screen','refine'} for stage in terminated_stages):
        return None, 'needs_interface_or_later_stage_diagnosis'
    # A stage name does not distinguish singular frames, off-patch geometry,
    # confidence failures or individual filters. Changing length/weight here
    # previously spent money without that causal evidence.
    return None, 'physical_failure_review_required'
