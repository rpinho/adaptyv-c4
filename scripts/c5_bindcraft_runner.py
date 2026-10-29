#!/usr/bin/env python3
"""
C5 BindCraft multi-template runner.

Runs BindCraft's multi-template hallucination with:
  - G12D pMHC (target): bind this
  - WT pMHC (detarget): avoid binding
  - Peptide-free HLA+b2m (detarget): avoid binding

Uses the Multi_BC_Original notebook's engine, driven headlessly.
"""

import copy
import json
import os
import sys
import time

import numpy as np

# BindCraft must be on the path
BINDCRAFT_DIR = "/root/BindCraft"
sys.path.insert(0, BINDCRAFT_DIR)

# Import BindCraft functions
from functions import *
from functions.colabdesign_utils import *
from functions.generic_utils import *
from functions.pyrosetta_utils import *

# Try to import multi-template utils
try:
    from functions.multi_colabdesign_utils import *
    HAS_MULTI = True
except ImportError:
    HAS_MULTI = False

SETTINGS_DIR = "/root/c5_design/settings"
OUTPUT_DIR = "/root/c5_design/output"

def load_settings():
    with open(os.path.join(SETTINGS_DIR, "c5_bindcraft_target.json")) as f:
        target_settings = json.load(f)
    with open(os.path.join(SETTINGS_DIR, "c5_bindcraft_advanced.json")) as f:
        advanced_settings = json.load(f)
    # Use default filters
    filters_path = os.path.join(BINDCRAFT_DIR, "settings_filters", "default_filters.json")
    with open(filters_path) as f:
        filters = json.load(f)
    return target_settings, advanced_settings, filters


def setup_templates(target_settings):
    """Build template_information dict from the target JSON."""
    templates = {}
    for t in target_settings["templates"]:
        name = t["template_name"]
        templates[name] = {
            "template_pdb": t["template_pdb"],
            "template_hostspot_residues": t.get("template_hostspot_residues", ""),
            "template_targeting_strategy": t.get("template_targeting_strategy", "target"),
            "template_mpnn_filters_path": t.get("template_mpnn_filters_path", ""),
            "loss_weighting": float(t.get("loss_weighting", 1.0)),
        }
    return templates


def main():
    print("=" * 60)
    print("C5 BINDCRAFT MULTI-TEMPLATE DESIGN RUN")
    print("=" * 60)
    print(f"Start time: {time.strftime('%Y-%m-%d %H:%M:%S UTC', time.gmtime())}")

    # Check GPU
    check_jax_gpu()

    # Load settings
    target_settings, advanced_settings, filters = load_settings()

    # General information (from target_settings)
    general_information = {
        "binder_name": target_settings["binder_name"],
        "chains": target_settings["chains"],
        "lengths": target_settings["lengths"],
        "number_of_final_designs": target_settings["number_of_final_designs"],
    }

    # Setup templates
    template_information = setup_templates(target_settings)

    print(f"\nTemplates:")
    for name, info in template_information.items():
        print(f"  {name}: {info['template_targeting_strategy']} "
              f"(weight={info['loss_weighting']}, "
              f"hotspots={info['template_hostspot_residues'] or 'none'})")

    print(f"\nBinder lengths: {general_information['lengths']}")
    print(f"Target designs: {general_information['number_of_final_designs']}")
    print(f"Max trajectories: {advanced_settings.get('max_trajectories', 'unlimited')}")

    # Check advanced settings
    advanced_settings = perform_advanced_settings_check(advanced_settings, BINDCRAFT_DIR)

    # Load AF2 models
    design_models, prediction_models, multimer_validation = load_af2_models(
        advanced_settings["use_multimer_design"]
    )

    # Generate directories
    design_path = target_settings["design_path"]
    design_paths = generate_directories(design_path)

    # Generate dataframes
    trajectory_labels, design_labels, final_labels = generate_dataframe_labels()

    trajectory_csv = os.path.join(design_path, 'trajectory_stats.csv')
    mpnn_csv = os.path.join(design_path, 'mpnn_design_stats.csv')
    final_csv = os.path.join(design_path, 'final_design_stats.csv')
    failure_csv = os.path.join(design_path, 'failure_csv.csv')

    create_dataframe(trajectory_csv, trajectory_labels)
    create_dataframe(mpnn_csv, design_labels)
    create_dataframe(final_csv, final_labels)

    # Initialize PyRosetta
    import pyrosetta as pr
    pr.init(
        f'-ignore_unrecognized_res -ignore_zero_occupancy -mute all '
        f'-holes:dalphaball {advanced_settings["dalphaball_path"]} '
        f'-corrections::beta_nov16 true -relax:default_repeats 1'
    )

    print(f"\nStarting design loop...")
    print(f"{'=' * 60}")

    script_start_time = time.time()
    trajectory_n = 1
    accepted_designs = 0
    chain = general_information["chains"]

    while True:
        # Check if we have enough designs
        if HAS_MULTI:
            if check_accepted_designs_multi(
                design_paths, template_information, general_information,
                final_labels, advanced_settings, target_settings, design_labels
            ):
                print("Target number of designs reached!")
                break
        else:
            if check_accepted_designs(
                design_paths, mpnn_csv, final_labels, final_csv,
                advanced_settings, target_settings, design_labels
            ):
                print("Target number of designs reached!")
                break

        if check_n_trajectories(design_paths, advanced_settings):
            print("Maximum trajectories reached.")
            break

        trajectory_start_time = time.time()
        seed = int(np.random.randint(0, high=999999, size=1, dtype=int)[0])
        samples = np.arange(
            min(general_information["lengths"]),
            max(general_information["lengths"]) + 1
        )
        length = int(np.random.choice(samples))
        helicity_value = load_helicity(advanced_settings)

        design_name = f'{general_information["binder_name"]}_l{length}_s{seed}'

        print(f"\n--- Trajectory {trajectory_n}: {design_name} (length={length}, seed={seed}) ---")

        # Create trajectory directories
        trajectory_dirs = [
            "Trajectory", "Trajectory/Relaxed", "Trajectory/LowConfidence",
            "Trajectory/Clashing", "Trajectory/Animation", "Trajectory/Plots"
        ]
        for d in trajectory_dirs:
            os.makedirs(os.path.join(design_path, d), exist_ok=True)

        try:
            if HAS_MULTI:
                # Multi-template hallucination
                trajectory_pdb, trajectory_data = binder_hallucination_multi(
                    template_information, design_paths, design_name,
                    advanced_settings, length, seed, chain,
                    helicity_value, design_models
                )
            else:
                # Fall back to single-template (G12D only)
                print("WARNING: multi_colabdesign_utils not available, using single-template mode")
                trajectory_pdb, trajectory_data = binder_hallucination(
                    template_information["g12d"]["template_pdb"],
                    template_information["g12d"]["template_hostspot_residues"],
                    design_paths, design_name,
                    advanced_settings, length, seed, chain,
                    helicity_value, design_models
                )

            elapsed = time.time() - trajectory_start_time
            print(f"  Trajectory completed in {elapsed:.0f}s")

            if trajectory_pdb is not None:
                # MPNN sequence design + validation
                print(f"  Running MPNN + validation...")

                mpnn_trajectories = mpnn_gen(
                    trajectory_pdb, design_paths, design_name,
                    advanced_settings, chain
                )

                if mpnn_trajectories:
                    for mpnn_design in mpnn_trajectories:
                        # Validate with AF2
                        pass_all = af2_initial_guess(
                            mpnn_design, design_paths, design_name,
                            advanced_settings, prediction_models,
                            multimer_validation, chain
                        )
                        if pass_all:
                            accepted_designs += 1

                total_elapsed = time.time() - script_start_time
                print(f"  Accepted so far: {accepted_designs} | "
                      f"Total time: {total_elapsed/3600:.1f}h")

        except Exception as e:
            print(f"  ERROR in trajectory {trajectory_n}: {e}")
            import traceback
            traceback.print_exc()

        trajectory_n += 1

    total_time = time.time() - script_start_time
    print(f"\n{'=' * 60}")
    print(f"DESIGN RUN COMPLETE")
    print(f"Trajectories: {trajectory_n - 1}")
    print(f"Accepted designs: {accepted_designs}")
    print(f"Total time: {total_time/3600:.1f} hours")
    print(f"End time: {time.strftime('%Y-%m-%d %H:%M:%S UTC', time.gmtime())}")

    # Save summary
    summary = {
        "trajectories": trajectory_n - 1,
        "accepted_designs": accepted_designs,
        "total_time_hours": round(total_time / 3600, 2),
        "settings": {
            "lengths": general_information["lengths"],
            "templates": list(template_information.keys()),
            "max_trajectories": advanced_settings.get("max_trajectories"),
        }
    }
    with open(os.path.join(design_path, "run_summary.json"), "w") as f:
        json.dump(summary, f, indent=2)

    return 0


if __name__ == "__main__":
    sys.exit(main())
