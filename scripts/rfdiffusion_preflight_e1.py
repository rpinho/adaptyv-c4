"""E1 (His707 patch) RFdiffusion CPU input/features/writer check."""
import hashlib
import json
import os
from pathlib import Path
import random
import sys

SOURCE = Path('/app/RFdiffusion')
TARGET = Path('/opt/c4/E1_his707_clean.pdb')
TARGET_SHA = 'fa71cc00c18ac38dd05aac4f0f6dcb6e0ad6cf43b5633e0e6f962f544cd4083c'
REVISION = '86507b6538f51fce57b5a72477165f03999ed7ae'
HOTSPOTS = (695, 704, 706, 707, 711, 714, 730)
SEED = 20261023

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def configuration():
    from hydra import compose, initialize_config_dir
    from hydra.core.hydra_config import HydraConfig
    from omegaconf import OmegaConf
    overrides = [
        f'inference.input_pdb={TARGET}',
        'inference.output_prefix=/tmp/c4-rfd-e1-preflight',
        'inference.model_directory_path=/opt/rfd-models',
        'inference.num_designs=1',
        f'inference.design_startnum={SEED}',
        'inference.deterministic=true',
        'inference.write_trajectory=false',
        'inference.final_step=1',
        'diffuser.T=50',
        'contigmap.contigs=["A260-268 A300-315 A322-330 A532-551 A616-630 A680-744/0 B317-326 B686-696/0 60-100"]',
        'ppi.hotspot_res=[' + ','.join('A' + str(i) for i in HOTSPOTS) + ']',
    ]
    with initialize_config_dir(version_base=None, config_dir=str(SOURCE / 'config/inference')):
        full = compose(config_name='base', overrides=overrides, return_hydra_config=True)
    HydraConfig.instance().set_config(full)
    return OmegaConf.masked_copy(full, [k for k in full if k != 'hydra'])


def main(destination):
    import subprocess
    import numpy as np
    import torch
    from Bio.PDB import PDBParser
    from rfdiffusion.inference import utils as iu
    from rfdiffusion.util import writepdb
    from omegaconf import OmegaConf
    assert not torch.cuda.is_available(), 'CPU image check must not allocate a GPU'
    assert sha(TARGET) == TARGET_SHA
    actual = subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=SOURCE, text=True).strip()
    assert actual == REVISION
    torch.manual_seed(SEED); np.random.seed(SEED); random.seed(SEED)
    sampler = iu.sampler_selector(configuration())
    assert sampler._conf.diffuser.T == 50 and sampler.t_step_input == 50
    assert sampler._conf.diffuser.partial_T is None and sampler._conf.inference.final_step == 1
    assert Path(sampler.ckpt_path).name == 'Complex_base_ckpt.pt'
    torch.manual_seed(SEED); np.random.seed(SEED); random.seed(SEED)
    xyz, seq = sampler.sample_init()
    refs = sampler.contig_map.ref
    target_indices_A = [i for i, x in enumerate(refs) if x[0] == 'A']
    target_indices_B = [i for i, x in enumerate(refs) if x[0] == 'B']
    binder_indices = [i for i, x in enumerate(refs) if x == ('_', '_')]
    assert len(target_indices_A) == 134, f'Expected 134 chain A residues, got {len(target_indices_A)}'
    assert len(target_indices_B) == 21, f'Expected 21 chain B residues, got {len(target_indices_B)}'
    assert 60 <= len(binder_indices) <= 100, f'Binder length {len(binder_indices)} outside 60-100'
    target_indices = target_indices_A + target_indices_B
    assert sampler.mask_str[0, target_indices].all() and sampler.mask_seq[0, target_indices].all()
    assert not sampler.mask_str[0, binder_indices].any() and not sampler.mask_seq[0, binder_indices].any()
    expected_hot = [i for i, x in enumerate(refs) if x in [('A', n) for n in HOTSPOTS]]
    features = sampler._preprocess(seq, xyz, sampler.t_step_input)
    actual_hotspots = torch.where(features[5][0, 0, :, -1] == 1)[0].tolist()
    assert actual_hotspots == expected_hot, (actual_hotspots, expected_hot)
    from diffusion_worker import tensor_mapping
    mapping_payload = json.dumps(tensor_mapping(sampler), sort_keys=True)
    checkpoint = Path(sampler.ckpt_path)
    result = dict(
        status='passed',
        scope='E1 (His707) CPU preparation/features/writer only; no neural forward',
        epitope='E1_his707',
        revision=actual, target_sha256=sha(TARGET), checkpoint=str(checkpoint),
        checkpoint_sha256=sha(checkpoint),
        target_residues_A=len(target_indices_A), target_residues_B=len(target_indices_B),
        binder_residues=len(binder_indices), seed=SEED,
        hotspot_author_ids=list(HOTSPOTS),
        actual_hotspot_tensor_indices=actual_hotspots,
        fixed_target_sequence_and_structure=True,
        diffusion_steps=50, final_step=1, partial_T=None,
        torch=torch.__version__, cuda_build=torch.version.cuda,
        worker_sha256=sha(__file__),
        config=OmegaConf.to_container(sampler._conf, resolve=True),
    )
    result['generation_worker_sha256'] = sha('/opt/c4/diffusion_worker.py')
    result['tensor_mapping_sha256'] = hashlib.sha256(mapping_payload.encode()).hexdigest()
    Path(destination).write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({k: v for k, v in result.items() if k != 'config'}), flush=True)


if __name__ == '__main__':
    main(sys.argv[1])
