"""Actual RFdiffusion CPU input/features/writer check; no neural forward."""
import hashlib
import json
import os
from pathlib import Path
import random
import sys

SOURCE = Path('/app/RFdiffusion')
TARGET = Path('/opt/c4/native-apical-search-only.pdb')
TARGET_SHA = '03abdb9bc5536396de3489e7a6cf3eec88066a1a5f36d5d1d78b7342ec1ec0b2'
REVISION = '86507b6538f51fce57b5a72477165f03999ed7ae'
HOTSPOTS = (205,206,208,209,210)
SEED = 20261022

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def configuration():
    from hydra import compose,initialize_config_dir
    from hydra.core.hydra_config import HydraConfig
    from omegaconf import OmegaConf
    overrides=[f'inference.input_pdb={TARGET}',
        'inference.output_prefix=/tmp/c4-rfd-preflight','inference.model_directory_path=/opt/rfd-models',
        'inference.num_designs=1',f'inference.design_startnum={SEED}',
        'inference.deterministic=true','inference.write_trajectory=false',
        'inference.final_step=1','diffuser.T=50',
        'contigmap.contigs=["A189-383/0 69-69"]',
        'ppi.hotspot_res=['+','.join('A'+str(i) for i in HOTSPOTS)+']']
    with initialize_config_dir(version_base=None,config_dir=str(SOURCE/'config/inference')):
        full=compose(config_name='base',overrides=overrides,return_hydra_config=True)
    HydraConfig.instance().set_config(full)
    return OmegaConf.masked_copy(full,[k for k in full if k!='hydra'])

def main(destination):
    import subprocess
    import numpy as np
    import torch
    from Bio.PDB import PDBParser
    from rfdiffusion.inference import utils as iu
    from rfdiffusion.util import writepdb, num2aa
    from omegaconf import OmegaConf
    assert not torch.cuda.is_available(), 'CPU image check must not allocate a GPU'
    assert sha(TARGET) == TARGET_SHA
    actual = subprocess.check_output(['git','rev-parse','HEAD'],cwd=SOURCE,text=True).strip()
    assert actual == REVISION
    torch.manual_seed(SEED);np.random.seed(SEED);random.seed(SEED)
    sampler = iu.sampler_selector(configuration())
    assert sampler._conf.diffuser.T == 50 and sampler.t_step_input == 50
    assert sampler._conf.diffuser.partial_T is None and sampler._conf.inference.final_step == 1
    assert Path(sampler.ckpt_path).name == 'Complex_base_ckpt.pt'
    assert sampler.preprocess_conf.d_t1d >= 24
    torch.manual_seed(SEED);np.random.seed(SEED);random.seed(SEED)
    xyz, seq = sampler.sample_init()
    assert sampler.binderlen == 69
    assert len(sampler.contig_map.ref) == 264
    refs = sampler.contig_map.ref
    target_indices = [i for i,x in enumerate(refs) if x[0] == 'A']
    binder_indices = [i for i,x in enumerate(refs) if x == ('_','_')]
    assert [refs[i] for i in target_indices] == [('A',n) for n in range(189,384)]
    assert binder_indices == list(range(69))
    assert sampler.mask_str[0,target_indices].all() and sampler.mask_seq[0,target_indices].all()
    assert not sampler.mask_str[0,binder_indices].any() and not sampler.mask_seq[0,binder_indices].any()
    expected = [i for i,x in enumerate(refs) if x in [('A',n) for n in HOTSPOTS]]
    receptor_local=[i for i,x in enumerate([refs[j] for j in target_indices]) if x in [('A',n) for n in HOTSPOTS]]
    assert list(sampler.hotspot_0idx) == receptor_local
    features = sampler._preprocess(seq,xyz,sampler.t_step_input)
    actual_hotspots = torch.where(features[5][0,0,:,-1] == 1)[0].tolist()
    assert actual_hotspots == expected, (actual_hotspots,expected)
    assert [sampler.idx_pdb[i] for i in target_indices] == list(range(189,384))
    assert all(sampler.chain_idx[i] == 'A' for i in target_indices)
    assert [sampler.idx_pdb[i] for i in binder_indices] == list(range(1,70))
    assert all(sampler.chain_idx[i] == 'B' for i in binder_indices)
    fixed = torch.tensor([sampler.target_feats['seq'][i].item() for i in range(195)])
    assert torch.equal(seq[target_indices].argmax(-1),fixed)
    final_seq = torch.where(seq.argmax(-1)==21,7,seq.argmax(-1))
    output = Path('/tmp/c4-rfd-writer-roundtrip.pdb')
    writepdb(str(output),xyz[:,:4],final_seq,sampler.binderlen,
             chain_idx=sampler.chain_idx,idx_pdb=sampler.idx_pdb)
    model = PDBParser(QUIET=True).get_structure('roundtrip',output)[0]
    original = PDBParser(QUIET=True).get_structure('input',TARGET)[0]['A']
    assert set(model.child_dict) == {'A','B'}
    assert [r.id[1] for r in model['A']] == list(range(189,384))
    assert [r.resname for r in model['A']] == [r.resname for r in original]
    assert [r.id[1] for r in model['B']] == list(range(1,70))
    assert all(r.resname == 'GLY' for r in model['B'])
    for r,o in zip(model['A'],original):
        for atom in ('N','CA','C','O'):
            assert atom in r and np.isfinite(r[atom].coord).all()
    checkpoint = Path(sampler.ckpt_path)
    result = dict(status='passed',scope='CPU preparation/features/writer only; no neural forward',
        revision=actual,target_sha256=sha(TARGET),checkpoint=str(checkpoint),checkpoint_sha256=sha(checkpoint),
        target_residues=195,binder_residues=69,seed=SEED,hotspot_author_ids=list(HOTSPOTS),
        actual_hotspot_tensor_indices=actual_hotspots,fixed_target_sequence_and_structure=True,
        hotspot_receptor_local_indices=receptor_local,diffusion_steps=50,final_step=1,partial_T=None,
        actual_output_target_chain='A',actual_output_target_author_range=[189,383],
        actual_output_binder_chain='B',actual_output_binder_author_range=[1,69],
        legacy_contig_output_mapping_not_used=True,torch=torch.__version__,cuda_build=torch.version.cuda,
        worker_sha256=sha(__file__),config=OmegaConf.to_container(sampler._conf,resolve=True))
    Path(destination).write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='config'}),flush=True)

if __name__ == '__main__':
    main(sys.argv[1])
