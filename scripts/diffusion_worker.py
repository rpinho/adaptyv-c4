"""One fixed target-conditioned backbone; no silent CPU fallback or successors."""
from pathlib import Path
import json
import sys
import time
from rfdiffusion_preflight import configuration,sha,TARGET,TARGET_SHA,REVISION

def main(destination):
    import torch
    from omegaconf import OmegaConf
    from rfdiffusion.inference import utils as iu
    from rfdiffusion.util import writepdb
    if not torch.cuda.is_available():raise RuntimeError('RFdiffusion CPU fallback refused')
    if 'A100' not in torch.cuda.get_device_name():raise RuntimeError('Declared A100 runtime required')
    root=Path(destination)
    proof=json.loads(Path('/opt/c4/rfdiffusion_preflight.json').read_text())
    if proof['status']!='passed' or sha(TARGET)!=TARGET_SHA:raise RuntimeError('Input CPU proof mismatch')
    if sha(proof['checkpoint'])!=proof['checkpoint_sha256']:raise RuntimeError('Checkpoint identity changed')
    conf=configuration();conf.inference.output_prefix=str(root/'backbone')
    # Upstream checkpoint assembly can overwrite defaults. Inspect actual resolved sampler.
    import random
    import numpy as np
    from rfdiffusion_preflight import SEED,HOTSPOTS
    torch.manual_seed(SEED);np.random.seed(SEED);random.seed(SEED)
    sampler=iu.sampler_selector(conf)
    if sampler._conf.diffuser.T!=50 or sampler.t_step_input!=50 or sampler._conf.diffuser.partial_T is not None or sampler._conf.inference.final_step!=1:
        raise RuntimeError('Resolved sampler differs from finite diffusion protocol')
    if sampler._conf.inference.num_designs!=1:raise RuntimeError('Backbone count differs')
    (root/'effective_config.json').write_text(json.dumps(OmegaConf.to_container(sampler._conf,resolve=True),indent=2)+'\n')
    torch.manual_seed(SEED);np.random.seed(SEED);random.seed(SEED)
    started=time.time();xyz,seq=sampler.sample_init()
    mapping=[dict(tensor_index=i,input_reference=list(ref),output_chain=sampler.chain_idx[i],
                  output_author_id=sampler.idx_pdb[i]) for i,ref in enumerate(sampler.contig_map.ref)]
    (root/'tensor_mapping.json').write_text(json.dumps(mapping,indent=2)+'\n')
    expected_hotspots=[i for i,x in enumerate(sampler.contig_map.ref) if x in [('A',n) for n in HOTSPOTS]]
    actual=torch.where(sampler._preprocess(seq,xyz,50)[5][0,0,:,-1]==1)[0].tolist()
    if actual!=expected_hotspots or len(actual)!=5:raise RuntimeError('Actual GPU-input hotspot features differ')
    # Same pinned sample_init/sample_step/writepdb interfaces as official run_inference.
    for t in range(50,0,-1):
        _,xyz,_,plddt=sampler.sample_step(t=t,x_t=xyz,seq_init=seq,final_step=1)
        with (root/'diffusion_steps.jsonl').open('a') as out:
            out.write(json.dumps({'step':51-t,'timestep':t,'unix':time.time(),'mean_pLDDT':float(plddt.mean())})+'\n')
    final_seq=torch.where(seq.argmax(-1)==21,7,seq.argmax(-1))
    expected=root/('backbone_'+str(conf.inference.design_startnum)+'.pdb')
    writepdb(str(expected),xyz[:,:4],final_seq,sampler.binderlen,
             chain_idx=sampler.chain_idx,idx_pdb=sampler.idx_pdb)
    if not expected.exists():raise RuntimeError('Expected single fresh backbone is missing')
    if len(list(root.glob('backbone_*.pdb')))!=1:raise RuntimeError('Unexpected backbone count')
    result=dict(status='completed',backbone_file=expected.name,backbone_sha256=sha(expected),
        backbone_count=1,diffusion_steps=50,revision=REVISION,gpu=torch.cuda.get_device_name(),
        elapsed_seconds=time.time()-started,qualified_candidates=0,
        interpretation='Poly-glycine backbone only; CPU native-context audit then sequence design required')
    (root/'generation_results.json').write_text(json.dumps(result,indent=2)+'\n')

if __name__=='__main__':main(sys.argv[1])
