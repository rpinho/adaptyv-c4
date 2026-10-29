"""Exactly one real AF2 forward/gradient call; software calibration, not a design campaign."""
import hashlib
import json
from pathlib import Path
import sys
import time

import jax
import numpy as np
from bindcraft.af2 import AlphaFoldDesignModel,MULTIMER_POOL
from bindcraft.model_weights import alphafold_parameters
from bindcraft.protein import write_structure
from bindcraft.protein_preparation import initialize_design_trajectory,receptor_chain_layouts
from bindcraft.sequence_optimization import LogitSequenceOptimizer
from bindcraft.settings import read_settings,build_design_settings


def main(out):
    out=Path(out);started=time.time()
    if any(device.platform!='gpu' for device in jax.devices()):
        raise RuntimeError('Refusing CPU fallback for the GPU diagnostic')
    settings=read_settings('/opt/c4/campaign.json',{
        'large_complex_initialization':True,
        'targets':[dict(name='hTfR1_full_monomer',target_path='/opt/c4/local-frame-target.pdb',chains='A',hotspots='A205,A206,A208,A209,A210',weight=1.)]})
    design=build_design_settings(settings)
    key=jax.random.fold_in(jax.random.PRNGKey(20261021),1)
    states,multi,losses=initialize_design_trajectory(design,jax.random.split(key)[0])
    binder_name=design.designed_binder_chain
    binder_length=len(next(iter(states.values()))[binder_name])
    if binder_length!=69:raise ValueError('Prepared sequence/length does not match the v5 first-trajectory control')
    model=AlphaFoldDesignModel(presets=MULTIMER_POOL,models=MULTIMER_POOL,
        data_dir=alphafold_parameters(download=False),key=key,num_recycle=settings['design_recycles'],
        large_complex_initialization=True,use_remat=True,dropout=True,
        multi_chain_binders=multi,amino_acid_bias=design.binder.amino_acid_bias)
    optimizer=LogitSequenceOptimizer(iterations=50,start_softmax_weight=0.,end_softmax_weight=.9,multi_chain_binders=multi)
    soft,hard,temp,scale=optimizer.sequence_parameters()
    chosen=model._sample_design_model()
    predictions,gradients,loss=model.sequence_gradients(states,losses,model=chosen,softmax_weight=soft,
        one_hot_weight=hard,temperature=temp,logit_scale=scale)
    jax.block_until_ready((gradients,loss))
    if not np.isfinite(float(loss)) or not all(np.isfinite(np.asarray(v)).all() for v in gradients.values()):
        raise ValueError('Nonfinite real forward loss or sequence gradient')
    prediction=next(iter(predictions.values()))
    binder=prediction.protein_complex[binder_name]
    target=next(p for name,p in prediction.protein_complex.items() if name!=binder_name)
    ca=np.asarray(binder.atoms[:,1],dtype=float)
    distances=np.linalg.norm(np.diff(ca,axis=0),axis=-1)
    patch=np.asarray(target.atoms[:,1])[np.isin(np.asarray(target.residue_index),[205,206,208,209,210])]
    if len(patch)!=5:raise ValueError('Physical hotspot mapping failed')
    unique=len(np.unique(ca,axis=0));finite=bool(np.isfinite(ca).all())
    noncollapsed=finite and unique>=.95*len(ca) and float(distances.min())>.1
    gradient_norm=float(np.linalg.norm(np.asarray(gradients[binder_name],dtype=np.float64)))
    gradient_norm_valid=bool(np.isfinite(gradient_norm) and gradient_norm>0)
    write_structure(prediction.protein_complex,str(out/'first_forward.cif'),receptor_chains=receptor_chain_layouts(design))
    result=dict(status='software_diagnostic_passed' if noncollapsed and gradient_norm_valid else 'software_diagnostic_failed',
        actual_forward_and_gradient=True,loss=float(loss),all_gradients_finite=True,binder_gradient_norm=gradient_norm,
        binder_CA_count=len(ca),unique_binder_CA_coordinates=unique,binder_CA_span_A=np.ptp(ca,axis=0).tolist(),
        adjacent_CA_min_A=float(distances.min()),adjacent_CA_median_A=float(np.median(distances)),adjacent_CA_max_A=float(distances.max()),
        minimum_hotspot_to_binder_CA_A=float(np.linalg.norm(patch[:,None,:]-ca[None,:,:],axis=-1).min()),
        seed=20261021,trajectory_number=1,trajectory_key=np.asarray(key).tolist(),binder_length=binder_length,
        selected_AF2_model=chosen,devices=[str(d) for d in jax.devices()],
        elapsed_seconds=time.time()-started,patch=json.loads(Path('/opt/bindcraft/frame_guard_manifest.json').read_text()),
        limitation='One unoptimized prediction and gradient. Not a qualified candidate, affinity, pH release or target-acquisition success.')
    (out/'diagnostic.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result),flush=True)
    if result['status']!='software_diagnostic_passed':raise RuntimeError('Software geometry/gradient diagnostic failed')


if __name__=='__main__':main(sys.argv[1])
