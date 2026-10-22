"""Observe real optimization steps and abort numerical collapse; leave design gates intact."""
import json
import math
from pathlib import Path
import sys


def geometry_metrics(ca,patch):
    if len(ca)<2 or len(patch)!=5:raise ValueError('Unexpected binder or requested hotspot mapping')
    if not all(math.isfinite(float(x)) for points in (ca,patch) for p in points for x in p):
        raise ValueError('Nonfinite predicted coordinates')
    distances=[math.dist(a,b) for a,b in zip(ca,ca[1:])]
    unique=len({tuple(float(x) for x in p) for p in ca})
    if unique<.95*len(ca) or min(distances)<=.1:
        raise ValueError('Predicted binder collapsed during real optimization')
    return dict(binder_CA_count=len(ca),unique_binder_CA_coordinates=unique,
                adjacent_CA_min_A=min(distances),adjacent_CA_max_A=max(distances),
                minimum_hotspot_to_binder_CA_A=min(math.dist(a,b) for a in ca for b in patch))


def preflight(config_path,output_path):
    """Prepare the REAL first trajectory and test actual frames without a network/GPU."""
    import hashlib
    import importlib.metadata
    import jax
    import numpy as np
    import bindcraft.af2 as af2
    from bindcraft.settings import read_settings,build_design_settings
    from bindcraft.loss import build_losses
    from bindcraft.protein_preparation import initialize_design_trajectory
    from bindcraft.protein import Protein
    from bindcraft.epitope_targeting import force_epitope_targeting
    from bindcraft.frame_guard import safe_initial_atom_positions
    from bindcraft.af.alphafold.common import residue_constants
    from bindcraft.af.alphafold.model.quat_affine import make_transform_from_reference
    settings=read_settings(config_path)
    raw_settings=json.loads(Path(config_path).read_text())
    losses=build_losses(settings,seed=20261021)
    effective_weight=float(losses['plddt_loss'].weight)
    if effective_weight!=raw_settings['weights_plddt_loss']:
        raise ValueError('Parsed/built pLDDT loss weight differs from declared input')
    if effective_weight==.5:
        baseline=build_losses(settings|{'weights_plddt_loss':.1},seed=20261021)
        if set(baseline)!=set(losses):raise ValueError('Confidence ablation changed active losses')
        for name,entry in losses.items():
            expected=.1 if name=='plddt_loss' else entry.weight
            if baseline[name].weight!=expected or baseline[name].required_states!=entry.required_states:
                raise ValueError('Confidence ablation changed another built objective')
            if baseline[name].function.keywords!=entry.function.keywords:
                raise ValueError('Confidence ablation changed bound loss parameters')
    design=build_design_settings(settings)
    reference=next(iter(Protein.from_structure(settings['targets'][0]['target_path'],chains='A').values()))
    reference_sequence=np.asarray(reference.sequence).argmax(-1)
    key=jax.random.split(jax.random.fold_in(jax.random.PRNGKey(20261021),1))[0]
    states,_,_=initialize_design_trajectory(design,key)
    # run_trajectory applies this AFTER initialization; include the real operation
    # before auditing any claim about the first search representation.
    if settings['forced_targeting']:states=force_epitope_targeting(states,settings)
    proteins=next(iter(states.values()))
    if len(proteins[design.designed_binder_chain])!=69:raise ValueError('Unexpected matched first-trajectory length')
    report=dict(status='passed',devices=[str(d) for d in jax.devices()],chains={})
    for name,protein in proteins.items():
        safe=safe_initial_atom_positions(protein.atoms,protein.atom_mask,residue_constants.atom_order)
        rotation,translation=make_transform_from_reference(*(safe[:,residue_constants.atom_order[a]] for a in ('N','CA','C')))
        np.testing.assert_allclose(np.linalg.det(np.asarray(rotation)),1,atol=1e-5)
        if not np.isfinite(np.asarray(translation)).all():raise ValueError('Nonfinite real prepared translation')
        report['chains'][name]=dict(residues=len(protein),fully_observed_backbone_residues=int(np.asarray(protein.atom_mask)[:,:3].all(-1).sum()))
        if name!=design.designed_binder_chain:
            changes=int((np.asarray(protein.sequence).argmax(-1)!=reference_sequence).sum())
            if bool(changes)!=bool(settings['forced_targeting']):raise ValueError('Prepared target representation does not match the declared native/surrogate mode')
            report['chains'][name]['sequence_differences_from_native_reference']=changes
    entries=importlib.metadata.entry_points(group='console_scripts',name='bindcraft')
    if len(entries)!=1:raise ValueError('Ambiguous CLI')
    next(iter(entries)).load()
    report['af2_sha256']=hashlib.sha256(Path(af2.__file__).read_bytes()).hexdigest()
    report['forced_targeting']=settings['forced_targeting']
    report['effective_plddt_loss_weight']=effective_weight
    report['target_scope']='Complete-target geometry; prepared sequence may be a forced-targeting surrogate. No forward, gradient or binding claim.'
    Path(output_path).write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report),flush=True)


def main(config_path,output_path):
    import hashlib
    import importlib.metadata
    import time
    import numpy as np
    import bindcraft.af2 as af2
    from bindcraft.protein import Protein
    from bindcraft.settings import read_settings,build_design_settings

    settings=read_settings(config_path)
    design=build_design_settings(settings)
    reference=next(iter(Protein.from_structure(settings['targets'][0]['target_path'],chains='A').values()))
    reference_sequence=np.asarray(reference.sequence).argmax(-1)
    binder_name=design.designed_binder_chain
    original=af2.AlphaFoldDesignModel.sequence_gradients
    calls=0
    output=Path(output_path)
    with output.open('x') as log:
        def monitored(self,*args,**kwargs):
            nonlocal calls
            predictions,gradients,loss=original(self,*args,**kwargs)
            if kwargs.get('compile_only',False):return predictions,gradients,loss
            calls+=1
            if not np.isfinite(float(loss)) or not all(np.isfinite(np.asarray(g)).all() for g in gradients.values()):
                raise ValueError('Nonfinite loss or gradient during optimization')
            prediction=next(p for p in predictions.values() if len(p.protein_complex)>1 and binder_name in p.protein_complex)
            binder=prediction.protein_complex[binder_name]
            others=[p for name,p in prediction.protein_complex.items() if name!=binder_name]
            if len(others)!=1:raise ValueError('This monitor requires the declared complete monomer target')
            target=others[0]
            ca=np.asarray(binder.atoms[:,1],dtype=float)
            if calls==1 and len(ca)!=69:raise ValueError('First trajectory length differs from matched v5 control')
            patch=np.asarray(target.atoms[:,1])[np.isin(np.asarray(target.residue_index),[205,206,208,209,210])]
            metrics=geometry_metrics(ca,patch)
            record=dict(step=calls,unix=time.time(),loss=float(loss),**metrics)
            record['target_sequence_differences_from_native']=int((np.asarray(target.sequence).argmax(-1)!=reference_sequence).sum())
            record['search_target_is_native']=record['target_sequence_differences_from_native']==0
            active_losses=args[1] if len(args)>1 else kwargs['losses']
            record['weighted_loss_components']={name:float(np.asarray(value)) for name,value in prediction.metrics.items()
                if name in active_losses and np.asarray(value).ndim==0}
            if calls==1:
                states=args[0] if args else kwargs['protein_states']
                logits=np.asarray(next(iter(states.values()))[binder_name].sequence)
                record['initial_sequence_array_sha256']=hashlib.sha256(logits.tobytes()).hexdigest()
                record['initial_sequence_array_shape']=list(logits.shape)
                record['initial_sequence_array_dtype']=str(logits.dtype)
            record['gradient_norm']=float(np.linalg.norm(np.asarray(gradients[binder_name],dtype=np.float64)))
            if not math.isfinite(record['gradient_norm']):raise ValueError('Nonfinite accumulated gradient norm')
            log.write(json.dumps(record)+'\n');log.flush()
            print('PHYSICAL_STEP '+json.dumps(record),flush=True)
            return predictions,gradients,loss
        af2.AlphaFoldDesignModel.sequence_gradients=monitored
        Path(output_path+'.provenance.json').write_text(json.dumps(dict(
            af2_import_path=af2.__file__,af2_sha256=hashlib.sha256(Path(af2.__file__).read_bytes()).hexdigest(),
            monitor_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            scope='Observation and numerical-abort wrapper; no loss/optimizer/acceptance changes'),indent=2)+'\n')
        sys.argv=['bindcraft','design',config_path]
        entries=importlib.metadata.entry_points(group='console_scripts',name='bindcraft')
        if len(entries)!=1:raise ValueError('Ambiguous installed BindCraft CLI entry point')
        next(iter(entries)).load()()


if __name__=='__main__':
    if sys.argv[1]=='--preflight':preflight(*sys.argv[2:])
    else:main(*sys.argv[1:])
