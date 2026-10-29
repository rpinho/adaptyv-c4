"""Fixed image-prepared sequences, both monomer models first; no optimization."""
import hashlib
import json
from pathlib import Path
import sys
import time
import numpy as np

MODELS=('model_1_ptm','model_2_ptm')
BB=('N','CA','C','O')
FIXED_SEQUENCES=(('v9','2699888c00c9c6af76b61c8547687cce87ccb0152754ae0d0a32c10e5f42b85d'),
                 ('seed23','8b0e2189b54820dd47b3d25e26969dd1a519f6c9a7267aee9d6f50a67fefa2b3'))

def validate_manifest(manifest):
    protocol=manifest['protocol'];expected_cap={1:8,2:4,3:8}.get(protocol)
    if expected_cap is None or manifest['models']!=list(MODELS) or manifest['recycles']!=3 or manifest['max_prediction_jobs']!=expected_cap:
        raise ValueError('Declared prediction scope differs')
    if manifest['designed_binder_coordinates_supplied'] is not False:raise ValueError('Coordinate-free declaration differs')
    entries=manifest['sequences']
    if protocol==1:
        if tuple((e['name'],e['sequence_sha256']) for e in entries)!=FIXED_SEQUENCES:
            raise ValueError('Exactly the two predeclared sequences are required')
    elif protocol==2 and (len(entries)!=1 or entries[0]['name']!='domain_primary' or entries[0]['length']!=69 or manifest.get('source_campaign')!='c4_native_domain_v16'):
        raise ValueError('Only the declared one-sequence domain follow-up is supported')
    if protocol==2:
        entry=entries[0]
        if not entry['source'].startswith('results/native_domain_v16/mirrored/campaign/') or len(entry['source_sha256'])!=64:
            raise ValueError('Domain source provenance differs')
    if protocol==3:
        if tuple(e['name'] for e in entries)!=('rfd_mpnn_1','rfd_mpnn_2') or any(e['length']!=69 for e in entries):
            raise ValueError('Exactly the two declared69-residue diffusion/MPNN samples are required')
        if manifest.get('source_campaign')!='c4_diffusion_v18' or manifest.get('mpnn_native_target_fixed_verified') is not True or manifest.get('independent_backbones')!=1:
            raise ValueError('Diffusion/MPNN provenance differs')
        if len(manifest.get('backbone_source_sha256',''))!=64 or any(not e['source'].startswith('results/diffusion_v18/mpnn/seqs/') or len(e['source_sha256'])!=64 for e in entries):
            raise ValueError('Diffusion source identity differs')
    for entry in entries:
        if len(entry['sequence'])!=entry['length'] or hashlib.sha256(entry['sequence'].encode()).hexdigest()!=entry['sequence_sha256']:
            raise ValueError('Sequence identity differs')

def dump(path,obj):
    temp=Path(str(path)+'.tmp');temp.write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n');temp.replace(path)

def geometry(points,mask):
    points=np.asarray(points,dtype=float);mask=np.asarray(mask,dtype=bool)
    if points.ndim!=3 or points.shape[1:]!=(4,3) or mask.shape!=points.shape[:2]:raise ValueError('Backbone shape differs')
    if not np.isfinite(points).all() or not mask.all():raise ValueError('Nonfinite/incomplete backbone')
    ranges={}
    for name,x,y in [('N_CA',points[:,0],points[:,1]),('CA_C',points[:,1],points[:,2]),
                     ('C_O',points[:,2],points[:,3]),('peptide_C_N',points[:-1,2],points[1:,0]),
                     ('adjacent_CA',points[:-1,1],points[1:,1])]:
        d=np.linalg.norm(x-y,axis=-1);ranges[name]=[float(d.min()),float(d.max())]
    flat=points.reshape(-1,3);ri=np.repeat(np.arange(len(points)),4)
    d=np.linalg.norm(flat[:,None]-flat[None],axis=-1)
    self_min=float(d[np.abs(ri[:,None]-ri[None,:])>1].min())
    bounds={'N_CA':(1.2,1.75),'CA_C':(1.2,1.8),'C_O':(1.05,1.5),
            'peptide_C_N':(1.1,1.7),'adjacent_CA':(2.7,4.2)}
    failed=[name for name,(low,high) in bounds.items() if ranges[name][0]<low or ranges[name][1]>high]
    if self_min<1.5:failed.append('nonadjacent_backbone_overlap')
    return dict(ranges_A=ranges,nonadjacent_backbone_min_A=self_min,passed=not failed,failed=failed)

def aligned_rmsd(x,y):
    x=np.asarray(x,dtype=float);y=np.asarray(y,dtype=float)
    x=x-x.mean(0);y=y-y.mean(0);u,_,vt=np.linalg.svd(x.T@y)
    correction=np.eye(3);correction[2,2]=np.linalg.det(u@vt)
    return float(np.sqrt(np.mean(np.sum((x@u@correction@vt-y)**2,axis=-1))))

def prepared(manifest):
    import jax
    import jax.numpy as jnp
    from bindcraft.protein import Protein,AMINO_ACIDS,ResidueFlags
    from bindcraft.settings import read_settings,build_design_settings
    from bindcraft.protein_preparation import initialize_design_trajectory
    validate_manifest(manifest)
    target_path=Path('/opt/c4/local-frame-target.pdb')
    if hashlib.sha256(target_path.read_bytes()).hexdigest()!=manifest['native_target_sha256']:raise ValueError('Native file identity changed')
    settings=read_settings('/opt/c4/campaign.json',{'forced_targeting':False,'targets':[
        dict(name='hTfR1_full_monomer',target_path=str(target_path),chains='A',hotspots='A205,A206,A208,A209,A210',weight=1.)]})
    design=build_design_settings(settings)
    states,_,_=initialize_design_trajectory(design,jax.random.PRNGKey(2026102200))
    binder=design.designed_binder_chain
    target=next(p for name,p in next(iter(states.values())).items() if name!=binder)
    if np.any(np.asarray(target.flags)&int(ResidueFlags.DESIGN)):raise ValueError('Target is designable')
    if not np.all(np.asarray(target.flags)&int(ResidueFlags.TEMPLATE)):raise ValueError('Native target template absent')
    seq=''.join(AMINO_ACIDS[i] for i in np.asarray(target.sequence).argmax(-1))
    if len(target)!=639 or hashlib.sha256(seq.encode()).hexdigest()!=manifest['native_sequence_sha256']:
        raise ValueError('Prepared native chemistry differs')
    np.testing.assert_array_equal(np.asarray(target.residue_index),np.arange(122,761))
    result=[]
    for entry in manifest['sequences']:
        sequence=entry['sequence']
        if len(sequence)!=entry['length'] or hashlib.sha256(sequence.encode()).hexdigest()!=entry['sequence_sha256']:raise ValueError('Sequence identity differs')
        atoms,mask,index=Protein.unresolved_atom_arrays(len(sequence))
        b=Protein(jax.nn.one_hot(jnp.array([AMINO_ACIDS.index(aa) for aa in sequence]),20,dtype=jnp.float16),
                  atoms,mask,jnp.full(len(sequence),int(ResidueFlags.SEQUENCE),dtype=jnp.uint8),index)
        if np.count_nonzero(b.atoms) or np.any(b.atom_mask):raise ValueError('Binder coordinates exposed')
        np.testing.assert_array_equal(np.asarray(b.flags),int(ResidueFlags.SEQUENCE))
        np.testing.assert_array_equal(np.asarray(b.sequence).sum(-1),1)
        np.testing.assert_array_equal(np.asarray(b.sequence).max(-1),1)
        complex_states={name:{key:(b if key==binder else value) for key,value in proteins.items()} for name,proteins in states.items()}
        result.append((entry,b,complex_states))
    return design,binder,result

def preflight(manifest_path,out):
    from bindcraft.MPNN_stage import predict_validation_ensemble
    manifest=json.loads(Path(manifest_path).read_text());design,binder,items=prepared(manifest)
    report=dict(status='passed',sequence_lengths=[len(b) for _,b,_ in items],
        binder_atoms_all_zero=True,binder_masks_all_false=True,binder_flags='SEQUENCE only',
        hard_one_hot_sequences=True,target_residues=639,target_native_identity_verified=True,
        prediction_models=list(MODELS),recycles=3,designed_coordinates_supplied=False,
        checkpoint_scope=manifest.get('checkpoint_scope','Held-out AF2 monomer family from original multimer design; same architecture lineage'),
        interface=predict_validation_ensemble.__module__+'.'+predict_validation_ensemble.__name__,
        scope=manifest['complex_scope'],monomer_output_chain='A',complex_output_chains='A target/B binder',
        worker_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        manifest_sha256=hashlib.sha256(Path(manifest_path).read_bytes()).hexdigest())
    dump(out,report);print(json.dumps(report),flush=True)

def main(out,manifest_path='/opt/c4/hard_refold_manifest.json'):
    import jax
    from bindcraft.af2 import AlphaFoldDesignModel
    from bindcraft.MPNN_stage import predict_validation_ensemble
    from bindcraft.model_weights import alphafold_parameters
    from bindcraft.protein import BINDER_ALONE,ATOM_INDEX,write_structure
    if any(d.platform!='gpu' for d in jax.devices()):raise ValueError('Refusing CPU inference fallback')
    out=Path(out);manifest=json.loads(Path(manifest_path).read_text())
    design,binder,items=prepared(manifest)
    model=AlphaFoldDesignModel(presets=MODELS,models=MODELS,data_dir=alphafold_parameters(download=False),
        num_recycle=3,dropout=False,large_complex_initialization=False,use_remat=True,max_cache_size=4)
    report=dict(status='running',started=time.time(),prediction_jobs=0,prediction_jobs_scheduled=0,models=list(MODELS),
        network_evaluations_per_prediction=4,scope=manifest['complex_scope'],sequences=[],
        qualification='Diagnostic only; no automatic generator/native admission or pH scoring',
        max_prediction_jobs=manifest['max_prediction_jobs'],checkpoint_scope=manifest.get('checkpoint_scope','Held-out AF2 monomer family, same architecture lineage'))
    dump(out/'refold_results.json',report)
    selected=[]
    def predict(states):
        if report['prediction_jobs_scheduled']+len(MODELS)>manifest['max_prediction_jobs']:raise ValueError('Prediction job cap would be exceeded')
        report['prediction_jobs_scheduled']+=len(MODELS);dump(out/'refold_results.json',report)
        return predict_validation_ensemble(model,states,MODELS)
    def observe(entry,phase,predictions):
        records=[];ca=[]
        for checkpoint,states in predictions.items():
            prediction=next(iter(states.values()));protein=prediction.protein_complex[binder]
            if sorted(prediction.protein_complex)[0]!=binder:raise ValueError('Unexpected metric chain order')
            conf=np.asarray(prediction.metrics['plddt'],dtype=float)[:len(protein)]
            if not np.isfinite(conf).all():raise ValueError('Nonfinite confidence')
            indices=[ATOM_INDEX[a] for a in BB]
            xyz=np.asarray(protein.atoms)[:,indices];mask=np.asarray(protein.atom_mask)[:,indices]
            physical=geometry(xyz,mask);ca.append(xyz[:,1])
            record=dict(checkpoint=checkpoint,binder_pLDDT=float(conf.mean()),geometry=physical)
            record['pTM']=float(prediction.metrics['ptm'])
            name=entry['name']+'_'+phase+'_'+checkpoint
            if phase=='monomer':
                write_structure({'A':protein},str(out/(name+'.cif')),plddt=conf)
            else:
                target=next(p for n,p in prediction.protein_complex.items() if n!=binder)
                all_conf=np.asarray(prediction.metrics['plddt'])
                write_structure({'A':target,'B':protein},str(out/(name+'.cif')),
                    plddt=np.concatenate([all_conf[len(protein):],all_conf[:len(protein)]]))
                record.update(i_pTM=float(prediction.metrics['iptm']),pTM=float(prediction.metrics['ptm']))
                patch=np.asarray(target.atoms)[np.isin(np.asarray(target.residue_index),[205,206,208,209,210]),ATOM_INDEX['CA']]
                if len(patch)!=5:raise ValueError('Hotspot mapping differs')
                record['minimum_hotspot_to_binder_CA_A']=float(np.linalg.norm(patch[:,None]-xyz[None,:,1],axis=-1).min())
            records.append(record);report['prediction_jobs']+=1
            dump(out/(name+'.json'),record)
            print('REFOLD_JOB '+json.dumps(dict(sequence=entry['name'],phase=phase,checkpoint=checkpoint,
                binder_pLDDT=record['binder_pLDDT'],geometry_passed=physical['passed'])),flush=True)
        agreement=aligned_rmsd(ca[0],ca[1])
        result=dict(models=records,mutual_internal_fold_CA_rmsd_A=agreement)
        if phase=='monomer':result['monomer_advance']=all(r['binder_pLDDT']>=.7 and r['geometry']['passed'] for r in records) and agreement<=3.5
        if phase=='complex':result['target_template_frame_binder_CA_rmsd_A']=float(np.sqrt(np.mean(np.sum((ca[0]-ca[1])**2,axis=-1))))
        return result
    # All four monomer jobs precede any large complex compilation/allocation.
    for index,(entry,b,complex_states) in enumerate(items):
        model.key=jax.random.PRNGKey(2026102200+index)
        ensemble=predict({BINDER_ALONE:{binder:b}})
        result=dict(name=entry['name'],sequence_sha256=entry['sequence_sha256'],monomer=observe(entry,'monomer',ensemble.model_predictions))
        report['sequences'].append(result)
        if result['monomer']['monomer_advance']:selected.append((index,entry,complex_states,result))
        dump(out/'refold_results.json',report)
    for index,entry,complex_states,result in selected:
        model.key=jax.random.PRNGKey(2026102300+index)
        ensemble=predict(complex_states)
        result['complex']=observe(entry,'complex',ensemble.model_predictions)
        dump(out/'refold_results.json',report)
    if report['prediction_jobs']>manifest['max_prediction_jobs']:raise ValueError('Prediction job cap exceeded')
    report.update(status='completed',finished=time.time(),plausible_monomer_sequences=len(selected))
    dump(out/'refold_results.json',report)

if __name__=='__main__':
    if sys.argv[1]=='--preflight':preflight(*sys.argv[2:])
    else:main(*sys.argv[1:])
