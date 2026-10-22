"""CPU numerical regression using BOTH actual pinned AF2 frame implementations."""
import argparse
import ast
import hashlib
import importlib
import json
from pathlib import Path
import sys


def validate(root):
    root=Path(root).resolve();sys.path.insert(0,str(root))
    import jax
    import jax.numpy as jnp
    import numpy as np
    from bindcraft.af.alphafold.model import geometry,quat_affine
    from bindcraft.frame_guard import safe_initial_atom_positions
    source=root/'bindcraft/af/alphafold/model/all_atom_multimer.py'
    # Execute the exact pinned upstream function, with its actual geometry classes.
    node=next(n for n in ast.parse(source.read_text()).body if isinstance(n,ast.FunctionDef) and n.name=='make_transform_from_reference')
    namespace={'geometry':geometry}
    exec(compile(ast.Module(body=[node],type_ignores=[]),str(source),'exec'),namespace)
    convert=namespace['make_transform_from_reference']
    order={'N':0,'CA':1,'C':2}
    atoms=np.zeros((6,37,3),np.float32);masks=np.zeros((6,37),np.float32)
    # unresolved; valid observed target; partial CA; collinear; nonfinite; padding.
    atoms[1,:3]=[[3,5,5],[3,4,5],[4,4,5]];masks[1,:3]=1
    atoms[2,1]=[7,8,9];masks[2,1]=1
    atoms[3,:3]=[[0,0,0],[1,0,0],[2,0,0]];masks[3,:3]=1
    atoms[4,0]=[np.nan,0,0];masks[4,:3]=1
    original=atoms.copy();original_masks=masks.copy()
    safe=np.asarray(safe_initial_atom_positions(jnp.asarray(atoms),jnp.asarray(masks),order))
    np.testing.assert_equal(atoms,original);np.testing.assert_equal(masks,original_masks)
    np.testing.assert_array_equal(safe[1],original[1])
    np.testing.assert_array_equal(safe[2,1],[7,8,9])
    assert np.isfinite(safe).all()
    def multi(x):
        return convert(*(geometry.Vec3Array.from_array(jnp.asarray(x[:,i])) for i in (0,1,2)))
    rigid=multi(safe)
    multimer=np.asarray(rigid.rotation.to_array())
    rot,translation=quat_affine.make_transform_from_reference(*(jnp.asarray(safe[:,i]) for i in (0,1,2)))
    monomer=np.asarray(rot)
    results={}
    for name,matrices in [('multimer',multimer),('monomer',monomer)]:
        np.testing.assert_allclose(matrices.transpose(0,2,1)@matrices,np.broadcast_to(np.eye(3),(6,3,3)),atol=1e-5)
        np.testing.assert_allclose(np.linalg.det(matrices),1,atol=1e-5)
        results[name]=dict(determinants=np.linalg.det(matrices).tolist(),finite=bool(np.isfinite(matrices).all()))
    original_valid=multi(original[1:2]).rotation.to_array()
    np.testing.assert_array_equal(multimer[1:2],np.asarray(original_valid))
    update=geometry.Rigid3Array(geometry.Rot3Array.identity((6,)),geometry.Vec3Array.from_array(jnp.ones((6,3))))
    moved=np.asarray((rigid@update).translation.to_array())-np.asarray(rigid.translation.to_array())
    np.testing.assert_allclose(np.linalg.norm(moved,axis=-1),np.sqrt(3),atol=1e-5)
    bad=np.asarray(multi(original[0:1]).rotation.to_array())
    assert np.linalg.det(bad)[0]==0
    bad_move=np.asarray((multi(original[0:1])@geometry.Rigid3Array(geometry.Rot3Array.identity((1,)),geometry.Vec3Array.from_array(jnp.ones((1,3))))).translation.to_array())
    assert np.all(bad_move==0)
    # A finite recycled prediction is retained when its sequence-derived mask exists.
    predicted=safe.copy();sequence_mask=np.ones_like(masks)
    np.testing.assert_array_equal(np.asarray(safe_initial_atom_positions(jnp.asarray(predicted),jnp.asarray(sequence_mask),order)),predicted)
    text=(root/'bindcraft/af2.py').read_text()
    assert text.count('safe_initial_atom_positions(')==2
    assert "safe_initial_atom_positions(previous_state['prev_pos'], model_inputs['atom37_atom_exists']" in text
    return dict(status='passed',jax=jax.__version__,devices=[str(d) for d in jax.devices()],
                cases=['unresolved','valid template','partial CA','collinear','nonfinite','padding','recycled valid prediction'],
                actual_upstream_frame_functions=results,old_zero_frame_determinant=float(np.linalg.det(bad)[0]),
                old_zero_frame_update_norm=float(np.linalg.norm(bad_move)),
                repaired_translation_update_norms=np.linalg.norm(moved,axis=-1).tolist(),
                actual_inputs_masks_unchanged=True,valid_template_frame_unchanged=True,
                patch=json.loads((root/'frame_guard_manifest.json').read_text()),
                source_sha256={str(p.relative_to(root)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [source,root/'bindcraft/af/alphafold/model/quat_affine.py',root/'bindcraft/af/alphafold/model/geometry/rotation_matrix.py']},
                limitation='Numerical frame regression, not an AlphaFold forward/gradient or affinity validation')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--upstream',required=True);p.add_argument('--out',required=True)
    args=p.parse_args();result=validate(args.upstream)
    dest=Path(args.out);dest.parent.mkdir(parents=True,exist_ok=True);dest.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
