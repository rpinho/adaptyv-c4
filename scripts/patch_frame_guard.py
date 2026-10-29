"""Apply a version-locked, initialization-only repair to pinned BindCraft2."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil

PINNED_AF2_SHA='aa185971676c5b6e3d2ebcef920e3b78acf16879b48b9073e74a527401ef79ae'


def apply(root,helper):
    path=Path(root)/'bindcraft/af2.py'
    original=path.read_bytes()
    if hashlib.sha256(original).hexdigest()!=PINNED_AF2_SHA:
        raise ValueError('Unexpected AF2 source revision or already patched source; refusing mutation')
    text=original.decode()
    replacements={
        'import jax.numpy as jnp': 'import jax.numpy as jnp\nfrom bindcraft.frame_guard import safe_initial_atom_positions',
        "'initial_atom_pos': atoms.astype(jnp.float32)": "'initial_atom_pos': safe_initial_atom_positions(atoms, atom_mask, residue_constants.atom_order)",
        "'initial_atom_pos': previous_state['prev_pos']": "'initial_atom_pos': safe_initial_atom_positions(previous_state['prev_pos'], model_inputs['atom37_atom_exists'], residue_constants.atom_order)",
    }
    for old,new in replacements.items():
        if text.count(old)!=1:raise ValueError(f'Expected exactly one source anchor: {old}')
        text=text.replace(old,new)
    compile(text,str(path),'exec')
    path.write_text(text)
    destination=path.parent/'frame_guard.py';shutil.copyfile(helper,destination)
    report=dict(original_af2_sha256=PINNED_AF2_SHA,patched_af2_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                helper_sha256=hashlib.sha256(destination.read_bytes()).hexdigest(),
                repaired_paths=['initial input','recycled previous prediction'],
                actual_atom_positions_masks_templates_modified=False)
    (Path(root)/'frame_guard_manifest.json').write_text(json.dumps(report,indent=2)+'\n')
    return report


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--upstream',required=True);p.add_argument('--helper',required=True)
    args=p.parse_args();print(json.dumps(apply(args.upstream,args.helper),indent=2))
