"""Initialization-only canonical frames for unobserved/degenerate AF2 residues."""
import jax.numpy as jnp


def safe_initial_atom_positions(atoms, atom_mask, atom_order):
    positions = atoms.astype(jnp.float32)
    n, ca, c = (atom_order[name] for name in ('N','CA','C'))
    nitrogen, alpha, carbon = positions[:,n], positions[:,ca], positions[:,c]
    first, second = carbon-alpha, nitrogen-alpha
    present = (atom_mask[:,n]>0)&(atom_mask[:,ca]>0)&(atom_mask[:,c]>0)
    finite = jnp.isfinite(jnp.stack((nitrogen,alpha,carbon),axis=1)).all(axis=(1,2))
    valid = present & finite & (jnp.linalg.norm(first,axis=-1)>1e-4) & (jnp.linalg.norm(second,axis=-1)>1e-4) & (jnp.linalg.norm(jnp.cross(first,second),axis=-1)>1e-4)
    origin = jnp.where(((atom_mask[:,ca]>0)&jnp.isfinite(alpha).all(-1))[:,None],alpha,0.0)
    canonical = jnp.broadcast_to(origin[:,None,:],positions.shape)
    canonical = canonical.at[:,n].set(origin+jnp.array([0.,1.,0.]))
    canonical = canonical.at[:,c].set(origin+jnp.array([1.,0.,0.]))
    return jnp.where(valid[:,None,None],positions,canonical)
