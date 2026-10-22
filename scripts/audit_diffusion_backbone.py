"""Post-release physical/context staging; never qualify a poly-Gly backbone."""
from pathlib import Path
import json
import numpy as np
from Bio.PDB import PDBParser
from hard_sequence_refold import geometry
from restore_apical_context import restore,sha,coordinates

ROOT=Path(__file__).resolve().parents[1]
STATE=ROOT/'results/diffusion_v18'

def main():
    supervisor=json.loads((STATE/'supervisor.json').read_text())
    if not supervisor.get('gpu_released_verified') or supervisor.get('result',{}).get('status')!='completed':
        raise RuntimeError('Completed generation and actual GPU release required')
    report=json.loads((STATE/'mirrored/generation_results.json').read_text())
    if report['backbone_count']!=1 or Path(report['backbone_file']).name!=report['backbone_file']:
        raise ValueError('Unexpected generation identity')
    source=STATE/'mirrored'/report['backbone_file']
    if sha(source)!=report['backbone_sha256']:raise ValueError('Backbone hash differs')
    out=dict(status='structural_rejection',source=str(source.relative_to(ROOT)),source_sha256=sha(source),
             qualified_candidates=0,sequence_design_allowed=False)
    try:
        model=PDBParser(QUIET=True,PERMISSIVE=False).get_structure('diffusion',source)
        if len(model)!=1 or set(model[0].child_dict)!=set('AB'):raise ValueError('Need exactly native A and generated B')
        binder=list(model[0]['B'])
        if [r.id for r in binder]!=[(' ',i,' ') for i in range(1,70)] or any(r.resname!='GLY' for r in binder):
            raise ValueError('Expected exactly69 generated poly-Gly residues')
        if any(a.is_disordered() for a in model.get_atoms()):raise ValueError('Alternate atom identity')
        points=np.array([[r[a].coord for a in ('N','CA','C','O')] for r in binder])
        check=geometry(points,np.ones(points.shape[:2],dtype=bool));out['binder_geometry']=check
        if not check['passed']:raise ValueError('Generated binder geometry failed: '+str(check['failed']))
        target=model[0]['A']
        patch=np.array([target[i]['CA'].coord for i in (205,206,208,209,210)])
        distances=np.linalg.norm(patch[:,None]-points[None,:,1],axis=-1).min(axis=1)
        out['hotspot_nearest_binder_CA_A']=distances.tolist()
        out['hotspot_CA_fraction_within8A']=float(np.mean(distances<=8))
        if out['hotspot_CA_fraction_within8A']<.5:raise ValueError('Insufficient declared-patch backbone contact')
        out['restoration']=restore(source,STATE/'restored_staging_only.pdb')
        out.update(status='native_context_backbone_staging_passed',sequence_design_allowed=True,
            interpretation='Physically staged backbone only. Generate two fixed MPNN samples, then coordinate-free native validation; no binding/pH/mouse evidence.')
    except (ValueError,KeyError) as exc:
        out['reason']=str(exc)
    (STATE/'BACKBONE_AUDIT.json').write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],flush=True)
    return out

if __name__=='__main__':main()
