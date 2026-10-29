"""CPU-only reduced-domain staging. Outputs require native full-target refolding.

This module is deliberately outside the ranked-design/pH admission path.
It never promotes a domain prediction to a native-complex prediction.
"""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from Bio.Data.PDBData import protein_letters_3to1 as AA
from Bio.PDB import MMCIFParser, PDBParser, PDBIO, Select

ROOT=Path(__file__).resolve().parents[1]
NATIVE_TEMPLATE=ROOT/'data/11qc.cif'
NATIVE_SHA256='4eaa5db67ed83a896713c4b42d1dc214c778bdac240fdf3a9a82e86a66f79495'
DOMAIN=tuple(range(189,384))
PATCH=(205,206,208,209,210)
BACKBONE=('N','CA','C','O')


def protein_residues(chain):
    return [r for r in chain if r.id[0]==' ' and r.resname in AA]


def native_model():
    if sha(NATIVE_TEMPLATE)!=NATIVE_SHA256:raise ValueError('Native template identity/hash changed')
    model=MMCIFParser(QUIET=True).get_structure('11qc',NATIVE_TEMPLATE)[0]
    for cid in 'AB':
        rr=protein_residues(model[cid])
        if [r.id[1] for r in rr]!=list(range(122,761)):
            raise ValueError('Native template identity/coverage changed')
    return model


def coordinates(rr,backbone=False):
    return np.array([a.coord for r in rr for a in r
                     if (a.name in BACKBONE if backbone else a.element!='H')],dtype=float)


def minimum_distance(x,y):
    if not len(x) or not len(y):raise ValueError('Missing clearance coordinates')
    return min(float(np.linalg.norm(chunk[:,None]-y[None],axis=-1).min())
               for chunk in np.array_split(x,max(1,len(x)//100)))


def audit_domain(chain,reference):
    if any(r.id[0]!=' ' or r.resname not in AA for r in chain):
        raise ValueError('Unknown/nonprotein domain residue must not be silently omitted')
    rr=protein_residues(chain)
    if [r.id[1] for r in rr]!=list(DOMAIN) or any(r.id[2]!=' ' for r in rr):
        raise ValueError('Need exactly contiguous native author residues189-383')
    if [r.resname for r in rr]!=[reference[n].resname for n in DOMAIN]:
        raise ValueError('Reduced target sequence differs from native')
    for r in rr:
        if any(name not in r for name in BACKBONE):raise ValueError('Incomplete domain backbone')
    if not np.isfinite(coordinates(rr)).all():raise ValueError('Nonfinite domain coordinates')
    bonds=[float(np.linalg.norm(a['C'].coord-b['N'].coord)) for a,b in zip(rr,rr[1:])]
    if min(bonds)<1 or max(bonds)>1.9:raise ValueError('Broken domain peptide geometry')
    return rr,dict(residues=len(rr),author_range=[189,383],native_sequence_verified=True,
                   minimum_peptide_CN_A=min(bonds),maximum_peptide_CN_A=max(bonds))


class ProteinOnly(Select):
    def accept_residue(self,r):return r.id[0]==' ' and r.resname in AA
    def accept_atom(self,a):return a.element!='H' and a.altloc in (' ','A')


def write(model,path):
    if any(not np.isfinite(a.coord).all() or min(a.coord)<-999.999 or max(a.coord)>9999.999
           for a in model.get_atoms()):raise ValueError('Coordinates exceed finite fixed-column PDB bounds')
    writer=PDBIO();writer.set_structure(model);writer.save(str(path),ProteinOnly())


def prepare(source,destination):
    model=PDBParser(QUIET=True).get_structure('source',str(source))[0]
    ref=native_model()
    source_aa=protein_residues(model['A'])
    if len(model['A'])!=639 or [r.id[1] for r in source_aa]!=list(range(122,761)) or [r.resname for r in source_aa]!=[r.resname for r in protein_residues(ref['A'])]:
        raise ValueError('Preparation requires full native639-residue author-numbered A')
    for chain in list(model):
        if chain.id!='A':model.detach_child(chain.id)
    for r in list(model['A']):
        if r.id[0]!=' ' or r.id[1] not in DOMAIN:model['A'].detach_child(r.id)
    _,audit=audit_domain(model['A'],ref['A'])
    write(model,destination)
    return audit|dict(stage='reduced_search_input',source_sha256=sha(source),
                      output_sha256=sha(destination),requires_native_full_context_refold=True)


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def restore(source,destination):
    model=PDBParser(QUIET=True).get_structure('domain_prediction',str(source))[0]
    if set(c.id for c in model)!=set('AB'):raise ValueError('Expected reduced target A and binder B')
    native=native_model()
    rr,audit=audit_domain(model['A'],native['A'])
    if any(r.id[0]!=' ' or r.resname not in AA for r in model['B']):
        raise ValueError('Unknown/nonprotein binder residue must not be silently omitted')
    binder=protein_residues(model['B'])
    if not 30<=len(binder)<=150:raise ValueError('Unexpected binder architecture')
    if any(any(name not in r for name in BACKBONE) for r in binder):raise ValueError('Incomplete binder backbone')
    if not np.isfinite(coordinates(binder)).all():raise ValueError('Nonfinite binder coordinates')
    peptide=[float(np.linalg.norm(a['C'].coord-b['N'].coord)) for a,b in zip(binder,binder[1:])]
    if min(peptide)<1 or max(peptide)>1.9:raise ValueError('Broken binder peptide geometry')
    # Fit ALL195 domain CAs, never only the favorable hotspot subset.
    x=np.array([r['CA'].coord for r in rr],dtype=float)
    y=np.array([native['A'][n]['CA'].coord for n in DOMAIN],dtype=float)
    xc,yc=x.mean(0),y.mean(0)
    u,_,vt=np.linalg.svd((x-xc).T@(y-yc))
    correction=np.eye(3);correction[2,2]=np.linalg.det(u@vt)
    rotation=u@correction@vt;translation=yc-xc@rotation
    deviations=np.linalg.norm(x@rotation+translation-y,axis=-1)
    rmsd=float(np.sqrt(np.mean(deviations**2)))
    patch_rmsd=float(np.sqrt(np.mean(deviations[[DOMAIN.index(n) for n in PATCH]]**2)))
    if rmsd>2 or patch_rmsd>2 or deviations.max()>5:
        raise ValueError('Reduced domain deformed relative to full native context')
    restored_binder=model['B'].copy();restored_binder.id='C'
    for atom in restored_binder.get_atoms():atom.coord=atom.coord@rotation+translation
    binder=protein_residues(restored_binder)
    target=[r for cid in 'AB' for r in protein_residues(native[cid])]
    clearance=minimum_distance(coordinates(binder,True),coordinates(target,True))
    if clearance<2:raise ValueError(f'Restored full-target backbone overlap:{clearance:.3f}A')
    glycans=np.array([a.coord for c in native for r in c if r.resname=='NAG' for a in r],dtype=float)
    glycan_clearance=minimum_distance(coordinates(binder),glycans)
    if glycan_clearance<2:raise ValueError('Restored binder overlaps observed NAG atoms')
    output=native.copy()
    for c in list(output):
        if c.id not in 'AB':output.detach_child(c.id)
    output.add(restored_binder)
    write(output,destination)
    return audit|dict(stage='restored_staging_only',source_sha256=sha(source),
          native_template_sha256=sha(NATIVE_TEMPLATE),output_sha256=sha(destination),
          alignment_residues=195,domain_CA_rmsd_A=rmsd,patch_CA_rmsd_A=patch_rmsd,
          maximum_domain_CA_displacement_A=float(deviations.max()),
          minimum_binder_full_target_backbone_distance_A=clearance,
          observed_NAG_min_distance_A=glycan_clearance,restored_target_residues=1278,
          rotation_determinant=float(np.linalg.det(rotation)),qualified_candidate=False,
          requires_native_full_context_refold=True,
          limitations=['Rigid context restoration is not native-complex prediction or binding validation',
                       'No pH scoring allowed before native refolding and existing admission gates',
                       'Observed NAG clearance checked; glycans are not retained in output PDB',
                       'Human only; residue121/unresolved glycans/assay tag absent'])


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('operation',choices=['prepare','restore'])
    parser.add_argument('--source',required=True);parser.add_argument('--out',required=True)
    parser.add_argument('--audit',required=True);args=parser.parse_args()
    for path in (args.out,args.audit):
        if Path(path).exists():raise FileExistsError('Use fresh staging output paths')
    result=(prepare if args.operation=='prepare' else restore)(args.source,args.out)
    Path(args.audit).write_text(json.dumps(result,indent=2)+'\n')
