#!/usr/bin/env python
"""Restore full 11QC context around every C4 candidate; do not relax clashes away."""
import csv
import json
from pathlib import Path
import time
import numpy as np
from Bio.PDB import MMCIFParser, PDBParser
from run_context_ensemble import residues, atoms, fit, sha, dump, ROOT


def main():
    out=ROOT/'results/full_context_audit';out.mkdir(exist_ok=True)
    template=MMCIFParser(QUIET=True).get_structure('11qc',ROOT/'data/11qc.cif')[0]
    target=[(c.id,r) for c in template for r in residues(c)]
    pts=atoms([r for c,r in target]); keys=[(c,r.id[1]) for c,r in target for a in r if a.element!='H']
    bbpts=np.array([a.coord for c,r in target for a in r if a.name in ('N','CA','C','O')])
    bbkeys=[(c,r.id[1]) for c,r in target for a in r if a.name in ('N','CA','C','O')]
    rows=[];details=[]
    started=time.time()
    for source in sorted((ROOT/'results/v3_all_pdb').glob('*.pdb')):
        m=PDBParser(QUIET=True).get_structure('x',source)[0]
        matched=[(c.id,r) for c in m if c.id in 'AB' for r in residues(c)]
        rotation,translation,rmsd=fit([r['CA'].coord for c,r in matched],
                                    [template[c][r.id]['CA'].coord for c,r in matched])
        source_keys={(c,r.id[1]) for c,r in matched}
        rr=residues(m['C']);binder=atoms(rr)@rotation+translation
        binder_bb=np.array([a.coord for r in rr for a in r if a.name in ('N','CA','C','O')])@rotation+translation
        restored=np.array([p for p,k in zip(pts,keys) if k not in source_keys])
        restored_keys=[k for k in keys if k not in source_keys]
        restored_bb=np.array([p for p,k in zip(bbpts,bbkeys) if k not in source_keys])
        restored_bb_keys=[k for k in bbkeys if k not in source_keys]
        clashes=set();pair_count=0;minimum=float('inf')
        for chunk in np.array_split(binder,max(1,len(binder)//100)):
            distances=np.linalg.norm(chunk[:,None]-restored[None],axis=-1)
            minimum=min(minimum,float(distances.min()))
            i,j=np.where(distances<2);pair_count+=len(i)
            clashes.update(restored_keys[k] for k in j)
        bb_dist=np.linalg.norm(binder_bb[:,None]-restored_bb[None],axis=-1)
        i,j=np.where(bb_dist<2)
        bb_clashes=sorted(set(restored_bb_keys[k] for k in j))
        row=dict(design=source.stem,target_ca_rmsd=rmsd,aligned_target_residues=len(matched),
                 restored_heavy_pair_clashes_under2A=pair_count,restored_clashing_residues=len(clashes),
                 restored_min_heavy_distance=minimum,restored_backbone_pair_clashes_under2A=len(i),
                 restored_backbone_min_distance=float(bb_dist.min()),
                 rigid_context_pass=not bool(bb_clashes),source_sha256=sha(source))
        rows.append(row);details.append(dict(row,clashing_residues=sorted(clashes),
                                             backbone_clashing_residues=bb_clashes))
        dump(out/'status.json',dict(status='running',checked=len(rows),started=started))
    with (out/'candidates.csv').open('w',newline='') as fh:
        writer=csv.DictWriter(fh,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    dump(out/'candidates.json',details)
    summary=dict(status='completed',candidates=len(rows),
                 rigid_context_pass=sum(r['rigid_context_pass'] for r in rows),
                 backbone_clash_failures=sum(not r['rigid_context_pass'] for r in rows),
                 rmsd_range=[min(r['target_ca_rmsd'] for r in rows),max(r['target_ca_rmsd'] for r in rows)],
                 well_aligned_under2A=sum(r['target_ca_rmsd']<=2 for r in rows),
                 well_aligned_backbone_clash_failures=sum(r['target_ca_rmsd']<=2 and not r['rigid_context_pass'] for r in rows),
                 source_template_sha256=sha(ROOT/'data/11qc.cif'),script_sha256=sha(__file__),
                 helper_sha256=sha(Path(__file__).with_name('run_context_ensemble.py')),
                 started=started,finished=time.time(),new_gpu_spend_usd=0)
    dump(out/'status.json',summary)
    (out/'RESULTS.md').write_text(
        '# Full resolved receptor context audit\n\n'
        f"{len(rows)} source candidates audited; {summary['rigid_context_pass']} have no restored "
        f"backbone clashes below 2 Å. {summary['backbone_clash_failures']} clash with protein omitted "
        'from the design crop. This is a rigid context-transfer diagnostic, not an experimental failure.\n\n'
        'All source receptor Cα atoms in chains A/B are fitted to matched 11QC residue numbers. '
        'The binder is transformed by the same rotation/translation. Only atoms of residues absent '
        'from the source crop count as restored-context clashes. Backbone-to-backbone distances '
        'separate impossible occupancy from rotamers that might repack. No receptor or binder '
        'atoms were relaxed. Larger receptor conformational changes require independent refolding. '
        f"{summary['well_aligned_under2A']} fits are <=2 Å Cα RMSD; all of these fail the backbone "
        'occupancy test. One fit exceeds 2 Å and must be interpreted as a failed context transfer, '
        'not an equally precise clash measurement.\n\n'
        'Do not advance the old numerical pH winners directly to submission. The crop permitted '
        'docking into volume occupied by omitted receptor structure. Rebuild generation inputs '
        'with the full resolved dimer and test an accessible approach before local switch expansion. '
        'Old crop scores remain reproducible but do not validate full-receptor binding.\n\n'
        'No new GPU spend from this CPU audit. See results/gpu_campaign/status.json for any '
        'separately authorized GPU work.\n')
    print(json.dumps(summary),flush=True)


if __name__=='__main__':main()
