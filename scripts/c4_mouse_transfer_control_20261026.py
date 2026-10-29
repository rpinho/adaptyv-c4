#!/usr/bin/env python3
"""Matched human-reference transfer control for the fixed-pose mouse audit."""
import argparse
import importlib.util
import json
import subprocess
from pathlib import Path

import numpy as np
from Bio.PDB import MMCIFParser, PDBParser

ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('mouse_audit',ROOT/'scripts/c4_mouse_compatibility_20261026.py')
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)


def run(out):
    out=out.resolve()
    stamp=subprocess.check_output(['date','-u','+%Y-%m-%dT%H:%M:%SZ'],text=True).strip()
    rawpath=out/'measurements.json'
    raw=json.loads(rawpath.read_text())
    hp=ROOT/'data/11qc.cif'
    hc=MMCIFParser(QUIET=True).get_structure('human_reference',hp)[0]['A']
    human=''
    from Bio import SwissProt
    with (ROOT/'data/P02786.txt').open() as f:
        human=SwissProt.read(f).sequence
    assert ''.join(m.AA[hc[h].resname] for h in m.DOMAIN)==human[188:383]
    target=m.atoms([hc[h] for h in m.DOMAIN])
    parser=PDBParser(QUIET=True)
    rows=[]
    for d in raw['designs']:
        p=ROOT/d['source_pdb']
        assert m.sha(p)==raw['inputs'][d['source_pdb']]
        model=parser.get_structure('design',p)[0]
        assert m.seqsha(m.sequence(model['A']))==d['sequence_sha256']
        binder=m.atoms(m.residues(model['A']))
        hr={r.id[1]+188:r for r in m.residues(model['B'])}
        fits={}
        for label,mf in d['fits'].items():
            hs=mf['fit_human_positions']
            x=np.array([hc[h]['CA'].coord for h in hs],dtype=float)
            y=np.array([hr[h]['CA'].coord for h in hs],dtype=float)
            rot,shift,rms=m.fit(x,y)
            transferred=[(n,a,xyz@rot+shift) for n,a,xyz in target]
            metrics=m.distances(binder,transferred)
            # pLDDT on the mouse AF residues in the conflict and contact region.
            fits[label]=dict(fit_CA_rmsd_A=rms,**metrics)
        bothmouse=all(f['target_backbone_pairs_lt2A']>0 for f in d['fits'].values())
        bothhuman=all(f['target_backbone_pairs_lt2A']>0 for f in fits.values())
        human_any=any(f['target_backbone_pairs_lt2A']>0 for f in fits.values())
        label=('mouse_conflict_with_human_transfer_confounded' if bothmouse and human_any else
               'mouse_conflict_human_reference_clear_both_fits' if bothmouse else
               'mouse_fit_sensitive_conflict' if any(f['target_backbone_pairs_lt2A']>0 for f in d['fits'].values()) else
               'mouse_backbone_clear_both_fits')
        rows.append(dict(sequence_sha256=d['sequence_sha256'],name=d['name'],
                         human_reference_fits=fits,mouse_conflict_both_fits=bothmouse,
                         human_reference_conflict_both_fits=bothhuman,
                         human_reference_conflict_any_fit=human_any,interpretation_label=label))
    record=dict(created_at_utc=stamp,inputs={str(rawpath.relative_to(ROOT)):m.sha(rawpath),
                 'data/11qc.cif':m.sha(hp),str(Path(__file__).relative_to(ROOT)):m.sha(Path(__file__))},
        method='Same 192 domain or local sequence-paired CA positions as mouse audit, transfer experimental human apical A189-383 onto each own predicted human target. No partner, glycan, or other-domain atoms used. A transfer control only, not full assay context.',
        interpretation='Human-reference clashes expose transfer/conformation sensitivity; even a mouse-only conflict is a fixed-pose observation, not exclusion of cross-reactivity.',
        counts={k:sum(d['interpretation_label']==k for d in rows) for k in sorted({d['interpretation_label'] for d in rows})},designs=rows)
    output=out/'human_transfer_control.json'
    assert not output.exists()
    output.write_text(json.dumps(record,indent=2)+'\n')
    print(json.dumps(record['counts'],indent=2))


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--out',type=Path,required=True)
    run(ap.parse_args().out)
