#!/usr/bin/env python
"""Write measured isolated-binder scout results and audit modeled coordinates."""
import collections
import csv
import json
from pathlib import Path
import statistics
from Bio.PDB import PDBParser
import numpy as np
from run_context_ensemble import ROOT, sha, dump


def main():
    out=ROOT/'results/free_networks_v3';data=json.loads((out/'status.json').read_text())
    queue=json.loads((out/'queue.json').read_text());groups=collections.defaultdict(list)
    parser=PDBParser(QUIET=True)
    reference=parser.get_structure('ref',queue['parent']['variant_path'])[0]['C']
    checks=[];maximum=0
    for row in data['models']:
        assert sha(row['pdb'])==row['sha256']
        model=parser.get_structure('new',row['pdb'])[0]['C']
        changes=dict(next(m for label,m in queue['variants'] if label==row['label']))
        for index,r in enumerate(reference,1):
            actual=model[r.id]
            expected=changes.get(index)
            if expected:
                from Bio.Data.PDBData import protein_letters_3to1 as AA
                assert AA[actual.resname]==expected
            else:assert actual.resname==r.resname
            for name in ['N','CA','C','O']:
                maximum=max(maximum,float(np.linalg.norm(actual[name].coord-r[name].coord)))
        groups[row['label']].append(row)
        checks.append(dict(seed=row['seed'],label=row['label'],sha256_ok=True,mutation_identity_ok=True))
    assert maximum<=.002
    rows=[]
    for label,replicas in groups.items():
        his={h['resnum'] for r in replicas for h in r['selected_his']}
        ranges={str(n):[min(h['pka'] for r in replicas for h in r['selected_his'] if h['resnum']==n),
                        max(h['pka'] for r in replicas for h in r['selected_his'] if h['resnum']==n)] for n in his}
        rows.append(dict(label=label,replicas=len(replicas),neutral_retained=sum(r['neutral_energy_filter'] for r in replicas),
                         ddE_min=min(r['ddE_vs_HH_REU'] for r in replicas),ddE_max=max(r['ddE_vs_HH_REU'] for r in replicas),
                         selected_his_pka_ranges=json.dumps(ranges),
                         both_selected_free_pka_ge7=len(his)==2 and all(lo>=7 for lo,hi in ranges.values())))
    with (out/'variants.csv').open('w',newline='') as fh:
        writer=csv.DictWriter(fh,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    dump(out/'verification.json',dict(models=len(checks),max_backbone_displacement_A=maximum,checks=checks,
                                    status_sha256=sha(out/'status.json'),summary_script_sha256=sha(__file__)))
    (out/'RESULTS.md').write_text(
        '# Isolated-binder pH-module scout\n\n'
        f"{len(data['models'])} models: 24 variants × 3 side-chain preparation seeds; all completed. "
        'All sequence identities/checksums verified; maximum backbone displacement '
        f'{maximum:.4f} Å. The native-code crash in the first attempt is preserved separately '
        'and fixed by using three-letter residue names with MutateResidue.\n\n'
        '**No variant produced two selected free-state His pKas >=7.** Independently '
        'repacking the isolated HH binder gives H4=6.33/H46=6.12 for all three seeds. '
        'The old bound-conformation extraction gave H46=7.25; its favorable value did '
        'not survive free-side-chain reorganization. These are fixed-backbone replicas, '
        'not an equilibrium ensemble or biological probabilities.\n\n'
        'R49E raises free H46 to6.40 reproducibly while remaining below the same-seed '
        'HH energy proxy; V8D raises H4 to6.42. L50E reaches H46=6.74 in two replicas '
        'but6.20 in the third. None is a strong two-His module. These are mechanism '
        'diagnostics, not improved binding candidates.\n\n'
        '**The source docking pose fails the full-receptor audit.** This scout only '
        'tests reusable isolated-binder sequence/conformation hypotheses. No receptor '
        'affinity, pH linkage, mouse cross-reactivity or experimental stability was '
        'established. Same-seed Rosetta ddE<=3 REU is an engineering filter; it is '
        'not folding free energy. Do not submit these variants on the basis of this run.\n\n'
        'Do not expand this nearest-neighbor carboxylate search: its proposed free-pKa '
        'benefit is too small on this scaffold. Advance fresh, physically accessible '
        'full-context binders, then install geometrically explicit networks.\n\n'
        'GPU spend from this CPU scout: $0. Separate active GPU spend and shutdown '
        'state: ../gpu_campaign/supervisor.json.\n')
    print(f'{len(checks)} models verified; backbone max {maximum:.4f} A')


if __name__=='__main__':main()
