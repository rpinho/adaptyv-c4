"""Extract verified full-backbone analysis inputs; never repair missing atoms."""
from collections import Counter
import argparse
import hashlib
import json
import math
from pathlib import Path
import sys

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
sys.path.insert(0,str(ROOT/'scripts'))
import novelty_level_predictor as n
from Bio.PDB.MMCIF2Dict import MMCIF2Dict


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--report',type=Path)
    ap.add_argument('--designs',type=Path)
    ap.add_argument('--out',type=Path,default=HERE)
    ap.add_argument('--coordinate-cache',type=Path)
    args=ap.parse_args()
    benchmark=args.report is None and args.designs is None
    output=args.out.resolve();output.mkdir(parents=True,exist_ok=True)
    if (output/'INPUTS.json').exists(): raise FileExistsError('Input inventory already frozen')
    source_report=args.report or ROOT/'results/novelty_predictor/exact_score/run/predictions.json'
    baseline=json.loads(source_report.read_text())
    designs=args.designs or ROOT/'results/novelty_predictor/analysis_inputs.csv'
    rows=n.load_designs(designs)
    assert len(rows)==len(baseline['designs'])
    folder=args.coordinate_cache or Path.home()/'.cache/adaptyv-c4/novelty_domain_source_audit/full_backbones_complete_ter'
    folder=folder.resolve();folder.mkdir(parents=True)
    inventory=[]
    for i,(row,result) in enumerate(zip(rows,baseline['designs']),1):
        assert row['name']==result['name']
        query=f'q{i:06d}'
        evidence=result['structure_evidence']
        record=dict(query=query,name=row['name'],sequence_sha256=result['sequence_sha256'],
            sequence_md5=hashlib.md5(row['sequence'].encode()).hexdigest(),length=len(row['sequence']))
        if evidence['status']!='coordinates_verified':
            record.update(status='coordinates_unavailable',reason=evidence.get('reason','Frozen coordinates unavailable'))
            inventory.append(record); continue
        spec=evidence['input']; source=Path(spec['path'])
        assert n.digest(source)==spec['sha256']
        residues=n.ca_residues(source,spec['chain'])
        if 'author_residue_range' in spec:
            lo,hi=spec['author_residue_range']; residues=[r for r in residues if lo<=r[0]<=hi]
        assert ''.join(r[1] for r in residues)==row['sequence']
        numbering={r[0]:j for j,r in enumerate(residues,1)}
        comps={r[0]:r[3] for r in residues}
        model_source=source
        if source.suffix.lower()!='.cif':
            from Bio.PDB import PDBParser,MMCIFIO
            structure=PDBParser(QUIET=True).get_structure(query,str(source))
            model_source=folder/(query+'_source.cif')
            writer=MMCIFIO();writer.set_structure(structure);writer.save(str(model_source))
            record['pdb_to_cif_format_conversion']=dict(path=str(model_source),sha256=n.digest(model_source))
        data=MMCIF2Dict(str(model_source))
        def column(auth,label): return data.get('_atom_site.'+auth,data.get('_atom_site.'+label))
        atom_names=column('auth_atom_id','label_atom_id')
        atom_chains=column('auth_asym_id','label_asym_id')
        atom_comps=column('auth_comp_id','label_comp_id')
        atom_nums=column('auth_seq_id','label_seq_id')
        count=len(atom_names)
        models=data.get('_atom_site.pdbx_PDB_model_num',['1']*count)
        alt=data.get('_atom_site.label_alt_id',['.']*count)
        elements=data.get('_atom_site.type_symbol',[name[0] for name in atom_names])
        occupancy=data.get('_atom_site.occupancy',['1']*count)
        bfactor=data.get('_atom_site.B_iso_or_equiv',['0']*count)
        xyz=list(zip(*(data['_atom_site.Cartn_'+d] for d in 'xyz')))
        lines=[]; observed=Counter(); copied_ca={}
        for atom,chain,comp,num,model,al,element,occ,bf,coords in zip(
                atom_names,atom_chains,atom_comps,atom_nums,models,alt,elements,occupancy,bfactor,xyz):
            if chain!=spec['chain'] or model!=models[0] or al not in ('.','?','A'): continue
            num=int(num)
            if num not in numbering: continue
            assert comp==comps[num], 'Nonprotein/ambiguous residue at selected author number'
            assert 1<=len(atom)<=4 and atom.isalnum()
            position=numbering[num]
            observed[position,atom]+=1
            assert observed[position,atom]==1, 'Duplicate selected atom'
            x,y,z=map(float,coords)
            assert all(math.isfinite(v) and -999.999<=v<=9999.999 for v in (x,y,z))
            if atom=='CA': copied_ca[position]=(x,y,z)
            atom_field=(' '+atom.ljust(3)) if len(atom)<4 else atom
            lines.append(f'ATOM  {len(lines)+1:5d} {atom_field} {comp:3s} A{position:4d}    '
                f'{x:8.3f}{y:8.3f}{z:8.3f}{float(occ):6.2f}{float(bf):6.2f}          {element:>2s}\n')
        for position in numbering.values():
            assert all(observed[position,atom]==1 for atom in ('N','CA','C','O')), 'Missing backbone atom; no repair allowed'
        path=folder/(query+'.pdb')
        footer=f'TER   {len(lines)+1:5d}      {residues[-1][3]:3s} A{len(residues):4d}\nEND\n'
        path.write_text(''.join(lines)+footer)
        reread=n.ca_residues(path,'A')
        assert ''.join(r[1] for r in reread)==row['sequence']
        text_ca={int(line[22:26]):tuple(float(line[a:b]) for a,b in ((30,38),(38,46),(46,54)))
                 for line in lines if line[12:16].strip()=='CA'}
        delta=max(abs(a-b) for position,coords in text_ca.items() for a,b in zip(coords,copied_ca[position]))
        assert delta<=0.00050000001, 'Coordinate change beyond PDB decimal rounding'
        parser_error=max(abs(a-b) for position,_,coords,_ in reread for a,b in zip(coords,text_ca[position]))
        max_magnitude=max(abs(v) for coords in text_ca.values() for v in coords)
        assert parser_error<=(2**-23)*max_magnitude, 'Unexpected parser quantization beyond float32 precision'
        record.update(status='full_backbone_verified',path=str(path),sha256=n.digest(path),
            source=spec,atoms=len(lines),backbone_atoms=4*len(residues),max_ca_rounding_angstrom=delta,
            max_biopython_float32_parser_error_angstrom=parser_error,
            numbering='Author-to-query-position bijection; sequence order unchanged; chain A',
            missing_occupancy_defaulted_to_one='_atom_site.occupancy' not in data,
            missing_bfactor_defaulted_to_zero='_atom_site.B_iso_or_equiv' not in data)
        inventory.append(record)
    available=[r for r in inventory if r['status']=='full_backbone_verified']
    assert len(available)==baseline['structural_search']['query_domains'], 'Parent must use one full domain per available protein'
    if benchmark:
        assert len(inventory)==74 and len(available)==68
        native=next(r for r in inventory if r['name']=='KNOWN_TFR1_POSITIVE_CONTROL')
        assert native['status']=='full_backbone_verified' and native['length']==195
        (output/'control_targets.txt').write_text(native['path']+'\n')
    (output/'all_targets.txt').write_text(''.join(r['path']+'\n' for r in available))
    n.write_json(output/'INPUTS.json',dict(source_report_sha256=n.digest(source_report),designs_csv_sha256=n.digest(designs),
        protocol_sha256=n.digest(HERE/'PROTOCOL.md'),coordinate_models=len(available),unavailable=len(inventory)-len(available),
        source_backbone_or_sequence_repaired=False,inputs=inventory))
    print(json.dumps(dict(verified_full_backbones=len(available),unavailable=len(inventory)-len(available),atoms_repaired=0)))


if __name__=='__main__': main()
