"""Select one predeclared domain output for diagnostics; never promote its scores."""
import csv
import hashlib
import json
from pathlib import Path
from Bio.Data.PDBData import protein_letters_3to1 as AA
from Bio.PDB import PDBParser
from cif2pdb import read_atom_site,convert
from restore_apical_context import restore

ROOT=Path(__file__).resolve().parents[1]
STATE=ROOT/'results/native_domain_v16'
OUT=ROOT/'results/domain_hard_refold_v17'

def sha(value):return hashlib.sha256(value).hexdigest()

def main():
    status=json.loads((STATE/'supervisor.json').read_text())
    if not status.get('gpu_released_verified') or status.get('result',{}).get('status')!='completed':
        raise RuntimeError('Completed domain probe and actual GPU release required')
    OUT.mkdir(exist_ok=True)
    campaign=STATE/'mirrored/campaign';ranked=campaign/'3_Ranked/!_Ranked.csv'
    rows=list(csv.DictReader(ranked.open())) if ranked.exists() else []
    if rows:
        name=rows[0]['design']
        if Path(name).name!=name:raise ValueError('Unsafe design identity')
        source=campaign/'3_Ranked'/(name+'.cif');selection='First ranked design by upstream ranking'
    else:
        sources=list((campaign/'1_Trajectories').glob('*/*_trajectory.cif'))
        if len(sources)!=1:raise ValueError('Exactly one declared trajectory snapshot required')
        source=sources[0];selection='Only saved first-trajectory snapshot, even if soft screen rejected'
    chains={'A':{},'B':{}}
    for atom in read_atom_site(source):
        if atom.get('pdbx_PDB_model_num','1')!='1' or atom.get('pdbx_PDB_ins_code','.') not in ('.','?',''):
            raise ValueError('Multiple models or insertion codes are not supported')
        if atom.get('label_alt_id','.') not in ('.','?',''):raise ValueError('Alternate conformers are not supported')
        cid=atom['auth_asym_id']
        if cid not in chains:raise ValueError('Unexpected domain output chain')
        number=int(atom['auth_seq_id']);aa=AA[atom['auth_comp_id']]
        if number in chains[cid] and chains[cid][number]!=aa:raise ValueError('Conflicting residue identities')
        chains[cid][number]=aa
    human=(ROOT/'data/P02786.txt').read_text()
    seq=''.join(human.split('\nSQ   ',1)[1].split('\n',1)[1].split('//')[0].split())
    if list(chains['A'])!=list(range(189,384)) or ''.join(chains['A'].values())!=seq[188:383]:
        raise ValueError('Domain output sequence/author mapping differs from native')
    binder=''.join(chains['B'].values())
    if list(chains['B'])!=list(range(1,70)) or len(binder)!=69:raise ValueError('Expected the complete declared69-residue binder')
    raw=OUT/'domain_source_raw.pdb';convert(str(source),str(raw))
    audit=dict(status='not_qualified',selection_rule=selection,source=str(source.relative_to(ROOT)),
               source_sha256=sha(source.read_bytes()),generator_accepted_rows=len(rows))
    try:audit['restoration']=restore(raw,OUT/'restored_staging_only.pdb');audit['status']='restored_staging_only'
    except ValueError as exc:audit['restoration_rejection']=str(exc)
    (OUT/'SOURCE_AUDIT.json').write_text(json.dumps(audit,indent=2)+'\n')
    native=ROOT/'results/modal_campaign_v5/human_11qc_local_frame.pdb'
    native_seq=''.join(AA[r.resname] for r in PDBParser(QUIET=True).get_structure('native',native)[0]['A'])
    manifest=dict(protocol=2,source_campaign='c4_native_domain_v16',models=['model_1_ptm','model_2_ptm'],
        recycles=3,max_prediction_jobs=4,designed_binder_coordinates_supplied=False,
        checkpoint_scope='AF2 monomer checkpoints may already have selected accepted domain designs; additional native-context test, not independent architecture',
        sequences=[dict(name='domain_primary',sequence=binder,length=len(binder),sequence_sha256=sha(binder.encode()),
            source=str(source.relative_to(ROOT)),source_sha256=sha(source.read_bytes()))],
        native_target_sha256=sha(native.read_bytes()),native_sequence_sha256=sha(native_seq.encode()),
        complex_scope='Complete639-residue native monomer prediction; restore/check dimer B afterward; not dimer-aware refolding')
    (OUT/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    (OUT/'input.fasta').write_text('>domain_primary\n'+binder+'\n')
    print('Selected one fixed domain output for coordinate-free validation; source scores confer no qualification.')

if __name__=='__main__':main()
