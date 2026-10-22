"""Extract exactly the two predeclared sequences; no coordinates enter the GPU input."""
import hashlib
import json
from pathlib import Path
from Bio.Data.PDBData import protein_letters_3to1 as AA
from Bio.PDB import PDBParser
from cif2pdb import read_atom_site

ROOT=Path(__file__).resolve().parents[1]
INPUTS=(
 ('v9','results/modal_repaired_native_v9/mirrored/campaign/1_Trajectories/c4_repaired_native_v9_denovo_l69_487e943cff8955d6/c4_repaired_native_v9_denovo_l69_487e943cff8955d6_trajectory.cif','8311b87eff4d30aa7c4323d4ed4051cf4ef7d50d6f631af9a3d609fd85949a5d','2699888c00c9c6af76b61c8547687cce87ccb0152754ae0d0a32c10e5f42b85d'),
 ('seed23','results/modal_native_replicates_v14/replicate_3/mirrored/campaign/1_Trajectories/c4_repaired_native_v9_replicate_3_denovo_l67_a207b66be0c44216/c4_repaired_native_v9_replicate_3_denovo_l67_a207b66be0c44216_trajectory.cif','a5e467a235cfa2237c36740e0d036c6bc8d6985332a72fe0ea4ddb2de5beaf69','8b0e2189b54820dd47b3d25e26969dd1a519f6c9a7267aee9d6f50a67fefa2b3'))

def digest(value):return hashlib.sha256(value).hexdigest()

def main():
    out=ROOT/'results/hard_sequence_refold_v15';out.mkdir(exist_ok=True)
    manifest=dict(protocol=1,models=['model_1_ptm','model_2_ptm'],recycles=3,
                  max_prediction_jobs=8,designed_binder_coordinates_supplied=False,sequences=[])
    for name,relative,source_hash,sequence_hash in INPUTS:
        source=ROOT/relative
        if digest(source.read_bytes())!=source_hash:raise ValueError('Declared artifact hash differs')
        residues={}
        for atom in read_atom_site(source):
            if atom['auth_asym_id']=='B':
                key=(atom['auth_seq_id'],atom.get('pdbx_PDB_ins_code','.'))
                value=AA[atom['auth_comp_id']]
                if key in residues and residues[key]!=value:raise ValueError('Conflicting residue identity')
                residues[key]=value
        sequence=''.join(residues.values())
        if digest(sequence.encode())!=sequence_hash:raise ValueError('Declared sequence hash differs')
        manifest['sequences'].append(dict(name=name,sequence=sequence,length=len(sequence),
            sequence_sha256=sequence_hash,source=relative,source_sha256=source_hash))
    target=ROOT/'results/modal_campaign_v5/human_11qc_local_frame.pdb'
    chain=PDBParser(QUIET=True).get_structure('native',target)[0]['A']
    seq=''.join(AA[r.resname] for r in chain)
    text=(ROOT/'data/P02786.txt').read_text()
    expected=''.join(text.split('\nSQ   ',1)[1].split('\n',1)[1].split('//')[0].split())[121:760]
    if seq!=expected or [r.id[1] for r in chain]!=list(range(122,761)):
        raise ValueError('Native target does not match UniProt122-760')
    manifest.update(native_target_sha256=digest(target.read_bytes()),native_sequence_sha256=digest(seq.encode()),
                    complex_scope='Complete639-residue native monomer prediction; restore/check dimer B afterward; not dimer-aware refolding')
    (out/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    (out/'inputs.fasta').write_text(''.join('>'+s['name']+'\n'+s['sequence']+'\n' for s in manifest['sequences']))
    print('Prepared two fixed sequences; native target identity verified; no binder coordinates exported.')

if __name__=='__main__':main()
