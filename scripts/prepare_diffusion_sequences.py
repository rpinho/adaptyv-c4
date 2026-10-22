"""Exactly two original MPNN samples after native-context backbone staging."""
from pathlib import Path
import hashlib
import json
import os
import re
import subprocess
import ast
import numpy as np
from Bio.PDB import PDBParser
from Bio.Data.PDBData import protein_letters_3to1 as AA

ROOT=Path(__file__).resolve().parents[1]
STATE=ROOT/'results/diffusion_v18'
OUT=ROOT/'results/diffusion_hard_refold_v19'
ALPHABET='ACDEFGHIKLMNPQRSTVWYX'

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def main():
    audit=json.loads((STATE/'BACKBONE_AUDIT.json').read_text())
    if audit.get('sequence_design_allowed') is not True:raise RuntimeError('Physically staged backbone required')
    source=STATE/'restored_staging_only.pdb'
    if sha(source)!=audit['restoration']['output_sha256']:raise ValueError('Staged backbone identity changed')
    preparation=json.loads((STATE/'MPNN_PREPARATION.json').read_text());repo=Path(preparation['source_directory'])
    if subprocess.check_output(['git','rev-parse','HEAD'],cwd=repo,text=True).strip()!=preparation['revision']:raise ValueError('MPNN revision changed')
    if any(sha(repo/name)!=expected for name,expected in preparation['files'].items()):raise ValueError('MPNN source/weight identity changed')
    folder=STATE/'mpnn'
    if folder.exists() or OUT.exists():raise FileExistsError('No resampling or manifest replacement')
    folder.mkdir();OUT.mkdir()
    command=[str(ROOT/'.venv-mpnn/bin/python'),str(repo/'protein_mpnn_run.py'),
        '--pdb_path',str(source),'--pdb_path_chains','C','--out_folder',str(folder),
        '--num_seq_per_target','2','--batch_size','1','--sampling_temp','0.1','--seed','20261022',
        '--omit_AAs','CX','--model_name','v_48_020','--save_probs','1']
    env=dict(os.environ,TORCH_FORCE_NO_WEIGHTS_ONLY_LOAD='1',OMP_NUM_THREADS='4')
    (folder/'COMMAND.json').write_text(json.dumps(command,indent=2)+'\n')
    with (folder/'worker.log').open('w') as log:
        subprocess.run(command,env=env,stdout=log,stderr=subprocess.STDOUT,check=True,timeout=600)
    fasta=folder/'seqs/restored_staging_only.fa';lines=fasta.read_text().splitlines()
    if len(lines)!=6 or not all(lines[i].startswith('>') for i in [0,2,4]):raise ValueError('Expected exactly two samples and one native reference')
    first=lines[0]
    fixed=ast.literal_eval(re.search(r'fixed_chains=(\[[^]]*\])',first)[1]);designed=ast.literal_eval(re.search(r'designed_chains=(\[[^]]*\])',first)[1])
    if set(fixed)!={'A','B'} or designed!=['C'] or lines[1]!='G'*69:raise ValueError('MPNN chain assignment differs')
    sequences=[lines[3],lines[5]]
    if any(len(s)!=69 or set(s)-set(ALPHABET[:-1]) or 'C' in s for s in sequences):raise ValueError('Unexpected designed sequence')
    model=PDBParser(QUIET=True).get_structure('fixed_native',source)[0]
    chains={c.id:''.join(AA[r.resname] for r in c) for c in model}
    with np.load(folder/'probs/restored_staging_only.npz',allow_pickle=False) as data:
        order=data['chain_order'].tolist()
        if len(order)!=1 or set(order[0])!={'A','B','C'}:raise ValueError('Unexpected saved chain ordering')
        native=''.join(chains[c] for c in order[0]);mask=data['mask'][0].astype(bool)
        expected_mask=np.array([c=='C' for c in order[0] for _ in chains[c]])
        if not np.array_equal(mask,expected_mask) or int(mask.sum())!=69:raise ValueError('MPNN did not design only binder')
        if data['S'].shape!=(2,len(native)):raise ValueError('Unexpected sampled full-sequence shape')
        for sample,seq in zip(data['S'],sequences):
            text=''.join(ALPHABET[i] for i in sample)
            if any(a!=b for a,b,m in zip(text,native,mask) if not m):raise ValueError('Native receptor sequence was changed')
            if ''.join(a for a,m in zip(text,mask) if m)!=seq:raise ValueError('FASTA/sample tensor disagree')
    native_target=ROOT/'results/modal_campaign_v5/human_11qc_local_frame.pdb'
    target_sequence=''.join(AA[r.resname] for r in PDBParser(QUIET=True).get_structure('native',native_target)[0]['A'])
    entries=[dict(name='rfd_mpnn_'+str(i),sequence=s,length=69,sequence_sha256=hashlib.sha256(s.encode()).hexdigest(),
                  source=str(fasta.relative_to(ROOT)),source_sha256=sha(fasta)) for i,s in enumerate(sequences,1)]
    manifest=dict(protocol=3,source_campaign='c4_diffusion_v18',models=['model_1_ptm','model_2_ptm'],recycles=3,
        max_prediction_jobs=8,designed_binder_coordinates_supplied=False,sequences=entries,
        checkpoint_scope='AF2 validation differs from RFdiffusion/RoseTTAFold generation; both validation checkpoints share AF2 architecture',
        native_target_sha256=sha(native_target),native_sequence_sha256=hashlib.sha256(target_sequence.encode()).hexdigest(),
        backbone_source_sha256=audit['source_sha256'],mpnn_native_target_fixed_verified=True,
        mpnn_samples=2,unique_sequences=len(set(sequences)),independent_backbones=1,
        complex_scope='Complete639-residue native monomer prediction; CPU dimer restoration afterward, not dimer-aware inference')
    (OUT/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    (OUT/'input.fasta').write_text(''.join('>'+e['name']+'\n'+e['sequence']+'\n' for e in entries))
    (folder/'RESULTS.json').write_text(json.dumps(dict(status='completed',samples=2,unique_sequences=len(set(sequences)),
        independent_backbones=1,target_A_B_fixed_verified=True,fasta_sha256=sha(fasta),qualified_candidates=0,
        interpretation='Sequence assignment only; MPNN scores are not affinity/novelty/pH evidence'),indent=2)+'\n')
    print('Two fixed MPNN samples saved; both native target chains verified unchanged',flush=True)

if __name__=='__main__':main()
