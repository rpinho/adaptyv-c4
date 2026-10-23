#!/usr/bin/env python3
"""One CPU workflow for name,sequence CSVs and optional verified coordinates.

Uses existing public reference snapshots and the locally calibrated domain tools.
No structure generation, GPU, cloud service, entry exporter or submission action.
Missing coordinates/consensus leave levels unestablished rather than inferred.
"""
import argparse
import copy
from datetime import datetime,timezone
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import uuid

import novelty_level_predictor as n

ROOT=Path(__file__).resolve().parents[1]
TOOLS=ROOT/'results/novelty_predictor/domain_consensus'
CACHE=Path.home()/'.cache/adaptyv-c4/novelty_domain_source_audit'


def selected_residues(path,chain,spec):
    residues=n.ca_residues(path,chain)
    if 'author_residue_range' in spec:
        lo,hi=spec['author_residue_range'];residues=[r for r in residues if lo<=r[0]<=hi]
    return residues


def normalize_structures(records,manifest_path):
    supplied=json.loads(manifest_path.read_text()) if manifest_path else dict(structures=[])
    specs={r['name']:r for r in supplied['structures']}
    if len(specs)!=len(supplied['structures']) or set(specs)-{r['name'] for r in records}:
        raise ValueError('Duplicate or unknown coordinate name')
    normalized=[]
    for record in records:
        spec=copy.deepcopy(specs.get(record['name'],dict(name=record['name'],status='unavailable',reason='No coordinates supplied')))
        if spec.get('status','available' if spec.get('path') else 'unavailable')!='available':
            spec['status']='unavailable';normalized.append(spec);continue
        path=Path(spec['path']).expanduser()
        if not path.is_absolute():path=manifest_path.resolve().parent/path
        path=path.resolve();actual=n.digest(path)
        if 'sha256' in spec and actual!=spec['sha256']:raise ValueError('Supplied coordinate fingerprint differs: '+record['name'])
        spec['source_hash_was_supplied']='sha256' in spec
        spec.update(path=str(path),sha256=actual,status='available')
        if 'chain' not in spec:
            if path.suffix.lower()=='.cif':
                from Bio.PDB.MMCIF2Dict import MMCIF2Dict
                data=MMCIF2Dict(str(path))
                chains=set(data.get('_atom_site.auth_asym_id',data.get('_atom_site.label_asym_id',[])))
            else:
                from Bio.PDB import PDBParser
                chains={c.id for c in PDBParser(QUIET=True).get_structure('input',str(path))[0]}
            matching=[]
            for chain in sorted(chains):
                try:sequence=''.join(r[1] for r in selected_residues(path,chain,spec))
                except (KeyError,ValueError):continue
                if sequence==record['sequence']:matching.append(chain)
            if len(matching)!=1:raise ValueError('Require one exact sequence-matched chain or an explicit chain: '+record['name'])
            spec.update(chain=matching[0],chain_inferred_by_unique_full_sequence=True)
        if ''.join(r[1] for r in selected_residues(path,spec['chain'],spec))!=record['sequence']:
            raise ValueError('Coordinate/sequence mismatch: '+record['name'])
        spec['domains']=[dict(start=1,end=len(record['sequence']),source='Whole-chain baseline used for parent-context candidate discovery')]
        normalized.append(spec)
    return dict(structures=normalized)


def verify_domain_installation():
    for filename in ('NATIVE_BUILD.json','WEIGHTS.json','INDEXING_CHECKS.json'):
        if not (TOOLS/filename).exists():raise FileNotFoundError('Calibrated CPU domain tools not installed: '+filename)
    native=json.loads((TOOLS/'NATIVE_BUILD.json').read_text())
    for path,fingerprint in native['binaries'].items():
        if n.digest(path)!=fingerprint:raise ValueError('Native domain executable changed')
    weights=json.loads((TOOLS/'WEIGHTS.json').read_text())
    for item in weights['weights']:
        if n.digest(item['path'])!=item['sha256']:raise ValueError('Merizo weight fingerprint differs')
    indexing=json.loads((TOOLS/'INDEXING_CHECKS.json').read_text())
    if not indexing['all_passed'] or n.digest(indexing['corrected_source'])!=indexing['corrected_sha256']:
        raise ValueError('Calibrated consensus implementation differs')
    source=json.loads((TOOLS/'TED_SOURCE_VERIFICATION.json').read_text())
    prefix='ted_consensus_1.0/'
    for item in source['files']:
        path=CACHE/'ted-pinned'/item['path'][len(prefix):]
        if n.digest(path)!=item['sha256']:raise ValueError('Pinned domain source/weight changed: '+item['path'])


def main(argv=None):
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('designs',type=Path)
    ap.add_argument('--structures',type=Path,help='JSON structures: name,path, optional chain/sha256; exact sequence match required')
    ap.add_argument('--databases',type=Path,default=ROOT/'results/novelty_predictor/patent_extension/DATABASES.json')
    ap.add_argument('--structure-database',type=Path,default=ROOT/'results/novelty_predictor/structural/DATABASE.json')
    ap.add_argument('--foldseek',default=str(Path.home()/'.local/opt/foldseek/foldseek/bin/foldseek'))
    ap.add_argument('--mmseqs',default='mmseqs')
    ap.add_argument('--timeout-structure',type=int,default=300,
        help='Seconds per structural search/alignment stage; increase for larger saved batches')
    ap.add_argument('--reuse-whole-chain-searches',type=Path,
        help='Reuse completed fingerprinted sequence/discovery checkpoints after an interrupted alignment')
    ap.add_argument('--positive-evidence-only',action='store_true',
        help='Bound alignment work and withhold conclusions not fixed by positive evidence')
    ap.add_argument('--out',type=Path,required=True)
    args=ap.parse_args(argv)
    if args.timeout_structure<=0:ap.error('--timeout-structure must be positive')
    records=n.load_designs(args.designs)
    if args.out.exists():raise FileExistsError('Fresh output directory required')
    specs=normalize_structures(records,args.structures)
    available=sum(r['status']=='available' for r in specs['structures'])
    if available:verify_domain_installation()
    out=args.out.resolve();out.mkdir(parents=True)
    state=dict(state='prepared',controller_pid=os.getpid(),inputs=len(records),coordinate_inputs=available,
        started_at=datetime.now(timezone.utc).isoformat(),steps=[],new_gpu_or_cloud_cpu_spend_usd=0,
        design_csv_sha256=n.digest(args.designs),workflow_sha256=n.digest(__file__))
    def save():n.write_json(out/'STATUS.json',state)
    def run(label,arguments):
        step=dict(name=label,command=[str(x) for x in arguments],started_at=datetime.now(timezone.utc).isoformat())
        state['steps'].append(step);state['state']='running';save()
        with (out/(label+'.log')).open('w') as log:
            process=subprocess.Popen(step['command'],stdout=log,stderr=subprocess.STDOUT,
                env={**os.environ,'CUDA_VISIBLE_DEVICES':'','OMP_NUM_THREADS':'4','OPENBLAS_NUM_THREADS':'4'})
            step['pid']=process.pid;save();code=process.wait()
            if code:raise subprocess.CalledProcessError(code,step['command'])
        step['completed']=True;save()
    try:
        parent_manifest=out/'whole_chain_structures.json';n.write_json(parent_manifest,specs)
        command=[sys.executable,ROOT/'scripts/novelty_level_predictor.py',args.designs.resolve(),
            '--databases',args.databases.resolve(),'--mmseqs',args.mmseqs,'--threads','4',
            '--timeout-per-database','600','--out',out/'whole_chain']
        if available:
            command += ['--structures',parent_manifest,'--structure-database',args.structure_database.resolve(),
                        '--foldseek',args.foldseek,'--timeout-structure',str(args.timeout_structure)]
        if args.reuse_whole_chain_searches:
            if not available:raise ValueError('Whole-chain checkpoint reuse requires coordinates')
            command += ['--cached-sequence-searches',args.reuse_whole_chain_searches.resolve(),
                '--cached-structural-discovery',args.reuse_whole_chain_searches.resolve()/'structure']
        if args.positive_evidence_only:command += ['--positive-evidence-only']
        run('whole_chain',command)
        final=out/'whole_chain/predictions.json'
        if available:
            # Safe generated paths keep the unchanged third-party shell wrappers
            # away from user-controlled file names or shell metacharacters.
            runtime=CACHE/'workflow_runtime'/uuid.uuid4().hex;runtime.mkdir(parents=True)
            state['domain_runtime']=str(runtime);save()
            run('prepare_backbones',[sys.executable,TOOLS/'prepare_inputs.py','--report',final,
                '--designs',args.designs.resolve(),'--out',out/'domain_inputs','--coordinate-cache',runtime/'full_backbones'])
            run('domain_models',[sys.executable,TOOLS/'run_models.py','all','--input-manifest',out/'domain_inputs/INPUTS.json',
                '--output-dir',runtime/'models'])
            copied=out/'domain_models';copied.mkdir()
            for path in (runtime/'models').iterdir():
                if path.is_file():shutil.copy2(path,copied/path.name)
            domain_manifest=out/'consensus_structures.json'
            run('make_domains',[sys.executable,TOOLS/'make_structure_manifest.py','--input-manifest',out/'domain_inputs/INPUTS.json',
                '--model-run',copied,'--designs',args.designs.resolve(),'--parent-structures',parent_manifest,'--out',domain_manifest])
            segmented=json.loads(domain_manifest.read_text())
            domains_available=sum(r['status']=='available' for r in segmented['structures'])
            if domains_available:
                run('domain_search',[sys.executable,ROOT/'scripts/novelty_level_predictor.py',args.designs.resolve(),
                    '--sequence-report',final,'--structures',domain_manifest,'--structure-database',args.structure_database.resolve(),
                    '--foldseek',args.foldseek,'--threads','4','--timeout-structure',str(args.timeout_structure),'--out',out/'domain_search']
                    +(['--positive-evidence-only'] if args.positive_evidence_only else []))
                run('parent_context',[sys.executable,ROOT/'scripts/augment_novelty_parent_context.py',
                    '--base-report',out/'domain_search/predictions.json','--parent-report',final,
                    '--foldseek',args.foldseek,'--threads','4','--timeout-stage',str(args.timeout_structure),'--out',out/'parent_context'])
                final=out/'parent_context/predictions.json'
            else:
                withheld=json.loads(final.read_text())
                for design in withheld['designs']:
                    design.update(local_predicted_level=None,predicted_level=None,
                        structure_status='no_accepted_consensus_domain',local_gate_assessment='unestablished',
                        structural_similarity=None,
                        structure_evidence=dict(status='no_accepted_consensus_domain',whole_chain_baseline_report=str(final)),
                        explanation='No high/medium consensus domain; numeric structural novelty withheld',organizer_level_established=False)
                folder=out/'consensus_unestablished';folder.mkdir();n.save_report(withheld,folder)
                final=folder/'predictions.json'
        report=json.loads(final.read_text())
        state.update(state='complete',final_report=str(final),final_csv=str(final.with_suffix('.csv')),
            final_report_sha256=n.digest(final),organizer_levels_established=0,
            local_levels_predicted=sum(d.get('local_predicted_level') is not None for d in report['designs']),
            missing_structure_inference='No structures are generated; supplied coordinates only')
        save();print(json.dumps({k:state[k] for k in ['state','inputs','coordinate_inputs','local_levels_predicted','final_report','organizer_levels_established','new_gpu_or_cloud_cpu_spend_usd']}))
    except Exception as error:
        state.update(state='failed',error=repr(error));save();raise


if __name__=='__main__':main()
