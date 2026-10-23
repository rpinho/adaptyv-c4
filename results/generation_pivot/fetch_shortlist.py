"""Read only ranked existing CIFs from the volume into a content-addressed cache."""
import argparse
import csv
import hashlib
import json
from pathlib import Path,PurePosixPath
import re
import tempfile

import modal


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--shortlist',type=Path,required=True)
    parser.add_argument('--candidates',type=Path,nargs='+',required=True)
    parser.add_argument('--out',type=Path,required=True)
    args=parser.parse_args()
    if args.out.exists():raise FileExistsError('Preserve existing coordinate receipt')
    candidates={};total=0
    for path in args.candidates:
        with path.open(newline='') as stream:
            for row in csv.DictReader(stream):candidates[row['name']]=row;total+=1
    if total!=len(candidates):raise ValueError('Duplicate candidate names')
    with args.shortlist.open(newline='') as stream:selected=list(csv.DictReader(stream))
    requests=[]
    for row in selected:
        record=candidates[row['source_id']]
        if row['sequence']!=record['sequence']:raise ValueError('Ranked sequence differs from source')
        arm,case=record['source_arm'],record['source_design_id']
        if not all(re.fullmatch(r'[A-Za-z0-9_]+',x) for x in (arm,case)):raise ValueError('Invalid source identifier')
        relative=PurePosixPath(record['designed_structure_path']).relative_to('/runs')
        if relative.parts!=(arm,'inference_samples',case,case+'_generated.cif'):raise ValueError('Unexpected structure path')
        requests.append((row['source_id'],str(relative),record))
    cache=Path.home()/'.cache/adaptyv-c4/rfd4_coordinates';cache.mkdir(parents=True,exist_ok=True)
    volume=modal.Volume.from_name('rfd4-runs');structures=[]
    for name,remote,record in requests:
        digest=hashlib.sha256()
        with tempfile.NamedTemporaryFile(dir=cache,delete=False,suffix='.partial') as stream:
            temporary=Path(stream.name)
            for chunk in volume.read_file(remote):stream.write(chunk);digest.update(chunk)
        if not temporary.stat().st_size:raise ValueError('Empty source CIF')
        fingerprint=digest.hexdigest();destination=cache/(fingerprint+'.cif')
        if destination.exists():
            if hashlib.sha256(destination.read_bytes()).hexdigest()!=fingerprint:raise ValueError('Cached content hash changed')
            temporary.unlink()
        else:temporary.rename(destination)
        structures.append(dict(name=name,path=str(destination),sha256=fingerprint,source_arm=record['source_arm'],
                               source_design_id=record['source_design_id'],remote_path=remote,status='available'))
    result=dict(structures=structures,scope='Read-only Volume API, no App or Function invoked',
                shortlist_sha256=hashlib.sha256(args.shortlist.read_bytes()).hexdigest(),new_gpu_or_remote_cpu_spend_usd=0)
    args.out.parent.mkdir(parents=True,exist_ok=True);args.out.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(fetched=len(structures),manifest=str(args.out),apps_started=0)))


if __name__=='__main__':main()
