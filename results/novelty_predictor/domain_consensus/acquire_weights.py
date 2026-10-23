"""Fetch pinned official Merizo shards; verify Git identities and local SHA256."""
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import urllib.request

HERE=Path(__file__).resolve().parent
CACHE=Path.home()/'.cache/adaptyv-c4/novelty_domain_source_audit'


def main():
    source=json.loads((CACHE/'MERIZO_WEIGHTS_SOURCE.json').read_text())
    folder=CACHE/'ted-pinned/programs/merizo/weights'
    folder.mkdir(exist_ok=True)
    receipts=[]
    for spec in source['files']:
        url='https://raw.githubusercontent.com/'+source['repository']+'/'+source['commit']+'/'+spec['path']
        path=folder/Path(spec['path']).name
        if not path.exists():
            partial=path.with_suffix('.partial')
            with urllib.request.urlopen(url,timeout=40) as src,partial.open('xb') as dst:
                while True:
                    block=src.read(4*1024*1024)
                    if not block: break
                    dst.write(block)
            partial.rename(path)
        assert path.stat().st_size==spec['size']
        git=hashlib.sha1(('blob '+str(spec['size'])+'\0').encode())
        sha=hashlib.sha256()
        with path.open('rb') as stream:
            for block in iter(lambda:stream.read(4*1024*1024),b''):
                git.update(block); sha.update(block)
        assert git.hexdigest()==spec['sha'], 'Weight content differs from pinned repository object'
        receipts.append(dict(url=url,path=str(path),bytes=spec['size'],git_blob_sha1=git.hexdigest(),sha256=sha.hexdigest()))
    receipt=dict(generated_at=datetime.now(timezone.utc).isoformat(),merizo_repository=source['repository'],
        merizo_weights_commit=source['commit'],weights=receipts,new_gpu_or_cloud_cpu_spend_usd=0)
    (HERE/'WEIGHTS.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps(dict(verified_weight_shards=len(receipts),total_bytes=sum(r['bytes'] for r in receipts))))


if __name__=='__main__': main()
