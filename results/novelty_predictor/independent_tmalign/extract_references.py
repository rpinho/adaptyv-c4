"""Export each numeric database key separately; reject filename collisions."""
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
from pathlib import Path
import subprocess

HERE = Path(__file__).resolve().parent
FOLDSEEK = '/Users/ricardo/.local/opt/foldseek/foldseek/bin/foldseek'


def main():
    inputs = json.loads((HERE/'PAIRS.json').read_text())
    prefix = inputs['database']['prefix']
    output = HERE/'references_by_key'; output.mkdir(exist_ok=True)
    work = Path.home()/'.cache/adaptyv-c4/standalone_tm_reference_export'
    work.mkdir(parents=True,exist_ok=True)

    def extract(pair):
        key, target = pair['database_key'], pair['hit']['target']
        folder = work/str(key); folder.mkdir(exist_ok=True)
        ids = folder/'id.txt'; ids.write_text(str(key)+'\n')
        commands = []
        with (folder/'export.log').open('w') as log:
            for suffix in ('','_h','_ca'):
                cmd=[FOLDSEEK,'createsubdb',str(ids),prefix+suffix,str(folder/('db'+suffix))]
                subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT,check=True,timeout=30)
                commands.append(cmd)
            header=(folder/'db_h').read_bytes().rstrip(b'\x00\n').decode()
            if header.split()[0] != target: raise ValueError('Header/key mapping differs from requested target')
            path=output/f'{key}.pdb'
            cmd=[FOLDSEEK,'convert2pdb',str(folder/'db'),str(path),'--pdb-output-mode','0','--threads','1']
            subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT,check=True,timeout=30)
            commands.append(cmd)
        count=sum(line.startswith('ATOM') and line[12:16].strip()=='CA' for line in path.read_text().splitlines())
        if count != pair['hit']['tlen']: raise ValueError('Exported reference length differs')
        return dict(database_key=key,target=target,path=str(path.resolve()),CA_count=count,
                    header=header,sha256=hashlib.sha256(path.read_bytes()).hexdigest(),commands=commands)

    with ThreadPoolExecutor(max_workers=4) as pool:
        records=list(pool.map(extract,inputs['pairs']))
    (HERE/'REFERENCE_MANIFEST.json').write_text(json.dumps(dict(status='verified',records=records,
        filename_policy='Numeric database keys, one entry per export; no interpretation of chain-letter filenames',
        rejected_batch_export='convert2pdb mode1 produced filename collisions/merged records; CA count guard rejected before scoring'),indent=2)+'\n')
    print('Verified27 numeric-key reference exports: exact headers and expected CA counts.')


if __name__ == '__main__':main()
