"""Compile unchanged STRIDE/UniDoc algorithms for this Mac; keep source receipts."""
import hashlib
import json
from pathlib import Path
import re
import shutil
import subprocess
import tarfile

HERE=Path(__file__).resolve().parent
CACHE=Path.home()/'.cache/adaptyv-c4/novelty_domain_source_audit'
SOURCE=CACHE/'ted-pinned'
BUILD=CACHE/'native_build'


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    if (HERE/'NATIVE_BUILD.json').exists(): raise FileExistsError('Build already frozen')
    BUILD.mkdir(exist_ok=True)
    compat=BUILD/'compat'; compat.mkdir(exist_ok=True)
    # Linux malloc.h is redundant with stdlib.h for these published sources.
    # An include-path shim preserves the algorithm sources byte-for-byte.
    (compat/'malloc.h').write_text('#include <stdlib.h>\n')
    stride=BUILD/'stride'; stride.mkdir(exist_ok=True)
    archive=SOURCE/'programs/chainsaw/stride/stride.tgz'
    with tarfile.open(archive,'r:gz') as tar:
        for member in tar:
            name=Path(member.name).name
            if not member.isfile() or not (name.endswith(('.c','.h')) or name=='Makefile'): continue
            path=stride/name
            with tar.extractfile(member) as src, path.open('xb') as dst: shutil.copyfileobj(src,dst)
    makefile=(stride/'Makefile').read_text()
    files=re.search(r'^SOURCE\s*=\s*(.+)$',makefile,re.M).group(1).split()
    if not all((stride/f).exists() for f in files): raise ValueError('Missing STRIDE source')
    # GNU89/fcommon match the historical C source conventions; no algorithm patch.
    stride_cmd=['/usr/bin/clang','-O2','-std=gnu89','-fcommon','-I',str(compat),
                *[str(stride/f) for f in files],'-lm','-o',str(stride/'stride')]
    unidoc=BUILD/'unidoc'; unidoc.mkdir(exist_ok=True)
    bindir=unidoc/'bin'; bindir.mkdir(exist_ok=True)
    unidoc_src=SOURCE/'programs/unidoc/src/UniDoc_struct.cpp'
    unidoc_cmd=['/usr/bin/clang++','-O3','-std=c++11','-I',str(compat),str(unidoc_src),
                '-o',str(bindir/'UniDoc_struct')]
    for name,cmd in [('stride',stride_cmd),('unidoc',unidoc_cmd)]:
        with (HERE/(name+'_build.log')).open('w') as log:
            subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT,check=True,timeout=180)
    shutil.copy2(stride/'stride',bindir/'stride')
    wrapper=SOURCE/'scripts/Run_UniDoc_from_scratch_structure_afdb.py'
    shutil.copy2(wrapper,unidoc/wrapper.name)
    receipt=dict(compiler=subprocess.check_output(['/usr/bin/clang','--version'],text=True).strip(),
        commands=[stride_cmd,unidoc_cmd],compatibility_shim=dict(path=str(compat/'malloc.h'),sha256=digest(compat/'malloc.h')),
        stride_source_archive_sha256=digest(archive),stride_source_files={f:digest(stride/f) for f in files},
        unidoc_source_files={p.name:digest(p) for p in unidoc_src.parent.iterdir() if p.suffix in ('.cpp','.h')},
        binaries={str(p):digest(p) for p in [stride/'stride',bindir/'stride',bindir/'UniDoc_struct']},
        wrapper=dict(path=str(unidoc/wrapper.name),sha256=digest(wrapper)),source_algorithms_modified=False,
        new_gpu_or_cloud_cpu_spend_usd=0)
    (HERE/'NATIVE_BUILD.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps(dict(native_binaries_built=2,algorithm_source_modified=False)))


if __name__=='__main__': main()
