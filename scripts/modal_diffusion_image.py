"""Separate pinned official-derived RFdiffusion runtime, built entirely on CPU."""
from pathlib import Path
import modal

ROOT=Path(__file__).resolve().parents[1]
DIGEST='sha256:fa8bbf2fd51291c861369adbf246c4909c5757763d1478cafacadcc76c4584c1'
REVISION='86507b6538f51fce57b5a72477165f03999ed7ae'
PYTHON='/app/RFdiffusion/.venv/bin/python'
image=(modal.Image.from_registry('rosettacommons/rfdiffusion@'+DIGEST,add_python='3.11')
    .entrypoint([])
    .env({'PATH':'/usr/local/bin:/usr/local/sbin:/usr/local/nvidia/bin:/usr/local/cuda/bin:/usr/sbin:/usr/bin:/sbin:/bin'})
    .run_commands("python -c 'import sys; assert sys.version_info[:2] == (3,11)' ")
    .apt_install('git','curl','ca-certificates')
    .run_commands(
        'git init /app/RFdiffusion && cd /app/RFdiffusion && git remote add c4 https://github.com/RosettaCommons/RFdiffusion.git && '
        f'git fetch --depth 1 c4 {REVISION} && git checkout --force --detach FETCH_HEAD',
        PYTHON+' -m ensurepip',
        PYTHON+' -m pip install numpy==1.26.4 biopython==1.85',
        "mkdir -p /opt/rfd-models && curl --fail --location --retry 2 https://files.ipd.uw.edu/pub/RFdiffusion/e29311f6f1bf1af907f9ef9f44b8328b/Complex_base_ckpt.pt -o /opt/rfd-models/Complex_base_ckpt.pt",
        "test $(stat -c %s /opt/rfd-models/Complex_base_ckpt.pt) -eq 483619179 && sha256sum /opt/rfd-models/Complex_base_ckpt.pt",
        PYTHON+' -m pip freeze > /opt/rfd-runtime-freeze.txt')
    .env({'DGLBACKEND':'pytorch','PYTHONUNBUFFERED':'1','OMP_NUM_THREADS':'4','PYTHONPATH':'/app/RFdiffusion'})
    .add_local_file(ROOT/'results/apical_context/native_189_383_search_only.pdb','/opt/c4/native-apical-search-only.pdb',copy=True)
    .add_local_file(ROOT/'scripts/rfdiffusion_preflight.py','/opt/c4/rfdiffusion_preflight.py',copy=True)
    .add_local_file(ROOT/'scripts/diffusion_worker.py','/opt/c4/diffusion_worker.py',copy=True)
    .run_commands('cd /app/RFdiffusion && '+PYTHON+' /opt/c4/rfdiffusion_preflight.py /opt/c4/rfdiffusion_preflight.json'))
