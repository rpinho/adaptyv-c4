"""Launch published domain model with explicit CPU/thread configuration."""
import json
import os
from pathlib import Path
import runpy
import sys

for variable in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[variable]='4'
os.environ['CUDA_VISIBLE_DEVICES']=''
import numpy as np
import torch
torch.set_num_threads(4)
torch.set_num_interop_threads(1)
torch.manual_seed(0)
np.random.seed(0)

cache=Path.home()/'.cache/adaptyv-c4/novelty_domain_source_audit'
method=sys.argv[1]
program=cache/'ted-pinned/programs'/method/('predict_afdb.py' if method=='merizo' else 'get_predictions.py')
sys.path.insert(0,str(program.parent))
if method=='chainsaw':
    from src import constants
    from src.utils.common import get_torch_device
    constants.STRIDE_EXE=str(cache/'native_build/stride/stride')
    if get_torch_device().type!='cpu': raise RuntimeError('CPU-only task; device differs')
elif method!='merizo':
    raise ValueError('Unknown method')
print(json.dumps(dict(runtime_device='cpu',method=method,torch=torch.__version__,threads=torch.get_num_threads(),seed=0)),flush=True)
sys.argv=[str(program),*sys.argv[2:]]
runpy.run_path(str(program),run_name='__main__')
