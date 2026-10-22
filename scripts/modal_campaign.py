"""Pinned C4 pilot on Modal; CPU image build precedes every GPU allocation.

Deploy: .venv-modal/bin/modal deploy scripts/modal_campaign.py
Launch/monitor series: .venv-modal/bin/python scripts/drive_modal_series.py
"""
from pathlib import Path

import modal
if modal.is_local():
    from modal_diffusion_image import image as diffusion_image
else:
    # Image construction is local; the deployed function already carries its image.
    diffusion_image=modal.Image.debian_slim(python_version='3.11')

ROOT = Path(__file__).resolve().parents[1]
REVISION = "b6fb3ae15ab6a702f1156cffc11bcc3323c2fd06"
APP_NAME = "adaptyv-c4-native-control"
VOLUME_NAME = "adaptyv-c4-results"
RUN_NAME = "c4_native_apical_ablation_v3"
LARGE_INIT_RUN = "c4_native_largeinit_v4"
app = modal.App(APP_NAME)
volume = modal.Volume.from_name(VOLUME_NAME, create_if_missing=True)

# Mirrors upstream's container recipe, retaining the previous CUDA 12 family.
# Source, dependencies and AF2 weights are fetched on CPU, never a rented GPU.
image = (
    modal.Image.debian_slim(python_version="3.12")
    .apt_install("git", "curl", "build-essential", "libgomp1")
    .run_commands(
        "git init /opt/bindcraft && cd /opt/bindcraft && "
        "git remote add origin https://github.com/PacesaLab/BindCraft2.git && "
        f"git fetch --depth 1 origin {REVISION} && git checkout --detach FETCH_HEAD"
    )
    .run_commands("python -m pip install -e '/opt/bindcraft[cuda12]'")
    .env({
        "NVIDIA_VISIBLE_DEVICES": "all",
        "NVIDIA_DRIVER_CAPABILITIES": "compute,utility",
        "BINDCRAFT_WEIGHTS": "/opt/bindcraft-weights",
        "PYTHONUNBUFFERED": "1",
        "XLA_PYTHON_CLIENT_PREALLOCATE": "false",
        "OMP_NUM_THREADS": "8",
    })
    .run_commands(
        "python -c \"import nvidia,pathlib; print('\\n'.join(sorted(str(p) "
        "for r in nvidia.__path__ for p in pathlib.Path(r).glob('*/lib'))))\" "
        "> /etc/ld.so.conf.d/bindcraft-cuda.conf && ldconfig",
        "bindcraft fetch-weights",
        "python -m bindcraft.selfcheck cuda12 --shipped-only",
        "python -m pip freeze > /opt/bindcraft-runtime-freeze.txt",
    )
    .add_local_file(ROOT / "results/gpu_campaign_v3/campaign.json", "/opt/c4/campaign.json", copy=True)
    .add_local_file(ROOT / "results/gpu_campaign_v2/human_11qc_full_monomer.pdb", "/opt/c4/target.pdb", copy=True)
    .add_local_file(ROOT / "results/gpu_campaign_v2/mirrored/campaign_metadata.json", "/opt/c4/reference_metadata.json", copy=True)
    .add_local_file(ROOT / "results/modal_campaign_v5/human_11qc_local_frame.pdb", "/opt/c4/local-frame-target.pdb", copy=True)
    .add_local_file(ROOT / "scripts/modal_search_policy.py", "/opt/c4/search_policy.py", copy=True)
    .add_local_file(ROOT / "scripts/frame_guard.py", "/opt/c4/frame_guard.py", copy=True)
    .add_local_file(ROOT / "scripts/patch_frame_guard.py", "/opt/c4/patch_frame_guard.py", copy=True)
    .add_local_file(ROOT / "scripts/validate_frame_guard.py", "/opt/c4/validate_frame_guard.py", copy=True)
    .add_local_file(ROOT / "scripts/frame_gpu_diagnostic.py", "/opt/c4/frame_gpu_diagnostic.py", copy=True)
    .add_local_file(ROOT / "scripts/monitored_design.py", "/opt/c4/monitored_design.py", copy=True)
    .run_commands(
        "python /opt/c4/patch_frame_guard.py --upstream /opt/bindcraft --helper /opt/c4/frame_guard.py",
        "JAX_PLATFORMS=cpu python /opt/c4/validate_frame_guard.py --upstream /opt/bindcraft --out /opt/c4/cpu_frame_validation.json",
        "python -c \"import sys,json;sys.path.insert(0,'/opt/c4');import search_policy as p;"
        "config=p.make_config(json.load(open('/opt/c4/campaign.json')),p.REPAIRED,{'campaign_seed':20261021},'/tmp/preflight-only');"
        "json.dump(config,open('/opt/c4/repaired_preflight.json','w'))\"",
        "JAX_PLATFORMS=cpu python /opt/c4/monitored_design.py --preflight /opt/c4/repaired_preflight.json /opt/c4/prepared_frame_preflight.json",
        "python -c \"import sys,json;sys.path.insert(0,'/opt/c4');import search_policy as p;"
        "config=p.make_config(json.load(open('/opt/c4/campaign.json')),p.FORCED,{'campaign_seed':20261021},'/tmp/forced-preflight-only');"
        "json.dump(config,open('/opt/c4/forced_preflight.json','w'))\"",
        "JAX_PLATFORMS=cpu python /opt/c4/monitored_design.py --preflight /opt/c4/forced_preflight.json /opt/c4/forced_prepared_frame_preflight.json",
        "python -c \"import sys,json;sys.path.insert(0,'/opt/c4');import search_policy as p;"
        "config=p.make_config(json.load(open('/opt/c4/campaign.json')),p.CONFIDENCE,{'campaign_seed':20261021},'/tmp/confidence-preflight-only');"
        "json.dump(config,open('/opt/c4/confidence_preflight.json','w'))\"",
        "JAX_PLATFORMS=cpu python /opt/c4/monitored_design.py --preflight /opt/c4/confidence_preflight.json /opt/c4/confidence_prepared_frame_preflight.json",
        "JAX_PLATFORMS=cpu python /opt/c4/monitored_design.py --seed-preflights /opt/c4/campaign.json /opt/c4",
    )
    .add_local_file(ROOT/'scripts/hard_sequence_refold.py','/opt/c4/hard_sequence_refold.py',copy=True)
    .add_local_file(ROOT/'results/hard_sequence_refold_v15/manifest.json','/opt/c4/hard_refold_manifest.json',copy=True)
    .run_commands('JAX_PLATFORMS=cpu python /opt/c4/hard_sequence_refold.py --preflight /opt/c4/hard_refold_manifest.json /opt/c4/hard_refold_preflight.json')
    .add_local_file(ROOT/'results/apical_context/native_189_383_search_only.pdb','/opt/c4/native-apical-search-only.pdb',copy=True)
    .add_local_file(ROOT/'results/apical_context/PREPARED_INPUT.json','/opt/c4/native-apical-preparation.json',copy=True)
    .add_local_file(ROOT/'scripts/cif2pdb.py','/opt/c4/cif2pdb.py',copy=True)
    .run_commands(
        "python -c \"import sys,json;sys.path.insert(0,'/opt/c4');import search_policy as p;"
        "config=p.make_config(json.load(open('/opt/c4/campaign.json')),p.DOMAIN,{'campaign_seed':20261021},'/tmp/domain-preflight-only');"
        "json.dump(config,open('/opt/c4/domain_preflight.json','w'))\"",
        'JAX_PLATFORMS=cpu python /opt/c4/monitored_design.py --preflight /opt/c4/domain_preflight.json /opt/c4/domain_prepared_frame_preflight.json')
)

if (ROOT/'results/domain_hard_refold_v17/manifest.json').exists():
    image=(image.add_local_file(ROOT/'results/domain_hard_refold_v17/manifest.json','/opt/c4/domain_refold_manifest.json',copy=True)
        .run_commands('JAX_PLATFORMS=cpu python /opt/c4/hard_sequence_refold.py --preflight /opt/c4/domain_refold_manifest.json /opt/c4/domain_refold_preflight.json'))

if (ROOT/'results/diffusion_hard_refold_v19/manifest.json').exists():
    image=(image.add_local_file(ROOT/'results/diffusion_hard_refold_v19/manifest.json','/opt/c4/diffusion_refold_manifest.json',copy=True)
        .run_commands('JAX_PLATFORMS=cpu python /opt/c4/hard_sequence_refold.py --preflight /opt/c4/diffusion_refold_manifest.json /opt/c4/diffusion_refold_preflight.json'))


@app.function(
    image=image,
    gpu=["H100", "A100-80GB"],
    cpu=(4, 4),
    memory=(65536, 65536),
    volumes={"/results": volume},
    timeout=7200,
    startup_timeout=300,
    retries=0,
    max_containers=1,
    min_containers=0,
    scaledown_window=2,
    single_use_containers=True,
)
def run_campaign(run_name: str = RUN_NAME, overrides: dict | None = None):
    import hashlib
    import importlib.util
    import json
    import os
    import re
    import shutil
    import subprocess
    import time

    import jax
    from bindcraft.model_weights import alphafold_parameter_file, alphafold_parameters

    spec = importlib.util.spec_from_file_location('search_policy','/opt/c4/search_policy.py')
    policy = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(policy)
    profile=policy.profile_for_name(run_name)
    monitored=profile in {policy.REPAIRED,policy.INTERFACE8,policy.FORCED,policy.CONFIDENCE,policy.DOMAIN}
    runtime_limit=1800 if monitored else 7200
    overrides = overrides or {}
    if set(overrides) - {"campaign_seed"}:
        raise ValueError("Only a new independent seed is allowed in this series")
    started = time.time()
    run = Path("/results") / run_name
    if run.exists():
        raise RuntimeError("Run directory already exists; inspect it before any duplicate run")
    run.mkdir()

    def save(name, data):
        path = run / name
        tmp = path.with_suffix(path.suffix + ".tmp")
        tmp.write_text(json.dumps(data, indent=2) + "\n")
        tmp.replace(path)

    proc = None
    status = {"status": "preflight", "started": started, "revision": REVISION}
    try:
        devices = jax.devices()
        if not devices or any(d.platform != "gpu" for d in devices):
            raise RuntimeError(f"Refusing CPU fallback: {devices}")
        gpu = subprocess.check_output([
            "nvidia-smi", "--query-gpu=name,memory.total,driver_version", "--format=csv,noheader"
        ], text=True).strip()
        if int(gpu.split(",")[1].strip().split()[0]) < 79000:
            raise RuntimeError(f"GPU has insufficient memory: {gpu}")

        reference = json.loads(Path("/opt/c4/reference_metadata.json").read_text())
        checksums = {}
        params = alphafold_parameters(download=False)
        for old_path, expected in reference["checkpoints"].items():
            name = Path(old_path).name
            if name.startswith("params_"):
                path = alphafold_parameter_file(params, name.removeprefix("params_").removesuffix(".npz"))
            else:
                path = "/opt/bindcraft/bindcraft/weights/proteinmpnn/weights_negative/" + name
            if not path:
                raise RuntimeError(f"Missing checkpoint: {name}")
            actual = hashlib.sha256(Path(path).read_bytes()).hexdigest()
            if actual != expected:
                raise RuntimeError(f"Checkpoint identity mismatch: {name}")
            checksums[name] = actual

        config = json.loads(Path("/opt/c4/campaign.json").read_text())
        config = policy.make_config(config,run_name,overrides,str(run / 'campaign'))
        (run / "campaign.json").write_text(json.dumps(config, indent=2) + "\n")
        shutil.copy(config['targets'][0]['target_path'], run / "target.pdb")
        shutil.copy("/opt/bindcraft-runtime-freeze.txt", run / "runtime-freeze.txt")
        if monitored:
            preflight=('confidence_prepared_frame_preflight.json' if profile==policy.CONFIDENCE else
                       'forced_prepared_frame_preflight.json' if profile==policy.FORCED else 'prepared_frame_preflight.json')
            if profile==policy.DOMAIN:preflight='domain_prepared_frame_preflight.json'
            if profile==policy.REPAIRED and run_name!=profile:
                seed=overrides.get('campaign_seed')
                if seed not in range(20261022,20261027) or run_name!=profile+f'_replicate_{seed-20261020}':
                    raise ValueError('Unknown native seed replicate')
                preflight=f'seed_preflight_{seed}.json'
            shutil.copy('/opt/c4/'+preflight,run/'prepared_frame_preflight.json')
        status.update(status="modeling", gpu=gpu, jax=jax.__version__, checkpoints=checksums,
                      container=os.environ.get("MODAL_TASK_ID"), max_runtime_seconds=runtime_limit)
        save("status.json", status)
        volume.commit()
        print(json.dumps(status), flush=True)
        with (run / "worker.log").open("w") as log:
            command=(["python","/opt/c4/monitored_design.py",str(run/'campaign.json'),str(run/'physical_steps.jsonl')]
                     if monitored else ["bindcraft","design",str(run/'campaign.json')])
            proc = subprocess.Popen(command,
                                    cwd="/opt/bindcraft", stdout=log, stderr=subprocess.STDOUT)
            deadline = started + runtime_limit-60
            while proc.poll() is None:
                if time.time() >= deadline:
                    proc.terminate()
                    try:
                        proc.wait(timeout=30)
                    except subprocess.TimeoutExpired:
                        proc.kill()
                        proc.wait(timeout=10)
                    raise TimeoutError("Finite GPU pilot runtime reached")
                status.update(elapsed_seconds=time.time() - started, worker_pid=proc.pid)
                save("status.json", status)
                volume.commit()
                time.sleep(30)
        status.update(status="completed" if proc.returncode == 0 else "worker_failed", returncode=proc.returncode)
    except BaseException as exc:
        status.update(status="failed", error=repr(exc))
        if proc is not None and proc.poll() is None:
            proc.kill()
            proc.wait(timeout=10)
        raise
    finally:
        status.update(finished=time.time(), elapsed_seconds=time.time() - started,
                      shutdown="function returning; single-use container must terminate")
        save("status.json", status)
        volume.commit()
        print(json.dumps(status), flush=True)
    return status


@app.function(image=image,gpu=['H100','A100-80GB'],cpu=(4,4),memory=(65536,65536),
              volumes={'/results':volume},timeout=1800,startup_timeout=300,retries=0,
              max_containers=1,min_containers=0,scaledown_window=2,single_use_containers=True)
def frame_diagnostic(run_name, overrides):
    import json
    import shutil
    import subprocess
    import time
    if run_name!='c4_frame_guard_diagnostic_v8' or overrides:
        raise ValueError('Only the declared one-call diagnostic is permitted')
    out=Path('/results')/run_name
    if out.exists():raise RuntimeError('Diagnostic output already exists; inspect before any repeated allocation')
    out.mkdir()
    status=dict(status='diagnostic_running',started=time.time())
    (out/'status.json').write_text(json.dumps(status,indent=2)+'\n')
    shutil.copy('/opt/c4/cpu_frame_validation.json',out/'cpu_frame_validation.json')
    shutil.copy('/opt/bindcraft/frame_guard_manifest.json',out/'frame_guard_manifest.json')
    volume.commit()
    try:
        with (out/'worker.log').open('w') as log:
            result=subprocess.run(['python','/opt/c4/frame_gpu_diagnostic.py',str(out)],
                cwd='/opt/bindcraft',stdout=log,stderr=subprocess.STDOUT,timeout=1500)
        status.update(status='completed' if result.returncode==0 else 'diagnostic_failed',returncode=result.returncode)
    except BaseException as exc:
        status.update(status='diagnostic_failed',error=repr(exc))
        raise
    finally:
        status.update(finished=time.time(),elapsed_seconds=time.time()-status['started'])
        (out/'status.json').write_text(json.dumps(status,indent=2)+'\n')
        volume.commit()
    return status


@app.function(image=image,gpu=['H100','A100-80GB'],cpu=(4,4),memory=(65536,65536),
              volumes={'/results':volume},timeout=2100,startup_timeout=300,retries=0,
              max_containers=1,min_containers=0,scaledown_window=2,single_use_containers=True)
def sequence_refold(run_name,overrides):
    import hashlib
    import json
    import shutil
    import subprocess
    import time
    import os
    import jax
    from bindcraft.model_weights import alphafold_parameters,alphafold_parameter_file
    routes={'c4_hard_sequence_refold_v15':('hard_refold_manifest.json','hard_refold_preflight.json')}
    if Path('/opt/c4/domain_refold_manifest.json').exists():
        routes['c4_domain_hard_refold_v17']=('domain_refold_manifest.json','domain_refold_preflight.json')
    if Path('/opt/c4/diffusion_refold_manifest.json').exists():
        routes['c4_diffusion_hard_refold_v19']=('diffusion_refold_manifest.json','diffusion_refold_preflight.json')
    if run_name not in routes or overrides:raise ValueError('Only an image-prepared fixed sequence diagnostic is allowed')
    manifest_name,preflight_name=routes[run_name]
    proof=json.loads(Path('/opt/c4',preflight_name).read_text())
    if proof.get('status')!='passed' or proof['manifest_sha256']!=hashlib.sha256(Path('/opt/c4',manifest_name).read_bytes()).hexdigest() or proof['worker_sha256']!=hashlib.sha256(Path('/opt/c4/hard_sequence_refold.py').read_bytes()).hexdigest():
        raise ValueError('Image proof/manifest/worker identity mismatch')
    if any(d.platform!='gpu' for d in jax.devices()):raise ValueError('CPU fallback refused')
    out=Path('/results')/run_name
    if out.exists():raise RuntimeError('Fresh diagnostic directory required')
    out.mkdir();started=time.time();proc=None
    gpu=subprocess.check_output(['nvidia-smi','--query-gpu=name,memory.total,driver_version','--format=csv,noheader'],text=True).strip()
    status=dict(status='preflight',started=started,gpu=gpu,container=os.environ.get('MODAL_TASK_ID'),worker_limit_seconds=1800)
    def save():
        (out/'status.json').write_text(json.dumps(status,indent=2)+'\n');volume.commit()
    try:
        reference=json.loads(Path('/opt/c4/reference_metadata.json').read_text())['checkpoints']
        params=alphafold_parameters(download=False);checks={}
        for name in ('model_1_ptm','model_2_ptm'):
            filename='params_'+name+'.npz';expected=next(v for k,v in reference.items() if Path(k).name==filename)
            actual=hashlib.sha256(Path(alphafold_parameter_file(params,name)).read_bytes()).hexdigest()
            if actual!=expected:raise ValueError('Held-out checkpoint identity differs')
            checks[filename]=actual
        status.update(checkpoints=checks,status='refolding')
        for source in (manifest_name,preflight_name):
            shutil.copy('/opt/c4/'+source,out/source)
        save()
        with (out/'worker.log').open('w') as log:
            proc=subprocess.Popen(['python','/opt/c4/hard_sequence_refold.py',str(out),'/opt/c4/'+manifest_name],cwd='/opt/bindcraft',stdout=log,stderr=subprocess.STDOUT)
            while proc.poll() is None:
                if time.time()-started>1740:
                    proc.terminate()
                    try:proc.wait(timeout=30)
                    except subprocess.TimeoutExpired:proc.kill();proc.wait(timeout=10)
                    raise TimeoutError('Finite refold allocation reached its stop point')
                status.update(elapsed_seconds=time.time()-started,worker_pid=proc.pid);save();time.sleep(15)
        status.update(status='completed' if proc.returncode==0 else 'worker_failed',returncode=proc.returncode)
    except BaseException as exc:
        status.update(status='failed',error=repr(exc));raise
    finally:
        if proc is not None and proc.poll() is None:proc.kill();proc.wait(timeout=10)
        status.update(finished=time.time(),elapsed_seconds=time.time()-started);save()
    return status


@app.function(image=diffusion_image,gpu=['A100-80GB','A100-40GB'],cpu=(4,4),memory=(65536,65536),
              volumes={'/results':volume},timeout=2100,startup_timeout=300,retries=0,
              max_containers=1,min_containers=0,scaledown_window=2,single_use_containers=True)
def diffusion_generate(run_name,overrides):
    import json
    import os
    import shutil
    import subprocess
    import time
    if run_name!='c4_diffusion_v18' or overrides:raise ValueError('Only the declared single-backbone pilot is allowed')
    out=Path('/results')/run_name
    if out.exists():raise RuntimeError('Fresh diffusion directory required')
    out.mkdir();started=time.time();proc=None
    status=dict(status='starting',started=started,container=os.environ.get('MODAL_TASK_ID'),
                worker_limit_seconds=1800,requested_runtime_upper_usd_per_hour=3.198528)
    def save():
        (out/'status.json').write_text(json.dumps(status,indent=2)+'\n');volume.commit()
    try:
        status['gpu']=subprocess.check_output(['nvidia-smi','--query-gpu=name,memory.total,driver_version','--format=csv,noheader'],text=True).strip()
        if 'A100' not in status['gpu']:raise RuntimeError('A100 compatibility envelope required')
        shutil.copy('/opt/c4/rfdiffusion_preflight.json',out/'CPU_PREFLIGHT.json')
        shutil.copy('/opt/rfd-runtime-freeze.txt',out/'runtime-freeze.txt')
        status['status']='diffusing';save()
        with (out/'worker.log').open('w') as log:
            proc=subprocess.Popen(['/app/RFdiffusion/.venv/bin/python','/opt/c4/diffusion_worker.py',str(out)],
                                  cwd='/app/RFdiffusion',stdout=log,stderr=subprocess.STDOUT)
            while proc.poll() is None:
                if time.time()-started>1740:
                    proc.terminate()
                    try:proc.wait(timeout=30)
                    except subprocess.TimeoutExpired:proc.kill();proc.wait(timeout=10)
                    raise TimeoutError('Finite diffusion allocation reached stop point')
                status.update(elapsed_seconds=time.time()-started,worker_pid=proc.pid);save();time.sleep(15)
        status.update(status='completed' if proc.returncode==0 else 'worker_failed',returncode=proc.returncode)
    except BaseException as exc:
        status.update(status='failed',error=repr(exc));raise
    finally:
        if proc is not None and proc.poll() is None:proc.kill();proc.wait(timeout=10)
        status.update(finished=time.time(),elapsed_seconds=time.time()-started);save()
    return status


@app.function(image=diffusion_image,cpu=1,memory=2048,timeout=120,startup_timeout=300,
              retries=0,max_containers=1,min_containers=0,scaledown_window=2,single_use_containers=True)
def diffusion_runtime_check():
    import sys
    import json
    import subprocess
    if sys.version_info<(3,10):raise RuntimeError('Modal runner Python is too old')
    child=json.loads(subprocess.check_output(['/app/RFdiffusion/.venv/bin/python','-c',
        'import sys,json,torch;print(json.dumps(dict(python=list(sys.version_info[:2]),torch=torch.__version__,cuda_available=torch.cuda.is_available())))'],text=True))
    if child['python']!=[3,9] or child['cuda_available']:raise RuntimeError('CPU bootstrap check differs')
    return dict(status='passed',modal_runner_python=list(sys.version_info[:2]),runner_executable=sys.executable,
                rf_runtime=child,scope='Actual serverless CPU bootstrap; no GPU requested')
