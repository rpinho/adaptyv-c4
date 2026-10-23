"""Run the registered native-context scan through the unmodified upstream metric.

Observational hooks save scored poses and check identities/backbones. They do
not change scores, packing, minimization, masks or thresholds. No cloud calls.
"""
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import subprocess
import sys
import time
import traceback

ROOT = Path(__file__).resolve().parents[2]
RUN = ROOT / 'results/rfd4_his_scan_20261023'
PYTHON = Path.home() / '.cache/adaptyv-c4/his_scan_cpu/venv/bin/python'


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def dump(path, value):
    path = Path(path)
    temporary = path.with_suffix(path.suffix + '.tmp')
    temporary.write_text(json.dumps(value, indent=2) + '\n')
    temporary.replace(path)


def worker(index):
    plan = json.loads((RUN / 'PLAN.json').read_text())
    job = plan['jobs'][index]
    out = RUN / 'sites' / ('site_%02d' % index)
    out.mkdir(parents=True, exist_ok=True)
    if (out / 'RESULTS.json').exists():
        raise FileExistsError('Preserve prior measured output')
    start = time.monotonic()
    receipt = dict(state='starting', job=job, plan_sha256=sha(RUN/'PLAN.json'),
                   script_sha256=sha(__file__), device='cpu', new_gpu_spend_usd=0)
    dump(out/'STATUS.json', receipt)
    try:
        os.environ['CLUSTER'] = 'local'
        os.environ['OMP_NUM_THREADS'] = str(plan['threads'])
        os.environ['MKL_NUM_THREADS'] = str(plan['threads'])
        import importlib.metadata as metadata
        import random
        import numpy as np
        import torch
        from atomworks.io import parse
        from atomworks.io.config import ParseConfig
        from biotite.structure import AtomArrayStack, get_residue_starts
        from tmol.io import atom_array_from_pose_stack, file_from_atom_array
        from rfproteina.metrics.physics import tmol_block_pair_metric as upstream

        torch.set_num_threads(plan['threads'])
        torch.set_num_interop_threads(1)
        random.seed(plan['seed'] + index)
        np.random.seed(plan['seed'] + index)
        torch.manual_seed(plan['seed'] + index)
        source_fingerprints = json.loads((RUN/'UPSTREAM_IMPLEMENTATION.json').read_text())
        assert sha(upstream.__file__) == source_fingerprints['installed_sources_byte_identical']['metrics/physics/tmol_block_pair_metric.py']
        context = json.loads((RUN/'context20/RESULTS.json').read_text())
        record = next(r for r in context['designs'] if r['source_id'] == job['source_id'])
        assert record['omitted_context_rigid_pass'] and record['full_native_staging_pass']
        source = RUN/'context20/restored_staging'/(record['fixture_name']+'.pdb')
        assert sha(source) == record['staging']['output_sha256']
        array = parse(str(source), config=ParseConfig(model=1, build_assembly=None,
                                                     add_missing_atoms=False))['asym_unit'][0]

        def atom_key(a, i):
            return str(a.chain_id[i]), int(a.res_id[i]), str(a.ins_code[i]), str(a.atom_name[i])

        def backbone(a):
            return {atom_key(a, i): a.coord[i].copy() for i in range(len(a))
                    if str(a.atom_name[i]) in ('N', 'CA', 'C', 'O')}

        def identities(a):
            return {(str(a.chain_id[i]), int(a.res_id[i]), str(a.ins_code[i])):str(a.res_name[i])
                    for i in get_residue_starts(a)}

        original_bb = backbone(array)
        original_identities = identities(array)
        assert original_identities[('C', job['position'], '')] == 'HIS'
        assert len([x for x in original_identities if x[0] in ('A', 'B')]) == 1278
        legs = []
        original_score = upstream.score_block_pair

        def observe_score(*args, **kwargs):
            pose = kwargs.get('pose_stack') if 'pose_stack' in kwargs else args[0]
            current = atom_array_from_pose_stack(pose)
            if isinstance(current, AtomArrayStack):
                current = current[0]
            bb = backbone(current)
            ids = identities(current)
            assert set(bb) == set(original_bb), 'Backbone atom inventory changed'
            maximum = max(float(np.linalg.norm(bb[k]-original_bb[k])) for k in bb)
            assert maximum <= .002, 'Backbone moved in a fixed-backbone scan'
            label = 'WT' if not legs else plan['mutate_to'][len(legs)-1]
            expected = dict(original_identities)
            expected[('C', job['position'], '')] = 'HIS' if label == 'WT' else label
            assert ids == expected, 'Unexpected mutation or lost residue'
            saved = Path.home()/'.cache/adaptyv-c4/his_scan_cpu/scored_poses'/('site_%02d'%index)
            saved.mkdir(parents=True, exist_ok=True)
            pdb = saved/(label+'.pdb')
            file_from_atom_array(current, path=pdb)
            result = original_score(*args, **kwargs)
            legs.append(dict(leg=label, result=result, saved_pdb=str(pdb), sha256=sha(pdb),
                             maximum_backbone_displacement_A=maximum, identities_verified=True,
                             retained_native_target_residues=1278))
            dump(out/'LEGS.json', legs)
            print('scored', label, 'site', index, flush=True)
            return result

        upstream.score_block_pair = observe_score
        metric = upstream.TMolMutationDdgMetric(
            mutation_selection=[job['mutation_selection']], mutate_to=plan['mutate_to'],
            region_b=job['region_b'], device='cpu', shell_radius=plan['shell_radius_A'],
            num_repeats=plan['num_repeats'], reuse_context=False,
            raise_on_dropped_residues=True)
        receipt.update(state='scoring', source_pdb_sha256=sha(source),
                       versions={x:metadata.version(x) for x in ('torch','tmol','atomworks','rfd4-proteina')},
                       cuda_available=torch.cuda.is_available(), native_target_residues=1278)
        dump(out/'STATUS.json', receipt)
        metric.setup()
        result = metric.compute(array)
        metric.teardown()
        assert len(legs) == 1 + len(plan['mutate_to'])
        assert all(math.isfinite(value) for value in result.values())
        for mutation, leg in zip(plan['mutate_to'], legs[1:], strict=True):
            prefix = 'C'+str(job['position'])
            assert abs(result[f'interface_ddg_delta__{prefix}__{mutation}'] -
                       (leg['result']['interface_ddg']-legs[0]['result']['interface_ddg'])) < 1e-8
        receipt.update(state='complete', results=result, scored_legs=legs,
                       elapsed_seconds=time.monotonic()-start,
                       pH_selectivity_measured=False, rigorous_binding_ddg_measured=False)
        dump(out/'RESULTS.json', receipt)
        dump(out/'STATUS.json', receipt)
        print(json.dumps(dict(state='complete', site=index, elapsed_seconds=receipt['elapsed_seconds'])), flush=True)
    except Exception as error:
        receipt.update(state='failed', error=repr(error), traceback=traceback.format_exc(),
                       elapsed_seconds=time.monotonic()-start)
        dump(out/'STATUS.json', receipt)
        raise


def execute():
    plan = json.loads((RUN/'PLAN.json').read_text())
    for path, fingerprint in plan['sources'].items():
        assert sha(ROOT/path) == fingerprint
    states = []
    for index, job in enumerate(plan['jobs']):
        out = RUN/'sites'/('site_%02d'%index)
        out.mkdir(parents=True, exist_ok=True)
        if (out/'RESULTS.json').exists():
            result = json.loads((out/'RESULTS.json').read_text())
            assert result['plan_sha256'] == sha(RUN/'PLAN.json') and result['job'] == job
            assert result['state'] == 'complete'
            states.append(dict(index=index, state='complete', reused=True))
            continue
        with (out/'worker.log').open('a') as log:
            try:
                command = subprocess.run([str(PYTHON), str(Path(__file__).resolve()), '--job-index', str(index)],
                                         cwd='/Users/ricardo/dev/RFD4-Proteina', stdout=log, stderr=subprocess.STDOUT,
                                         timeout=plan['site_timeout_seconds'])
                state = 'complete' if command.returncode == 0 and (out/'RESULTS.json').exists() else 'failed'
            except subprocess.TimeoutExpired:
                state = 'timed_out'
        states.append(dict(index=index, state=state, reused=False))
        dump(RUN/'SCAN_STATUS.json', dict(state='running', planned_sites=len(plan['jobs']),
                                        attempted=states, completed=sum(s['state']=='complete' for s in states),
                                        new_gpu_spend_usd=0))
        if index == 0 and state != 'complete':
            raise RuntimeError('First native-context site failed; diagnose before scaling')
    dump(RUN/'SCAN_STATUS.json', dict(state='finished', planned_sites=len(plan['jobs']),
                                    attempted=states, completed=sum(s['state']=='complete' for s in states),
                                    new_gpu_spend_usd=0))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--job-index', type=int)
    parser.add_argument('--execute', action='store_true')
    args = parser.parse_args()
    if args.job_index is not None:
        worker(args.job_index)
    elif args.execute:
        execute()
    else:
        print('CPU-only; registered plan:', RUN/'PLAN.json')
