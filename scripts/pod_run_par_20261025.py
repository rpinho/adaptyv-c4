#!/usr/bin/env python3
"""N-capable AF2M screen runner for the RunPod lane (Astra ruling 2026-10-25 00:54).

Same fold settings as /workspace/pod_run.sh (colabfold alphafold2_multimer_v3,
num-recycle 3, num-models 5, rank multimer, msa-mode mmseqs2_uniref_env) and the
same per-design skip-on-valid-scores idempotency. Adds:
  * a global concurrency cap enforced by atomic claim, so a design is folded once
  * the provenance fields Astra requires on every attempt, serial included
  * a batch record (concurrency, input hash, timestamps, exit codes)
  * a hard no-new-admission boundary at 03:00 ET, checked before each claim
  * terminal error.json is never retried; failed attempts kept in failures/

Claims are keyed on the DESIGN ONLY (never on the batch), so two batches on the
same pod can never fold the same sequence twice. A claim whose owner is dead and
whose record is older than STALE_S is taken over and the takeover is recorded.

usage: pod_run_par_20261025.py <csv> <outdir> <concurrency> [--limit N] [--dry-run]
"""
import argparse
import csv
import glob
import hashlib
import json
import os
import pathlib
import shutil
import subprocess
import sys
import threading
import time

CUTOFF_ET = 1792911600   # 2026-10-25 03:00:00 America/New_York
STALE_S = 600


def admit(now, cutoff=CUTOFF_ET):
    """No design may be admitted at or after the cutoff epoch."""
    return now < cutoff


def claim_key(name):
    """Design-only key: identical across batches and pods for the same design."""
    return hashlib.sha256(f'design|{name}'.encode()).hexdigest()[:16]


def pid_alive(pid):
    try:
        os.kill(int(pid), 0)
        return True
    except (ProcessLookupError, ValueError):
        return False
    except PermissionError:
        return True


def terminal(d):
    """A design is terminal if it has a usable score or a non-empty error.json."""
    d = pathlib.Path(d)
    sp = d / 'scores.json'
    if sp.exists():
        try:
            if isinstance(json.load(open(sp)).get('iptm'), (int, float)):
                return 'scored'
        except Exception:
            pass
    ep = d / 'error.json'
    if ep.exists() and ep.stat().st_size > 0:
        return 'error'
    return None


def sha_file(p):
    try:
        return hashlib.sha256(pathlib.Path(p).read_bytes()).hexdigest()
    except Exception:
        return None


def settings_block():
    return {
        'predictor': 'colabfold_alphafold2_multimer_v3',
        'model_type': 'alphafold2_multimer_v3',
        'num_recycle': 3,
        'num_models': 5,
        'rank': 'multimer',
        'msa_mode': 'mmseqs2_uniref_env',
        'colabfold': os.environ.get('COLABFOLD_PIN', '1.5.5'),
        'python': os.environ.get('AF2_PYTHON_PIN', '3.10'),
        'jax': os.environ.get('JAX_PIN', '0.4.23'),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('csv')
    ap.add_argument('outdir')
    ap.add_argument('concurrency', nargs='?', type=int, default=1)
    ap.add_argument('--limit', type=int, default=None)
    ap.add_argument('--dry-run', action='store_true')
    a = ap.parse_args()

    csv_path = pathlib.Path(a.csv)
    out_root = pathlib.Path(a.outdir)
    cap = a.concurrency
    if cap < 1:
        print('ERROR concurrency must be >= 1', file=sys.stderr)
        return 2

    pod_id = os.environ.get('RUNPOD_POD_ID', 'unknown')
    batch_id = os.environ.get('R3_BATCH_ID',
                              time.strftime('b%Y%m%dT%H%M%SZ', time.gmtime()))
    run_mode = os.environ.get('R3_RUN_MODE', 'screen')
    validation = os.environ.get('R3_VALIDATION', '0') == '1'
    val_ref = os.environ.get('R3_VALIDATION_REF', '')
    cutoff = int(os.environ.get('R3_CUTOFF_EPOCH', CUTOFF_ET))
    settings = settings_block()
    settings_sha = hashlib.sha256(
        json.dumps(settings, sort_keys=True).encode()).hexdigest()

    rows = list(csv.DictReader(open(csv_path)))
    if a.limit is not None:
        rows = rows[:a.limit]
    input_sha = sha_file(csv_path)
    claims = out_root / '.claims'
    recs = out_root / 'records'
    failures = out_root / 'failures'
    qd = out_root / 'queue'
    if not a.dry_run:
        for d in (claims, recs, failures, qd):
            d.mkdir(parents=True, exist_ok=True)

    batch_rec = {
        'batch_id': batch_id, 'pod_id': pod_id, 'run_mode': run_mode,
        'input_csv': str(csv_path), 'input_sha256': input_sha,
        'cohort_sha256': input_sha,
        'configured_concurrency': cap, 'validation_control': validation,
        'validation_reference': val_ref, 'settings': settings,
        'settings_sha256': settings_sha, 'cutoff_epoch': cutoff,
        'queue_index': {r['name']: i for i, r in enumerate(rows)},
        'pass_no': {r['name']: r.get('pass_no') for r in rows},
        'designs': [r['name'] for r in rows],
        'started_epoch': round(time.time(), 1),
    }
    if not a.dry_run:
        (qd / f'queue_{batch_id}.json').write_text(json.dumps(batch_rec, indent=2))

    lock = threading.Lock()
    print_lock = threading.Lock()
    claimed = set()
    done = []
    observed = {'peak': 0, 'active': 0}
    started = time.time()

    def log(msg):
        with print_lock:
            print(msg, flush=True)

    def claim_next():
        with lock:
            now = time.time()
            for r in rows:
                name = r['name']
                if name in claimed:
                    continue
                state = terminal(out_root / name)
                if state:
                    claimed.add(name)
                    return ('skip_' + state, r)
                key = claim_key(name)
                cf = claims / key
                if cf.exists():
                    try:
                        c = json.loads(cf.read_text())
                    except Exception:
                        c = {}
                    age = now - float(c.get('t', 0))
                    if pid_alive(c.get('pid')) or age < STALE_S:
                        claimed.add(name)
                        return ('claimed_elsewhere', r)
                    cf.write_text(json.dumps(
                        {'name': name, 'batch_id': batch_id, 'pid': os.getpid(),
                         't': now, 'takeover_from': c, 'takeover': True}, indent=2))
                    claimed.add(name)
                    return ('work', r)
                if not admit(now, cutoff):
                    return ('cutoff', r)
                cf.write_text(json.dumps({'name': name, 'batch_id': batch_id,
                                          'pid': os.getpid(), 't': now}, indent=2))
                claimed.add(name)
                return ('work', r)
            return (None, None)

    def fold(r, slot):
        name = r['name']
        od = out_root / name
        od.mkdir(parents=True, exist_ok=True)
        attempt_id = f'{batch_id}-{name}-slot{slot}'
        work = pathlib.Path(f'/tmp/{attempt_id}')
        work.mkdir(parents=True, exist_ok=True)
        fa = work / f'{name}.fasta'
        fa.write_text(f'>{name}\n{r["binder"]}:{r["target"]}\n')
        cf = work / 'out'
        cf.mkdir(exist_ok=True)
        cmd = ['colabfold_batch', str(fa), str(cf),
               '--model-type', 'alphafold2_multimer_v3', '--num-recycle', '3',
               '--num-models', '5', '--rank', 'multimer',
               '--msa-mode', 'mmseqs2_uniref_env']
        with lock:
            observed['active'] += 1
            observed['peak'] = max(observed['peak'], observed['active'])
        t0 = time.time()
        try:
            p = subprocess.run(cmd, capture_output=True, text=True)
            rc = p.returncode
            err = (p.stderr or '')[-400:]
        except Exception as e:
            rc = -1
            err = f'launch_error {e!r}'
        t1 = time.time()
        with lock:
            observed['active'] -= 1
        sc = sorted(glob.glob(str(cf / '*_scores_rank_001_*.json')))
        pdb = sorted(glob.glob(str(cf / '*_rank_001_*.pdb')))
        a3m = sorted(glob.glob(str(cf / '*.a3m')))
        base = {
            'design_id': name, 'name': name, 'pod_id': pod_id,
            'attempt_id': attempt_id, 'batch_id': batch_id, 'slot': slot,
            'worker_id': f'{pod_id}-w{slot}', 'run_mode': run_mode,
            'configured_concurrency': cap,
            'observed_concurrency_peak': observed['peak'],
            'queue_index': batch_rec['queue_index'].get(name),
            'pass_no': r.get('pass_no'),
            'input_sha256': input_sha, 'source': r.get('source'),
            'settings_sha256': settings_sha,
            'input_fasta_sha256': sha_file(fa),
            'msa_sha256': sha_file(a3m[0]) if a3m else None,
            'validation_control': validation, 'validation_reference': val_ref,
            'started_epoch': round(t0, 1), 'ended_epoch': round(t1, 1),
            'wall_s': round(t1 - t0, 1), 'rc': rc,
            'seed': r.get('seed'), 'predictor': settings['predictor'],
            'msa_mode': settings['msa_mode'], 'num_recycle': 3, 'num_models': 5,
        }
        if rc != 0 or not sc:
            json.dump({**base, 'status': 'INCOMPLETE', 'error': err},
                      open(recs / f'{attempt_id}.json', 'w'), indent=2)
            json.dump({**base, 'status': 'INCOMPLETE',
                       'error': err or 'no rank_001 scores'},
                      open(od / 'error.json', 'w'), indent=2)
            fp = failures / f'{name}.json'
            prev = []
            if fp.exists():
                try:
                    prev = json.loads(fp.read_text()).get('attempts', [])
                except Exception:
                    prev = []
            prev.append({'attempt_id': attempt_id, 'batch_id': batch_id,
                         'rc': rc, 'error': err, 'wall_s': round(t1 - t0, 1),
                         'status': 'INCOMPLETE'})
            fp.write_text(json.dumps({'name': name, 'terminal': True,
                                      'attempts': prev}, indent=2))
            log(f'[{name}] FAILED rc={rc} {err[-160:]}')
            return
        d = json.load(open(sc[0]))
        pl = d.get('plddt', [])
        res = {'name': name, 'source': r.get('source'),
               'binder_len': len(r['binder']), 'target_len': len(r['target']),
               'iptm': d.get('iptm'), 'ptm': d.get('ptm'),
               'mean_plddt': sum(pl) / len(pl) if pl else None,
               'model': os.path.basename(sc[0]),
               'predictor': settings['predictor'],
               'num_recycle': 3, 'num_models': 5,
               'subprocess_wall_s': round(t1 - t0, 1),
               'gpu': os.environ.get('RUNPOD_GPU_NAME', 'A40'),
               'run_host': 'runpod', 'batch_id': batch_id,
               'attempt_id': attempt_id, 'pod_id': pod_id,
               'worker_id': base['worker_id'], 'run_mode': run_mode,
               'configured_concurrency': cap,
               'observed_concurrency_peak': observed['peak'],
               'queue_index': base['queue_index'], 'pass_no': base['pass_no'],
               'input_sha256': input_sha, 'settings_sha256': settings_sha,
               'validation_control': validation,
               'validation_reference': val_ref}
        shutil.copy2(sc[0], od / 'scores_rank_001.json')
        if pdb:
            shutil.copy2(pdb[0], od / 'complex_rank_001.pdb')
        models = {}
        for f in sorted(glob.glob(str(cf / '*_scores_rank_*_alphafold2_multimer_v3_model_*.json'))):
            b = os.path.basename(f)
            shutil.copy2(f, od / b)
            try:
                md = json.load(open(f))
                models[b] = {'iptm': md.get('iptm'), 'ptm': md.get('ptm')}
            except Exception:
                models[b] = {'error': 'unreadable'}
        for f in sorted(glob.glob(str(cf / '*_rank_*_alphafold2_multimer_v3_model_*.pdb'))):
            shutil.copy2(f, od / os.path.basename(f))
        res['models_complete'] = len(models) == 5
        res['models'] = models
        json.dump(res, open(od / 'scores.json', 'w'), indent=2)
        art = {}
        for f in sorted(od.iterdir()):
            if f.is_file():
                art[f.name] = sha_file(f)
        json.dump({**base, 'status': 'COMPLETE', 'iptm': res['iptm'],
                   'models_complete': res['models_complete'], 'artifacts': art},
                  open(recs / f'{attempt_id}.json', 'w'), indent=2)
        log(f'[{name}] iptm={res["iptm"]:.3f} plddt={res["mean_plddt"]:.1f} '
            f'slot={slot} wall={t1 - t0:.0f}s')

    if a.dry_run:
        print(json.dumps({'batch_id': batch_id, 'pod_id': pod_id, 'designs': len(rows),
                          'configured_concurrency': cap, 'cutoff_epoch': cutoff,
                          'admitted_now': admit(time.time(), cutoff),
                          'input_sha256': input_sha,
                          'settings_sha256': settings_sha,
                          'claim_key_example': claim_key(rows[0]['name']) if rows else None},
                         indent=2))
        return 0

    def worker(slot):
        while True:
            kind, r = claim_next()
            if kind is None:
                return
            if kind == 'cutoff':
                log(f'CUTOFF_NO_NEW_ADMISSION {r["name"]}')
                return
            if kind.startswith('skip_'):
                done.append(r['name'])
                continue
            if kind == 'claimed_elsewhere':
                done.append(r['name'])
                continue
            fold(r, slot)
            done.append(r['name'])

    log(f'running {len(rows)} design(s) cap={cap} batch={batch_id} pod={pod_id} '
        f'mode={run_mode} cutoff={cutoff} settings={settings_sha[:12]}')
    threads = [threading.Thread(target=worker, args=(i,), daemon=True) for i in range(cap)]
    for t in threads:
        t.start()
    for t in threads:
        t.join()
    end = time.time()
    summary = {**{k: v for k, v in batch_rec.items() if k not in ('queue_index', 'pass_no')},
               'processed': len(done),
               'observed_concurrency_peak': observed['peak'],
               'ended_epoch': round(end, 1),
               'wall_s': round(end - started, 1),
               'cutoff_reached': not admit(end, cutoff)}
    json.dump(summary, open(recs / f'batch_{batch_id}.json', 'w'), indent=2)
    print('ALL_DONE', flush=True)
    return 0


if __name__ == '__main__':
    sys.exit(main())
