"""Finite CPU-only annotation of existing sequences; no design or official typing."""
import collections
import csv
from datetime import datetime, timezone
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import random
import sys
import time

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
CACHE = Path.home() / '.cache/adaptyv-c4/novelty_antibody_scope'
sys.path.insert(0, str(CACHE / 'runtime'))
os.environ.update(CUDA_VISIBLE_DEVICES='', OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1',
    MKL_NUM_THREADS='1', VECLIB_MAXIMUM_THREADS='1')
import numpy as np
import torch
from anarcii import Anarcii
from Bio import SeqIO


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def save(path, value):
    path.write_text(json.dumps(value, indent=2) + '\n')


def annotated_sequence(record):
    return ''.join(aa for _, aa in (record.get('numbering') or []) if aa != '-')


def main():
    if (HERE / 'STATUS.json').exists():
        raise FileExistsError('This finite audit has already been started')
    state = dict(state='prepared', pid=os.getpid(), started_at=datetime.now(timezone.utc).isoformat(),
        cpu=True, numerical_threads=1, mode='accuracy', seq_type='antibody', seed=0,
        generator=False, organizer_typing_established=False, new_gpu_or_cloud_cpu_spend_usd=0)
    save(HERE / 'STATUS.json', state)
    try:
        random.seed(0); np.random.seed(0); torch.manual_seed(0)
        torch.set_num_threads(1); torch.set_num_interop_threads(1)
        runtime = CACHE / 'runtime/anarcii'
        receipt = dict(versions={m: importlib.metadata.version(m) for m in
            ('anarcii', 'torch', 'numpy', 'biopython', 'gemmi', 'msgpack', 'psutil', 'sortedcontainers')},
            wheel_sha256=sha(CACHE / 'anarcii-2.0.8-py3-none-any.whl'),
            weights={str(p.relative_to(runtime)): sha(p) for p in runtime.rglob('*.pt')},
            source_files={str(p.relative_to(runtime)): sha(p) for p in runtime.rglob('*.py')},
            classic_anarci_pinned_source='79f6c575056dedef86cb8f405ebb039197923eec',
            classic_anarci_not_run='Does not ship HMM/germlines; dynamic source build not executed',
            organizer_equivalence_established=False)
        save(HERE / 'RUNTIME.json', receipt)
        model = Anarcii(seq_type='antibody', mode='accuracy', batch_size=32, cpu=True, ncpu=1)
        assert str(model.device) == 'cpu'
        fasta = CACHE / 'ANARCI/Example_scripts_and_sequences/12e8.fasta'
        supplied = {r.id: str(r.seq) for r in SeqIO.parse(str(fasta), 'fasta')}
        controls = dict(CONTROL_12E8_H=supplied['12E8:H|PDBID|CHAIN|SEQUENCE'],
            CONTROL_12E8_L=supplied['12E8:L|PDBID|CHAIN|SEQUENCE'])
        old = list(csv.DictReader((ROOT / 'results/novelty_predictor/patent_extension/analysis_inputs.csv').open()))
        native = next(r for r in old if r['name'] == 'KNOWN_TFR1_POSITIVE_CONTROL')
        controls['CONTROL_NATIVE_RECEPTOR'] = native['sequence']
        save(HERE / 'CONTROLS.json', dict(source_fasta_sha256=sha(fasta),
            expected=dict(CONTROL_12E8_H=['H'], CONTROL_12E8_L=['K', 'L'], CONTROL_NATIVE_RECEPTOR=[None]),
            sequences_sha256={k: hashlib.sha256(v.encode()).hexdigest() for k, v in controls.items()},
            selection='Publisher-labeled H/L controls only; other author chain IDs M/P not assumed to encode class'))
        state['state'] = 'running_controls'; save(HERE / 'STATUS.json', state)
        start = time.monotonic(); raw = model.number(controls)
        save(HERE / 'control_raw.json', raw)
        assert raw['CONTROL_12E8_H']['chain_type'] == 'H'
        assert raw['CONTROL_12E8_L']['chain_type'] in ('K', 'L')
        assert raw['CONTROL_NATIVE_RECEPTOR']['chain_type'] is None
        for name, result in raw.items():
            if result['chain_type']:
                text = annotated_sequence(result)
                assert text and text in controls[name], 'Control annotation changes sequence'
        state.update(controls_passed=True, control_seconds=time.monotonic() - start)
        state['state'] = 'running_batch'; save(HERE / 'STATUS.json', state)
        inputs = ROOT / 'results/novelty_predictor/incoming_rfd4_20261022/analysis_inputs.csv'
        rows = list(csv.DictReader(inputs.open()))
        assert len(rows) == 400
        sequences = {r['name']: r['sequence'] for r in rows}
        start = time.monotonic(); annotated = model.number(sequences)
        save(HERE / 'batch_raw.json', annotated)
        assert set(annotated) == set(sequences), 'Missing/extra model output'
        evidence = []; classes = collections.Counter(); errors = collections.Counter()
        for row in rows:
            name = row['name']; result = annotated[name]
            kind = result['chain_type']; classes[kind or 'unassigned'] += 1
            if result.get('error'): errors[str(result['error'])] += 1
            text = annotated_sequence(result)
            if kind:
                assert text and text in row['sequence'], 'Annotation changes source sequence: ' + name
            cdr3 = ''.join(aa for (pos, insertion), aa in (result.get('numbering') or [])
                if 105 <= pos <= 117 and aa != '-') if kind == 'H' else None
            evidence.append(dict(name=name, source_arm=row['source_arm'], annotation_chain_type=kind,
                annotation_score=result.get('score'), error=result.get('error'),
                annotation_sequence_sha256=hashlib.sha256(text.encode()).hexdigest() if text else None,
                imgt_cdrh3=cdr3, cdrh3_length=len(cdr3) if cdr3 else None,
                official_typing_established=False, antibody_novelty_level=None))
        save(HERE / 'ANNOTATIONS.json', dict(input_sha256=sha(inputs), designs=evidence,
            warning='Sequence annotation only. Not structural ProteinTyper, not proof of official class, and no CDRH3 reference novelty comparison.'))
        state.update(state='complete', inputs=400, classes=classes, errors=errors,
            batch_seconds=time.monotonic() - start, all_assigned_residues_map_to_input=True,
            antibody_novelty_levels_established=0, finished_at=datetime.now(timezone.utc).isoformat())
        save(HERE / 'STATUS.json', state)
        print(json.dumps(state))
    except Exception as error:
        state.update(state='failed', error=repr(error)); save(HERE / 'STATUS.json', state); raise


if __name__ == '__main__':
    main()
