#!/usr/bin/env python3
"""Local C4 evidence-integrity check. No network, model calls, GPU, or submissions."""
import argparse
import hashlib
import json
import subprocess
import tempfile
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results/c4_mouse_compatibility_20261026'


def verify(entries):
    failures=[]
    for name,expected in entries.items():
        p=Path(name)
        if not p.is_absolute():
            p=ROOT/p
        actual=hashlib.sha256(p.read_bytes()).hexdigest() if p.is_file() else None
        if actual!=expected:
            failures.append(dict(path=name,expected=expected,actual=actual))
    return failures


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--self-test',action='store_true')
    args=ap.parse_args()
    if args.self_test:
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/'fixture';p.write_bytes(b'c4 control')
            expected=hashlib.sha256(p.read_bytes()).hexdigest()
            assert not verify({str(p):expected})
            assert len(verify({str(p):'0'*64}))==1
            p.unlink()
            assert len(verify({str(p):expected}))==1
        print('PASS: valid digest accepted; changed and missing controls detected')
        return 0
    stamp=subprocess.check_output(['date','-u','+%Y-%m-%dT%H:%M:%SZ'],text=True).strip()
    manifest=json.loads((OUT/'monitor_manifest.json').read_text())
    failures=verify(manifest['files'])
    state=dict(checked_at_utc=stamp,checked_files=len(manifest['files']),
               status='FAIL' if failures else 'PASS',failures=failures,
               scope='C4 local artifact integrity only; no scientific revalidation or C5 polling')
    (OUT/'monitor_status.json').write_text(json.dumps(state,indent=2)+'\n')
    if failures:
        (OUT/'ALARM_C4_CPU_CHECK.txt').write_text(json.dumps(state,indent=2)+'\n')
    with (OUT/'monitor_checks.jsonl').open('a') as f:
        f.write(json.dumps(state)+'\n')
    print(json.dumps(state))
    return int(bool(failures))


if __name__=='__main__':
    raise SystemExit(main())
