#!/usr/bin/env python3
"""Report whether a pod is busy with an arm. Counts argv[0]/argv[1] only, never a substring of a
script text, and never itself - a detector that greps for its own search string counts itself (this
bug made p3 look busy for ~10 minutes while it was actually finished)."""
import os, sys

ARM = 0
FOLD = 0
me = os.getpid()
for pid in filter(str.isdigit, os.listdir('/proc')):
    if int(pid) == me:
        continue
    try:
        raw = open(f'/proc/{pid}/cmdline', 'rb').read()
    except Exception:
        continue
    argv = [a for a in raw.split(b'\0') if a]
    if not argv:
        continue
    prog = os.path.basename(argv[0].decode('utf-8', 'replace'))
    # a `python -c "..."` one-liner is never a real worker
    if prog.startswith('python') and len(argv) > 1 and argv[1] == b'-c':
        continue
    if prog == 'pod_run_par.py' or (len(argv) > 1 and argv[1].endswith(b'pod_run_par.py')):
        ARM += 1
    if prog == 'colabfold_batch' or any(a.endswith(b'colabfold_batch') for a in argv[:2]):
        FOLD += 1
print(f'arm={ARM} folds={FOLD}')
