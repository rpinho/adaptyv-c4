"""Negative software integration control using actual400 baseline records.

The empty native export is deliberately simulated; this is not a scientific run.
"""
import contextlib
import io
import json
from pathlib import Path
import tempfile

import finalize as f


def main():
    original = f.HERE
    with tempfile.TemporaryDirectory() as folder:
        fixture = Path(folder)
        raw = fixture/'empty_export.tsv'
        raw.touch()
        union = json.loads((original/'CANDIDATE_UNION.json').read_text())
        (fixture/'CANDIDATE_UNION.json').write_text(json.dumps(union))
        (fixture/'PLAN.md').write_bytes((original/'PLAN.md').read_bytes())
        state = dict(state='alignments_complete', candidate_union_sha256=f.n.digest(fixture/'CANDIDATE_UNION.json'),
                     raw_hits=str(raw), raw_hits_sha256=f.n.digest(raw), rule_compatibility={'software_fixture': True})
        (fixture/'STATUS.json').write_text(json.dumps(state))
        f.HERE = fixture
        with contextlib.redirect_stdout(io.StringIO()):
            f.main()
        result = json.loads((fixture/'VALIDATION.json').read_text())
        assert result['local_level_counts'] == {'None': 257, '1': 73, '2': 70}
        assert result['preserved_points'] == 143 and result['exported_native_rows'] == 0
        assert len(result['incomplete_domains']) == 259
        result.update(test='Actual400 base records, simulated empty remaining export; not a scientific run',
                      state='software_control_passed', source_sha256=f.n.digest(f.__file__), control_script_sha256=f.n.digest(__file__))
        (original/'FINALIZER_CONTROL.json').write_text(json.dumps(result, indent=2)+'\n')
    f.HERE = original
    print('Negative finalizer integration control passed:143 points preserved,257 remain unestablished.')


if __name__ == '__main__':
    main()
