"""A failed export may reuse only a demonstrably completed native alignment."""
import copy
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'results/novelty_predictor/adaptive_completion'))
from run_completion import EXPORT_TIMEOUT,export_resume,stage_timeout


class ExportResume(unittest.TestCase):
    def test_preserves_completed_steps_but_rejects_changed_inputs_or_partial_alignment(self):
        with tempfile.TemporaryDirectory() as folder:
            aligned=Path(folder)/'aligned'
            aligned.write_bytes(b'1\t0\n\0')
            Path(str(aligned)+'.dbtype').write_bytes((5).to_bytes(4,'little'))
            Path(str(aligned)+'.index').write_text('0\t0\t5\n')
            commands=[['tool','prepare'],['tool','align'],['tool','convertalis','query','target',str(aligned),'new.tsv','--exact-tmscore','1']]
            steps=[dict(command=c,completed=True) for c in commands]
            steps[2]=dict(command=[*commands[2][:5],'old.tsv',*commands[2][6:]])
            previous=dict(state='failed',current_stage=2,candidate_union_sha256='same',steps=steps)
            completed,receipt=export_resume(previous,commands,'same',aligned)
            self.assertEqual(len(completed),2)
            self.assertIn(str(aligned),receipt['files'])
            completed[0]['completed']=False
            self.assertTrue(previous['steps'][0]['completed'])
            for field in ('input','completion','union'):
                changed=copy.deepcopy(previous)
                if field=='input':changed['steps'][1]['command'].append('different')
                if field=='completion':changed['steps'][1]['completed']=False
                if field=='union':changed['candidate_union_sha256']='different'
                with self.assertRaises(ValueError):export_resume(changed,commands,'same',aligned)
            self.assertGreaterEqual(stage_timeout(2),3600)
            self.assertEqual(stage_timeout(2),EXPORT_TIMEOUT)


if __name__=='__main__':unittest.main()
