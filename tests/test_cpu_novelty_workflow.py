"""Input-integrity guards for the generic CPU novelty workflow."""
import importlib.util
import json
from pathlib import Path
import sys
import tempfile
import unittest

SCRIPTS=Path(__file__).resolve().parents[1]/'scripts'
sys.path.insert(0,str(SCRIPTS))
SPEC=importlib.util.spec_from_file_location('cpu_workflow',SCRIPTS/'run_cpu_novelty.py')
workflow=importlib.util.module_from_spec(SPEC);SPEC.loader.exec_module(workflow)


class WorkflowInputTests(unittest.TestCase):
    def fixture(self,root,second_comp):
        path=root/'two_chains.pdb'
        lines=[]
        for chain,comp in [('A','ALA'),('B',second_comp)]:
            for i in range(1,7):
                lines.append(f'ATOM  {len(lines)+1:5d}  CA  {comp} {chain}{i:4d}    {i*3.8:8.3f}{0.:8.3f}{0.:8.3f}  1.00  0.00           C\n')
        path.write_text(''.join(lines)+'END\n')
        manifest=root/'structures.json'
        manifest.write_text(json.dumps(dict(structures=[dict(name='fixture',path=path.name)])))
        return path,manifest

    def test_chain_is_identified_by_unique_full_sequence(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);path,manifest=self.fixture(root,'CYS')
            record=workflow.normalize_structures([dict(name='fixture',sequence='C'*6)],manifest)['structures'][0]
            self.assertEqual(record['chain'],'B')
            self.assertEqual(record['sha256'],workflow.n.digest(path))
            self.assertFalse(record['source_hash_was_supplied'])

    def test_ambiguous_chains_and_changed_expected_hash_are_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);path,manifest=self.fixture(root,'ALA')
            with self.assertRaisesRegex(ValueError,'one exact'):
                workflow.normalize_structures([dict(name='fixture',sequence='A'*6)],manifest)
            data=json.loads(manifest.read_text());data['structures'][0].update(chain='A',sha256='0'*64)
            manifest.write_text(json.dumps(data))
            with self.assertRaisesRegex(ValueError,'fingerprint differs'):
                workflow.normalize_structures([dict(name='fixture',sequence='A'*6)],manifest)


if __name__=='__main__':unittest.main()
