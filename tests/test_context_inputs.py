"""Reject the representations that made earlier pH numbers misleading."""
import json
from pathlib import Path
import sys
import tempfile
import unittest
from Bio.PDB import PDBParser

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from postprocess_gpu_campaign import canonicalize
from run_context_ensemble import ROOT, write_protein


class ContextAdmissionTests(unittest.TestCase):
    def test_crop_cannot_enter_full_context_linkage(self):
        source=json.loads((ROOT/'results/next_campaign/parent_priorities.json').read_text())[0]['source_path']
        with tempfile.TemporaryDirectory() as tmp:
            with self.assertRaisesRegex(ValueError,'target sequence differs'):
                canonicalize(source,Path(tmp)/'canonical.pdb')

    def test_lysinated_target_is_rejected(self):
        full=ROOT/'results/gpu_campaign/human_11qc_full_dimer.pdb'
        source=json.loads((ROOT/'results/next_campaign/parent_priorities.json').read_text())[0]['source_path']
        parser=PDBParser(QUIET=True)
        model=parser.get_structure('full',full)[0]
        binder=parser.get_structure('source',source)[0]['C'].copy()
        for atom in binder.get_atoms():atom.coord=atom.coord+2000
        model.add(binder)
        model['A'][122].resname='LYS'
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/'lysinated.pdb';write_protein(model,path)
            with self.assertRaisesRegex(ValueError,'target sequence differs'):
                canonicalize(path,Path(tmp)/'canonical.pdb')

    def test_native_full_target_is_preserved(self):
        full=ROOT/'results/gpu_campaign/human_11qc_full_dimer.pdb'
        source=json.loads((ROOT/'results/next_campaign/parent_priorities.json').read_text())[0]['source_path']
        parser=PDBParser(QUIET=True);model=parser.get_structure('full',full)[0]
        binder=parser.get_structure('source',source)[0]['C'].copy()
        # Deliberately separated test fixture exercises sequence/topology admission;
        # it is not a binding candidate and is never written to results.
        for atom in binder.get_atoms():atom.coord=atom.coord+2000
        model.add(binder)
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/'native.pdb';write_protein(model,path)
            result=canonicalize(path,Path(tmp)/'canonical.pdb')
            self.assertTrue(result['target_sequence_verified'])
            self.assertGreater(result['minimum_binder_target_backbone_distance'],2)

    def test_monomer_fallback_restores_the_opposite_subunit(self):
        full=ROOT/'results/gpu_campaign/human_11qc_full_dimer.pdb'
        source=json.loads((ROOT/'results/next_campaign/parent_priorities.json').read_text())[0]['source_path']
        parser=PDBParser(QUIET=True);model=parser.get_structure('full',full)[0]
        model.detach_child('B')
        binder=parser.get_structure('source',source)[0]['C'].copy();binder.id='B'
        for atom in binder.get_atoms():atom.coord=atom.coord+2000
        model.add(binder)
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/'monomer.pdb';output=Path(tmp)/'canonical.pdb';write_protein(model,path)
            result=canonicalize(path,output)
            self.assertTrue(result['opposite_subunit_restored'])
            restored=parser.get_structure('restored',output)[0]
            self.assertEqual(len(list(restored['B'])),639)
            self.assertEqual(set(c.id for c in restored),set('ABC'))


if __name__=='__main__':unittest.main()
