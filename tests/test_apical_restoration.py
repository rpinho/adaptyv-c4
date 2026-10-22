"""Rigid restoration fixtures; none is a binding candidate."""
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

import numpy as np
from Bio.PDB import PDBParser, MMCIFIO, PDBIO

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
import restore_apical_context as domain


class ApicalRestorationTests(unittest.TestCase):
    def fixture(self,clashing=False):
        model=domain.native_model()
        model.detach_child('B')
        for c in list(model):
            if c.id!='A':model.detach_child(c.id)
        binder=model['A'].copy();binder.id='B'
        keep=set(range(400,430) if clashing else range(122,152))
        for r in list(binder):
            if r.id[0]!=' ' or r.id[1] not in keep:binder.detach_child(r.id)
        if not clashing:
            for a in binder.get_atoms():a.coord=a.coord+500
        for r in list(model['A']):
            if r.id[0]!=' ' or r.id[1] not in domain.DOMAIN:model['A'].detach_child(r.id)
        model.add(binder)
        return model

    def test_all_domain_residues_restore_under_rigid_transform(self):
        model=self.fixture()
        expected_binder={(r.id,a.name):a.coord.copy() for r in model['B'] for a in r}
        rotation=np.array([[0,-1,0],[1,0,0],[0,0,1]])
        for a in model.get_atoms():a.coord=a.coord@rotation+np.array([70,-30,20])
        with tempfile.TemporaryDirectory() as tmp:
            source=Path(tmp)/'domain.pdb';out=Path(tmp)/'staging.pdb';domain.write(model,source)
            audit=domain.restore(source,out)
            result=PDBParser(QUIET=True).get_structure('restored',out)[0]
            self.assertEqual([len(domain.protein_residues(result[c])) for c in 'ABC'],[639,639,30])
            self.assertLess(audit['domain_CA_rmsd_A'],.002)
            self.assertLess(audit['patch_CA_rmsd_A'],.002)
            self.assertAlmostEqual(audit['rotation_determinant'],1)
            self.assertTrue(audit['requires_native_full_context_refold'])
            self.assertFalse(audit['qualified_candidate'])
            native=domain.native_model()
            np.testing.assert_allclose(result['C'][122]['CA'].coord,native['A'][122]['CA'].coord+500,atol=.002)
            for cid in 'AB':
                for r in domain.protein_residues(native[cid]):
                    for a in r:np.testing.assert_allclose(result[cid][r.id][a.name].coord,a.coord,atol=.001)
            for (rid,name),coord in expected_binder.items():
                np.testing.assert_allclose(result['C'][rid][name].coord,coord,atol=.002)

    def test_omitted_native_protein_overlap_is_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            source=Path(tmp)/'bad.pdb';out=Path(tmp)/'never.pdb';domain.write(self.fixture(clashing=True),source)
            with self.assertRaisesRegex(ValueError,'full-target backbone overlap'):domain.restore(source,out)
            self.assertFalse(out.exists())

    def test_surrogate_sequence_and_disconnected_crop_are_rejected(self):
        for mutate in (True,False):
            model=self.fixture()
            if mutate:model['A'][205].resname='ALA'
            else:model['A'].detach_child((' ',200,' '))
            with self.subTest(mutate=mutate),tempfile.TemporaryDirectory() as tmp:
                source=Path(tmp)/'bad.pdb';domain.write(model,source)
                with self.assertRaisesRegex(ValueError,'sequence differs|contiguous native'):
                    domain.restore(source,Path(tmp)/'never.pdb')

    def test_patch_deformation_is_not_hidden_by_global_fit(self):
        model=self.fixture()
        for n in domain.PATCH:model['A'][n]['CA'].coord+=np.array([6,0,0])
        with tempfile.TemporaryDirectory() as tmp:
            source=Path(tmp)/'bad.pdb';domain.write(model,source)
            with self.assertRaisesRegex(ValueError,'domain deformed'):domain.restore(source,Path(tmp)/'never.pdb')

    def test_prepare_preserves_native_contiguous_sequence_and_coordinates(self):
        source=domain.ROOT/'results/gpu_campaign_v2/human_11qc_full_monomer.pdb'
        with tempfile.TemporaryDirectory() as tmp:
            out=Path(tmp)/'domain.pdb';audit=domain.prepare(source,out)
            model=PDBParser(QUIET=True).get_structure('prepared',out)[0]
            self.assertEqual([r.id[1] for r in model['A']],list(domain.DOMAIN))
            self.assertTrue(audit['native_sequence_verified'])
            self.assertTrue(audit['requires_native_full_context_refold'])
            original=PDBParser(QUIET=True).get_structure('full',source)[0]
            for r in model['A']:
                for a in r:np.testing.assert_array_equal(a.coord,original['A'][r.id][a.name].coord)

    def test_unknown_residue_is_rejected_instead_of_silently_removed(self):
        for cid in 'AB':
            model=self.fixture();next(iter(model[cid])).resname='UNK'
            with self.subTest(chain=cid),tempfile.TemporaryDirectory() as tmp:
                source=Path(tmp)/'unknown.pdb'
                writer=PDBIO();writer.set_structure(model);writer.save(str(source))
                with self.assertRaisesRegex(ValueError,'Unknown/nonprotein'):
                    domain.restore(source,Path(tmp)/'never.pdb')

    def test_modified_native_template_fails_identity_check(self):
        model=domain.native_model();model['A'][142].resname='ALA'
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/'changed.cif';writer=MMCIFIO();writer.set_structure(model);writer.save(str(path))
            with patch.object(domain,'NATIVE_TEMPLATE',path),self.assertRaisesRegex(ValueError,'identity/hash changed'):
                domain.native_model()

    def test_broken_binder_peptide_is_rejected_before_staging(self):
        model=self.fixture();model['B'][123]['N'].coord+=np.array([5,0,0])
        with tempfile.TemporaryDirectory() as tmp:
            source=Path(tmp)/'broken.pdb';domain.write(model,source)
            with self.assertRaisesRegex(ValueError,'Broken binder peptide'):
                domain.restore(source,Path(tmp)/'never.pdb')

    def test_glycan_overlap_blocks_otherwise_clear_staging(self):
        template=domain.native_model();model=self.fixture()
        # Synthetic glycan displacement isolates this gate from protein clashes.
        nag=next(a for c in template for r in c if r.resname=='NAG' for a in r)
        nag.coord=model['B'][122]['CA'].coord.copy()
        with tempfile.TemporaryDirectory() as tmp:
            source=Path(tmp)/'clear_protein.pdb';domain.write(model,source)
            with patch.object(domain,'native_model',return_value=template),self.assertRaisesRegex(ValueError,'NAG'):
                domain.restore(source,Path(tmp)/'never.pdb')


if __name__=='__main__':unittest.main()
