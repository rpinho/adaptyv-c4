"""Coordinate-free refolding gates must distinguish geometry, confidence and pose."""
import sys
from pathlib import Path
import unittest
import copy
import json
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from hard_sequence_refold import geometry,aligned_rmsd,BB,validate_manifest
from restore_apical_context import native_model

class HardRefoldTests(unittest.TestCase):
    def points(self):
        native=native_model()['A']
        return np.array([[native[n][a].coord for a in BB] for n in range(122,152)])

    def test_native_peptide_passes_but_collapsed_backbone_does_not(self):
        xyz=self.points();mask=np.ones(xyz.shape[:2],dtype=bool)
        self.assertTrue(geometry(xyz,mask)['passed'])
        bad=xyz.copy();bad[:,0]=bad[:,1]
        self.assertFalse(geometry(bad,mask)['passed'])

    def test_incomplete_or_nonfinite_output_fails_closed(self):
        xyz=self.points();mask=np.ones(xyz.shape[:2],dtype=bool);mask[0,0]=False
        with self.assertRaises(ValueError):geometry(xyz,mask)
        xyz[0,0,0]=np.nan
        with self.assertRaises(ValueError):geometry(xyz,np.ones_like(mask))

    def test_self_overlap_is_not_hidden_by_normal_local_geometry(self):
        xyz=self.points();xyz[10:20]=xyz[:10]
        result=geometry(xyz,np.ones(xyz.shape[:2],dtype=bool))
        self.assertFalse(result['passed']);self.assertIn('nonadjacent_backbone_overlap',result['failed'])

    def test_internal_fold_agreement_removes_rigid_motion_only(self):
        ca=self.points()[:,1];rotation=np.array([[0,-1,0],[1,0,0],[0,0,1]])
        self.assertLess(aligned_rmsd(ca,ca@rotation+50),.0001)
        distorted=ca.copy();distorted[:10]+=np.array([15,0,0])
        self.assertGreater(aligned_rmsd(ca,distorted),3.5)

    def test_different_sequences_models_or_coordinates_cannot_expand_this_diagnostic(self):
        manifest=json.loads((Path(__file__).resolve().parents[1]/'results/hard_sequence_refold_v15/manifest.json').read_text())
        validate_manifest(manifest)
        for field,value in [('models',['model_5_multimer_v3']),('max_prediction_jobs',10),('designed_binder_coordinates_supplied',True)]:
            with self.subTest(field=field),self.assertRaises(ValueError):validate_manifest(manifest|{field:value})
        changed=copy.deepcopy(manifest);changed['sequences'][0]['sequence']='A'*69
        with self.assertRaises(ValueError):validate_manifest(changed)

if __name__=='__main__':unittest.main()
