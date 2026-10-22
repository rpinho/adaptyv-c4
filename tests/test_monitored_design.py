"""Observation must catch collapse without changing optimizer outputs or treating proximity as binding."""
from pathlib import Path
import sys
import unittest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from monitored_design import geometry_metrics


class PhysicalMonitorTests(unittest.TestCase):
    def test_off_patch_geometry_is_reported_without_fabricating_acquisition(self):
        ca=[[i*3.8,0.,0.] for i in range(69)]
        patch=[[0.,50.,float(i)] for i in range(5)]
        result=geometry_metrics(ca,patch)
        self.assertEqual(result['unique_binder_CA_coordinates'],69)
        self.assertEqual(result['minimum_hotspot_to_binder_CA_A'],50.)
        self.assertNotIn('qualified',result)

    def test_exact_and_near_collapse_stop_the_worker(self):
        for ca in ([[0.,0.,0.]]*69,[[i*.001,0.,0.] for i in range(69)]):
            with self.assertRaisesRegex(ValueError,'collapsed'):
                geometry_metrics(ca,[[1.,2.,3.]]*5)

    def test_nonfinite_coordinates_stop_the_worker(self):
        ca=[[i*3.8,0.,0.] for i in range(69)];ca[20][0]=float('nan')
        with self.assertRaisesRegex(ValueError,'Nonfinite'):
            geometry_metrics(ca,[[1.,2.,3.]]*5)

    def test_wrong_hotspot_mapping_stops_instead_of_logging_another_patch(self):
        with self.assertRaisesRegex(ValueError,'hotspot'):
            geometry_metrics([[0.,0.,0.],[3.8,0.,0.]],[[1.,2.,3.]])


if __name__=='__main__':unittest.main()
