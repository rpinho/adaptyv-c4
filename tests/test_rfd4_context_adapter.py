"""Artificial duplicate-native-fragment control, not a generated design."""
import contextlib
import csv
import io
import json
from pathlib import Path
import sys
import tempfile
import unittest

from Bio.PDB import Chain,Model,MMCIFIO

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
import audit_rfd4_context as adapter


class ContextAdapter(unittest.TestCase):
    def test_scrambled_chain_letters_and_a_known_native_occupancy_collision(self):
        with tempfile.TemporaryDirectory() as folder:
            folder=Path(folder);native=adapter.staging.native_model();model=Model.Model(0)
            for cid,reference,numbers in [('Q','A',range(189,384)),('Z','B',range(122,182))]:
                chain=Chain.Chain(cid)
                for n in numbers:chain.add(native[reference][n].copy())
                model.add(chain)
            source=folder/'duplicate_native_fragment.cif';writer=MMCIFIO();writer.set_structure(model);writer.save(str(source))
            binder_sequence=adapter.sequence(model['Z'])
            manifest=folder/'structures.json';manifest.write_text(json.dumps(dict(structures=[dict(name='software_control',path=str(source),sha256=adapter.sha(source))])))
            shortlist=folder/'shortlist.csv'
            with shortlist.open('w',newline='') as stream:
                writer=csv.DictWriter(stream,fieldnames=['name','source_id','sequence']);writer.writeheader()
                writer.writerow(dict(name='artificial_duplicate_native_fragment',source_id='software_control',sequence=binder_sequence))
            with contextlib.redirect_stdout(io.StringIO()):result=adapter.run(shortlist,manifest,folder/'audit')
            record=result['designs'][0]
            self.assertEqual(record['preparation']['target_source_chain'],'Q')
            self.assertEqual(record['preparation']['binder_source_chain'],'Z')
            self.assertLess(record['unchanged_gate']['target_ca_rmsd'],.01)
            self.assertGreater(record['unchanged_gate']['restored_backbone_pair_clashes_under2A'],100)
            self.assertLess(record['unchanged_gate']['restored_backbone_min_distance'],.01)
            self.assertFalse(record['omitted_context_rigid_pass'])
            self.assertFalse(record['full_native_staging_pass'])
            self.assertFalse(result['organizer_acceptance_measured'])
            with self.assertRaises(ValueError):adapter.prepare_crop(source,binder_sequence+'A',folder/'wrong.pdb')


if __name__=='__main__':unittest.main()
