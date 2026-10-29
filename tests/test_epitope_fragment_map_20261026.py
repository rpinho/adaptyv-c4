"""Regression for the 6X18 missing-residue boundary and submitted C3 immunity.

No existing results are written. Actual complex evidence joins by extracted sequence hash.
"""
import csv
import hashlib
import importlib.util
import itertools
from pathlib import Path
import unittest

from Bio.Data.PDBData import protein_letters_3to1 as AA
from Bio.PDB import PDBParser

ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('epitope_audit', ROOT/'results/c3_sprint_8h_20261025/epitope_agonist_audit.py')
audit=importlib.util.module_from_spec(spec);spec.loader.exec_module(audit)


class FragmentMappingRegression(unittest.TestCase):
    def test_gap_boundary_and_c_terminal_binder_control(self):
        mapping=audit.fragment_map()
        self.assertEqual(len(mapping),111)
        self.assertEqual([mapping[i] for i in (1,101,102,107,110,111)], [38,138,145,150,153,154])
        self.assertTrue(set(range(139,145)).isdisjoint(mapping.values()))
        # Synthetic future binder contacting the C-terminal footprint: old map mislabels it.
        patch={107,110,111}
        site=audit.glp1_site()
        corrected={audit.map_fragment_residue(i,mapping) for i in patch}&site
        old={i+37 for i in patch}&site
        self.assertEqual(corrected,{150,153,154})
        self.assertNotEqual(corrected,old)
        self.assertEqual(old,{147})

    def test_assessable_denominator_derived_from_actual_image(self):
        site=audit.glp1_site()
        image=set(audit.fragment_map().values())
        self.assertEqual((len(site),len(site&image),len(site-image)),(52,22,30))
        self.assertEqual(len(site & set(range(38,149))),19)  # Negative control: old bug.

    def test_insertion_codes_fail_explicitly(self):
        mapping=audit.fragment_map()
        self.assertEqual(audit.map_fragment_residue((102,''),mapping),145)
        with self.assertRaisesRegex(ValueError,'Insertion-coded'):
            audit.map_fragment_residue((102,'A'),mapping)

    def test_real_submitted_structures_sequence_identity_and_patch_immunity(self):
        # Numbering reference is read independently from mmCIF label_seq_id fields.
        reference={}
        for line in (ROOT/'data/6x18.cif').read_text().splitlines():
            if line.startswith('ATOM'):
                f=line.split()
                if f[6]=='F' and f[8].isdigit():
                    reference[int(f[8])]=AA[f[5]]
        csvpath=ROOT/'results/c3_sprint_8h_20261025/c3_candidate_D_swap18_20261025.csv'
        rows=list(csv.DictReader(csvpath.read_text().splitlines()))
        digest=lambda s:hashlib.sha256(s.encode()).hexdigest()
        expected={digest(r['sequence']) for r in rows}
        found={}
        parser=PDBParser(QUIET=True)
        site=audit.glp1_site()
        files=sorted((ROOT/'results/c3_sampling_depth_v1_20261025/out').glob('*/*/*rank_001*.pdb'))
        self.assertGreaterEqual(len(files),18)
        for p in files:
            chains=parser.get_structure('c',p)[0]
            seqs={c.id:''.join(AA[r.resname] for r in c if r.id[0]==' ' and r.resname in AA) for c in chains}
            matches=[c for c,s in seqs.items() if digest(s) in expected]
            if not matches:
                continue
            self.assertEqual(len(matches),1)
            bc=matches[0]; sid=digest(seqs[bc])
            tc=next(c for c,s in seqs.items() if c!=bc and len(s)==111)
            target=seqs[tc]
            mapping=audit.fragment_map()
            self.assertEqual(sum(target[i-1]==reference.get(i+37) for i in range(1,112)),101)
            self.assertEqual(sum(target[i-1]==reference[mapping[i]] for i in range(1,112)),111)
            profile=audit.interface(p)
            self.assertIsNotNone(profile)
            patch=profile['patch']
            old={i+37 for i in patch};new={mapping[i] for i in patch}
            self.assertEqual(len(old&site),len(new&site))
            self.assertGreater(len(new&site),0)  # Working contact control.
            found[sid]=(patch,new)
        self.assertEqual(set(found),expected)
        for (pa,ma),(pb,mb) in itertools.combinations(found.values(),2):
            self.assertEqual(audit.jac(pa,pb),audit.jac(ma,mb))

    def test_submitted_bytes_and_claim_scope_unchanged(self):
        files={
            'c3_candidate_D_swap18_20261025.csv':'ba5e3789c0dfdf8cf44e087f0d281b37335fe71fd370d0b64d11ea09d6d1f4d2',
            'independent/METHODOLOGY_D_swap18_20261025.txt':'6da47973814530c4aae581fe4a5c41edb5c1b159990e87a4598ee954e12d49db'}
        for name,expected in files.items():
            p=ROOT/'results/c3_sprint_8h_20261025'/name
            self.assertEqual(hashlib.sha256(p.read_bytes()).hexdigest(),expected)
        method=(ROOT/'results/c3_sprint_8h_20261025/independent/METHODOLOGY_D_swap18_20261025.txt').read_text().lower()
        self.assertNotIn('assessable',method)
        self.assertNotIn('6x18',method)


if __name__=='__main__':
    unittest.main(verbosity=2)
