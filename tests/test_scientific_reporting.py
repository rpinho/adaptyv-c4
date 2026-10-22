import sys
from pathlib import Path
import unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from campaign_admission import generator_stages
from diversity_audit import audit

class ReportingTests(unittest.TestCase):
    def test_search_rejection_is_not_lost_or_called_cpu_error(self):
        log='rejected at screen design stage i_pTM=0.19 pLDDT=0.31 due to [pLDDT]\n'
        self.assertEqual(generator_stages(log),[dict(outcome='rejected',stage='screen',iptm=.19,plddt=.31,detail='due to [pLDDT]')])
    def test_related_variants_are_one_sequence_cluster(self):
        rows=[dict(name='one',sequence='ACDEFGHIKLMNPQRSTVWY'),dict(name='sibling',sequence='ACDEFGHIKLMNPQRSTVWA'),dict(name='other',sequence='YYYYYYYYYYYYYYYYYYYY')]
        r=audit(rows);self.assertEqual(r['cluster_count'],2);self.assertEqual(r['largest_cluster'],2)
        self.assertFalse(r['external_known_protein_novelty_assessed'])

if __name__=='__main__':unittest.main()
