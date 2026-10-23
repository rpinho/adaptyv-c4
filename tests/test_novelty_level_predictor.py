"""Boundary and missing-evidence tests, independent of GPU/generator work."""
import csv
import importlib.util
from pathlib import Path
import unittest

SPEC = importlib.util.spec_from_file_location('novelty_predictor', Path(__file__).resolve().parents[1] / 'scripts/novelty_level_predictor.py')
novelty = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(novelty)


def domain(end, tm):
    return dict(start=1, end=end, tm_score=tm, normalization='query_domain_length')


class NoveltyTests(unittest.TestCase):
    def test_sequence_boundaries(self):
        for identity, expected in [(0.30,'low'), (0.300001,'moderate'), (0.70,'moderate'), (0.700001,'high')]:
            with self.subTest(identity=identity): self.assertEqual(novelty.sequence_class(identity), expected)

    def test_published_general_scale(self):
        for sequence, expected in [('low',[4,3,2]), ('moderate',[3,2,2]), ('high',[2,1,1])]:
            with self.subTest(sequence=sequence):
                self.assertEqual([novelty.general_level(sequence,s) for s in ('low','moderate','high')], expected)

    def test_tm_and_coverage_boundaries(self):
        for end,tm,expected in [(70,.5,'low'), (71,.5,'moderate'), (71,.499999,'low'), (70,.8,'high'), (70,.799999,'low'), (69,1.,'low')]:
            with self.subTest(end=end,tm=tm): self.assertEqual(novelty.structural_class(100,[domain(end,tm)])['category'], expected)

    def test_no_double_counting_redundant_domain_hits(self):
        with self.assertRaisesRegex(ValueError,'Overlapping'):
            novelty.structural_class(100,[domain(40,.9),domain(40,.9)])

    def test_tm_normalization_must_be_query_domain(self):
        with self.assertRaisesRegex(ValueError,'normalized'):
            novelty.structural_class(100,[dict(start=1,end=80,tm_score=.9,normalization='reference_length')])

    def test_missing_database_can_only_increase_similarity(self):
        self.assertEqual(novelty.possible_levels('low',['patents']),[1,2,3,4])
        self.assertEqual(novelty.possible_levels('moderate',['patents']),[1,2,3])
        self.assertEqual(novelty.possible_levels('high',['patents']),[1,2])
        self.assertEqual(novelty.possible_levels('low',[]),[2,3,4])
        self.assertEqual(novelty.possible_levels('low',[],structure='low'),[4])

    def test_reference_hit_identity_is_a_fraction_not_percent(self):
        import tempfile
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/'hits.tsv'
            path.write_text('q000001\tref\t80\t40\t50\t50\t50\t1\t50\t1\t50\t0\t100\t1\t1\n')
            with self.assertRaisesRegex(ValueError,'fraction'):
                novelty.parse_hits(path,'swissprot',{'q000001':dict(sequence='A'*50)})

    def test_identity_uses_integer_counts_before_boundary_comparison(self):
        import tempfile
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/'hits.tsv'
            path.write_text('q000001\tref\t0.7\t701\t1000\t1000\t1000\t1\t1000\t1\t1000\t0\t100\t1\t1\n')
            hit=novelty.parse_hits(path,'swissprot',{'q000001':dict(sequence='A'*1000)})['q000001'][0]
            self.assertEqual(hit['fident'],.701)
            self.assertEqual(novelty.sequence_class(hit['fident']),'high')

    def test_highest_identity_is_not_highest_score(self):
        hits=[dict(fident=.4,bits=200,database='pdb',target='long'),dict(fident=.8,bits=100,database='pdb',target='closer')]
        self.assertEqual(novelty.identity_ranked_hits(hits)[0]['target'],'closer')

    def test_input_does_not_silently_strip_colons_or_fix_case(self):
        import tempfile
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/'designs.csv'
            with path.open('w',newline='') as f:
                w=csv.writer(f);w.writerow(['name','sequence']);w.writerow(['bad','ACD:efg'])
            with self.assertRaisesRegex(ValueError,'silent repairs'): novelty.load_designs(path)

    def test_domains_matching_different_targets_aggregate_without_averaging_tm(self):
        matches = [dict(start=1,end=40,tm_score=.51,normalization='query_domain_length'),
                   dict(start=41,end=80,tm_score=.81,normalization='query_domain_length')]
        result = novelty.structural_class(100,matches)
        self.assertEqual(result['category'],'moderate')
        self.assertEqual(result['moderate_coverage'],.8)
        self.assertEqual(result['high_coverage'],.4)

    def test_discontinuous_domain_coverage_counts_its_residues_once(self):
        d=dict(segments=[[1,35],[66,100]],tm_score=.8,normalization='query_domain_length')
        result=novelty.structural_class(100,[d])
        self.assertEqual(result['category'],'high')
        self.assertEqual(result['high_coverage'],.7)
        d['tm_score']=.5
        self.assertEqual(novelty.structural_class(100,[d])['category'],'low')
        d['segments'][0][1]=36
        self.assertEqual(novelty.structural_class(100,[d])['category'],'moderate')

    def test_another_domain_can_occupy_a_discontinuous_domains_gap(self):
        matches=[dict(segments=[[1,20],[81,100]],tm_score=.8,normalization='query_domain_length'),
                 dict(start=21,end=80,tm_score=.5,normalization='query_domain_length')]
        result=novelty.structural_class(100,matches)
        self.assertEqual(result['category'],'moderate')
        self.assertEqual(result['moderate_coverage'],1.)
        self.assertEqual(result['high_coverage'],.4)
        matches[1]['start']=20
        with self.assertRaisesRegex(ValueError,'Overlapping'): novelty.structural_class(100,matches)

    def test_discontinuous_domains_reject_ambiguous_or_reordered_ranges(self):
        for d in [dict(segments=[]),dict(segments=[[1,20],[20,30]]),dict(segments=[[21,30],[1,20]]),
                  dict(segments=[[1,101]]),dict(segments=[[True,10]]),dict(segments=[[1,20]],start=1,end=20)]:
            with self.subTest(domain=d),self.assertRaises(ValueError):
                novelty.structural_class(100,[dict(**d,tm_score=.8,normalization='query_domain_length')])

    def test_discontinuous_query_keeps_gaps_and_coordinates(self):
        import tempfile
        residues=[(100+i,'A',(float(i),float(i+1),float(i+2)),'ALA') for i in range(1,11)]
        d=dict(segments=[[1,3],[8,10]])
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/'domain.pdb';path.write_text(novelty.domain_query_pdb(residues,d))
            rr=novelty.ca_residues(path,'A')
        self.assertEqual([r[0] for r in rr],[1,2,3,8,9,10])
        self.assertEqual([r[2] for r in rr],[(float(i),float(i+1),float(i+2)) for i in [1,2,3,8,9,10]])

    def test_discontinuous_query_tm_denominator_excludes_gap(self):
        import tempfile
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/'hits.tsv'
            path.write_text('q\tref\t40\t40\t1\t40\t1\t40\t40\t0.8\t0.8\tAAAA\tAAAA\n')
            d=dict(segments=[[1,20],[81,100]])
            self.assertEqual(novelty.parse_structure_hits(path,{'q':d})['q'][0]['qlen'],40)
            path.write_text(path.read_text().replace('\t40\t40\t1','\t100\t40\t1',1))
            with self.assertRaisesRegex(ValueError,'domain length'): novelty.parse_structure_hits(path,{'q':d})

    def test_structure_hit_ranking_uses_query_normalized_not_target_normalized_tm(self):
        import tempfile
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/'hits.tsv'
            path.write_text('q\tshort\t100\t40\t1\t40\t1\t40\t40\t0.4\t1.0\tAAAA\tAAAA\n'
                            'q\tfull\t100\t100\t1\t100\t1\t100\t100\t0.81\t0.81\tAAAA\tAAAA\n')
            hits=novelty.parse_structure_hits(path,{'q':dict(start=1,end=100)})['q']
            self.assertEqual(hits[0]['target'],'full')
            self.assertEqual(hits[0]['qtmscore'],.81)

    def test_mismatched_domain_length_rejected_before_classification(self):
        import tempfile
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/'hits.tsv'
            path.write_text('q\tref\t40\t40\t1\t40\t1\t40\t40\t1.0\t1.0\tAAAA\tAAAA\n')
            with self.assertRaisesRegex(ValueError,'domain length'):
                novelty.parse_structure_hits(path,{'q':dict(start=1,end=100)})

    def test_discovery_union_repeated_targets_do_not_add_evidence(self):
        import tempfile
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/'hits.tsv'
            path.write_text('q\tsame\t100\t100\t1\t100\t1\t100\t100\t0.7\t0.7\tAAAA\tAAAA\n'
                            'q\tsame\t100\t100\t1\t100\t1\t100\t100\t0.81\t0.81\tAAAA\tAAAA\n')
            hits=novelty.parse_structure_hits(path,{'q':dict(start=1,end=100)})['q']
            self.assertEqual(len(hits),1)
            self.assertEqual(hits[0]['qtmscore'],.81)

    def test_ca_extraction_preserves_chain_identity(self):
        path=Path(__file__).resolve().parents[1]/'data/11qc.cif'
        rr=novelty.ca_residues(path,'A')
        by_number={r[0]:r for r in rr}
        self.assertEqual([(p,by_number[p][1]) for p in (205,206,208,209,210)],
                         [(205,'K'),(206,'N'),(208,'R'),(209,'L'),(210,'V')])
        self.assertEqual(len(rr),639)

    def test_exact_export_contract_changes_the_observed_novelty_level(self):
        """Upstream-like exporter defaults approximate; missing exact flag changes level.

        This exercises the full structural command path through a local tool stub,
        rather than inspecting whether a command happens to contain a string.
        The numbers/coordinates are synthetic fixtures, not scientific evidence.
        """
        import argparse
        import json
        import tempfile
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);out=root/'out';out.mkdir()
            pdb=root/'fixture.pdb'
            pdb.write_text(''.join(f'ATOM  {i:5d}  CA  ALA A{i:4d}    {3.8*i:8.3f}{0.:8.3f}{0.:8.3f}  1.00  0.00           C\n'
                                   for i in range(1,11))+'TER\nEND\n')
            db=root/'db';db.write_bytes(b'test-only database fingerprint')
            db_manifest=root/'db.json'
            db_manifest.write_text(json.dumps(dict(prefix=str(db),files={'':novelty.digest(db)})))
            structures=root/'structures.json'
            structures.write_text(json.dumps(dict(structures=[dict(name='fixture',status='available',
                path=str(pdb),sha256=novelty.digest(pdb),chain='A',
                domains=[dict(start=1,end=10,source='synthetic test-only domain')])])) )
            binary=root/'foldseek_fixture'
            binary.write_text('''#!/usr/bin/env python3
import sys
from pathlib import Path
if sys.argv[1] == 'version':
    print('test-only external export contract')
elif sys.argv[1] == 'convertalis':
    exact = '--exact-tmscore' in sys.argv and sys.argv[sys.argv.index('--exact-tmscore')+1]=='1'
    score = '0.79' if exact else '0.81'
    Path(sys.argv[5]).write_text('q000001d001\\tref\\t10\\t10\\t1\\t10\\t1\\t10\\t10\\t'+score+'\\t'+score+'\\tAAAAAAAAAA\\tAAAAAAAAAA\\n')
''')
            binary.chmod(0o755)
            args=argparse.Namespace(out=out,structures=structures,structure_database=db_manifest,
                foldseek=str(binary),threads=1,timeout_structure=10,
                structure_prefilter_mode='union',structure_max_seqs=50000)
            records=[dict(name='fixture',sequence='A'*10),dict(name='missing',sequence='A'*10)]
            report=dict(designs=[dict(name='fixture',sequence_evidence_class='no_significant_hit',
                                     sequence_database_gaps=['patents']),
                                dict(name='missing',sequence_evidence_class='no_significant_hit',
                                     sequence_database_gaps=['patents'],local_predicted_level=2,predicted_level=2,
                                     structural_similarity=dict(category='high'),
                                     local_level_interval_with_unfinished_alignments=[2],
                                     possible_levels_if_structure_policy_accepted=[2])])
            novelty.run_structural_half(records,report,args)
            self.assertEqual(report['designs'][0]['local_predicted_level'],3)
            self.assertEqual(report['designs'][0]['structural_similarity']['category'],'moderate')
            self.assertIsNone(report['designs'][1]['local_predicted_level'])
            self.assertIsNone(report['designs'][1]['predicted_level'])
            self.assertIsNone(report['designs'][1]['structural_similarity'])
            self.assertNotIn('local_level_interval_with_unfinished_alignments',report['designs'][1])
            self.assertNotIn('possible_levels_if_structure_policy_accepted',report['designs'][1])


if __name__ == '__main__': unittest.main()
