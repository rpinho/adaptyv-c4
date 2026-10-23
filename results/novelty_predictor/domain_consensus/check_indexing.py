"""Demonstrate upstream shift; verify one explicit correction and real filters."""
import hashlib
import importlib.util
import json
from pathlib import Path
import subprocess
import sys

HERE=Path(__file__).resolve().parent
CACHE=Path.home()/'.cache/adaptyv-c4/novelty_domain_source_audit'
SOURCE=CACHE/'ted-pinned/scripts'
COMPAT=CACHE/'consensus_indexing_compat'


def load(path,name):
    spec=importlib.util.spec_from_file_location(name,path)
    module=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def main():
    COMPAT.mkdir(exist_ok=True)
    original=(SOURCE/'domain_consensus.py').read_text()
    bad='assign = domstr_to_assignment(c, nres, offset=1)'
    good='assign = domstr_to_assignment(c, nres, offset=0)'
    assert original.count(bad)==1
    corrected=original.replace(bad,good)
    for name in ('get_consensus.py','score_utils.py','filter_domains_consensus.py'):
        (COMPAT/name).write_bytes((SOURCE/name).read_bytes())
    (COMPAT/'domain_consensus.py').write_text(corrected)
    sys.path.insert(0,str(SOURCE))
    old=load(SOURCE/'domain_consensus.py','ted_original')
    new=load(COMPAT/'domain_consensus.py','ted_corrected')
    cases=[
        (['1-100']*3,('1-100','na','na')),
        (['25-74']*3,('25-74','na','na')),
        (['1-100']*2,('na','1-100','na')),
        (['1-40_61-100']*3,('1-40_61-100','na','na')),
        (['1-50','51-100','1-50'],('na','1-50','51-100')),
        ([],('na','na','na'))]
    comparisons=[]
    for inputs,expected in cases:
        _,before,_=old.calculate_domain_consensus(inputs,100,consensus_levels=[3,2,1],iou_threshold=.7)
        _,after,_=new.calculate_domain_consensus(inputs,100,consensus_levels=[3,2,1],iou_threshold=.7)
        assert tuple(after)==expected,(inputs,after,expected)
        comparisons.append(dict(inputs=inputs,expected=expected,original=before,corrected=after))
    assert any(tuple(c['original'])!=tuple(c['expected']) for c in comparisons)
    fixtures=HERE/'indexing_fixtures'; fixtures.mkdir(exist_ok=True)
    filters=[]
    for i,(domain,offset,expected) in enumerate([
        ('1-24',0,'NULL'),('1-25',0,'1-25'),('1-4_10-34',0,'10-34'),
        ('1-5_10-29',0,'1-5_10-29'),('0-99',1,'1-100')]):
        source=fixtures/f'filter_{i}.input.tsv'; output=fixtures/f'filter_{i}.output.tsv'
        source.write_text(f'q000001\tdummy_md5\t100\t1\t{domain}\t1.0\n')
        subprocess.run([sys.executable,str(SOURCE/'filter_domains.py'),str(source),'-o',str(output),
            '--offset_resi',str(offset)],check=True,capture_output=True,timeout=30)
        measured=output.read_text().strip().split('\t')[4]
        assert measured==expected,(domain,measured,expected)
        filters.append(dict(input=domain,offset=offset,expected=expected,measured=measured))
    receipt=dict(all_passed=True,unmodified_consensus_negative_control_failed=True,
        correction='Only domstr_to_assignment input offset 1->0; input strings already use one-based inclusive residue IDs',
        exact_old_line=bad,exact_new_line=good,
        original_sha256=hashlib.sha256(original.encode()).hexdigest(),
        corrected_sha256=hashlib.sha256(corrected.encode()).hexdigest(),
        original_source=str(SOURCE/'domain_consensus.py'),corrected_source=str(COMPAT/'domain_consensus.py'),
        consensus_cases=comparisons,published_filter_cases=filters,
        iou_or_size_or_novelty_thresholds_changed=False)
    (HERE/'INDEXING_CHECKS.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps(dict(corrected_consensus_cases=len(comparisons),published_filter_cases=len(filters),
                         upstream_indexing_bug_reproduced=True,all_passed=True)))


if __name__=='__main__': main()
