"""Build explicit five-source inputs, preserving every prior analysis sequence."""
import csv
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT/'scripts'))
from novelty_level_predictor import digest, write_json


def main():
    if (HERE/'DATABASES.json').exists():
        raise FileExistsError('Frozen comparison inputs already exist')
    acquisition = json.loads((HERE/'ACQUISITION.json').read_text())
    source = ROOT/'results/novelty_predictor/DATABASES.json'
    manifest = json.loads(source.read_text())
    patent = dict(id='patents', status='available', provider='NCBI pataa',
        fasta_path=acquisition['fasta_path'], fasta_sha256=acquisition['fasta_sha256'],
        scope_complete=False, scope=acquisition['scope'],
        release=acquisition['source_metadata']['last-updated'],
        source_url='https://ftp.ncbi.nlm.nih.gov/blast/db/pataa.tar.gz',
        exported_records=acquisition['exported_records'],
        exported_residues=acquisition['exported_residues'],
        transformations=acquisition['sequence_transformations'],
        acquisition_receipt_sha256=digest(HERE/'ACQUISITION.json'),
        reason_scope_incomplete='Public Patent-division proteins only; organizer corpus/provider/version unknown')
    manifest['databases'] = [patent if d['id']=='patents' else d for d in manifest['databases']]
    manifest['prior_manifest_sha256'] = digest(source)
    manifest['omitted_scope'] = 'No silent substitute: public NCBI pataa is identified explicitly; organizer patent corpus remains unverified'
    rows = list(csv.DictReader((ROOT/'results/novelty_predictor/analysis_inputs.csv').open()))
    control = list(csv.DictReader((HERE/'PATENT_POSITIVE_CONTROL.csv').open()))
    assert len(rows)==74 and len(control)==1
    assert control[0]['name'] not in {r['name'] for r in rows}
    with (HERE/'analysis_inputs.csv').open('w', newline='') as stream:
        writer=csv.DictWriter(stream, fieldnames=['name','sequence'])
        writer.writeheader(); writer.writerows(rows+control)
    baseline=json.loads((ROOT/'results/novelty_predictor/exact_score/run/predictions.json').read_text())
    structures=dict(structures=[])
    for design in baseline['designs']:
        evidence=design['structure_evidence']
        if evidence['status']=='coordinates_verified':
            structures['structures'].append(evidence['input'])
        else:
            structures['structures'].append(dict(name=design['name'],status='unavailable',
                reason=evidence.get('reason','Frozen baseline coordinates unavailable')))
    available=[r for r in structures['structures'] if r['status']=='available']
    assert sum(len(r['domains']) for r in available)==baseline['structural_search']['query_domains']==68
    assert {r['name'] for r in structures['structures']}=={r['name'] for r in rows}
    structures['structures'].append(dict(name=control[0]['name'], status='unavailable',
        reason='Patent-source sequence-only database control; no structure inferred or borrowed'))
    write_json(HERE/'DATABASES.json', manifest)
    write_json(HERE/'STRUCTURES.json', structures)
    write_json(HERE/'INPUT_SCOPE.json', dict(
        inherited_analysis_inputs_sha256=digest(ROOT/'results/novelty_predictor/analysis_inputs.csv'),
        inherited_count=74, inherited_designs='73 retired designs + native receptor positive control, unchanged',
        added_control=control[0]['name'], added_control_provenance=acquisition['positive_control'],
        total_rows=75, submission_candidates_created=0,
        inherited_structure_report_sha256=digest(ROOT/'results/novelty_predictor/exact_score/run/predictions.json'),
        expected_coordinate_domains=68,
        database_manifest_sha256=digest(HERE/'DATABASES.json'),
        analysis_inputs_sha256=digest(HERE/'analysis_inputs.csv'),
        structure_manifest_sha256=digest(HERE/'STRUCTURES.json'),
        protocol_sha256=digest(HERE/'PROTOCOL.md')))
    print(json.dumps(dict(inputs=75, public_sequence_sources=5, patent_scope_complete=False)))


if __name__=='__main__': main()
