"""Map existing RFD4 outputs into the unchanged native-context audit.

Sequence identifies both chains. Only chain/residue labels change in the crop;
the real audit performs its existing all-target-CA fit and clash measurements.
"""
import argparse
import csv
import json
from pathlib import Path

from Bio.PDB import MMCIFParser,Model,Chain
from Bio.Data.PDBData import protein_letters_3to1 as AA

import audit_full_context as gate
import restore_apical_context as staging
from run_context_ensemble import sha,write_protein,dump

ROOT=Path(__file__).resolve().parents[1]


def sequence(chain):return ''.join(AA[r.resname] for r in staging.protein_residues(chain))


def prepare_crop(source,binder_sequence,destination):
    native=staging.native_model()
    expected=''.join(AA[native['A'][n].resname] for n in staging.DOMAIN)
    source_model=MMCIFParser(QUIET=True).get_structure('source',str(source))[0]
    chains=[c for c in source_model if staging.protein_residues(c)]
    targets=[c for c in chains if sequence(c)==expected]
    binders=[c for c in chains if sequence(c)==binder_sequence]
    if len(chains)!=2 or len(targets)!=1 or len(binders)!=1 or targets[0] is binders[0]:
        raise ValueError('Need exactly one exact195-residue native target and one exact binder sequence')
    model=Model.Model(0)
    for chain,cid,numbers in ((targets[0],'A',staging.DOMAIN),(binders[0],'B',range(1,len(binder_sequence)+1))):
        rr=staging.protein_residues(chain)
        if len(rr)!=len(chain):raise ValueError('Nonprotein source residues cannot be silently omitted')
        if len(rr)!=len(numbers):raise ValueError('Incomplete source residue mapping')
        out=Chain.Chain(cid)
        for original,number in zip(rr,numbers):
            if any(name not in original for name in staging.BACKBONE):raise ValueError('Incomplete backbone')
            residue=original.copy();residue.id=(' ',number,' ');out.add(residue)
        model.add(out)
    staging.audit_domain(model['A'],native['A'])
    staging.write(model,destination)
    return model,dict(source_sha256=sha(source),crop_sha256=sha(destination),
                      target_source_chain=targets[0].id,binder_source_chain=binders[0].id,
                      mapped_target_residues=195,author_range=[189,383],
                      chain_identity_by_exact_sequence=True,source_coordinates_not_fitted_or_relaxed=True)


def run(shortlist,structures,out):
    out=Path(out);out.mkdir(parents=True)
    with Path(shortlist).open(newline='') as stream:selected=list(csv.DictReader(stream))
    specs_list=json.loads(Path(structures).read_text())['structures']
    specs={s['name']:s for s in specs_list}
    if len(specs)!=len(specs_list):raise ValueError('Duplicate structure identity')
    protected={str(ROOT/'scripts'/name):sha(ROOT/'scripts'/name) for name in
               ('audit_full_context.py','run_context_ensemble.py','restore_apical_context.py','postprocess_gpu_campaign.py')}
    fixture=out/'audit_fixture';inputs=fixture/'results/v3_all_pdb';inputs.mkdir(parents=True)
    (fixture/'data').mkdir();(fixture/'data/11qc.cif').symlink_to(ROOT/'data/11qc.cif')
    crops=out/'mapped_crops';crops.mkdir();restored=out/'restored_staging';restored.mkdir()
    records=[]
    for index,row in enumerate(selected,1):
        source_id=row.get('source_id') or row['name']
        record=dict(source_id=source_id,shortlist_name=row['name'],fixture_name='candidate_%03d'%index)
        records.append(record)
        try:
            spec=specs[source_id];source=Path(spec['path']).expanduser()
            if not source.is_absolute():source=Path(structures).resolve().parent/source
            if sha(source)!=spec['sha256']:raise ValueError('Source CIF fingerprint changed')
            crop=crops/(record['fixture_name']+'.pdb')
            model,receipt=prepare_crop(source,row['sequence'],crop)
            record.update(status='prepared',preparation=receipt)
            # The audit contract is source receptor A/B and binder C. Preserve
            # the original generated coordinates; the unchanged gate does fit().
            binder=model['B'];model.detach_child('B');binder.id='C';model.add(binder)
            write_protein(model,inputs/(record['fixture_name']+'.pdb'))
        except Exception as error:record.update(status='preparation_unestablished',error=repr(error))
    prepared=[r for r in records if r['status']=='prepared']
    if prepared:
        original_root=gate.ROOT
        try:gate.ROOT=fixture;gate.main()
        finally:gate.ROOT=original_root
        measured={r['design']:r for r in json.loads((fixture/'results/full_context_audit/candidates.json').read_text())}
        for record in prepared:
            actual=measured[record['fixture_name']];record['unchanged_gate']=actual
            if actual['target_ca_rmsd']>2:
                record['status']='inconclusive_context_transfer';continue
            record['omitted_context_rigid_pass']=actual['rigid_context_pass']
            # Reuse the existing stricter staging helper; never write invented
            # coordinates or relax the generated binder to clear a collision.
            try:
                receipt=staging.restore(crops/(record['fixture_name']+'.pdb'),restored/(record['fixture_name']+'.pdb'))
                record.update(status='measured',full_native_staging_pass=True,staging=receipt)
            except Exception as error:
                record.update(status='measured',full_native_staging_pass=False,staging_error=repr(error))
    assert all(sha(path)==fingerprint for path,fingerprint in protected.items())
    result=dict(state='complete',inputs=len(selected),measured=len(prepared),designs=records,
                omitted_context_rigid_pass=sum(bool(r.get('omitted_context_rigid_pass')) for r in records),
                full_native_staging_pass=sum(bool(r.get('full_native_staging_pass')) for r in records),
                protected_source_sha256=protected,shortlist_sha256=sha(shortlist),structures_sha256=sha(structures),
                organizer_acceptance_measured=False,independent_refolding_performed=False,new_gpu_spend_usd=0,
                limitations=['Rigid context transfer is a diagnostic, not binding/pH or full-complex prediction',
                    'The unchanged audit human-text template includes historical-pool claims; use its numeric JSON and this report',
                    'Staging peptide/backbone/NAG checks remain unchanged; missing/inconclusive inputs are retained, not passed'])
    dump(out/'RESULTS.json',result)
    (out/'RESULTS.md').write_text('# RFD4 restored-context diagnostic\n\n'
        +str(len(selected))+' shortlisted inputs; '+str(len(prepared))+' mapped into the unchanged audit.\n\n'
        +str(result['omitted_context_rigid_pass'])+' rigid omitted-context passes at target CA fit<=2A; '
        +str(result['full_native_staging_pass'])+' pass the existing stricter full-native staging checks.\n\n'
        +'See RESULTS.json for every mapping, actual gate row, fit and staging error. The raw gate prose contains historical-pool boilerplate; its numerical JSON is retained unchanged. No binding, pH selectivity, independent refolding or organizer acceptance is established. No new GPU spend from this local CPU diagnostic.\n')
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--shortlist',type=Path,required=True);parser.add_argument('--structures',type=Path,required=True)
    parser.add_argument('--out',type=Path,required=True);args=parser.parse_args()
    result=run(args.shortlist,args.structures,args.out)
    print(json.dumps({k:v for k,v in result.items() if k not in ('designs','protected_source_sha256','limitations')}))


if __name__=='__main__':main()
