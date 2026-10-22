"""After verified GPU release only: structural staging, never automatic qualification."""
import json
import sys
from pathlib import Path
from cif2pdb import convert
from postprocess_gpu_campaign import canonicalize
from hard_sequence_refold import dump

ROOT=Path(__file__).resolve().parents[1]

def main(name='hard_sequence_refold_v15'):
    if name not in ('hard_sequence_refold_v15','domain_hard_refold_v17'):raise ValueError('Unknown diagnostic output')
    out=ROOT/'results'/name
    status=json.loads((out/'supervisor.json').read_text())
    if not status.get('gpu_released_verified'):raise RuntimeError('Actual GPU release required first')
    predictions=json.loads((out/'mirrored/refold_results.json').read_text())
    if predictions.get('status')!='completed':raise RuntimeError('Completed diagnostic evidence required')
    expected={s['name']+'_complex_'+r['checkpoint']+'.cif' for s in predictions['sequences']
              for r in s.get('complex',{}).get('models',[])}
    sources=sorted((out/'mirrored').glob('*_complex_model_*.cif'))
    if {p.name for p in sources}!=expected:raise RuntimeError('Complex exports differ from completed prediction jobs')
    audits=[]
    for source in sources:
        raw=out/(source.stem+'_raw.pdb');canonical=out/(source.stem+'_canonical.pdb')
        record=dict(source=str(source),submission_recommended=False)
        try:
            convert(str(source),str(raw));record['audit']=canonicalize(raw,canonical)
            record['status']='native_context_staging_passed'
        except (ValueError,KeyError) as exc:record.update(status='structural_rejection',reason=str(exc))
        audits.append(record)
    dump(out/'native_context_audits.json',dict(status='completed',models=audits,expected_complex_jobs=len(expected),
        qualified_candidates=0,new_gpu_spend_usd=0,
        interpretation='Structural diagnostic only; all final quality/pose/novelty criteria require review before pH/submission'))

if __name__=='__main__':main(*sys.argv[1:])
