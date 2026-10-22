"""Conservative native-context accounting; a structural proxy is not a binder."""
import csv
import hashlib
import json
import math
from pathlib import Path
import re


def trajectory_id(row):
    """Sequence redesigns share a trajectory identity; unknown identity fails closed."""
    design=row.get('design','')
    match=re.search(r'_([0-9a-f]{16})(?:_seq\d+)?$',design)
    explicit=row.get('hash','')
    if explicit and not re.fullmatch(r'[0-9a-f]{16}',explicit):
        raise ValueError('Invalid trajectory hash')
    if explicit and match and explicit!=match.group(1):
        raise ValueError('Conflicting trajectory hashes')
    identity=explicit or (match.group(1) if match else None)
    if not identity:raise ValueError('Missing trajectory identity')
    if Path(design).name!=design or not design:raise ValueError('Unsafe design name')
    return identity


def unique_parents(rows,limit=2):
    result=[];seen=set()
    for row in rows:
        identity=trajectory_id(row)
        if identity in seen:continue
        seen.add(identity);result.append(row)
        if len(result)==limit:break
    return result


def followup_path(root,run_name):
    suffix='v3' if run_name=='c4_native_apical_ablation_v3' else run_name
    return Path(root)/'results'/('modal_followup_'+suffix)


def read_outcome(campaign,followup):
    campaign=Path(campaign);followup=Path(followup)
    table=campaign/'3_Ranked/!_Ranked.csv'
    rows=[]
    if table.exists():
        with table.open() as fh:rows=list(csv.DictReader(fh))
    return _evaluate(rows,followup)


def _evaluate(rows,followup):
    report=dict(generator_accepted_rows=len(rows),native_context_admitted_trajectory_ids=[],
                engineering_advance_trajectory_ids=[],followup_completed_without_errors=False,
                rejected_records=[],interpretation='Native-context structural and pH proxies only; not binding, affinity, diversity or experimental validation')
    status_path=followup/'status.json'
    if not status_path.exists():
        report['error']='Missing scientific CPU follow-up status';return report
    try:status=json.loads(status_path.read_text())
    except (ValueError,OSError) as exc:
        report['error']=repr(exc);return report
    if not rows:
        report['followup_completed_without_errors']=status.get('status')=='no_accepted_scaffolds' and not status.get('errors')
        if not report['followup_completed_without_errors']:report['error']='Empty generator result lacks matching completed CPU status'
        return report
    if status.get('status')!='completed' or status.get('errors'):
        report['error']='CPU follow-up incomplete or contains scientific errors';return report
    by_design={r.get('design'):r for r in rows}
    admitted=set();advances=set()
    for record in status.get('designs',[]):
        name=record.get('design')
        try:
            if name not in by_design:raise ValueError('CPU record not present in ranked generator output')
            identity=trajectory_id(by_design[name])
            if record.get('status') not in {'modeled','no_eligible_pair'} or record.get('error'):
                raise ValueError('CPU record did not complete native-context admission')
            audit=record.get('audit',{})
            if audit.get('target_sequence_verified') is not True:raise ValueError('Native target sequence unverified')
            for key,minimum in [('minimum_binder_target_backbone_distance',2.0),('target_subunit_backbone_min_distance',1.8),('observed_NAG_min_distance',2.0)]:
                value=float(audit[key])
                if not math.isfinite(value) or value<minimum:raise ValueError(f'Failed native-context clearance: {key}')
            canonical=followup/(name+'_canonical.pdb')
            if hashlib.sha256(canonical.read_bytes()).hexdigest()!=audit.get('canonical_sha256'):
                raise ValueError('Canonical structure hash mismatch')
            admitted.add(identity)
            if record.get('engineering_advance') is True:advances.add(identity)
        except (ValueError,TypeError,KeyError,OSError) as exc:
            report['rejected_records'].append(dict(design=name,error=str(exc)))
    # Partial/invalid scientific evidence cannot trigger another paid batch.
    report['followup_completed_without_errors']=not report['rejected_records'] and bool(admitted)
    report['native_context_admitted_trajectory_ids']=sorted(admitted)
    report['engineering_advance_trajectory_ids']=sorted(advances)
    return report
