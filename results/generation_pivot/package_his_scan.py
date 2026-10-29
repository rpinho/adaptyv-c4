"""Attach measured scan evidence and invoke the existing submission emitter.

Preserves the registered shortlist order and sequences. No qualification model,
mutation design, cloud execution or submission is implemented here.
"""
import csv
import hashlib
import html
import json
import os
from pathlib import Path
import subprocess
import sys
from datetime import datetime, timezone

ROOT = Path(__file__).resolve().parents[2]
RUN = ROOT/'results/rfd4_his_scan_20261023'
CANON = Path('/Users/ricardo/.worktrees/adaptyv-c4/claude-submission')


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read(path):
    return json.loads(Path(path).read_text())


def csvrows(path):
    with Path(path).open(newline='') as stream:
        return list(csv.DictReader(stream))


def write_csv(path, rows):
    with Path(path).open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)


def main():
    plan = read(RUN/'PLAN.json')
    status = read(RUN/'SCAN_STATUS.json')
    assert status['state'] == 'finished', 'Wait for the registered scan to finish'
    rows = csvrows(ROOT/'results/rfd4_final_review_20261023/review/annotated20.csv')
    ordered_ids = [r['source_id'] for r in rows]
    contexts = read(RUN/'SHORTLIST20_CONTEXT_VALIDATION.json')
    assert [r['source_id'] for r in contexts['designs']] == ordered_ids
    results, observations, failures = [], [], []
    for index, job in enumerate(plan['jobs']):
        path = RUN/'sites'/('site_%02d'%index)/'RESULTS.json'
        if not path.exists():
            failures.append(dict(index=index, job=job, status=status['attempted'][index]))
            continue
        result = read(path)
        assert result['state'] == 'complete' and result['plan_sha256'] == sha(RUN/'PLAN.json')
        assert result['job'] == job
        for leg in result['scored_legs']:
            assert sha(leg['saved_pdb']) == leg['sha256']
            assert leg['identities_verified'] and leg['maximum_backbone_displacement_A'] <= .002
            assert leg['retained_native_target_residues'] == 1278
        results.append(result)
        values = result['results']
        site = 'C'+str(job['position'])
        observations.append(dict(source_id=job['source_id'], rank=job['rank'], histidine=site,
                                 wt_interaction_score=values[f'interface_ddg_wt__{site}'],
                                 delta_H_to_ALA=values[f'interface_ddg_delta__{site}__ALA'],
                                 delta_H_to_GLN=values[f'interface_ddg_delta__{site}__GLN'],
                                 delta_H_to_HIS_control=values[f'interface_ddg_delta__{site}__HIS'],
                                 maximum_backbone_displacement_A=max(x['maximum_backbone_displacement_A'] for x in result['scored_legs']),
                                 units='TMol beta2016 score units', pH_selectivity='unmeasured'))
    assert len(results) == status['completed']
    if observations:
        write_csv(RUN/'MUTATION_ENERGIES.csv', observations)
    planned = {j['source_id'] for j in plan['jobs']}
    for row in rows:
        evidence = [r for r in observations if r['source_id'] == row['source_id']]
        expected = sum(j['source_id'] == row['source_id'] for j in plan['jobs'])
        row.update(his_scan_status=('complete' if len(evidence)==expected else 'incomplete') if row['source_id'] in planned else 'not_scanned',
                   his_scan_sites_measured=str(len(evidence)), his_scan_sites_planned=str(expected),
                   mutation_delta_units='TMol beta2016 score units',
                   his_to_ALA_delta_min='', his_to_ALA_delta_max='',
                   his_to_GLN_delta_min='', his_to_GLN_delta_max='',
                   his_self_control_max_absolute_delta='',
                   his_sites_with_positive_ALA_and_GLN_deltas=str(sum(x['delta_H_to_ALA']>0 and x['delta_H_to_GLN']>0 for x in evidence)) if evidence else '',
                   pH_selectivity_measurement='unmeasured',
                   mutation_scan_scope='Matched backbone-fixed shell relaxation; site/receptor interaction energy, not binding or pH measurement')
        if evidence:
            for aa in ('ALA','GLN'):
                values = [x['delta_H_to_'+aa] for x in evidence]
                row['his_to_'+aa+'_delta_min'] = str(min(values))
                row['his_to_'+aa+'_delta_max'] = str(max(values))
            row['his_self_control_max_absolute_delta'] = str(max(abs(x['delta_H_to_HIS_control']) for x in evidence))
    annotated = RUN/'submission_annotated.csv'
    destination = RUN/'submission_for_ricardo_review.csv'
    if destination.exists():
        raise FileExistsError('Preserve the emitted review CSV')
    write_csv(annotated, rows)
    emitter = CANON/'scripts/build_submission.py'
    expected_emitter = read(ROOT/'results/rfd4_final_review_20261023/VALIDATION.json')['canonical_emitter_sha256']
    assert sha(emitter) == expected_emitter
    extras = ','.join(k for k in rows[0] if k not in ('name','sequence','molecule_class'))
    with (RUN/'emitter.log').open('w') as log:
        subprocess.run([sys.executable, str(emitter), str(annotated), '--out', str(destination),
                        '--twin-strep','excluded','--extra-columns',extras], check=True,
                       stdout=log, stderr=subprocess.STDOUT)
    assert csvrows(destination) == rows
    assert sha(emitter) == expected_emitter
    original = csvrows(ROOT/'results/rfd4_final_review_20261023/shortlist20.csv')
    assert [(r['name'],r['sequence']) for r in rows] == [(r['name'],r['sequence']) for r in original]
    report = dict(state='complete', recorded_at=datetime.now(timezone.utc).isoformat(),
                  task4=dict(measured=20, omitted_pass=contexts['omitted_context_pass'],
                             native_staging_pass=contexts['full_native_staging_pass'],
                             omitted_backbone_pairs_under2A=contexts['total_omitted_backbone_pairs_under2A'],
                             target_CA_RMSD_range_A=contexts['target_CA_RMSD_range_A'],
                             unchanged_validated_gate=True),
                  task5=dict(planned_designs=5, planned_sites=len(plan['jobs']), completed_sites=len(results),
                             completed_designs=sum(r['his_scan_status']=='complete' for r in rows),
                             mutation_deltas=len(results)*len(plan['mutate_to']), failures=failures,
                             scored_legs=sum(len(r['scored_legs']) for r in results),
                             maximum_backbone_displacement_A=max((r['maximum_backbone_displacement_A'] for r in observations),default=None),
                             maximum_absolute_self_substitution_delta=max((abs(r['delta_H_to_HIS_control']) for r in observations),default=None),
                             sites_with_positive_ALA_and_GLN_deltas=sum(r['delta_H_to_ALA']>0 and r['delta_H_to_GLN']>0 for r in observations),
                             positive_delta_count_is_descriptive_not_a_pH_gate=True,
                             results=observations),
                  format_ready_for_ricardo_review=True, submission_rows=len(rows),
                  unchanged_sequence_and_rank_order=True, submission_csv_sha256=sha(destination),
                  emitter_sha256=expected_emitter, emitter_manifest_sha256=sha(destination.with_suffix('.manifest.json')),
                  biological_qualification_established=False, pH_selectivity_measured=False,
                  rigorous_binding_ddg_measured=False, organizer_novelty_measured=False,
                  independent_refold_performed=False, submission_performed=False,
                  new_gpu_or_cloud_spend_usd=0, spend=read(ROOT/'results/rfd4_his_sweep_20261023/SPEND_AFTER_RELEASE.json'),
                  limitations=['Mutation deltas compare relaxed site-to-receptor interaction scores; not a rigorous binding ddG',
                               'Neither pH nor protonation states are varied by this scan',
                               'Self-substitution controls quantify repacking/minimization variability; do not assume they are zero',
                               'Only contacting His in top5 are scanned; remaining15 designs are explicitly not_scanned',
                               'Native receptor is rigidly restored, not independently refolded; observed glycans and assay tag are absent from scored protein poses',
                               'No organizer novelty/type acceptance, mouse-binding or expression/developability validation',
                               'Twin-Strep exclusion and upload replacement versus accumulation remain for Ricardo to resolve'])
    (RUN/'RESULTS.json').write_text(json.dumps(report, indent=2)+'\n')
    text = (
        'RFD4-PROTEINA / HUMAN TfR1 — WORKFLOW FOR RICARDO REVIEW\n\n'
        'Four 200-design generation arms were geometrically ranked. The exact final shortlist20 '
        'was rerun through the unchanged 6W3H-validated omitted-context audit and the stricter '
        'full-native staging helper. All20 passed, with no omitted-context backbone pairs below2A.\n\n'
        f"The upstream TMolMutationDdgMetric ran on CPU for {len(results)}/{len(plan['jobs'])} registered contacting-His sites "
        'from the top5 designs. Each site used an8A repacking/minimization shell, two relaxation '
        'repeats and matched WT/H-to-ALA/H-to-GLN/H-to-HIS control legs. The complete1278-residue '
        'observed receptor dimer was retained. Backbone and residue identities were verified for '
        'every scored leg. Positive deltas mean that a substitution worsened that modeled '
        'site/receptor interaction. Units are TMol beta2016 score units.\n\n'
        'These are mutation-energy hypotheses, not measured pH sensitivity, a rigorous binding '
        'free energy, pKa or KD. Self-substitution deviations are reported without assuming zero. '
        'No energy threshold was added and the original20 sequences/rank order are preserved. '
        'Other15 designs have no mutation-scan evidence; their CSV fields explicitly say not_scanned. '
        'The canonical build_submission.py validates format only. Organizer novelty/typing, '
        'biological qualification, mouse affinity, independent refolding and expression remain '
        'unmeasured. All sequences exclude the Twin-Strep tag, an inference requiring confirmation. '
        'Ricardo resolves tag and replace-versus-accumulate questions and performs any actual '
        'upload/submission. This packet has not been submitted.\n'
    )
    (RUN/'WORKFLOW_FOR_REVIEW.txt').write_text(text)
    out = ROOT/'outbox/adaptyv-c4-his-scan-review.html'
    def link(path, title):
        return '<a href="'+html.escape(os.path.relpath(path,out.parent),quote=True)+'">'+html.escape(title)+'</a>'
    table = ''.join('<tr>'+''.join('<td>'+html.escape(str(v))+'</td>' for v in (
        r['rank'],r['histidine'],f"{r['wt_interaction_score']:.4f}",f"{r['delta_H_to_ALA']:.4f}",
        f"{r['delta_H_to_GLN']:.4f}",f"{r['delta_H_to_HIS_control']:.4f}"))+'</tr>' for r in observations)
    out.write_text('''<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>Native-context and His-scan review</title><style>body{font:16px/1.55 system-ui,sans-serif;background:#f6f8fb;color:#192a3c;margin:0}main{max-width:1060px;padding:32px;margin:auto}section{background:white;border:1px solid #ccd6e1;padding:22px;margin:20px 0;border-radius:10px}table{border-collapse:collapse;width:100%}th,td{padding:10px;border-bottom:1px solid #ddd;text-align:left}a{color:#1456a0}pre{white-space:pre-wrap;font:inherit}</style><main>
<h1>Restored-context audit and CPU His-scan</h1><p>Ricardo review packet · 23 October 2026 · No submission performed</p>
<p>'''+link(destination,'CSV for Ricardo review')+' · '+link(destination.with_suffix('.manifest.json'),'Format manifest')+' · '+link(RUN/'RESULTS.json','Full measured results')+' · '+link(RUN/'PLAN.json','Preregistered scan')+'''</p>
<section><h2>Exact shortlist20 context audit</h2><p>20/20 passed the unchanged omitted-context gate and native staging checks. Zero omitted-context backbone pairs under2Å; target CA fit '''+
        f"{contexts['target_CA_RMSD_range_A'][0]:.3f}–{contexts['target_CA_RMSD_range_A'][1]:.3f}Å."+'''</p><p>Rigid context clearance does not establish binding or independent refolding.</p></section>
<section><h2>Measured mutation energies</h2><p>'''+str(len(results))+'/'+str(len(plan['jobs']))+''' sites completed. Positive H→X deltas mean the substitution worsened the modeled site/receptor interaction. Units are TMol beta2016 score units. H→H is the self-substitution/repacking control.</p>
<table><thead><tr><th>Design rank</th><th>Site</th><th>WT interaction</th><th>Δ H→A</th><th>Δ H→Q</th><th>Δ H→H control</th></tr></thead><tbody>'''+table+'''</tbody></table>
<p>'''+str(report['task5']['sites_with_positive_ALA_and_GLN_deltas'])+''' sites have positive deltas for both His-removal substitutions; this describes modeled contact support, not a pH gate. Maximum absolute self-substitution delta: '''+f"{report['task5']['maximum_absolute_self_substitution_delta']:.4f}"+''' score units. Small mutation deltas should be read alongside these controls.</p>
<p><strong>No pH-dependent calculation was performed.</strong> Mutation energies do not establish protonation, pKa, KD or the sign/magnitude of a pH response. The15 unscanned designs remain labeled not_scanned.</p></section>
<section><h2>Workflow and limits</h2><pre>'''+html.escape(text)+'''</pre></section>
<section><h2>Spend and control</h2><p>New GPU/cloud spend for this task: $0. Recorded Runpod: $56.59. Own tracked prior estimate: $73.31 plus unreconciled Modal host charges; excludes the other lane and is not an account-wide bill. No GPU app was started.</p></section></main></html>''')
    print(json.dumps(dict(state='complete', sites=len(results), csv=str(destination), report=str(out), new_gpu_spend_usd=0)))


if __name__ == '__main__':
    main()
