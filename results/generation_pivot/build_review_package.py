"""Package the completed CPU review; invoke the unchanged canonical emitter.

No generation, ranking, cloud mutation, biological qualification or submission.
"""
import csv
import hashlib
import html
import json
import os
from pathlib import Path
import subprocess
import sys
from collections import Counter
from datetime import datetime, timezone

ROOT = Path(__file__).resolve().parents[2]
RUN = ROOT / 'results/rfd4_final_review_20261023'
CANON = Path('/Users/ricardo/.worktrees/adaptyv-c4/claude-submission')


def read(path):
    return json.loads(path.read_text())


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def rows(path):
    with path.open(newline='') as stream:
        return list(csv.DictReader(stream))


def main():
    destination = RUN / 'review'
    destination.mkdir(exist_ok=True)
    annotated_path = destination / 'annotated20.csv'
    submission = destination / 'submission_review.csv'
    if annotated_path.exists() or submission.exists():
        raise FileExistsError('Preserve the existing review package')
    pool = {row['source_id']: row for row in rows(RUN / 'analysis_pool/annotated80.csv')}
    ranked = rows(RUN / 'shortlist20.csv')
    provenance = read(RUN / 'shortlist20.provenance.json')
    assert len(pool) == 80 and len(ranked) == 20
    selected = []
    for rank, evidence in zip(ranked, provenance['shortlist'], strict=True):
        item = dict(pool[rank['source_id']])
        assert item['sequence'] == rank['sequence'] == evidence['sequence']
        assert rank['source_id'] == evidence['source_id']
        item.update(rank)
        item['hypothesis_basis'] = evidence['hypothesis_basis']
        item['his_fraction'] = str(item['sequence'].count('H') / len(item['sequence']))
        selected.append(item)
    with annotated_path.open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(selected[0]), lineterminator='\n')
        writer.writeheader()
        writer.writerows(selected)
    extras = [key for key in selected[0] if key not in ('name', 'sequence', 'molecule_class')]
    emitter = CANON / 'scripts/build_submission.py'
    emitter_hash = sha(emitter)
    with (destination / 'emitter.log').open('w') as log:
        subprocess.run([sys.executable, str(emitter), str(annotated_path),
                        '--out', str(submission), '--twin-strep', 'excluded',
                        '--extra-columns', ','.join(extras)], check=True,
                       stdout=log, stderr=subprocess.STDOUT)
    assert sha(emitter) == emitter_hash
    emitted = rows(submission)
    assert emitted == selected  # Including rank, extra evidence and literal sequence.
    assert len({row['sequence'] for row in emitted}) == 20
    assert max(Counter(row['cluster'] for row in selected).values()) <= 3
    assert all(row['native_context_staging_pass'] == 'True' for row in selected)

    context = read(RUN / 'context80/RESULTS.json')
    comparison = read(RUN / 'ARM_COMPARISON.json')
    spend = read(ROOT / 'results/rfd4_his_sweep_20261023/SPEND_AFTER_RELEASE.json')
    apps = read(RUN / 'final_modal_apps.json')
    own_ids = ['ap-MLsDay51dMvQnlyUqCMuwc', 'ap-3oixAY66KO1YAveVPzHW1T']
    owned = [app for app in apps if app['app_id'] in own_ids]
    assert len(owned) == 2 and all(a['state'] == 'stopped' and str(a['tasks']) == '0' for a in owned)
    for source, expected in context['protected_source_sha256'].items():
        assert sha(Path(source)) == expected
    old = read(ROOT / 'results/rfd4_geom_20261023/shortlist_context/RESULTS.json')
    new_by_id = {row['source_id']: row for row in context['designs']}
    overlap = [row for row in old['designs'] if row['source_id'] in new_by_id]
    errors = []
    for before in overlap:
        after = new_by_id[before['source_id']]
        assert before['full_native_staging_pass'] == after['full_native_staging_pass']
        assert before['omitted_context_rigid_pass'] == after['omitted_context_rigid_pass']
        for key in ('target_ca_rmsd', 'restored_backbone_pair_clashes_under2A', 'restored_backbone_min_distance'):
            errors.append(abs(before['unchanged_gate'][key] - after['unchanged_gate'][key]))
    assert errors and max(errors) < 1e-5
    fits = [r['unchanged_gate']['target_ca_rmsd'] for r in context['designs']]
    validation = dict(
        state='complete', recorded_at=datetime.now(timezone.utc).isoformat(),
        generated_rows=comparison['generated_rows'], unique_generated_sequences=comparison['unique_sequences'],
        reviewed_rows=len(pool), ranked_rows=len(selected), unique_review_sequences=20,
        selected_arms=dict(Counter(r['source_arm'] for r in selected)),
        selected_hypotheses=dict(Counter(r['hypothesis_basis'] for r in selected)),
        selected_sequence_families=len({r['cluster'] for r in selected}),
        maximum_selected_per_family=max(Counter(r['cluster'] for r in selected).values()),
        omitted_context_pass=context['omitted_context_rigid_pass'], full_native_staging_pass=context['full_native_staging_pass'],
        target_fit_rmsd_range_A=[min(fits), max(fits)],
        selected_his_fraction_at_least_0p30=sum(float(r['his_fraction']) >= .30 for r in selected),
        fraction_0p30_is_descriptive_review_flag_not_a_gate=True,
        overlap_runtime_parity=dict(designs=len(overlap), maximum_numeric_difference=max(errors), passed=True),
        canonical_emitter_sha256=emitter_hash, canonical_emitter_unchanged=True,
        submission_csv_sha256=sha(submission), canonical_manifest_sha256=sha(submission.with_suffix('.manifest.json')),
        annotated_input_sha256=sha(annotated_path), exact_annotation_rank_and_sequence_parity=True,
        source_context_sha256=sha(RUN/'context80/RESULTS.json'),
        source_rank_provenance_sha256=sha(RUN/'shortlist20.provenance.json'),
        protected_gate_sources_unchanged=True, owned_generation_apps=owned,
        spend=spend, submission_performed=False, pH_selectivity_measured=False,
        affinity_measured=False, independent_refold_performed=False, organizer_novelty_measured=False,
        new_gpu_or_remote_cpu_spend_for_packaging_usd=0)
    (RUN/'VALIDATION.json').write_text(json.dumps(validation, indent=2)+'\n')

    workflow = (
        'REVIEW DRAFT FOR RICARDO — NOT SUBMITTED\n\n'
        'Open-weight RFD4-Proteina generated four arms of 200 human TfR1 apical-domain '
        'binder poses each, with inline geometric metrics. One unconditioned arm and '
        'His-prior 0.10, 0.20 and 0.30 arms used the same target/epitope, seed and generation settings. '
        'The per-arm top-20 geometry pool was rigidly restored to the complete observed human '
        'receptor dimer using unchanged native-context checks. All 80 passed those geometric checks. '
        'The final 20 were selected by explicit His-N/basic-receptor proximity hypotheses, then '
        'hotspot/clash/bond geometry, with at most three sequences per declared family.\n\n'
        'His contacts are descriptive distances in generated poses against observed native sidechain '
        'rotamers; they do not measure protonation, affinity, its direction of change or pH selectivity. '
        'The 4.5 A contact convention is not a validated functional threshold. Priors are soft, '
        'not achieved sequence frequencies. Seven selected sequences contain at least 30% His; '
        'this is a descriptive composition flag, not a validated developability cutoff. '
        'Many geometry tuples tie, so lexical source-ID ordering can determine tied ranks. '
        'The family rule is sequence-based and does not establish structural independence.\n\n'
        'Outstanding: no independent complex refolding, interface ddG/His scan, biochemical pH '
        'binding measurement, expression/developability validation or mouse binding validation. '
        'The new 800 outputs have not received a fresh complete novelty screen; earlier local '
        'novelty labels are not transferred by design index. Any exact-sequence historical evidence '
        'is identified separately in the CSV and is not an organizer score. The organizer novelty '
        'gate remains unmeasured. The CSV is format-valid, not submission-qualified.\n\n'
        'C-terminal Twin-Strep is excluded from design sequences in this artifact. Tag accessibility '
        'and unresolved glycans were not modeled. Ricardo must resolve the tag convention and '
        'whether subsequent uploads replace or accumulate against the first-20 rule before '
        'choosing an actual entry. Ricardo alone controls accounts, uploads and submission. '
        'No account, upload or entry was created. Do not treat this hypothesis set as a final '
        'recommendation to consume 20 competition slots.\n'
    )
    (destination/'WORKFLOW_REVIEW.txt').write_text(workflow)

    outbox = ROOT/'outbox'
    outbox.mkdir(exist_ok=True)
    output = outbox/'adaptyv-c4-generation-review.html'
    def link(path, label):
        return '<a href="'+html.escape(os.path.relpath(path, outbox), quote=True)+'">'+html.escape(label)+'</a>'
    arm_rows = ''.join('<tr>'+''.join('<td>'+html.escape(str(value))+'</td>' for value in (
        a['arm'], a['design_rows'], a['his_at_least3'], round(a['mean_his_count'], 3),
        round(a['median_his_percent'], 2), a['cys_present']))+'</tr>' for a in comparison['arms'])
    candidates = ''.join('<tr>'+''.join('<td>'+html.escape(str(value))+'</td>' for value in (
        row['rank'], row['name'], len(row['sequence']), row['n_his'], f"{100*float(row['his_fraction']):.1f}",
        row['n_cys'], row['hotspot_obedience'], row['his_nitrogens_near_basic_4p5A'],
        f"{float(row['target_fit_CA_RMSD_A']):.3f}"))+'</tr>' for row in selected)
    output.write_text('''<!doctype html><html lang="en"><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1"><title>Oct 23 generation review</title>
<style>body{font:16px/1.5 system-ui,sans-serif;color:#142235;background:#f5f7fa;margin:0}main{max-width:1200px;margin:auto;padding:32px}h1{font-size:30px}h2{margin-top:32px}a{color:#1456a0}.card{background:white;border:1px solid #ccd6e1;padding:20px;border-radius:10px;margin:16px 0}.notice{border-left:5px solid #b76400}table{border-collapse:collapse;width:100%;font-size:14px}td,th{text-align:left;padding:9px;border-bottom:1px solid #ddd;vertical-align:top}th{background:#edf1f7}td:nth-child(2){overflow-wrap:anywhere}.scroll{overflow:auto}pre{white-space:pre-wrap;font:inherit}small{color:#425169}</style>
<main><h1>800 generated outputs. 20 ranked hypotheses for review.</h1>
<p>23 October 2026 · CPU review complete · Both owned generation apps stopped, zero tasks.</p>
<div class="card notice"><strong>Format-valid review set. Functional and organizer qualification unmeasured.</strong>
<p>798 unique sequences; finite inline ranking geometry for all 800 outputs. All 80 audited poses passed unchanged omitted-context and full-native staging checks. These results establish geometric clearance only.</p></div>
<p>'''+link(submission,'Review CSV')+' · '+link(destination/'submission_review.manifest.json','Canonical format manifest')+' · '+link(RUN/'VALIDATION.json','Validation and spend')+' · '+link(destination/'WORKFLOW_REVIEW.txt','Workflow draft')+'''</p>
<h2>Measured arm comparison</h2><div class="card scroll"><table><thead><tr><th>Arm</th><th>Outputs</th><th>At least 3 His</th><th>Mean His count</th><th>Median %His</th><th>Any Cys</th></tr></thead><tbody>'''+arm_rows+'''</tbody></table>
<p>Higher priors increased His content. This does not establish better binding or pH selectivity. “Any Cys” is sequence composition, not a free-cysteine determination.</p></div>
<h2>What ranking measures</h2><p>Twenty geometry-ranked outputs per arm formed an 80-pose pool. The same canonical ranker selected the final review set using the registered adaptive contact-hypothesis policy and a global family cap of three. All 20 have a His nitrogen within 4.5 Å of a native Lys/Arg group. This proximity rule is a hypothesis, not an energy or pH prediction. Many geometry ties remain; source-ID ordering can decide tied ranks.</p>
<p>Seven selected sequences have at least 30% His. Review their composition before interpreting high-prior enrichment as progress. No validated developability cutoff was applied.</p>
<h2>Ranked review set</h2><div class="card scroll"><table><thead><tr><th>Rank</th><th>Source name</th><th>AA</th><th>His</th><th>%His</th><th>Cys</th><th>Hotspot</th><th>His near basic</th><th>Fit Å</th></tr></thead><tbody>'''+candidates+'''</tbody></table></div>
<h2>Evidence and limits</h2><div class="card"><pre>'''+html.escape(workflow)+'''</pre></div>
<p>'''+link(RUN/'analysis_pool/annotated80.csv','All 80 annotations')+' · '+link(RUN/'analysis_pool/annotated80.evidence.json','Per-His distances')+' · '+link(RUN/'context80/RESULTS.json','Unchanged gate and native staging results')+' · '+link(RUN/'shortlist20.provenance.json','Complete final ranking provenance')+'''</p>
<h2>Spend after release</h2><div class="card"><p>Recorded Runpod: $56.59. Prior own Modal estimate: $10.72. New four-arm GPU-time upper estimate: $6.0029. Prior tracked total plus new GPU-time estimate: $73.3118, <strong>plus unreconciled host CPU/RAM/build/storage charges</strong>. This excludes Claude’s separate lane and is not an account-wide invoice.</p>
<p>First pair: $3.1297; higher-prior pair: $2.8732. Both owned generation apps are verified stopped with zero tasks. No successor GPU run was launched. Packaging and context review used local CPU only.</p></div>
<small>Novelty recovery is complete and committed as 126a675. Cached alignment was reused; the exact export returned 3,523,855 rows. Its older 400-sequence cohort is separate from these outputs.</small></main></html>''')
    print(json.dumps(dict(state='complete', review_csv=str(submission), html=str(output), rows=20,
                          owned_gpu_apps_stopped=True, new_gpu_spend_usd=0)))


if __name__ == '__main__':
    main()
