#!/usr/bin/env python3
"""Summarise mouse audit by sequence SHA; no reordering or nomination."""
import argparse
import csv
import json
import subprocess
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]


def run(out):
    now=subprocess.check_output(['date','-u','+%Y-%m-%dT%H:%M:%SZ'],text=True).strip()
    raw=json.loads((out/'measurements.json').read_text())
    control=json.loads((out/'human_transfer_control.json').read_text())
    ctr={d['sequence_sha256']:d for d in control['designs']}
    tables={r['sequence_sha256']:r for r in csv.DictReader((out/'mouse_compatibility.csv').read_text().splitlines())}
    assert len(ctr)==len(tables)==len(raw['designs'])==20
    assessments=[];sweep=[];md=[]
    labels={
        'mouse_backbone_clear_both_fits':'Backbone transfer compatible in both fits; side-chain conflicts unresolved; cross-reactivity remains plausible but unestablished.',
        'mouse_conflict_human_reference_clear_both_fits':'Current static mouse pose incompatible at the 2.0 A backbone screen in both fits; matched human backbone clear. Modest differential warning, not biological exclusion.',
        'mouse_conflict_with_human_transfer_confounded':'Current static mouse pose has backbone conflicts in both fits; human-reference transfer also conflicts. Species attribution unresolved; no supported cross-reactivity exclusion.',
        'mouse_fit_sensitive_conflict':'Backbone compatibility depends on fit; side-chain conflicts unresolved. No robust cross-reactivity verdict.'}
    for d in raw['designs']:
        sid=d['sequence_sha256'];c=ctr[sid];t=tables[sid]
        changed=[str(r['human_aa'])+str(r['human_uniprot'])+' -> '+('/'.join('gap' if n is None else a+str(n) for n,a in zip(r['mouse_options'],r['mouse_aa_options'])))
                 for r in d['contact_residue_evidence'] if not r['invariant']]
        r=dict(t)
        r['contact_substitutions']='; '.join(changed)
        r['human_reference_domain_backbone_pairs_lt2A']=c['human_reference_fits']['domain']['target_backbone_pairs_lt2A']
        r['human_reference_local_backbone_pairs_lt2A']=c['human_reference_fits']['local']['target_backbone_pairs_lt2A']
        r['assessment']=labels[c['interpretation_label']]
        r['unchanged_pose_passes_full_atom_screen_both_fits']=all(f['pairs_lt2A']==0 for f in d['fits'].values())
        r['unchanged_pose_passes_full_atom_screen_any_fit']=any(f['pairs_lt2A']==0 for f in d['fits'].values())
        assert any(f['pairs_lt2A']>0 for f in d['fits'].values())
        r['biological_crossreactivity_excluded_by_this_test']=False
        r['mouse_binding_measured']=False
        assessments.append(r)
        for cutoff in (4.,4.5,5.):
            hs=[i+188 for i in d['human_baseline']['contact_residues'][str(cutoff)]]
            # Use all-optima invariance, not a preferred mouse alignment.
            # All positions outside the central cutoff derive identity from the canonical sequences.
            from Bio import SwissProt
            with (ROOT/'data/P02786.txt').open() as f:
                human=SwissProt.read(f).sequence
            with (ROOT/'data/Q62351.txt').open() as f:
                mouse=SwissProt.read(f).sequence
            inv=sum(all(n is not None and human[h-1]==mouse[n-1] for n in raw['mapping_options'][str(h)]) for h in hs)
            for threshold in (1.5,1.8,2.):
                mb=[d['fits'][f]['min_target_backbone_A']<threshold for f in ('domain','local')]
                hb=[c['human_reference_fits'][f]['min_target_backbone_A']<threshold for f in ('domain','local')]
                sweep.append(dict(sequence_sha256=sid,name=d['name'],human_contact_cutoff_A=cutoff,
                    contact_count=len(hs),invariant_contact_count=inv,invariant_fraction=inv/len(hs),
                    backbone_overlap_threshold_A=threshold,mouse_conflict_domain=mb[0],mouse_conflict_local=mb[1],
                    human_reference_conflict_domain=hb[0],human_reference_conflict_local=hb[1],
                    persistent_mouse_conflict=all(mb),human_transfer_confounded=all(mb) and any(hb)))
        md.append(f"| {d['name'].replace('tfr1_c4_acid_','')} | `{sid[:12]}` | {t['invariant_contact_count']}/{t['human_contact_count']} | {t['domain_target_backbone_pairs_lt2A']}/{t['local_target_backbone_pairs_lt2A']} | {r['human_reference_domain_backbone_pairs_lt2A']}/{r['human_reference_local_backbone_pairs_lt2A']} | {c['interpretation_label'].replace('_',' ')} |")
    for name,rows in (('per_design_assessment.csv',assessments),('threshold_sensitivity.csv',sweep)):
        p=out/name
        assert not p.exists(),p
        with p.open('w',newline='') as f:
            w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    contents=f'''# C4 mouse compatibility: fixed-pose audit

Written {now}. Evidence only. No ranking or submission was produced.

Twenty of twenty v4 sequences match unique actual binder chains by SHA256; every
target chain equals human P02786 residues 189-383. Mouse Q62351 is independently
sequence-verified. Retain all six full-length BLOSUM62 optimal alignments (gap open
-10, extension -0.5). Domain fits use 192 unambiguous CA pairs; local fits use the
unambiguous CA pairs within 12 A of any contacted human target CA, selected before
examining fit results. Only aligned apical mouse residues and human A189-383 enter
the transfer test. These are human-pose transfers onto a mouse monomer prediction,
not new mouse complex predictions or full native-context validation.

Human contacts use any heavy atom <=4.5 A. Conservation spans 60.87%-95.24% across
designs. All twenty have noninvariant contacts; fourteen contact K261 -> mouse E263
(a formal charge reversal). This does not establish a lost salt bridge, energetic
penalty, or pH response. F297 -> mouse L299 is frequent. site2_152 contacts the
ambiguous L212 -> gap/P214 region; every alignment is retained, never resolved to
improve the result. The old small hotspot summary does not describe every actual
v4 contact surface. See the full contact-level JSON and substitution lists in CSV.

At the prespecified <2.0 A threshold, 13/20 have mouse target-backbone overlaps in
both fits. Twelve of those also conflict with experimental human-reference apical
transfer in at least one fit: they cannot isolate a species effect. site1_185 has
the remaining differential backbone warning: mouse minima 1.915/1.784 A versus
human-reference 2.081/2.065 A, with only 2/1 mouse overlapping pairs. The human
reference still has side-chain overlaps, so this is not a human-pass/mouse-fail
binding test. Fits have different source conformations and are not independent.

Three designs are fit-sensitive. Four avoid target-backbone overlap under both fits:
site1_145, site2_172, site1_40 and site1_56. All four still have <2.0 A heavy-atom
overlaps in both fits, including local minima 0.826, 1.316, 0.661 and 0.452 A.
They remain backbone-compatible hypotheses with unresolved side-chain conflicts,
not clash-free designs. No design passes the full atom screen under BOTH fits.
site1_123 alone passes the local full-atom screen by 0.0015 A (minimum 2.0015 A),
while its domain fit has three backbone overlaps and minimum 1.565 A. This is a
fragile fit-dependent clearance, not robust compatibility. All other designs have
full-atom overlaps under both fits.

There is no evidence here that any sequence is biologically incapable of mouse
cross-reactivity. The current poses with persistent backbone overlaps cannot be
retained unchanged against this mouse model under this test. Resolving them would
require a pose or conformation change; neither rescue nor impossibility was tested.
Contact retention is proximity only: an overlapped atom can count as a contact.

Controls: twenty exact human self-transfers preserve contacts and minima; twenty
deliberately injected mouse backbone overlaps fire; twenty far-away controls have
no contacts; K205 -> S208 charge-loss mapping control fires. The matched experimental
human transfer is a separate scientific control, and exposes the transfer artifact.
These controls validate arithmetic and scope, not predictive accuracy.

Sensitivity: persistent mouse backbone conflicts are 9/20, 11/20 and 13/20 at
thresholds <1.5, <1.8 and <2.0 A. Of those, 8, 9 and 12 respectively also fail a human
transfer control. The mouse-only warning set changes with threshold; site1_185 is
not a cutoff-independent rejection. The companion table reports all 20 designs
at contact cutoffs 4.0/4.5/5.0 A and backbone thresholds 1.5/1.8/2.0 A (180 rows).
Do not choose a new cutoff based on which candidate it admits.

The table below preserves the input roster order; it is NOT a new ranking. M/H
are mouse/human-reference backbone overlap pair counts for domain/local fits.
Full molecule hashes and all changed contact residues are in per_design_assessment.csv.

| Design shorthand | Sequence SHA256 prefix | Invariant contacts | M domain/local | H domain/local | Interpretation |
|---|---|---:|---:|---:|---|
'''+ '\n'.join(md)+'''

Binding, KD, mouse cross-reactivity, protonation-dependent interface rearrangement,
pH7.4 binding/pH6.0 release, mouse pose relaxation, homodimer context, complete assay
glycans, and partner-chain accessibility remain unmeasured. Mouse/human reference
coordinates have different provenance; local fitting is a sensitivity check, not
independent validation. No new GPU job was launched.

Consultations: exact question and evidence digest are recorded under
results/decisions/20261026-c4-mouse-consult-*. Raw Astra and Kimi reviews are retained
separately. The evidence packet contains the full 20-row metrics and frozen v1 rule.

Recommendation: hold current C4 submissions. Keep v4 as the human-surrogate research
incumbent; do not promote the four backbone-clear rows into a new file. Geometry
and conservation do not resolve the required mouse/pH objectives. This is an
evidence limitation, not a claim of twenty nonbinders. Ricardo submits personally.

Reproduce in a fresh results directory, CPU only:

    OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 python3 scripts/c4_mouse_compatibility_20261026.py --out results/c4_mouse_reproduce_20261026
    OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 python3 scripts/c4_mouse_transfer_control_20261026.py --out results/c4_mouse_reproduce_20261026
    python3 scripts/c4_mouse_summary_20261026.py --out results/c4_mouse_reproduce_20261026

The scripts refuse to overwrite measurements or summary files. JSON input digests
bind all 20 original complex files, sequences, reference structures and the rule.
Biopython 1.85 and NumPy 2.0.2 were used with the installed Python 3.9; .venv does
not exist on this checkout. No historical results directory was rerun or overwritten.
'''
    assert not (out/'RESULTS.md').exists()
    (out/'RESULTS.md').write_text(contents)
    print('Wrote 20 assessments and 180 prespecified sensitivity rows; original order preserved')


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--out',type=Path,required=True)
    run(ap.parse_args().out)
