#!/usr/bin/env python3
"""Does the C4 apical preview's ranking criterion predict acid-ON linkage? (CPU only)

Kimi 05:49 / Astra 05:49 both ruled: freeze C4 GPU work; audit whether the pKa-based
selection establishes pH-conditional binding. The 06:04 linkage audit answered "what is L
for the 20 preview designs" but NOT "does the selection rule enrich for L>0". This script
answers the second question against the whole 343-design pool.

Fail-closed:
  * sign controls must pass first (same arithmetic as c4_linkage_audit_20261025)
  * a design with a missing/unreadable complex is reported as MISSING, never scored as 0
  * no GPU, no network, read-only on inputs
Writes only results/c4_selection_audit_20261025/.
"""
import csv, json, math, statistics, sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / 'scripts'))
import c4_linkage_audit_20261025 as A  # noqa: E402

RESULTS = ROOT / 'results'
OUT = RESULTS / 'c4_selection_audit_20261025'
PREVIEW = RESULTS / 'c4_apical_preview_20261024.csv'
CONS = RESULTS / 'c4_apical_consolidated_results.json'
LINKS = RESULTS / 'c4_linkage_audit_20261025' / 'audit.json'
PDBS = RESULTS / 'c4_rebuild_20261024' / 'remote_pdbs'
IPTM_GATE = 0.60


def spearman(x, y):
    """Average-rank Spearman; ties handled. Returns rho and two-sided normal p."""
    n = len(x)
    if n < 3:
        return None, None

    def ranks(v):
        order = sorted(range(n), key=lambda i: v[i])
        r = [0.0] * n
        i = 0
        while i < n:
            j = i
            while j + 1 < n and v[order[j + 1]] == v[order[i]]:
                j += 1
            avg = (i + j) / 2.0 + 1
            for k in range(i, j + 1):
                r[order[k]] = avg
            i = j + 1
        return r

    rx, ry = ranks(x), ranks(y)
    mx, my = statistics.mean(rx), statistics.mean(ry)
    num = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    den = math.sqrt(sum((a - mx) ** 2 for a in rx) * sum((b - my) ** 2 for b in ry))
    if den == 0:
        return 0.0, 1.0
    rho = num / den
    # Fisher z approximation; adequate for n=116, reported as indicative only.
    z = math.sqrt(n - 3) * math.atanh(max(-0.999999, min(0.999999, rho)))
    p = 2 * (1 - 0.5 * (1 + math.erf(abs(z) / math.sqrt(2))))
    return rho, p


def fisher_two_sided(a, b, c, d):
    """2x2 table [[a,b],[c,d]]; exact two-sided Fisher p."""
    from math import comb
    n = a + b + c + d
    r1, c1 = a + b, a + c
    def prob(k):
        return comb(r1, k) * comb(n - r1, c1 - k) / comb(n, c1)
    p0 = prob(a)
    lo, hi = max(0, c1 - (n - r1)), min(r1, c1)
    return min(1.0, sum(prob(k) for k in range(lo, hi + 1) if prob(k) <= p0 + 1e-12))


def main():
    ctl = A.controls()
    print('sign controls:', json.dumps(ctl))
    if not ctl.get('all_ok'):
        raise SystemExit('SIGN CONTROLS FAILED - audit invalid')

    cons = json.loads(CONS.read_text())
    by_name = {r['name']: r for r in cons['results']}
    preview = list(csv.DictReader(open(PREVIEW)))
    preview_names = [r['name'] for r in preview]
    prior = {r['name']: r for r in json.loads(LINKS.read_text())['records']}

    records, missing = [], []
    for name, meta in by_name.items():
        pdb = PDBS / (name + '_complex.pdb')
        row = dict(name=name,
                   iptm=meta.get('iptm'),
                   avg_complex_his_pka=meta.get('avg_complex_his_pka'),
                   n_elevated_complex=meta.get('n_elevated_complex'),
                   n_his_sequence=str(meta.get('sequence', '')).count('H'),
                   in_preview=name in preview_names)
        if not pdb.exists():
            row['status'] = 'MISSING'; row['L6'] = None
            missing.append(name)
        else:
            s = A.score(str(pdb), 'B', 'A', 6.0, 7.4)
            if s is None:
                row['status'] = 'PROPKA_FAILED'; row['L6'] = None
            else:
                row['status'] = 'scored'; row['L6'] = round(s['L'], 4)
                row['n_groups'] = s['n_groups']
        # Prefer the 06:04 audit's value for preview rows so the two reports cannot disagree.
        if name in prior and prior[name].get('L_6.0_7.4') is not None:
            row['L6'] = prior[name]['L_6.0_7.4']
        records.append(row)

    scored = [r for r in records if r.get('L6') is not None]
    elig = [r for r in scored if (r.get('iptm') or 0) >= IPTM_GATE]
    pv = [r for r in scored if r['in_preview']]
    non_pv_elig = [r for r in elig if not r['in_preview']]
    all_non_pv = [r for r in scored if not r['in_preview']]

    def rate(rs):
        return (sum(1 for r in rs if r['L6'] > 0), len(rs))

    pv_pos, pv_n = rate(pv)
    en_pos, en_n = rate(non_pv_elig)
    al_pos, al_n = rate(all_non_pv)

    rho_his, p_his = spearman([r['avg_complex_his_pka'] for r in elig],
                              [r['L6'] for r in elig])
    rho_iptm, p_iptm = spearman([r['iptm'] for r in elig],
                                [r['L6'] for r in elig])

    # Deterministic replay of the documented selection rule. The original builder lived in
    # /tmp and is not frozen, so this tests reproducibility; it does not overwrite the preview.
    def kmer_set(seq, k=4):
        return {seq[i:i + k] for i in range(max(0, len(seq) - k + 1))}
    def jac(a, b):
        if not a or not b:
            return 0.0
        return len(a & b) / len(a | b)

    ordered = sorted(elig, key=lambda r: -(r['avg_complex_his_pka'] or 0))
    picked, picked_sets = [], []
    site1 = site2 = 0
    for r in ordered:
        if len(picked) >= 20:
            break
        is1 = '_site1_' in r['name']
        if is1 and site1 >= 8:
            continue
        if (not is1) and site2 >= 12:
            continue
        ks = kmer_set(by_name[r['name']]['sequence'])
        if any(jac(ks, s) >= 0.60 for s in picked_sets):
            continue
        picked.append(r['name']); picked_sets.append(ks)
        site1 += is1; site2 += (not is1)
    replay = set(picked)
    overlap = len(replay & set(preview_names))

    summary = dict(
        audit='c4_selection_audit_20261025',
        question='Does the preview pKa ranking criterion enrich for acid-ON linkage (L6>0)?',
        sign_controls=ctl,
        n_consolidated=len(records), n_scored=len(scored), n_missing=len(missing),
        missing=missing, iptm_gate=IPTM_GATE, n_eligible=len(elig),
        preview_n=len(pv), preview_positive=pv_pos, preview_rate=round(pv_pos / pv_n, 4) if pv_n else None,
        eligible_nonpreview_n=len(non_pv_elig), eligible_nonpreview_positive=en_pos,
        eligible_nonpreview_rate=round(en_pos / en_n, 4) if en_n else None,
        all_nonpreview_n=len(all_non_pv), all_nonpreview_positive=al_pos,
        all_nonpreview_rate=round(al_pos / al_n, 4) if al_n else None,
        fisher_preview_vs_eligible_nonpreview=fisher_two_sided(pv_pos, pv_n - pv_pos, en_pos, en_n - en_pos),
        fisher_preview_vs_all_nonpreview=fisher_two_sided(pv_pos, pv_n - pv_pos, al_pos, al_n - al_pos),
        spearman_his_pka_vs_L6_eligible=[round(rho_his, 4) if rho_his is not None else None,
                                         round(p_his, 4) if p_his is not None else None],
        spearman_iptm_vs_L6_eligible=[round(rho_iptm, 4) if rho_iptm is not None else None,
                                      round(p_iptm, 4) if p_iptm is not None else None],
        replay_rule_overlap_of_20=overlap,
        replay_rule_names=picked,
        preview_mean_L6=round(statistics.mean([r['L6'] for r in pv]), 4) if pv else None,
        eligible_mean_L6=round(statistics.mean([r['L6'] for r in elig]), 4) if elig else None,
        interpretation=None,
    )
    # Plain-language verdict, computed from the numbers, not asserted in advance.
    pr, er = summary['preview_rate'], summary['eligible_nonpreview_rate']
    if pr is None or er is None:
        verdict = 'INCONCLUSIVE: insufficient scored coverage.'
    elif summary['fisher_preview_vs_eligible_nonpreview'] >= 0.05:
        verdict = ('NO MEASURED ENRICHMENT: the preview selects acid-ON designs at %.1f%% vs %.1f%% '
                   'for other ipTM>=0.60 designs (Fisher p=%.3f). The pKa ranking does not '
                   'demonstrably select the conditional-binding direction.' % (pr * 100, er * 100,
                   summary['fisher_preview_vs_eligible_nonpreview']))
    elif pr > er:
        verdict = ('ENRICHED: preview %.1f%% vs %.1f%% (Fisher p=%.3f). Model-level association only; '
                   'not experimental evidence.' % (pr * 100, er * 100,
                   summary['fisher_preview_vs_eligible_nonpreview']))
    else:
        verdict = ('ANTI-ENRICHED: preview %.1f%% vs %.1f%% (Fisher p=%.3f). The criterion is '
                   'selecting against the stated objective.' % (pr * 100, er * 100,
                   summary['fisher_preview_vs_eligible_nonpreview']))
    summary['interpretation'] = verdict + ' PROPKA linkage on model structures; no experimental claim.'

    OUT.mkdir(parents=True, exist_ok=True)
    json.dump(dict(summary=summary, records=records), open(OUT / 'audit.json', 'w'), indent=2)
    with open(OUT / 'audit.csv', 'w', newline='') as fh:
        w = csv.writer(fh)
        w.writerow(['name', 'in_preview', 'iptm', 'avg_complex_his_pka', 'n_his_sequence', 'L6', 'status'])
        for r in sorted(scored, key=lambda z: -(z.get('L6') or -999)):
            w.writerow([r['name'], r['in_preview'], r.get('iptm'), r.get('avg_complex_his_pka'),
                        r.get('n_his_sequence'), r.get('L6'), r.get('status')])
    md = ['# C4 selection audit — does the pKa ranking predict acid-ON linkage?', '',
          '**Status: evidence only. Nothing promoted, no candidate replaced, nothing submitted.**', '',
          '## Controls', '```', json.dumps(ctl), '```', '',
          'CORRECTION 2026-10-26: positive L6 = tighter binding at pH 6.0 than 7.4.',
          'That is ACID-ON. The C4 objective is ACID-OFF (KD at 7.4, no binding at 6.0),',
          'so the C4 objective is L6 < 0, not L6 > 0. The original sentence in this file',
          'read "... = the C4 objective" and inverted the gate; see',
          'results/c4_independent_check_20261026/FINDING_L6_SIGN.md. Everything below is',
          'an audit of the acid-ON ranking as it was actually built, not a C4 endorsement.', '',
          '## Coverage', '| metric | value |', '|---|---|',
          '| consolidated designs | %d |' % len(records),
          '| scored | %d |' % len(scored),
          '| missing/unreadable | %d |' % len(missing),
          '| eligible (ipTM >= %.2f) | %d |' % (IPTM_GATE, len(elig)), '',
          '## Result', '| comparison | L6>0 | n | rate |', '|---|---|---|---|',
          '| preview (20) | %d | %d | %s |' % (pv_pos, pv_n, summary['preview_rate']),
          '| eligible non-preview | %d | %d | %s |' % (en_pos, en_n, summary['eligible_nonpreview_rate']),
          '| all non-preview | %d | %d | %s |' % (al_pos, al_n, summary['all_nonpreview_rate']), '',
          '- Fisher exact, preview vs eligible non-preview: **p=%.4f**' % summary['fisher_preview_vs_eligible_nonpreview'],
          '- Fisher exact, preview vs all non-preview: **p=%.4f**' % summary['fisher_preview_vs_all_nonpreview'],
          '- Spearman avg_complex_his_pka vs L6 (eligible): **rho=%s, p=%s**' % tuple(summary['spearman_his_pka_vs_L6_eligible']),
          '- Spearman ipTM vs L6 (eligible): **rho=%s, p=%s**' % tuple(summary['spearman_iptm_vs_L6_eligible']),
          '- Preview mean L6: **%s**; eligible mean L6: **%s**' % (summary['preview_mean_L6'], summary['eligible_mean_L6']),
          '- Deterministic replay of the documented rule reproduces **%d/20** preview names.' % overlap, '',
          '## Verdict', '', verdict, '']
    (OUT / 'audit.md').write_text('\n'.join(md))
    print('\n'.join(md))
    print('wrote', OUT)


if __name__ == '__main__':
    main()
