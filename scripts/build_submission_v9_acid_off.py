"""Build C4 submission v9 — ACID-OFF scoring (corrected pH direction).

Challenge 4 spec: bind at pH 7.4, NO binding at pH 6.0.
Rank by: acid-OFF switch score (His near cationic, buried His) + ipTM.
Penalize His near anionic (wrong-direction stabilization).

Tiers:
  T1: switch_score > 0, ipTM >= 0.75, no anionic contacts
  T2: switch_score > 0, ipTM >= 0.5
  T3: high ipTM with buried His (some acid-OFF potential)
  T4: highest ipTM controls (no clear switching)
"""
import argparse
import csv
import json


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--scored', required=True)
    p.add_argument('--max', type=int, default=20)
    p.add_argument('--out', default='results/submission_v9.csv')
    args = p.parse_args()

    with open(args.scored) as f:
        designs = json.load(f)

    # T1: strong acid-OFF + high ipTM + clean (no anionic)
    t1 = [d for d in designs
          if d['total_switch_score'] > 2
          and d['iptm'] >= 0.75
          and d['n_anionic_contacts'] == 0]
    t1.sort(key=lambda x: (-x['total_switch_score'], -x['iptm']))

    # T2: acid-OFF + decent ipTM (allows some anionic)
    t1_names = set(d['name'] for d in t1)
    t2 = [d for d in designs
          if d['name'] not in t1_names
          and d['total_switch_score'] > 0
          and d['iptm'] >= 0.5]
    t2.sort(key=lambda x: (-x['total_switch_score'], -x['iptm']))

    # T3: high ipTM + buried His
    used = t1_names | set(d['name'] for d in t2)
    t3 = [d for d in designs
          if d['name'] not in used
          and d['iptm'] >= 0.8
          and d['n_buried_his'] >= 2]
    t3.sort(key=lambda x: (-x['n_buried_his'], -x['iptm']))

    # T4: highest ipTM controls
    used2 = used | set(d['name'] for d in t3)
    t4 = [d for d in designs
          if d['name'] not in used2
          and d['iptm'] >= 0.9]
    t4.sort(key=lambda x: -x['iptm'])

    selected = []
    tier_labels = [
        ('T1:acid-off+highIPTM', t1),
        ('T2:acid-off+midIPTM', t2),
        ('T3:buried+highIPTM', t3),
        ('T4:control', t4),
    ]
    for label, tier in tier_labels:
        for d in tier:
            if len(selected) >= args.max:
                break
            d['tier'] = label
            selected.append(d)

    fmt = '{:<40s} {:>6s} {:>4s} {:>4s} {:>4s} {:>5s} {:>6s}  {}'
    print(fmt.format('Name', 'ipTM', 'Cat', 'Ani', 'Bur', 'Sw', 'pLDDT', 'Tier'))
    print('-' * 100)
    for s in selected:
        print(fmt.format(
            s['name'],
            '{:.3f}'.format(s['iptm']),
            str(s['n_cationic_contacts']),
            str(s['n_anionic_contacts']),
            str(s['n_buried_his']),
            '{:.1f}'.format(s['total_switch_score']),
            '{:.1f}'.format(s['mean_plddt']),
            s.get('tier', '')))

    counts = {}
    for s in selected:
        t = s.get('tier', 'unknown').split(':')[0]
        counts[t] = counts.get(t, 0) + 1
    print('\nTier counts: {}'.format(
        ', '.join('{}={}'.format(k, v) for k, v in sorted(counts.items()))))

    rows = []
    for s in selected:
        rows.append({
            'name': s['name'],
            'sequence': s['sequence'],
            'molecule_class': 'single_chain',
        })

    with open(args.out, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=['name', 'sequence', 'molecule_class'])
        w.writeheader()
        w.writerows(rows)
    print('\nSaved {} designs to {}'.format(len(rows), args.out))


if __name__ == '__main__':
    main()
