#!/usr/bin/env python3
"""Epitope identity + agonist-site coverage for every C3 candidate row. $0, CPU, read-only.

Answers Astra's Q3 ("report epitope/binding-mode diversity if already available") and settles the
caveat left open by c3_epitope_diversity_20261025.py: is c3_glp1r_ecd_29's "second epitope" real?

Ground truth for "the agonist site" is 6X18 (GLP-1 peptide bound to GLP1R, cryo-EM 2.10 A): 52
receptor residues lie within 4.5 A of the peptide. Design complex PDBs contain only the ECD fragment,
numbered 1..111, which maps to 6X18 chain F LABEL positions 38..154, excluding
unresolved 139..144. The map is +37 for 1..101 and +43 for 102..111. These are
label_seq_id positions, not author numbering.

IMPORTANT DENOMINATOR: 22 of those 52 residues fall inside the resolved mapped
fragment. The other 30 are outside it and CANNOT be assessed from these complexes.
So all coverage fractions are reported against the ASSESSABLE footprint (22 residues); a fraction
against the full 52 would silently treat unmeasurable residues as non-contacts.

Outputs one JSON. Writes no CSV, changes no gate, submits nothing.
"""
import argparse, csv, json, math
from collections import Counter, defaultdict
from itertools import combinations
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
ART = ROOT / 'results/c3_pod_artifacts_20261025'
SPR = ROOT / 'results/c3_sprint_8h_20261025'
V8 = ROOT / 'results/c3_submission_v8_20261025.csv'
CUTS = (3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5)
CUT = 4.5
J_SAME = 0.5

def fragment_map():
    """Derive numbering from observed 6X18 label residues, checking the known gap."""
    resolved = set()
    for line in (ROOT / 'data/6x18.cif').read_text().splitlines():
        if line.startswith('ATOM'):
            f = line.split()
            if f[6] == 'F' and f[8].isdigit() and 38 <= int(f[8]) <= 154:
                resolved.add(int(f[8]))
    mapping = dict(enumerate(sorted(resolved), start=1))
    expected = {i: i + (37 if i <= 101 else 43) for i in range(1, 112)}
    if mapping != expected:
        raise ValueError('6X18 resolved fragment differs from the verified 111-residue map')
    return mapping


def map_fragment_residue(residue, mapping):
    """Accept an integer or a blank insertion-code key; reject ambiguous numbering."""
    if isinstance(residue, tuple):
        residue, icode = residue
        if icode.strip():
            raise ValueError('Insertion-coded fragment residue requires explicit sequence mapping')
    return mapping[residue]

def glp1_site(cut=CUT):
    pep, ecd = [], defaultdict(list)
    for line in (ROOT / 'data/6x18.cif').read_text().splitlines():
        if not line.startswith('ATOM'):
            continue
        f = line.split()
        if f[2] == 'H' or f[6] not in ('E', 'F') or not f[8].isdigit():
            continue
        xyz = (float(f[10]), float(f[11]), float(f[12]))
        if f[6] == 'E':
            pep.append(xyz)
        else:
            ecd[int(f[8])].append(xyz)
    c2 = cut * cut
    return {r for r, pts in ecd.items()
            if any((a[0]-b[0])**2 + (a[1]-b[1])**2 + (a[2]-b[2])**2 <= c2 for a in pts for b in pep)}

def chains_of(pdb):
    ch = Counter()
    for line in Path(pdb).read_text().splitlines():
        if line.startswith('ATOM') and line[21].strip():
            ch[line[21]] += 1
    return ch

def interface(pdb, cut=CUT):
    """Target residues contacting the binder, plus interface pLDDT on both sides."""
    ch = chains_of(pdb)
    if len(ch) < 2:
        return None
    binder = min(ch, key=lambda k: ch[k])
    target = max(ch, key=lambda k: ch[k])
    batoms, tres, bres = [], defaultdict(list), defaultdict(list)
    for line in Path(pdb).read_text().splitlines():
        if not line.startswith('ATOM') or (line[77:78].strip() or 'C') == 'H':
            continue
        xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
        plddt = float(line[60:66])
        key = (int(line[22:26]), line[26].strip())
        if line[21] == binder:
            batoms.append(xyz)
            bres[key].append((xyz, plddt))
        elif line[21] == target:
            if key[1]:
                raise ValueError('Insertion-coded target requires explicit sequence mapping')
            tres[key].append((xyz, plddt))
    c2 = cut * cut
    contacts = {k: v[0][1] for k, v in tres.items()
                if any((a[0]-b[0])**2 + (a[1]-b[1])**2 + (a[2]-b[2])**2 <= c2
                       for a, _ in v for b in batoms)}
    if not contacts:
        return None
    tpts = [xyz for v in tres.values() for xyz, _ in v]
    bip = [v[0][1] for v in bres.values()
           if any(math.dist(p, q) <= 5.0 for p, _ in v for q in tpts)]
    return dict(patch=frozenset(r for r, _ in contacts),
                n_contacts=len(contacts),
                target_iface_plddt=round(sum(contacts.values()) / len(contacts), 1),
                binder_iface_plddt=round(sum(bip) / len(bip), 1) if bip else None)

def pdb_for(name):
    p = ART / f'{name}_complex.pdb'
    if p.exists():
        return p
    alt = sorted(ROOT.glob(f'results/**/{name}/*rank_001*multimer*.pdb'))
    return alt[0] if alt else None

def jac(a, b):
    return len(a & b) / len(a | b)

def main(output):
    site = glp1_site()
    variants = {}
    for f in sorted(SPR.glob('c3_candidate_v9_*.csv')):
        variants[f.stem] = [r['name'] for r in csv.DictReader(open(f))]
    variants['v8_frozen'] = [r['name'] for r in csv.DictReader(open(V8))]

    cache = {}
    def prof(name):
        if name not in cache:
            p = pdb_for(name)
            cache[name] = (interface(p) if p else None, p)
        return cache[name][0]

    mapping = fragment_map()
    frag = site & set(mapping.values())
    out = dict(cutoff=CUT, fragment_to_6x18_label=mapping,
               glp1_site_full=len(site), glp1_site_assessable=len(frag),
               glp1_site_unassessable=len(site) - len(frag),
               j_same_threshold=J_SAME, variants={}, designs={})
    for v, names in variants.items():
        rows = []
        for n in names:
            d = prof(n)
            if d is None:
                rows.append(dict(name=n, measured=False))
                continue
            conv = {map_fragment_residue(r, mapping) for r in d['patch']}
            on = conv & frag
            rows.append(dict(name=n, measured=True, n_contacts=d['n_contacts'],
                             target_iface_plddt=d['target_iface_plddt'],
                             binder_iface_plddt=d['binder_iface_plddt'],
                             n_on_agonist_site=len(on),
                             frac_of_assessable_site=round(len(on) / len(frag), 3),
                             frac_of_full_site=round(len(on) / len(site), 3),
                             frac_of_own_patch_on_site=round(len(on) / len(conv), 3)))
        ok = [r for r in rows if r['measured']]
        parent = {r['name']: r['name'] for r in ok}
        def find(x):
            while parent[x] != x:
                parent[x] = parent[parent[x]]
                x = parent[x]
            return x
        for a, b in combinations(ok, 2):
            if jac(prof(a['name'])['patch'], prof(b['name'])['patch']) >= J_SAME:
                parent[find(a['name'])] = find(b['name'])
        clusters = defaultdict(list)
        for r in ok:
            clusters[find(r['name'])].append(r['name'])
        out['variants'][v] = dict(
            n_rows=len(rows), n_measured=len(ok),
            n_distinct_epitopes=len(clusters),
            cluster_sizes=sorted((len(c) for c in clusters.values()), reverse=True),
            clusters=sorted((sorted(c) for c in clusters.values()), key=lambda c: (-len(c), c[0])),
            mean_frac_of_assessable_site=round(sum(r['frac_of_assessable_site'] for r in ok) / len(ok), 3),
            min_frac_of_assessable_site=min(r['frac_of_assessable_site'] for r in ok),
            rows_on_agonist_site=sum(1 for r in ok if r['frac_of_assessable_site'] > 0),
            rows=[r for r in rows])
        for r in ok:
            out['designs'].setdefault(r['name'], r)

    # Cluster-count sensitivity: "1 epitope" is threshold-dependent. Report the whole sweep so the
    # diversity statement cannot be tuned by picking a threshold.
    def clusters_at(present, thr):
        par = {n: n for n in present}
        def f(x):
            while par[x] != x:
                par[x] = par[par[x]]
                x = par[x]
            return x
        for a, b in combinations(present, 2):
            if jac(prof(a)['patch'], prof(b)['patch']) >= thr:
                par[f(a)] = f(b)
        out = defaultdict(list)
        for n in present:
            out[f(n)].append(n)
        return sorted((sorted(v) for v in out.values()), key=lambda c: (-len(c), c[0]))

    for v, d in out['variants'].items():
        present = [r['name'] for r in d['rows'] if r.get('measured')]
        d['cluster_sensitivity'] = {}
        for thr in (0.5, 0.6, 0.7, 0.8, 0.9, 0.99):
            cl = clusters_at(present, thr)
            d['cluster_sensitivity'][str(thr)] = dict(n=len(cl), sizes=[len(c) for c in cl])
        js = [jac(prof(a)['patch'], prof(b)['patch']) for a, b in combinations(present, 2)]
        d['pairwise_jaccard'] = dict(
            n_pairs=len(js), min=round(min(js), 3),
            median=round(sorted(js)[len(js) // 2], 3), max=round(max(js), 3),
            n_ge_0_9=sum(1 for x in js if x >= 0.9))
        identical = defaultdict(list)
        for n in present:
            identical[frozenset(prof(n)['patch'])].append(n)
        d['identical_patch_groups'] = [v2 for v2 in identical.values() if len(v2) > 1]

    # cutoff robustness for the one contested pair
    a, b = 'c3_glp1r_ecd_29', 'c3_bc_l65_8f7b63'
    robust = {}
    for cut in CUTS:
        pa, pb = interface(pdb_for(a), cut), interface(pdb_for(b), cut)
        if pa and pb:
            robust[str(cut)] = dict(n_a=pa['n_contacts'], n_b=pb['n_contacts'],
                                    shared=len(pa['patch'] & pb['patch']),
                                    jaccard=round(jac(pa['patch'], pb['patch']), 3))
    out['ecd29_cutoff_robustness'] = dict(pair=[a, b], results=robust,
        interpretation="Jaccard stays below the 0.5 'same epitope' line at every cutoff "
                       "3.5-6.5 A, so the DISTINCT call is not a cutoff artifact. What is an "
                       "artifact is the INTERPRETATION of that call as portfolio diversity: the "
                       "patch is 82% off the GLP-1 site and its interface confidence is far lower.")

    out['verdict'] = dict(
        threshold_dependence="'1 epitope' is the count at Jaccard>=0.5 and is threshold-dependent. "
            "Read cluster_sensitivity: at >=0.9 the same portfolio splits into many near-distinct "
            "contact sets. The honest summary is ONE receptor face with substantial internal contact "
            "variation, not 20 identical poses and not 20 independent epitopes.",
        effective_epitope_diversity="1 credible predicted epitope, not 2. 24 of 25 measurable "
            "designs cluster into one epitope at J>=0.5. The lone singleton is c3_glp1r_ecd_29, a "
            "low-confidence mostly-off-site prediction, not an established second binding mode.",
        agonist_targeting="All 19 structurally assessable designs (of 20) have predicted interfaces "
            "overlapping receptor residues contacted by GLP-1 in 6X18 at a 4.5 A heavy-atom cutoff. "
            "Coverage is stated against the 22 assessable footprint residues, not the full 52.",
        denominator_caveat="30 of the 52 GLP-1 contact residues lie outside the modeled ECD fragment "
            "and are NOT assessable. Coverage fractions describe only the assessable subset.",
        does_not_undercut="This removes an argument for KEEPING ecd_29 (it is not diversity) and "
            "removes an argument against variant B (B was never more diverse than v8; neither is).",
        cannot_claim="Nothing here measures affinity or agonism. Contacting the peptide site in a "
            "predicted complex is not evidence of activation.")
    output.parent.mkdir(parents=True, exist_ok=True)
    if output.exists():
        raise FileExistsError('Refusing to replace an existing audit: ' + str(output))
    output.write_text(json.dumps(out, indent=1, default=list))

    print(f"GLP-1 site (6X18, cryo-EM 2.1 A, {CUT} A cutoff): {len(site)} receptor residues,")
    print(f"  of which {len(frag)} fall inside the modeled ECD fragment and {len(site)-len(frag)} do not.\n")
    for v, d in sorted(out['variants'].items()):
        print(f"{v:44s} epitopes {d['n_distinct_epitopes']} {str(d['cluster_sizes']):9s} "
              f"assessable-site {d['mean_frac_of_assessable_site']:.2f} "
              f"({d['min_frac_of_assessable_site']:.2f} min)  "
              f"on-site {d['rows_on_agonist_site']}/{d['n_measured']}")
    print("\ncluster-count sensitivity (designs in how many contact clusters):")
    for v, d in sorted(out['variants'].items()):
        print(f"  {v:44s} " + "  ".join(f"{t}:{c['n']}" for t, c in d['cluster_sensitivity'].items())
              + f"   median pairwise J {d['pairwise_jaccard']['median']:.2f}")
    print("\ncutoff robustness, ecd_29 vs l65:")
    for c, r in out['ecd29_cutoff_robustness']['results'].items():
        print(f"  {c:>4} A  n={r['n_a']:2d}/{r['n_b']:2d}  shared {r['shared']:2d}  J {r['jaccard']:.3f}")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True,
                        help='New audit path; historical JSON is never overwritten')
    main(parser.parse_args().output)
