#!/usr/bin/env python
"""Install ranked two-state histidine sites and matched glutamine controls.

Candidate order and exclusions come exclusively from switch_sites.score_design.
The default models one pair from an eight-site shortlist. Existing histidines can
occupy a site: HH retains them and QQ removes their titration. Pair estimates only
order modelling work; relaxed PROPKA linkage, including acids, decides the pH gate.

--plan-only needs only Python's standard library. Modelling needs PyRosetta on CPU.
The complex is first repacked/minimised with a fixed backbone. The reference and
both twins then undergo the same fixed-backbone local repack/minimise.
Rosetta ddG is a structural filter in REU, not a pH-specific binding measurement.
"""
import argparse
import hashlib
import itertools
import json
import math
from pathlib import Path

from switch_sites import score_design

AA1 = dict(zip(
    'ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL'.split(),
    'ARNDCQEGHILKMFPSTWYV'))


def make_plan(pdb, binder_chain='B', target_chain='A', max_positions=8,
              sites_per_variant=2, max_pairs=1):
    """Use the scorer's order unchanged; never substitute a burial sort."""
    if len(binder_chain) != 1 or not target_chain or binder_chain in target_chain:
        raise ValueError('one binder chain and disjoint target chain(s) are required')
    if min(max_positions, sites_per_variant, max_pairs) < 1:
        raise ValueError('position, site and pair limits must be positive')
    if sites_per_variant > max_positions:
        raise ValueError('sites per variant exceeds the position shortlist')
    # The scorer keys residues by (chain, integer number); refuse ambiguous PDBs.
    identifiers = set()
    with open(pdb) as fh:
        for line in fh:
            if line.startswith('ATOM'):
                if line[26:27].strip() or line[16:17].strip() not in ('', 'A'):
                    raise ValueError('insertion codes/alternate conformers need canonicalisation')
                identifiers.add(line[21])
    if not set(target_chain + binder_chain) <= identifiers:
        raise ValueError('requested chains are absent from the input PDB')
    unassigned = identifiers - set(target_chain + binder_chain)
    if unassigned:
        raise ValueError('assign all protein chains to a partner: %s' % sorted(unassigned))
    ranked = score_design(str(pdb), binder_chain, target_chain)
    eligible = [r for r in ranked if not r['note'] and r['resname'] in AA1
                and r['est_contribution'] < 0 and r['d_burial'] > 0]
    chosen = eligible[:max_positions]
    sets = list(itertools.combinations(chosen, sites_per_variant))
    sets.sort(key=lambda pair: sum(r['est_contribution'] for r in pair))
    pairs = []
    for pair in sets[:max_pairs]:
        sites = [dict(r, wt=AA1[r['resname']]) for r in pair]
        pairs.append(dict(sites=sites,
                          estimated_site_sum=round(sum(r['est_contribution'] for r in pair), 3)))
    return dict(source_pdb=str(Path(pdb).resolve()),
                source_sha256=hashlib.sha256(Path(pdb).read_bytes()).hexdigest(),
                code_sha256={name: hashlib.sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
                             for name in ('install_switch.py', 'switch_sites.py',
                                          'ph_linkage.py', 'run_switch_campaign.py')},
                binder_chain=binder_chain, target_chain=target_chain,
                ranking_method='switch_sites.score_design',
                estimate_caveat='Coordinate heuristic, not measured pKa or total complex linkage; '
                                'site sums assume additivity only to order modelling.',
                max_positions=max_positions, sites_per_variant=sites_per_variant,
                max_pairs=max_pairs, ranking=ranked, shortlist=chosen, pairs=pairs)


def init(seed=20261021):
    global pyrosetta
    import pyrosetta
    if not pyrosetta.rosetta.basic.was_init_called():
        pyrosetta.init(' '.join([
            '-mute all', '-ex1', '-ex2aro', '-use_input_sc',
            '-ignore_unrecognized_res', '-ignore_zero_occupancy false',
            '-constant_seed', '-jran %d' % seed,
        ]))
    pyrosetta.rosetta.numeric.random.rg().set_seed(seed)


def map_sites(pose, sites):
    """PDB numbers are not pose indices (C24 follows the fragmented target)."""
    mapped = []
    info = pose.pdb_info()
    for site in sites:
        i = info.pdb2pose(site['chain'], site['resnum'])
        if not i or info.chain(i) != site['chain'] or info.number(i) != site['resnum']:
            raise ValueError('PDB site not present in pose: %s%d' % (site['chain'], site['resnum']))
        if pose.residue(i).name1() != site['wt']:
            raise ValueError('residue identity mismatch at %s%d' % (site['chain'], site['resnum']))
        mapped.append(i)
    return mapped


def dG_separated(pose, target_chain='A', binder_chain='B', scorefxn=None):
    from pyrosetta.rosetta.protocols.analysis import InterfaceAnalyzerMover
    partners = '%s_%s' % (target_chain, binder_chain)
    # Rosetta changed this setter from str to DockingPartners in newer releases.
    partner_type = getattr(pyrosetta.rosetta.core.pose, 'DockingPartners', None)
    if partner_type is not None:
        partners = partner_type.docking_partners_from_string(partners)
    iam = InterfaceAnalyzerMover()
    iam.set_interface(partners)
    iam.set_scorefunction(scorefxn or pyrosetta.get_fa_scorefxn())
    iam.set_compute_interface_energy(True)
    iam.set_compute_packstat(False)
    iam.set_calc_dSASA(True)
    iam.set_pack_separated(True)
    iam.apply(pose.clone())
    dg = float(iam.get_separated_interface_energy())
    sasa = float(iam.get_interface_delta_sasa())
    if not math.isfinite(dg) or not math.isfinite(sasa) or sasa <= 0:
        raise ValueError('invalid interface energy/area')
    return dg, sasa


def mutate_and_relax(pose, positions, aas, scorefxn, radius=8.0):
    from pyrosetta.toolbox import mutate_residue
    if len(positions) != len(aas) or not positions:
        raise ValueError('one mutation identity is required per selected site')
    work = pose.clone()
    for resi, aa in zip(positions, aas):
        # Even a native H or Q gets the same repacking as its matched twin.
        mutate_residue(work, resi, aa, pack_radius=radius, pack_scorefxn=scorefxn)
    mm = pyrosetta.rosetta.core.kinematics.MoveMap()
    mm.set_bb(False)
    mm.set_chi(False)
    centres = [work.residue(i).nbr_atom_xyz() for i in positions]
    for j in range(1, work.total_residue() + 1):
        if any(c.distance(work.residue(j).nbr_atom_xyz()) < radius for c in centres):
            mm.set_chi(j, True)
    mover = pyrosetta.rosetta.protocols.minimization_packing.MinMover(
        mm, scorefxn, 'lbfgs_armijo_nonmonotone', 0.001, True)
    mover.apply(work)
    return work


def prepare_complex(pose, scorefxn):
    """Remove input rotamer clashes before comparing bound and repacked free states.

    InterfaceAnalyzer packs separated partners globally. A locally repacked input
    against globally repacked partners otherwise charges pre-existing core strain
    to binding. Keep the prepared coordinates for downstream PROPKA and auditing.
    """
    work = pose.clone()
    task = pyrosetta.standard_packer_task(work)
    task.restrict_to_repacking()
    task.or_include_current(True)
    pyrosetta.rosetta.protocols.minimization_packing.PackRotamersMover(scorefxn, task).apply(work)
    mm = pyrosetta.rosetta.core.kinematics.MoveMap()
    mm.set_bb(False)
    mm.set_chi(True)
    pyrosetta.rosetta.protocols.minimization_packing.MinMover(
        mm, scorefxn, 'lbfgs_armijo_nonmonotone', 0.001, True).apply(work)
    return work


def save_json(path, data):
    tmp = path.with_suffix(path.suffix + '.tmp')
    tmp.write_text(json.dumps(data, indent=2, allow_nan=False) + '\n')
    tmp.replace(path)


def run_plan(plan, out_dir, max_ddg=3.0, seed=20261021):
    """Model each pair and retain diagnostic failures in an atomic manifest."""
    if not math.isfinite(max_ddg):
        raise ValueError('max_ddg must be finite')
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    stem = Path(plan['source_pdb']).stem
    manifest = out / (stem + '_variants.json')
    if manifest.exists():
        raise ValueError('manifest already exists; use a fresh output directory')
    result = dict(plan, status='running', seed=seed, max_ddg=max_ddg,
                  variants=[], pair_references=[], errors=[])
    save_json(manifest, result)
    if not plan['pairs']:
        result['status'] = 'no_eligible_pair'
        save_json(manifest, result)
        return result
    init(seed)
    result['pyrosetta_version'] = pyrosetta.version()
    scorefxn = pyrosetta.get_fa_scorefxn()
    result['scorefunction'] = scorefxn.get_name()
    pose = pyrosetta.pose_from_pdb(plan['source_pdb'])
    pose = prepare_complex(pose, scorefxn)
    prepared_name = stem + '_prepared.pdb'
    pose.dump_pdb(str(out / prepared_name))
    result['preparation'] = dict(
        protocol='whole-complex repack and side-chain minimisation; backbone fixed',
        total_energy=float(scorefxn(pose)), pdb=prepared_name,
        sha256=hashlib.sha256((out / prepared_name).read_bytes()).hexdigest(),
        ranking_geometry='original source PDB, before side-chain preparation')
    save_json(manifest, result)
    target, binder = plan['target_chain'], plan['binder_chain']
    for pair_n, pair in enumerate(plan['pairs'], 1):
        sites = pair['sites']
        positions = map_sites(pose, sites)
        tag = '__'.join('%s%s%d' % (s['chain'], s['wt'], s['resnum']) for s in sites)
        reference_name = '%s_%s_reference.pdb' % (stem, tag)
        try:
            pyrosetta.rosetta.numeric.random.rg().set_seed(seed)
            reference = mutate_and_relax(pose, positions, [s['wt'] for s in sites], scorefxn)
            base_dg, base_sasa = dG_separated(reference, target, binder, scorefxn)
            reference.dump_pdb(str(out / reference_name))
        except Exception as exc:
            result['errors'].append(dict(pair=pair_n, stage='reference', error=repr(exc)))
            save_json(manifest, result)
            continue
        result['pair_references'].append(dict(
            pair=pair_n, dG=base_dg, dSASA=base_sasa, pdb=reference_name,
            sha256=hashlib.sha256((out / reference_name).read_bytes()).hexdigest()))
        for aa, kind in (('H', 'switch'), ('Q', 'control')):
            name = '%s_%s_%s_%s.pdb' % (stem, tag, aa * len(sites), kind)
            row = dict(pair=pair_n, sites=sites, mut=aa * len(sites), kind=kind,
                       reference=reference_name, base_dG=base_dg, kept=False, pdb=None)
            try:
                # All three start from the same prepared pose and operation order.
                pyrosetta.rosetta.numeric.random.rg().set_seed(seed)
                variant = mutate_and_relax(pose, positions, [aa] * len(sites), scorefxn)
                dg, sasa = dG_separated(variant, target, binder, scorefxn)
                ddg = dg - base_dg
                keep = base_dg < 0 and dg < 0 and ddg <= max_ddg
                row.update(dG=dg, dSASA=sasa, ddG=ddg, kept=keep,
                           reference_binding_ok=base_dg < 0,
                           status='kept' if keep else 'discarded')
                if not keep:
                    row['rejection_reason'] = ('nonnegative reference dG' if base_dg >= 0 else
                                               'nonnegative variant dG' if dg >= 0 else 'ddG exceeds limit')
                if keep:
                    variant.dump_pdb(str(out / name))
                    row['pdb'] = name
                    row['sha256'] = hashlib.sha256((out / name).read_bytes()).hexdigest()
            except Exception as exc:
                row.update(status='failed', error=repr(exc))
                result['errors'].append(dict(pair=pair_n, stage=kind, error=repr(exc)))
            result['variants'].append(row)
            save_json(manifest, result)
        print('%s pair %d: %d/%d twins kept' % (
            stem, pair_n, sum(r['kept'] for r in result['variants'] if r['pair'] == pair_n), 2), flush=True)
    result['status'] = 'completed_with_errors' if result['errors'] else 'completed'
    save_json(manifest, result)
    return result


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('pdb')
    ap.add_argument('--out-dir', required=True)
    ap.add_argument('--target-chain', default='A', help='all target chains, e.g. AB')
    ap.add_argument('--binder-chain', default='B')
    ap.add_argument('--max-positions', type=int, default=8)
    ap.add_argument('--sites-per-variant', type=int, default=2)
    ap.add_argument('--max-pairs', type=int, default=1)
    ap.add_argument('--max-ddg', type=float, default=3.0)
    ap.add_argument('--seed', type=int, default=20261021)
    ap.add_argument('--plan-only', action='store_true')
    a = ap.parse_args()
    plan = make_plan(a.pdb, a.binder_chain, a.target_chain, a.max_positions,
                     a.sites_per_variant, a.max_pairs)
    out = Path(a.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    if a.plan_only:
        save_json(out / (Path(a.pdb).stem + '_plan.json'), dict(plan, status='ranking_only'))
    else:
        result = run_plan(plan, out, a.max_ddg, a.seed)
        if result['status'] != 'completed':
            raise SystemExit(result['status'])


if __name__ == '__main__':
    main()
