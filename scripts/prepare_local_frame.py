"""Rigidly translate the complete receptor so an unresolved binder starts near its patch."""
import argparse
import hashlib
import json
import math
from pathlib import Path


def records(path):
    result = []
    for line in Path(path).read_text().splitlines(keepends=True):
        if line.startswith(("ATOM  ", "HETATM")):
            result.append((line, tuple(float(line[start:start+8]) for start in (30,38,46))))
    return result


def mean(points):
    return tuple(sum(p[i] for p in points)/len(points) for i in range(3))


def prepare(source, context, destination, offset=15.0):
    source, context, destination = map(Path, (source, context, destination))
    raw = records(source)
    ca = [xyz for line,xyz in raw if line.startswith('ATOM  ') and line[12:16].strip() == 'CA']
    patch = [xyz for line,xyz in raw if line[21]=='A' and int(line[22:26]) in {205,206,208,209,210} and line[12:16].strip()=='CA']
    if len(ca)!=639 or len(patch)!=5:
        raise ValueError('Requires the complete639-residue A target and all five declared hotspots')
    centroid, hotspot = mean(ca), mean(patch)
    outward = tuple(a-b for a,b in zip(hotspot,centroid))
    norm = math.sqrt(sum(x*x for x in outward))
    shift = tuple(a+offset*b/norm for a,b in zip(hotspot,outward))
    context_atoms = [xyz for line,xyz in records(context) if line[76:78].strip() not in {'H','D'}]
    clearance = min(math.dist(shift,xyz) for xyz in context_atoms)
    if clearance < 6:
        raise ValueError(f'Binder origin would be inside the restored context: clearance{clearance:.3f}A')
    output = []
    rounding = []
    for line in source.read_text().splitlines(keepends=True):
        if line.startswith(('ATOM  ','HETATM')):
            old = tuple(float(line[start:start+8]) for start in (30,38,46))
            new = tuple(a-b for a,b in zip(old,shift))
            text = ''.join(f'{v:8.3f}' for v in new)
            rounding.append(math.dist(new,tuple(float(text[start:start+8]) for start in (0,8,16))))
            line = line[:30]+text+line[54:]
        output.append(line)
    destination.parent.mkdir(parents=True,exist_ok=True)
    destination.write_text(''.join(output))
    written = records(destination)
    if len(written)!=len(raw) or any(a[0][:30]!=b[0][:30] or a[0][54:]!=b[0][54:] for a,b in zip(raw,written)):
        raise ValueError('Atom identity or non-coordinate PDB columns changed')
    after_ca = [xyz for line,xyz in written if line.startswith('ATOM  ') and line[12:16].strip()=='CA']
    max_distance_error = max(abs(math.dist(ca[i],ca[j])-math.dist(after_ca[i],after_ca[j]))
                             for i in range(len(ca)) for j in range(i))
    if max(rounding)>0.001 or max_distance_error>0.002:
        raise ValueError('Rigid translation changed internal geometry beyond PDB precision')
    return dict(source=str(source),context=str(context),destination=str(destination),
                source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
                destination_sha256=hashlib.sha256(destination.read_bytes()).hexdigest(),
                source_hotspot_centroid=hotspot,source_receptor_centroid=centroid,
                old_origin_to_hotspot_A=math.dist((0,0,0),hotspot),
                translation=shift,new_origin_to_hotspot_A=offset,
                new_origin_full_context_atom_clearance_A=clearance,
                atoms=len(raw),residues=len(ca),non_coordinate_columns_identical=True,
                max_coordinate_rounding_A=max(rounding),max_CA_pair_distance_change_A=max_distance_error,
                limitation='Placement is an initialization hypothesis; not evidence of affinity or a corrected v3 failure cause')


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--source',default='results/gpu_campaign_v2/human_11qc_full_monomer.pdb')
    parser.add_argument('--context',default='results/gpu_campaign/human_11qc_full_dimer.pdb')
    parser.add_argument('--out',default='results/modal_campaign_v5/human_11qc_local_frame.pdb')
    args=parser.parse_args()
    result=prepare(args.source,args.context,args.out)
    Path(args.out).with_suffix('.audit.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
