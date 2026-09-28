"""Offline pilot: conservative ring identity strata and capacity-bounded torsion blocks."""
import argparse
from collections import Counter,defaultdict
import hashlib
import json
from pathlib import Path
import sqlite3

import numpy as np

VERSION='macrocycle-torsion-blocks-v1'


def digest(value):
    return hashlib.sha256(json.dumps(value,sort_keys=True,separators=(',',':')).encode()).hexdigest()


def angle(points):
    p=np.asarray(points,float);axis=p[2]-p[1]
    norm=np.linalg.norm(axis)
    if norm<1e-8:raise ValueError('degenerate_ring_geometry')
    axis/=norm;left=p[0]-p[1];right=p[3]-p[2]
    left-=left.dot(axis)*axis;right-=right.dot(axis)*axis
    if min(np.linalg.norm(left),np.linalg.norm(right))<1e-8:raise ValueError('degenerate_ring_geometry')
    return float(np.arctan2(np.dot(np.cross(axis,left),right),np.dot(left,right)))


def amide_state(degrees):
    absolute=abs(degrees)
    return 'cis' if absolute<=30 else ('trans' if absolute>=150 else 'boundary')


def describe(chem,xyz,ring=None,minimum_ring=12):
    from rdkit import Chem
    xyz=np.asarray(xyz,float);n=len(chem.atomic_numbers)
    if xyz.shape!=(n,3) or not np.isfinite(xyz).all():raise ValueError('invalid_coordinates_or_atom_count')
    edges={};neighbors=defaultdict(list);mol=Chem.RWMol()
    for z in chem.atomic_numbers:mol.AddAtom(Chem.Atom(int(z)))
    for bond in chem.bonds:
        a,b,order=int(bond['begin']),int(bond['end']),int(bond['order'])
        if not 0<=a<n or not 0<=b<n or a==b or tuple(sorted((a,b))) in edges:raise ValueError('invalid_bond_graph')
        edges[tuple(sorted((a,b)))]=order;neighbors[a].append(b);neighbors[b].append(a)
        # Ring perception needs connectivity only, not reconstruction/sanitization.
        mol.AddBond(a,b,Chem.BondType.SINGLE)
    rings=[list(map(int,r)) for r in Chem.GetSymmSSSR(mol)]
    mapping='supplied' if ring is not None else 'unique_isolated_macrocycle'
    if ring is None:
        large=[r for r in rings if len(r)>=minimum_ring]
        if len(large)!=1:raise ValueError('no_unique_macrocycle')
        ring=large[0]
        if any(r!=ring and len(set(r)&set(ring))>1 for r in rings):raise ValueError('fused_or_bridged_requires_mapping')
    if not isinstance(ring,list) or any(type(i) is not int for i in ring) or len(set(ring))!=len(ring) or len(ring)<minimum_ring:
        raise ValueError('invalid_ordered_ring_map')
    if any(i<0 or i>=n for i in ring):raise ValueError('ring_atom_out_of_range')
    size=len(ring)
    def order(a,b):return edges.get(tuple(sorted((a,b))),0)
    if any(not order(ring[i],ring[(i+1)%size]) for i in range(size)):raise ValueError('ring_map_not_closed_bond_path')
    def amide(a,b):
        if {int(chem.atomic_numbers[a]),int(chem.atomic_numbers[b])}!={6,7} or order(a,b)!=1:return False
        carbon=a if chem.atomic_numbers[a]==6 else b
        return any(chem.atomic_numbers[c]==8 and order(carbon,c)==2 for c in neighbors[carbon])
    variants=[]
    for direction in (ring,list(reversed(ring))):
        direction_angles=[angle(xyz[[direction[(i-1)%size],direction[i],direction[(i+1)%size],direction[(i+2)%size]]]) for i in range(size)]
        for offset in range(size):
            ids=direction[offset:]+direction[:offset]
            tokens=tuple((int(chem.atomic_numbers[a]),int(chem.formal_charges[a]),order(a,ids[(i+1)%size])) for i,a in enumerate(ids))
            angles=direction_angles[offset:]+direction_angles[:offset]
            states=tuple(amide_state(np.degrees(angles[i])) if amide(a,ids[(i+1)%size]) else '-' for i,a in enumerate(ids))
            variants.append((tokens,states,tuple(np.round(angles,6)),ids,angles))
    tokens,states,_,ids,angles=min(variants,key=lambda v:v[:3])
    values=np.asarray(angles)
    units=[];ring_set=set(ring)
    for nitrogen in ids:
        if chem.atomic_numbers[nitrogen]!=7:continue
        for ca in neighbors[nitrogen]:
            if ca not in ring_set or chem.atomic_numbers[ca]!=6:continue
            for carbon in neighbors[ca]:
                if carbon==nitrogen or carbon not in ring_set or chem.atomic_numbers[carbon]!=6:continue
                following=[b for b in neighbors[carbon] if b in ring_set and amide(carbon,b)]
                if len(following)!=1:continue
                branches=set();pending=[b for a in (nitrogen,ca,carbon) for b in neighbors[a] if b not in ring_set]
                while pending:
                    atom=pending.pop()
                    if atom in branches:continue
                    branches.add(atom);pending.extend(b for b in neighbors[atom] if b not in ring_set and b not in branches)
                units.append(dict(nitrogen=nitrogen,alpha_carbon=ca,carbonyl_carbon=carbon,next_nitrogen=following[0],
                    attached_heavy_atoms=sorted(branches)))
    return dict(ring_atoms=ids,ring_size=size,ring_mapping=mapping,ring_tokens=tokens,
        peptide_units=units,peptide_unit_scope='Graph-derived N-C-C(=O) units; no name-token or synthesis building-block assignment',
        amide_states=states,amide_count=sum(s!='-' for s in states),
        eligible_backbone_single_bonds=sum(order(a,ids[(i+1)%size])==1 and states[i]=='-' for i,a in enumerate(ids)),
        torsions_radians=angles,descriptor=np.concatenate((np.sin(values),np.cos(values))).tolist(),
        hard_group=digest([size,tokens,states]),
        limitations=['Backbone-only baseline; side chains/building blocks and independent closed-ring degrees of freedom are not modeled',
                    'Symmetric cyclic correspondence uses a lexicographic torsion tie-break; boundary stability needs benchmarking'])


def partition(records,capacity):
    if type(capacity) is not int or capacity<1:raise ValueError('capacity must be positive')
    if len({r['conformer_id'] for r in records})!=len(records):raise ValueError('Duplicate conformer identity')
    groups=defaultdict(list)
    for row in records:groups[row['hard_group']].append(row)
    memberships=[];blocks=[]
    for group,members in sorted(groups.items()):
        def split(rows,path):
            x=np.asarray([r['descriptor'] for r in rows],float)
            if x.ndim!=2 or not np.isfinite(x).all():raise ValueError('Invalid descriptor matrix')
            if len(rows)>capacity:
                axis=int(np.argmax(x.var(axis=0)))
                ordered=sorted(rows,key=lambda r:(r['descriptor'][axis],r['conformer_id']))
                midpoint=len(rows)//2
                split(ordered[:midpoint],path+'0');split(ordered[midpoint:],path+'1');return
            block='BLK-'+digest([VERSION,capacity,group,path,[r['conformer_id'] for r in sorted(rows,key=lambda r:r['conformer_id'])]])[:20]
            center=x.mean(axis=0)
            blocks.append(dict(block_id=block,hard_group=group,split_path=path,conformers=len(rows),
                molecules=len({r['molecule_id'] for r in rows}),centroid=center.tolist(),
                maximum_distance=float(np.linalg.norm(x-center,axis=1).max())))
            memberships.extend(dict(conformer_id=r['conformer_id'],molecule_id=r['molecule_id'],global_id=r['global_id'],block_id=block) for r in rows)
        split(members,'')
    return sorted(memberships,key=lambda r:r['global_id']),blocks


def run(batch,output,limit=10000,capacity=10000,seed=20260925,mapping=None):
    from .gaussian_batch import ArtifactCatalogReader,_sha256
    from .chemical_companion import ChemicalCompanionReader
    if not 1<=limit<=100000:raise ValueError('Pilot limit must be 1..100000 conformers')
    if capacity<1:raise ValueError('capacity must be positive')
    batch=Path(batch).resolve();output=Path(output).resolve()
    if output.exists():raise ValueError('Use a new output directory; frozen memberships are not overwritten')
    catalog=batch/'artifacts/catalog.json';chemical=batch/'chemical/catalog.json'
    shape=ArtifactCatalogReader(catalog);chem=ChemicalCompanionReader(chemical)
    chemical_metadata=json.loads(chemical.read_text())
    if chemical_metadata.get('library_id')!=shape.catalog['library_id']:raise ValueError('Cross-library catalogs')
    if chemical_metadata.get('artifact_v1_catalog_sha256')!=_sha256(catalog):raise ValueError('Chemical companion refers to a different artifact catalog')
    maps={}
    if mapping:
        for line in Path(mapping).read_text().splitlines():
            row=json.loads(line)
            if row['conformer_id'] in maps:raise ValueError('Duplicate mapping conformer ID')
            maps[row['conformer_id']]=row
    total=int(shape.catalog['conformers']);count=min(limit,total)
    # Floyd sampling uses O(sample size) memory even for a billion-row catalog.
    rng=np.random.default_rng(seed);chosen=set()
    for j in range(total-count,total):
        k=int(rng.integers(0,j+1));chosen.add(j if k in chosen else k)
    selected=sorted(chosen);records=[];unassigned=[]
    output.mkdir(parents=True)
    with sqlite3.connect((batch/'registry.sqlite3').as_uri()+'?mode=ro',uri=True) as db:
        db.row_factory=sqlite3.Row
        with (output/'descriptors.jsonl').open('w',encoding='utf-8') as stream:
            for gid in selected:
                s=shape.get(gid);c=chem.get(gid)
                if (s.conformer_id,s.molecule_id)!=(c.conformer_id,c.molecule_id):raise ValueError('Artifact/chemical identity mismatch')
                identity=db.execute('SELECT c.*,m.source_name,m.library_id FROM conformer c JOIN molecule m ON m.id=c.molecule_id WHERE c.id=?',(s.conformer_id,)).fetchone()
                if identity is None or identity['molecule_id']!=s.molecule_id or identity['library_id']!=shape.catalog['library_id']:raise ValueError('Registry identity mismatch')
                provenance={k:identity[k] for k in ('source_name','source_record_name','source_path','source_record_index','content_sha256','topology_sha256','conformer_index')}
                base=dict(global_id=gid,molecule_id=s.molecule_id,conformer_id=s.conformer_id,**provenance)
                try:
                    supplied=maps.get(s.conformer_id)
                    if supplied and supplied.get('content_sha256')!=identity['content_sha256']:raise ValueError('ring_mapping_source_hash_mismatch')
                    descriptor=describe(c,s.shape_points,supplied['ring_atoms'] if supplied else None)
                    row=dict(base,**descriptor);stream.write(json.dumps(row)+'\n')
                    records.append({k:row[k] for k in ('global_id','molecule_id','conformer_id','hard_group','descriptor')})
                except ValueError as exc:unassigned.append(dict(base,reason=str(exc)))
    memberships,blocks=partition(records,capacity)
    for name,value in [('memberships',memberships),('blocks',blocks),('unassigned',unassigned)]:
        with (output/(name+'.jsonl')).open('w',encoding='utf-8') as stream:
            for row in value:stream.write(json.dumps(row)+'\n')
    report=dict(status='complete',version=VERSION,scope='Offline conformer pilot; no screening or block rejection',
        library_conformers=total,sampled_conformers=count,assigned_conformers=len(records),unassigned_conformers=len(unassigned),
        blocks=len(blocks),capacity=capacity,seed=seed,global_ids=selected,
        unassigned_reasons=dict(Counter(r['reason'] for r in unassigned)),
        input_hashes={str(p):_sha256(p) for p in [catalog,chemical]+([Path(mapping)] if mapping else [])},
        code_sha256=_sha256(Path(__file__)),output_hashes={p.name:_sha256(p) for p in output.glob('*.jsonl')},
        audit=dict(building_block_mapping='not_available_in_standard_catalogs',atom_index_basis='zero-based heavy-atom artifact order',
                   chemistry_hydrogens='not_reconstructed',integrity='Catalog hashes and sampled ID joins, not full artifact/source byte verification'),
        limitations=['Not production full-library clustering or validated screening acceleration',
                    'Incremental assignment, typed side-chain geometry and energy-basin partitioning are not implemented',
                    'Target capacity is an upper bound; small hard groups and residual blocks remain small'])
    (output/'report.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
    return report


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--batch',type=Path,required=True);parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--limit',type=int,default=10000);parser.add_argument('--capacity',type=int,default=10000)
    parser.add_argument('--seed',type=int,default=20260925);parser.add_argument('--ring-mapping',type=Path)
    args=parser.parse_args();print(json.dumps(run(args.batch,args.output,args.limit,args.capacity,args.seed,args.ring_mapping),indent=2))


if __name__=='__main__':main()
