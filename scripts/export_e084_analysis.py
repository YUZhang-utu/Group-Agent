"""Export an immutable pocket-state result and its complete analysis inputs."""
import argparse
from collections import Counter
import csv
import hashlib
import html
import itertools
import json
from pathlib import Path
import shutil
import zipfile

import numpy as np
from scipy.spatial import cKDTree
from aidd_agent.consensus_admission import read_structure, ca_map
from aidd_agent.screening_selection import check_hashes


def write_csv(path, rows, fields=None):
    rows=list(rows)
    fields=fields or list(rows[0])
    with Path(path).open('w',newline='',encoding='utf-8-sig') as stream:
        writer=csv.DictWriter(stream,fieldnames=fields);writer.writeheader();writer.writerows(rows)


def digest(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for part in iter(lambda:stream.read(1024*1024),b''):h.update(part)
    return h.hexdigest()


def export(source, output):
    source=Path(source).resolve();output=Path(output).resolve();r=json.loads(source.read_text())
    check_hashes(r['sources'])
    output.mkdir(parents=True,exist_ok=False)
    data=np.load(r['artifact']);grid=data['grid'];masks=data['masks'];features=data['chemical_fields']
    distance=data['distance'];overlap=data['overlap'];chem=data['chemical_distance'];local=data['local_difference']
    rows=r['structures'];ids=[row['id'] for row in rows];lookup={x:i for i,x in enumerate(ids)}
    assert ids==r['pairwise']['ids'] and np.allclose(distance,r['pairwise']['distance'])
    groups=r['clusters'];labels={g['id']:f'C{i+1}' for i,g in enumerate(groups)}
    owner={sid:g for g in groups for sid in g['members']}
    assert len(owner)==len(rows) and set(owner)==set(ids)
    pdb_counts=Counter(row['pdb_id'] for row in rows);weights=np.array([1/pdb_counts[row['pdb_id']] for row in rows])
    volume=r['policy']['spacing']**3;summary=[];members=[]
    for group in groups:
        indices=[lookup[sid] for sid in group['members']];rep=lookup[group['representative']]
        values=[rows[i]['cavity_volume'] for i in indices]
        summary.append(dict(cluster=labels[group['id']],cluster_id=group['id'],representative=group['representative'],
            chains=len(indices),distinct_PDBs=group['distinct_pdb_support'],
            representative_volume_A3=rows[rep]['cavity_volume'],volume_min_A3=min(values),volume_max_A3=max(values),
            max_member_distance_to_representative=group['representative_max_distance'],
            min_member_IoU_to_representative=float(overlap[rep,indices].min()),
            maximum_within_cluster_chemical_difference=group['maximum_chemical_difference'],
            representative_resolution_A=rows[rep]['resolution'],
            small_support=group['small_support'],review_flags='; '.join(group['review_flags'])))
    for i,row in enumerate(rows):
        group=owner[row['id']];rep=lookup[group['representative']];indices=[lookup[x] for x in group['members']]
        other=[]
        for g in groups:
            if g['id']==group['id']:continue
            indices_other=[lookup[x] for x in g['members']]
            other.append((float(np.average(distance[i,indices_other],weights=weights[indices_other])),labels[g['id']]))
        nearest=min(other) if other else (None,None)
        entry=dict(cluster=labels[group['id']],cluster_id=group['id'],structure=row['id'],PDB=row['pdb_id'],
            target_chain=row['target_chain'],is_representative=i==rep,
            representative=group['representative'],volume_A3=row['cavity_volume'],
            distance_to_representative=float(distance[i,rep]),IoU_to_representative=float(overlap[i,rep]),
            chemical_difference_to_representative=float(chem[i,rep]),
            local_difference_to_representative=float(local[i,rep]),
            mean_distance_to_own_cluster=float(np.average(distance[i,indices],weights=weights[indices])),
            closest_other_cluster=nearest[1],mean_distance_to_closest_other_cluster=nearest[0],
            core_CA_RMSD_A=row['core_rmsd'],alignment_residues=';'.join(map(str,row['alignment_residues'])),
            resolution_A=row['resolution'],PDB_chain_weight=float(weights[i]),
            ligand_queries=';'.join(q['query_id'] for q in row['queries']),
            PDB_link='https://www.rcsb.org/structure/'+row['pdb_id'])
        entry.update({channel+'_field_volume_A3':float(features[i,k].sum()*volume) for k,channel in enumerate(r['channels'])})
        members.append(entry)
    pairs=[]
    for i,j in itertools.combinations(range(len(rows)),2):
        pairs.append(dict(first=ids[i],second=ids[j],first_cluster=labels[owner[ids[i]]['id']],
            second_cluster=labels[owner[ids[j]]['id']],same_cluster=owner[ids[i]]['id']==owner[ids[j]]['id'],
            distance=float(distance[i,j]),cavity_IoU=float(overlap[i,j]),chemical_difference=float(chem[i,j]),
            local_difference=float(local[i,j]),intersection_volume_A3=float((masks[i]&masks[j]).sum()*volume),
            first_only_volume_A3=float((masks[i]&~masks[j]).sum()*volume),
            second_only_volume_A3=float((masks[j]&~masks[i]).sum()*volume)))
    group_pairs=[]
    for a,b in itertools.combinations(groups,2):
        ia=[lookup[x] for x in a['members']];ib=[lookup[x] for x in b['members']]
        value=float(np.average(distance[np.ix_(ia,ib)],weights=np.outer(weights[ia],weights[ib])))
        group_pairs.append(dict(first_cluster=labels[a['id']],second_cluster=labels[b['id']],
            weighted_average_linkage_distance=value,cutoff=r['policy']['cluster_distance'],
            separation_margin=value-r['policy']['cluster_distance'],
            closest_member_distance=float(distance[np.ix_(ia,ib)].min()),
            representative_distance=float(distance[lookup[a['representative']],lookup[b['representative']]])))
    assert all(p['separation_margin']>0 for p in group_pairs)
    # Reconstruct the full weighted average-linkage tree, recording the cutoff.
    work=distance.copy();np.fill_diagonal(work,np.inf);active=list(range(len(rows)))
    group_members={i:[i] for i in active};mass={i:weights[i] for i in active};node={i:ids[i] for i in active};trace=[]
    at_cutoff=None
    while len(active)>1:
        sub=work[np.ix_(active,active)];a,b=np.unravel_index(np.argmin(sub),sub.shape)
        i,j=active[a],active[b];d=float(sub[a,b])
        if d>r['policy']['cluster_distance'] and at_cutoff is None:
            at_cutoff={frozenset(ids[k] for k in group_members[v]) for v in active}
        name=f'M{len(trace)+1:03d}'
        trace.append(dict(merge=name,left=node[i],right=node[j],distance=d,
            accepted_at_default_cutoff=d<=r['policy']['cluster_distance'],
            left_members=';'.join(ids[k] for k in group_members[i]),right_members=';'.join(ids[k] for k in group_members[j]),
            merged_weight=float(mass[i]+mass[j])))
        total=mass[i]+mass[j]
        for k in active:
            if k not in {i,j}:work[i,k]=work[k,i]=(mass[i]*work[i,k]+mass[j]*work[j,k])/total
        mass[i]=total;group_members[i]+=group_members.pop(j);node[i]=name;active.remove(j)
    at_cutoff=at_cutoff or {frozenset(ids[k] for k in group_members[v]) for v in active}
    assert at_cutoff=={frozenset(g['members']) for g in groups}
    held=[]
    for item in r['held_for_review']:
        for reason in item['reasons']:
            held.append(dict(candidate=item['id'],candidate_type='chain' if ':' in item['id'] else 'entry',
                reason_code=reason.split(':')[0],detail=reason,additional_detail=item.get('detail','')))
    write_csv(output/'01_cluster_summary.csv',summary)
    write_csv(output/'02_all_68_members.csv',members)
    write_csv(output/'03_all_2278_pairs.csv',pairs)
    write_csv(output/'04_between_cluster_distances.csv',group_pairs)
    write_csv(output/'05_weighted_merge_history.csv',trace)
    write_csv(output/'06_held_reason_records.csv',held)
    write_csv(output/'07_threshold_sensitivity.csv',r['sensitivity'])
    for name,array in [('distance',distance),('cavity_IoU',overlap),('chemical_difference',chem),('local_difference',local)]:
        write_csv(output/('matrix_'+name+'.csv'),[dict(structure=sid,**{ids[j]:float(array[i,j]) for j in range(len(ids))}) for i,sid in enumerate(ids)])
    # Atom-name matched displacement is descriptive, not causal attribution.
    aligned={}
    for group in groups:
        row=rows[lookup[group['representative']]]
        s=read_structure(row['structure_path'],r['target']);matrix=np.array(row['transform'])
        atoms=[dict(a,xyz=np.asarray(a['xyz'])@matrix[:3,:3].T+matrix[:3,3]) for a in s['chains'][row['target_chain']]
               if not a['label_alt_id']]
        aligned[group['id']]=atoms
    residue_rows=[]
    for ga,gb in itertools.combinations(groups,2):
        aa=aligned[ga['id']];bb=aligned[gb['id']]
        da={(a['canonical_residue'],a['auth_atom_id']):a for a in aa}
        db={(a['canonical_residue'],a['auth_atom_id']):a for a in bb}
        for residue in r['pocket_residues']:
            keys=sorted(k for k in set(da)&set(db) if k[0]==residue)
            if not keys:continue
            shifts=[float(np.linalg.norm(da[k]['xyz']-db[k]['xyz'])) for k in keys]
            largest=int(np.argmax(shifts));ca=(residue,'CA')
            residue_rows.append(dict(first_cluster=labels[ga['id']],second_cluster=labels[gb['id']],
                first_representative=ga['representative'],second_representative=gb['representative'],
                canonical_residue=residue,residue_name=da[keys[0]]['auth_comp_id'],matched_heavy_atoms=len(keys),
                CA_displacement_A=float(np.linalg.norm(da[ca]['xyz']-db[ca]['xyz'])) if ca in da and ca in db else '',
                matched_heavy_atom_RMSD_A=float(np.sqrt(np.mean(np.square(shifts)))),
                largest_atom_displacement_A=shifts[largest],largest_displacement_atom=keys[largest][1]))
    write_csv(output/'08_representative_residue_displacements.csv',residue_rows)
    reference_id=r['reference']['id'];reference_row=rows[lookup[reference_id]]
    reference_structure=read_structure(reference_row['structure_path'],r['target'])
    fixed=ca_map(reference_structure['chains'][reference_row['target_chain']])
    audits=[];ca_details=[]
    for row in rows:
        structure=read_structure(row['structure_path'],r['target']);moving=ca_map(structure['chains'][row['target_chain']])
        matrix=np.array(row['transform']);common=sorted(set(fixed)&set(moving)&set(r['pocket_residues']))
        residuals=[float(np.linalg.norm(np.asarray(moving[k])@matrix[:3,:3].T+matrix[:3,3]-fixed[k])) for k in common]
        for residue,residual in zip(common,residuals):
            ca_details.append(dict(structure=row['id'],cluster=labels[owner[row['id']]['id']],
                canonical_residue=residue,CA_displacement_from_reference_A=residual,
                retained_in_alignment_core=residue in row['alignment_residues']))
        displaced=sum(v>5 for v in residuals);fraction=displaced/max(len(common),1)
        critical=max(residuals,default=0)>15 or fraction>.25
        audits.append(dict(structure=row['id'],cluster=labels[owner[row['id']]['id']],
            fitted_core_RMSD_A=row['core_rmsd'],all_observed_pocket_CA_RMSD_A=float(np.sqrt(np.mean(np.square(residuals)))),
            max_pocket_CA_displacement_A=max(residuals,default=0),pocket_CA_count=len(common),
            pocket_CA_over_5A=displaced,fraction_over_5A=fraction,
            export_review_flag='critical_geometry_review' if critical else 'localized_difference_review' if displaced else 'no_large_CA_displacement',
            original_membership_preserved=True))
    write_csv(output/'10_independent_pocket_alignment_audit.csv',audits)
    write_csv(output/'11_all_pocket_CA_residuals.csv',ca_details)
    # Package all inspected CIFs, including held entries, not just selected medoids.
    inputs=output/'inputs';(inputs/'structures').mkdir(parents=True)
    source_map=[]
    for path,expected in r['sources'].items():
        p=Path(path)
        if p.suffix.lower()=='.cif':
            target=inputs/'structures'/p.name;shutil.copy2(p,target)
            source_map.append(dict(original_path=path,package_path=target.relative_to(output).as_posix(),sha256=expected))
    from Bio.PDB.MMCIF2Dict import MMCIF2Dict
    metadata=[]
    for item in source_map:
        path=output/item['package_path']; d=MMCIF2Dict(str(path))
        def values(key):
            value=d.get(key,[])
            return '; '.join(value if isinstance(value,list) else [str(value)])
        metadata.append(dict(PDB=path.stem,experimental_method=values('_exptl.method'),
            deposited_crystal_resolution_A=values('_refine.ls_d_res_high'),
            deposited_EM_resolution_A=values('_em_3d_reconstruction.resolution'),
            deposited_R_free=values('_refine.ls_R_factor_R_free'),title=values('_struct.title'),
            citation_DOIs=values('_citation.pdbx_database_id_DOI'),package_CIF=item['package_path']))
    write_csv(output/'09_all_133_PDB_metadata.csv',metadata)
    for name in ('protein.json','ligand-instances.json','report.json'):
        p=Path(r['diversity_report']).parent/name;shutil.copy2(p,inputs/('diversity-'+name))
    shutil.copy2(source,output/'original_pocket_report.json')
    shutil.copy2(r['artifact'],output/'pocket-grids.npz')
    shutil.copy2(source.with_name('report.html'),output/'interactive_pocket_comparison.html')
    (output/'input_path_map.json').write_text(json.dumps(source_map,indent=2))
    (output/'settings.json').write_text(json.dumps(dict(policy=r['policy'],channels=r['channels'],
        reference=r['reference'],pocket_residues=r['pocket_residues'],quality_atom_scope=r['quality_atom_scope'],
        structures=[{k:row[k] for k in ('id','pdb_id','target_chain','transform','alignment_residues')} for row in rows]),indent=2))
    description='''# MDM2 pocket-state analysis package

This package exports an existing experimental-structure analysis without changing
its assignments. Six clusters contain 55, 8, 2, 1, 1 and 1 receptor chains, from
60 distinct PDB entries. C1..C6 are readable aliases, not stability rankings.
180 held candidates mix chain-level and entry-level records; they are not 180 PDBs.

IMPORTANT INDEPENDENT REVIEW FINDING: C5 / 7BJ0:A is not currently an accepted
ordinary pocket state. The original trimmed-core alignment RMSD is about 1.82 A,
but the full 24 observed pocket CA atoms have RMSD about 23.77 A, maximum 45.80 A,
and 15 of 24 exceed 5 A. The trimming criterion did not protect against severe
displacement of the rest of the pocket. This could involve alignment, assembly
or a genuinely distinct structure; the cause has not been established. Retain it
for inspection and do not use it as a validated consensus/docking receptor.
The original assignments are preserved explicitly, not silently repaired.

## Read first

- 01_cluster_summary.csv: cluster sizes, representatives, volumes and review flags.
- 02_all_68_members.csv: exact membership, representative distances and input IDs.
- 03_all_2278_pairs.csv: every unique pair, overlap and free-space differences.
- 04_between_cluster_distances.csv: why final groups remain separate at the cutoff.
- 05_weighted_merge_history.csv: every merge, including above-cutoff hypothetical
  merges. Only rows marked accepted_at_default_cutoff contribute to the result.
- 06_held_reason_records.csv: one row per reason; several rows can refer to one
  candidate. Do not sum reason rows as a number of excluded structures.
- 07_threshold_sensitivity.csv: 12 / 6 / 3 clusters at distances 0.30 / 0.40 / 0.50.
- 08_representative_residue_displacements.csv: coordinate differences between all
  representative pairs, using canonical residue correspondence and atom names.
- 09_all_133_PDB_metadata.csv: experimental method, deposited resolution, R-free,
  title and citation DOI strings extracted directly from the original CIF files.
  Resolution in the member/cluster tables comes from the earlier eligible-ligand
  record; a blank there does not imply missing deposited resolution. Use table 09
  for the original metadata. Multiple deposited values remain semicolon-separated.
- 10_independent_pocket_alignment_audit.csv: all-pocket CA residuals, independent
  of the subset retained for fitting. Critical flags are export-time triage
  (>15 A maximum or >25% of pocket CA atoms displaced over 5 A), not a calibrated
  scientific acceptance rule or a change to the original clustering.
- 11_all_pocket_CA_residuals.csv: every observed pocket CA displacement and whether
  the residue was retained in the original trimmed alignment core.
- matrix_*.csv: full symmetric matrices in the same 68-structure order.
- interactive_pocket_comparison.html: offline XY/XZ free-space comparison.
- pocket-grids.npz: shared grid, Boolean cavity masks, six chemical feature fields,
  distance/overlap/chemical/local matrices. Array order is in original_pocket_report.
- inputs/structures: all 133 inspected original experimental CIF files, including
  held structures. input_path_map.json maps old absolute paths to package paths.
- settings.json: parameters, reference, local quality scope and alignment matrices.

## Exact distance and clustering definition

IoU = intersection / union of free-space voxels in one shared reference region.
ShapeDistance = max(1 - IoU, 0.5 * LocalDifference).
D = 0.85 * ShapeDistance + 0.15 * ChemicalDifference.
LocalDifference is the largest local IoU loss in a 4-A sphere centered on a
reference pocket CA, only if changed space in that sphere is at least 30 A^3.
ChemicalDifference averages channel IoU loss over active channels in the common
free space. Field volumes are voxel support, not atom counts or probabilities.
The reference-ligand region is expanded by 4 A, sampled at 1 A, and protein VDW
spheres plus a 0.5-A probe are excluded. This is not the entire protein cavity.

Weighted average linkage merges the closest clusters while their distance is at
most 0.40. Every PDB has total weight 1 distributed across its admitted chains.
The final cluster-to-cluster weighted distances exceed 0.40; individual pairs
across clusters can still be closer. Greedy merge history, not a nearest-medoid
classifier, determines membership. A member can exceed 0.40 from its medoid.
Representative selection minimizes weighted within-cluster distance, allowing a
0.02 cost tolerance for ligand-reference availability, resolution and stable ID.
PDB weights control duplicate chains, not publication or ligand-series bias.

## Interpretation boundaries

These are geometric/chemical states, not ligand scaffold families or equilibrium
populations. A large cluster can include meaningful local differences. Review
the representative-coverage and chemical-difference flags, especially for C1.
Chemical fields are residue templates and approximate heavy-atom directions,
not electrostatic energies or complete protonation/hydrogen-bond assignments.
Residue displacement is not an energetic contribution or proof of causality;
symmetric equivalent atom names can inflate atom-matched RMSD. Grid free-space
volume does not prove that a specific ligand fits, binds or can enter kinetically.
Small clusters are retained. Missing atoms, alternate occupancies, target mapping
and unsupported assembly contexts are separate reasons for holding candidates.
Default settings were chosen during development on this MDM2 set. No independent
docking recall, PLANTS run, N-E rescoring or user adoption is claimed here.

## Suggested review priorities

C1 contains most chains and has both representative-coverage and chemical-difference
flags. C2/C3 also have chemical-difference flags. C6 is a borderline separate
group: C1-C6 average linkage is about 0.415 versus cutoff 0.40, although the two
representatives are closer than the cutoff. This illustrates why centroid/medoid
distance alone does not explain average-linkage membership. Inspect local space
and contacts rather than interpreting six clusters as six established biological
families. Candidate chains/entries held for quality are separate from rare clusters.

## Verification

All recorded input hashes were checked. Exported pair distances match the saved
matrices. The full weighted merge trace reproduces the exact six original member
sets. All final cluster distances exceed the selected cutoff. SHA256SUMS.json
records every packaged file. Original analysis source code hashes are in the
unaltered original_pocket_report.json; implementation commit was 33800c5.
'''
    (output/'README.md').write_text(description,encoding='utf-8')
    table='<table><tr>'+''.join('<th>'+html.escape(k)+'</th>' for k in summary[0])+'</tr>'
    for row in summary:table+='<tr>'+''.join('<td>'+html.escape(str(v))+'</td>' for v in row.values())+'</tr>'
    table+='</table>'
    (output/'START_HERE.html').write_text('<!doctype html><meta charset="utf-8"><title>MDM2 analysis</title><style>body{font:15px system-ui;margin:28px}table{border-collapse:collapse}td,th{padding:8px;border:1px solid #ccc}th{background:#e7eef5}pre{white-space:pre-wrap;max-width:1100px}</style><h1>MDM2: six experimental pocket clusters</h1><p><a href="interactive_pocket_comparison.html">Interactive pocket comparison</a></p>'+table+'<pre>'+html.escape(description)+'</pre>',encoding='utf-8')
    files=sorted(p for p in output.rglob('*') if p.is_file())
    (output/'SHA256SUMS.json').write_text(json.dumps({p.relative_to(output).as_posix():digest(p) for p in files},indent=2))
    archive=output.with_suffix('.zip')
    if archive.exists():raise FileExistsError(archive)
    with zipfile.ZipFile(archive,'w',zipfile.ZIP_DEFLATED) as z:
        for p in sorted(output.rglob('*')):
            if p.is_file():z.write(p,p.relative_to(output.parent))
    with zipfile.ZipFile(archive) as z:assert z.testzip() is None
    print(json.dumps(dict(output=str(output),archive=str(archive),archive_bytes=archive.stat().st_size,
                         clusters=len(groups),members=len(rows),pairs=len(pairs),input_CIFs=len(source_map)),indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source',type=Path,required=True);parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();export(args.source,args.output)
