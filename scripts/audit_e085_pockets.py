"""Reproducible E085 controlled distance audit; writes fresh CSV/NPZ artifacts."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import shutil

import numpy as np
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform
from scipy.spatial import cKDTree
from Bio.PDB.MMCIF2Dict import MMCIF2Dict

from aidd_agent import pocket_states as ps


def dump(path, value):
    path.write_text(json.dumps(value,indent=2),encoding='utf-8')


def table(path, rows):
    rows=list(rows)
    if not rows: return
    with path.open('w',newline='',encoding='utf-8-sig') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)


def matrix(path, values, ids):
    table(path,[dict(id=sid,**dict(zip(ids,map(float,row)))) for sid,row in zip(ids,values)])


def run(old_path, new_path, output):
    out=Path(output);out.mkdir(parents=True,exist_ok=False)
    old=json.loads(Path(old_path).read_text());new=json.loads(Path(new_path).read_text())
    ps.check_hashes(old['sources']);ps.check_hashes(new['sources'])
    oldz=np.load(old['artifact']);z=np.load(new['artifact'])
    rows=new['structures'];ids=[r['id'] for r in rows];oldids=[r['id'] for r in old['structures']]
    assert len(ids)==len(oldids)-sum(r['pdb_id']=='7BJ0' for r in old['structures'])
    assert not any(r['pdb_id']=='7BJ0' for r in rows)
    index=[oldids.index(sid) for sid in ids]
    assert np.array_equal(z['grid'],oldz['grid'])
    assert np.array_equal(z['masks'],oldz['masks'][index])
    grid=z['grid'];masks=z['masks'];p=new['policy']
    source=Path(new['diversity_report']).parent
    reference=ps.read_structure(source/'structures/5C5A.cif',new['target'])
    rm=ps.ca_map(reference['chains']['A']);centers=[rm[r] for r in new['pocket_residues']]
    diagnostics=ps.local_diagnostics(masks,grid,centers,p)
    assert np.allclose(diagnostics['legacy'],oldz['local_difference'][np.ix_(index,index)])
    metadata=MMCIF2Dict(str(source/'structures/22IZ.cif'))
    cache=source/'raw/7f451a28613a21011ff0b189a0f3805fc9035f484fc76282a6c1cbb30746575f.json'
    assert any(r['identifier']=='22IZ' for r in json.loads(cache.read_text())['result_set'])
    assert metadata['_entry.id']==['22IZ']
    dump(out/'01_22IZ_provenance.json',dict(id='22IZ',entry_id=metadata['_entry.id'],
        title=metadata['_struct.title'],source_query=json.loads((source/'entries.json').read_text())['query'],
        cached_search=str(cache),cif=str(source/'structures/22IZ.cif'),
        hashes=ps.fingerprint([cache,source/'entries.json',source/'structures/22IZ.cif']),
        authoritative_url='https://www.rcsb.org/structure/22IZ',
        official_release_date='2026-08-26',official_resolution_angstrom=1.52))
    dump(out/'02_exclusions.json',dict(excluded_entries=new['excluded_entries'],
        records=[r for r in new['held_for_review'] if r['id'].split(':')[0]=='7BJ0'],
        retained_chains=len(rows),retained_entries=len({r['pdb_id'] for r in rows})))
    table(out/'03_field_volumes.csv',[dict(id=r['id'],channel=k,
        source_features=r['feature_source_counts'][k],
        old_volume=float(oldz['chemical_fields'][oi,ki].sum()*p['spacing']**3),
        corrected_volume=r['field_volumes'][k])
        for r,oi in zip(rows,index) for ki,k in enumerate(ps.CHANNELS)])
    oldc1=next(c for c in old['clusters'] if c['representative']=='6Y4Q:A')
    ci=[ids.index(sid) for sid in oldc1['members']];rep=ids.index('6Y4Q:A')
    # Fixed diagnostic cuts: 50% raw local mismatch, or 25% support-weighted mismatch.
    # Complete linkage limits every pair; these are candidate subdivisions, not adopted states.
    candidates={}
    for mode,cut in [('raw',.5),('smooth',.25)]:
        d=diagnostics[mode][np.ix_(ci,ci)]
        candidates[mode]=fcluster(linkage(squareform(d,checks=True),method='complete'),cut,criterion='distance')
    c1rows=[dict(id=ids[i],legacy_stratum='zero' if diagnostics['legacy'][i,rep]==0 else 'nonzero',
        legacy_local_to_old_rep=float(diagnostics['legacy'][i,rep]),
        raw_local_to_old_rep=float(diagnostics['raw'][i,rep]),
        smooth_local_to_old_rep=float(diagnostics['smooth'][i,rep]),
        max_changed_region_volume=float(diagnostics['maximum_changed_volume'][i,rep]),
        corrected_distance_to_old_rep=float(z['distance'][i,rep]),
        raw_local_candidate=int(candidates['raw'][j]),smooth_local_candidate=int(candidates['smooth'][j]))
        for j,i in enumerate(ci)]
    table(out/'04_C1_member_reanalysis.csv',c1rows)
    # Test the user's two-mode hypothesis directly, without calling strata states.
    strata={name:[ci[j] for j,r in enumerate(c1rows) if r['legacy_stratum']==name]
            for name in ['zero','nonzero']}
    dump(out/'04b_C1_requested_two_strata.json',
         {name:[ids[i] for i in group] for name,group in strata.items()})
    pair_summary=[]
    for a,b in [('zero','zero'),('nonzero','nonzero'),('zero','nonzero')]:
        pairs=[(i,j) for i in strata[a] for j in strata[b] if a!=b or i>j]
        for mode,values in {**diagnostics,'corrected_distance':z['distance']}.items():
            v=np.array([values[i,j] for i,j in pairs])
            pair_summary.append(dict(first=a,second=b,metric=mode,pairs=len(v),
                mean=float(v.mean()),median=float(np.median(v)),maximum=float(v.max())))
    table(out/'04c_C1_stratum_pair_evidence.csv',pair_summary)
    local_splits=[]
    for cutoff in np.round(np.arange(.20,.601,.05),2):
        for g in ps.clusters(diagnostics['smooth'][np.ix_(ci,ci)],[rows[i] for i in ci],float(cutoff)):
            local_splits.extend(dict(cutoff=float(cutoff),id=sid,cluster=g['id'],
                representative=g['representative'],size=len(g['members']),
                max_member_distance=g['representative_max_distance'],
                radius_exceeds_cutoff=g['representative_max_distance']>cutoff) for sid in g['members'])
    table(out/'04d_C1_local_reclustering.csv',local_splits)
    table(out/'05_all_pairs.csv',[dict(first=ids[i],second=ids[j],
        old_distance=float(oldz['distance'][index[i],index[j]]),
        corrected_distance=float(z['distance'][i,j]),overlap=float(z['overlap'][i,j]),
        old_chemical=float(oldz['chemical_distance'][index[i],index[j]]),
        corrected_chemical=float(z['chemical_distance'][i,j]),
        **{f'local_{k}':float(v[i,j]) for k,v in diagnostics.items()})
        for i in range(len(ids)) for j in range(i)])
    for key in ['distance','chemical_distance','overlap']:
        matrix(out/f'matrix_{key}.csv',z[key],ids)
    for key in ['raw','smooth','legacy']:
        matrix(out/f'matrix_local_{key}.csv',diagnostics[key],ids)
    # Recover exactly the same aligned atoms; vary only descriptor shell width.
    aligned=[]
    for r in rows:
        atoms=ps.read_structure(r['structure_path'],new['target'])['chains'][r['target_chain']]
        transform=np.array(r['transform'])
        atoms=[dict(a,xyz=np.array(a['xyz'])@transform[:3,:3].T+transform[:3,3]) for a in atoms]
        near=cKDTree(grid).query(np.array([a['xyz'] for a in atoms]))[0]<8
        aligned.append([a for a,keep in zip(atoms,near) if keep])
    fields_by_width={2:z['chemical_fields']}
    support=[]
    for ri,(r,atoms) in enumerate(zip(rows,aligned)):
        free_tree=cKDTree(grid[masks[ri]])
        for channel,features in ps.atom_features(atoms,with_elements=True).items():
            gap=min((float(free_tree.query(point)[0])-p['feature_radius']-
                     (ps.RADII[element]+p['probe_radius'] if channel not in {'donor','acceptor'} else 0)
                     for point,element in features),default=float('inf'))
            volume=r['field_volumes'][channel]
            assert (gap<=0)==(volume>0)
            support.append(dict(id=r['id'],channel=channel,source_features=len(features),
                volume=volume,minimum_free_grid_distance_minus_field_radius=gap,
                explanation='intersects_accessible_grid' if volume else 'no_accessible_grid_in_field_support'))
    table(out/'10_field_support_audit.csv',support)
    for width in [1,3]:
        q=dict(p,feature_radius=width)
        fields_by_width[width]=np.array([ps.describe(grid,a,q)[1] for a in aligned])
        # Isolate nondirectional shell width; keep directional proxies unchanged.
        fields_by_width[width][:,:2]=fields_by_width[2][:,:2]
    sensitivity=[];membership=[];scenarios={}
    for width,scale in [(2,30),(1,30),(3,30),(2,10),(2,50)]:
        q=dict(p,feature_radius=width,local_change_volume=scale)
        distance,overlap,chemical,local=ps.compare(masks,fields_by_width[width],grid,centers,q)
        scenarios[f'shell{width}_scale{scale}_smooth']=(distance,chemical)
        if width==2 and scale==30:
            for mode in ['legacy','raw']:
                d=(1-q['chemical_weight'])*np.maximum(1-overlap,.5*diagnostics[mode])+q['chemical_weight']*chemical
                scenarios[f'shell2_scale30_{mode}']=(d,chemical)
    scenarios['old_fields_legacy_local_without_7BJ0']=(oldz['distance'][np.ix_(index,index)],oldz['chemical_distance'][np.ix_(index,index)])
    c6='8GCG:A'
    for scenario,(distance,chemical) in scenarios.items():
        np.savez_compressed(out/(scenario+'.npz'),ids=np.array(ids),distance=distance,chemical=chemical)
        for cutoff in np.round(np.arange(.20,.601,.05),2):
            groups=ps.clusters(distance,rows,float(cutoff))
            c6group=next(g for g in groups if c6 in g['members'])
            sensitivity.append(dict(scenario=scenario,cutoff=float(cutoff),clusters=len(groups),
                sizes=';'.join(str(n) for n in sorted([len(g['members']) for g in groups],reverse=True)),
                representative_radius_violations=int(sum(g['representative_max_distance']>cutoff for g in groups)),
                C6_cluster_size=len(c6group['members']),C6_representative=c6group['representative'],
                C6_members=';'.join(c6group['members'])))
            membership.extend(dict(scenario=scenario,cutoff=float(cutoff),id=sid,cluster=g['id'],
                representative=g['representative'],representative_max_distance=g['representative_max_distance'])
                for g in groups for sid in g['members'])
    table(out/'06_sensitivity.csv',sensitivity);table(out/'07_sensitivity_memberships.csv',membership)
    c1summary={mode:[dict(candidate=int(label),members=[ids[i] for i,l in zip(ci,labels) if l==label])
                    for label in sorted(set(labels))] for mode,labels in candidates.items()}
    dump(out/'08_C1_candidate_groups.json',c1summary)
    summary=dict(status='exploratory_not_adopted',n_chains=len(ids),
        dead_channels_before=[k for k in ps.CHANNELS if not any(oldz['chemical_fields'][:,ps.CHANNELS.index(k),:].ravel())],
        zero_volume_rows_after={k:sum(r['field_volumes'][k]==0 for r in rows) for k in ps.CHANNELS},
        source_absent_rows={k:sum(r['feature_source_counts'][k]==0 for r in rows) for k in ps.CHANNELS},
        C1_legacy_zero_to_rep=sum(r['legacy_stratum']=='zero' for r in c1rows),
        C1_raw_zero_to_rep=sum(r['raw_local_to_old_rep']==0 for r in c1rows),
        C1_candidate_sizes={k:[len(g['members']) for g in v] for k,v in c1summary.items()},
        corrected_default_clusters=len(new['clusters']),
        baseline_sensitivity=[r for r in sensitivity if r['scenario']=='shell2_scale30_smooth'],
        code_hashes=ps.fingerprint([Path(__file__),Path(ps.__file__)]),
        input_hashes=ps.fingerprint([old_path,new_path]))
    dump(out/'09_summary.json',summary)
    shutil.copy2(new_path,out/'corrected_report.json');shutil.copy2(new['artifact'],out/'corrected_grids.npz')
    shutil.copy2(Path(new_path).parent/'report.html',out/'pocket_comparison.html')
    dump(out/'SHA256SUMS.json',{str(f.relative_to(out)):hashlib.sha256(f.read_bytes()).hexdigest()
                              for f in sorted(out.rglob('*')) if f.is_file()})
    print(json.dumps(summary,indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--old',required=True);p.add_argument('--new',required=True)
    p.add_argument('--output',required=True);args=p.parse_args();run(args.old,args.new,args.output)
