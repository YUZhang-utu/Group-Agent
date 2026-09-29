"""Source-backed boundary routing with disjoint reference calibration and explicit residuals."""
import argparse
from collections import Counter,defaultdict
import hashlib
import json
from pathlib import Path
import sqlite3
import time

import numpy as np

from .boundary_pair_review import sha,readonly,prepare_record
from .boundary_block_candidates import WIDTHS
from .mol2 import iter_mol2_blocks,parse_mol2_block


def priority(mid):return hashlib.sha256(('E091-v1:'+mid).encode()).hexdigest()


def feature_matrix(value,row,variant):
    steric=np.asarray([[r['side_heavy_atoms']/10,r['proximal_branch_excess']/3,
        *np.asarray(r['local_lower_vdw'])/5,*np.asarray(r['local_upper_vdw'])/5,
        *np.sqrt(np.maximum(r['spatial_eigenvalues'],0))/5,r['maximum_reach']/5] for r in value['sterics']])
    n=len(value['units']);backbone=np.asarray(row['descriptor']).reshape(n,WIDTHS[variant])[:,:6]
    vector=np.concatenate([np.asarray(value['properties']),steric,backbone],axis=1)
    if vector.shape!=(n,29) or not np.isfinite(vector).all():raise ValueError('Invalid continuous feature vector')
    return vector


def distances(query,prototypes):
    delta=(prototypes-query)**2
    return np.stack([np.sqrt(delta[:,:,a:b].mean(axis=(1,2))) for a,b in [(0,11),(11,23),(23,29)]],axis=1)


def diverse_indices(vectors,limit):
    selected=[0];nearest=distances(vectors[0],vectors).mean(1)
    while len(selected)<min(limit,len(vectors)):
        nearest[selected]=-1;index=int(nearest.argmax());selected.append(index)
        nearest=np.minimum(nearest,distances(vectors[index],vectors).mean(1))
    return selected


def fit_class(records,minimum_references=32,prototype_count=12):
    """records are (CID, MID, vector); molecular splits precede prototype selection."""
    records=sorted(records,key=lambda r:(priority(r[1]),r[0]))
    if len({r[1] for r in records})!=len(records):raise ValueError('Reference molecules must be distinct')
    if len(records)<minimum_references:return dict(status='insufficient_distinct_references',references=len(records))
    n=len(records);a=n//2;b=a+(n-a)//2
    fit,cal,check=records[:a],records[a:b],records[b:]
    if min(len(fit),len(cal),len(check))<2:raise ValueError('Reference split too small')
    ids=diverse_indices(np.array([r[2] for r in fit]),prototype_count)
    prototypes=np.array([fit[i][2] for i in ids])
    d=np.array([distances(r[2],prototypes) for r in cal])
    scales=np.maximum(np.median(d,axis=(0,1)),1e-8)
    cal_scores=(d/scales).max(2).min(1)
    cutoff=max(float(np.quantile(cal_scores,.95)),1e-8)
    check_scores=[float((distances(r[2],prototypes)/scales).max(1).min()/cutoff) for r in check]
    retention=sum(v<=1 for v in check_scores)/len(check_scores)
    return dict(status='routable' if retention>=.8 else 'independent_check_failed',references=n,
        prototype_ids=[fit[i][0] for i in ids],prototype_vectors=prototypes.tolist(),
        scales=scales.tolist(),cutoff=cutoff,fit_ids=[r[0] for r in fit],
        calibration_ids=[r[0] for r in cal],check_ids=[r[0] for r in check],
        check_scores=check_scores,check_retention=retention)


def score_candidate(vector,model,rotation):
    if model['status']!='routable':return None
    if type(rotation) is not int or not 0<=rotation<len(vector):raise ValueError('Invalid candidate rotation')
    query=np.roll(vector,-rotation,axis=0)
    prototypes=model['_prototypes'] if '_prototypes' in model else np.asarray(model['prototype_vectors'])
    if query.shape!=prototypes.shape[1:]:raise ValueError('Candidate/reference dimensions disagree')
    scales=model['_scales'] if '_scales' in model else np.asarray(model['scales'])
    d=distances(query,prototypes);scores=(d/scales).max(1)/model['cutoff']
    i=int(scores.argmin())
    normalized=d[i]/scales/model['cutoff']
    return dict(score=float(scores[i]),prototype_id=model['prototype_ids'][i],
                property_rms=float(d[i,0]),steric_rms=float(d[i,1]),backbone_circular_rms=float(d[i,2]),
                normalized_channels=normalized.tolist(),limiting_channel=('property','steric','backbone')[int(normalized.argmax())])


def run(build,candidates,output,reference_limit=64):
    started=time.monotonic();timings={}
    if reference_limit<32:raise ValueError('At least 32 reference slots required')
    build=Path(build).resolve();candidates=Path(candidates).resolve();output=Path(output).resolve()
    if any(output==p or p in output.parents for p in (build,candidates)):raise ValueError('Fresh output outside inputs required')
    receipt=json.loads((build/'report.json').read_text());up=json.loads((candidates/'report.json').read_text())
    if receipt['status']!='complete' or up['status']!='complete':raise ValueError('Completed inputs required')
    if up['sources'].get(str(build/'report.json'))!=sha(build/'report.json'):raise ValueError('Build provenance mismatch')
    if sha(candidates/'candidates.sqlite')!=up['output_hashes']['candidates.sqlite']:raise ValueError('Candidate evidence changed')
    diag_path=next(Path(p) for p in up['sources'] if p!=str(build/'report.json') and Path(p).name=='report.json')
    if sha(diag_path)!=up['sources'][str(diag_path)]:raise ValueError('Diagnostic changed')
    diag=json.loads(diag_path.read_text());model_path=Path(diag['model'])
    if sha(model_path/'report.json')!=diag['provenance']['model_report_sha256']:raise ValueError('Original model changed')
    variant=receipt['schema']['variant'];output.mkdir(parents=True,exist_ok=False)
    cache_path=output/'features.sqlite';target_ids=set();references={};counts=Counter()
    with readonly(candidates/'candidates.sqlite') as incoming,readonly(model_path/'blocks.sqlite') as blocks, \
         readonly(build/'descriptors.sqlite') as source,sqlite3.connect(cache_path) as cache:
        cache.executescript('''CREATE TABLE wanted(cid TEXT PRIMARY KEY,role TEXT,group_id TEXT,path TEXT,idx INTEGER);
            CREATE TABLE feature(cid TEXT PRIMARY KEY,n INTEGER,vector BLOB,source_sha256 TEXT);''')
        def want(cid,role,expected_group):
            value=source.execute('SELECT payload FROM descriptor WHERE cid=?',(cid,)).fetchone()
            if value is None:raise ValueError('Missing source descriptor')
            row=json.loads(value[0]);p=row['provenance']
            if row['hard_group']!=expected_group:raise ValueError('Source class mismatch')
            cache.execute('INSERT INTO wanted VALUES(?,?,?,?,?)',(cid,role,expected_group,p['source_path'],p['source_record_index']))
        for cid,g,deviation,encoded in incoming.execute('SELECT cid,original_group,maximum_deviation,candidates FROM candidate'):
            matches=json.loads(encoded)
            if matches and deviation is not None and deviation<=45:
                want(cid,'boundary',g);target_ids.update(m['hard_group'] for m in matches);counts['boundary_features']+=1
                if counts['boundary_features']%10000==0:
                    cache.commit();print(f"Queued {counts['boundary_features']} boundary feature records",flush=True)
        # DISTINCT MID query is indexed by group; no chemistry is recomputed here.
        for number,g in enumerate(sorted(target_ids),1):
            chosen=[]
            for mid,cid in blocks.execute('SELECT mid,min(cid) FROM point WHERE group_id=? GROUP BY mid',(g,)):
                chosen.append((priority(mid),cid,mid))
                if len(chosen)>2*reference_limit:chosen=sorted(chosen)[:reference_limit]
            chosen=sorted(chosen)[:reference_limit];references[g]=[(cid,mid) for _,cid,mid in chosen]
            for _,cid,mid in chosen:want(cid,'reference',g)
            if number%25==0:cache.commit();print(f'Reference selection {number}/{len(target_ids)} classes',flush=True)
        cache.execute('CREATE UNIQUE INDEX source_record ON wanted(path,idx)');cache.commit()
        (output/'reference_manifest.json').write_text(json.dumps(references),encoding='utf-8')
        timings['selection_seconds']=time.monotonic()-started;phase=time.monotonic()
        paths=[r[0] for r in cache.execute('SELECT DISTINCT path FROM wanted ORDER BY path')]
        extracted=0
        for number,path in enumerate(paths,1):
            print(f'Feature extraction file {number}/{len(paths)}: {path}',flush=True)
            records=iter(cache.execute('SELECT idx,cid FROM wanted WHERE path=? ORDER BY idx',(path,)))
            wanted=next(records,None)
            for index,text in iter_mol2_blocks(Path(path)):
                if index and index%50000==0:print(f'Scanned {index} source records in current file',flush=True)
                if wanted is None:break
                if index!=wanted[0]:continue
                cid=wanted[1];row=json.loads(source.execute('SELECT payload FROM descriptor WHERE cid=?',(cid,)).fetchone()[0])
                record=parse_mol2_block(Path(path),index,text);p=prepare_record(record,row,variant)
                v=feature_matrix(p,row,variant)
                cache.execute('INSERT INTO feature VALUES(?,?,?,?)',(cid,len(v),v.astype('<f8').tobytes(),record.content_sha256))
                extracted+=1;wanted=next(records,None)
                if extracted%1000==0:cache.commit();print(f'Computed {extracted} source-verified profiles',flush=True)
            if wanted is not None:raise ValueError('Requested source records missing')
        cache.commit()
        if cache.execute('SELECT count(*) FROM feature').fetchone()[0]!=cache.execute('SELECT count(*) FROM wanted').fetchone()[0]:
            raise ValueError('Feature extraction incomplete')
        timings['source_feature_seconds']=time.monotonic()-phase;phase=time.monotonic()
        def vector(cid):
            row=cache.execute('SELECT n,vector FROM feature WHERE cid=?',(cid,)).fetchone()
            if row is None:raise ValueError('Missing computed feature')
            return np.frombuffer(row[1],dtype='<f8').reshape(row[0],29)
        models={g:fit_class([(cid,mid,vector(cid)) for cid,mid in ids]) for g,ids in references.items()}
        (output/'routing_models.json').write_text(json.dumps(models),encoding='utf-8')
        for m in models.values():
            if m['status']=='routable':
                m['_prototypes']=np.asarray(m['prototype_vectors']);m['_scales']=np.asarray(m['scales'])
        timings['model_fit_seconds']=time.monotonic()-phase;phase=time.monotonic()
        routed=Counter();reasons=Counter();by_band=defaultdict(Counter);out_path=output/'boundary_assignments.sqlite'
        with sqlite3.connect(out_path) as dest:
            dest.execute('CREATE TABLE assignment(cid TEXT PRIMARY KEY,original_group TEXT,state TEXT,target_group TEXT,rotation INTEGER,score REAL,evidence TEXT)')
            for i,(cid,g,deviation,encoded) in enumerate(incoming.execute('SELECT cid,original_group,maximum_deviation,candidates FROM candidate'),1):
                matches=json.loads(encoded);accepted=[];evaluated=[]
                category='up_to_35' if deviation is not None and deviation<=35 else '35_to_45' if deviation is not None and deviation<=45 else 'above_45_or_unresolved'
                if deviation is None or deviation>45:reason='outside_45_degree_policy_or_unresolved'
                elif not matches:reason='no_compatible_definite_class'
                else:
                    v=vector(cid)
                    for m in matches:
                        score=score_candidate(v,models[m['hard_group']],m['rotation'])
                        if score is not None:
                            evidence=dict(**m,**score);evaluated.append(evidence)
                            if score['score']<=1:accepted.append(evidence)
                    statuses={models[m['hard_group']]['status'] for m in matches}
                    reason=('all_target_references_sparse' if statuses=={'insufficient_distinct_references'}
                            else 'all_target_checks_failed' if statuses=={'independent_check_failed'}
                            else 'reference_models_failed_or_sparse') if not evaluated else 'outside_joint_calibrated_envelope'
                if accepted:
                    accepted.sort(key=lambda r:(r['score'],r['hard_group'],r['rotation']))
                    best=accepted[0];state='assigned_provisional';reason='joint_compatible_prototype'
                    dest.execute('INSERT INTO assignment VALUES(?,?,?,?,?,?,?)',(cid,g,state,best['hard_group'],best['rotation'],best['score'],json.dumps(dict(accepted=accepted,reason=reason))))
                else:
                    state='special';best=min(evaluated,key=lambda r:r['score']) if evaluated else None
                    if best is not None:reason='outside_calibrated_'+best['limiting_channel']
                    dest.execute('INSERT INTO assignment VALUES(?,?,?,?,?,?,?)',(cid,g,state,None,None,None,json.dumps(dict(reason=reason,best_rejected=best))))
                routed[state]+=1;reasons[reason]+=1;by_band[category][state]+=1
                if i%10000==0:dest.commit();print(f'Routed {i} boundary records',flush=True)
            dest.commit()
            if dest.execute('SELECT count(*) FROM assignment').fetchone()[0]!=up['boundary_conformers']:raise ValueError('Boundary membership accounting failed')
        timings['routing_seconds']=time.monotonic()-phase
    total=up['total_conformers'];regular=up['definite_conformers']+routed['assigned_provisional']
    result=dict(status='complete',readiness='provisional_routing_not_validated_search_dispatch',
        total_conformers=total,unchanged_definite_conformers=up['definite_conformers'],boundary_counts=dict(routed),
        provisional_regular_conformers=regular,special_conformers=routed['special'],
        special_fraction=routed['special']/total,engineering_target_special_fraction=.01,
        maximum_special_for_target=total//100,additional_assignments_needed_for_target=max(0,routed['special']-total//100),
        engineering_target_met=routed['special']/total<=.01,reasons=dict(reasons),bands={k:dict(v) for k,v in by_band.items()},
        class_model_statuses=dict(Counter(m['status'] for m in models.values())),
        wall_seconds=time.monotonic()-started,phase_timings=timings,special_search_benchmark='not_run',
        selected_reference_conformers=sum(map(len,references.values())),computed_profiles=extracted,
        logical_capacity_limit=None,sources={str(build/'report.json'):sha(build/'report.json'),str(candidates/'report.json'):sha(candidates/'report.json')},
        output_hashes={p:sha(output/p) for p in ('features.sqlite','reference_manifest.json','routing_models.json','boundary_assignments.sqlite')},
        code_hashes={p:sha(Path(__file__).with_name(p)) for p in ('calibrated_boundary_routing.py','boundary_pair_review.py','macrocycle_property_profiles.py','macrocycle_descriptors.py')},
        limitations=['Provisional candidate-specific placement; complete reference-search recall remains unvalidated',
            'Independent reference checks gate models; they are not an unbiased final evaluation after model selection',
            'Backbone routing uses circular torsions, not pairwise atom RMSD; original sin/cos features retained',
            'Sparse/failed reference models do not force acceptance; their candidates may remain special',
            'Ordinary source memberships remain unchanged; join boundary overrides by CID to the frozen model',
            'Original artifacts unchanged; no incremental growth or production dispatch adapter claimed',
            'Special fraction is measured, not forced to meet the engineering target'])
    if regular+routed['special']!=total:raise ValueError('Full population accounting failed')
    (output/'report.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    return result


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('build','candidates','output'):p.add_argument('--'+name,type=Path,required=True)
    p.add_argument('--reference-limit',type=int,default=64)
    a=p.parse_args();print(json.dumps(run(a.build,a.candidates,a.output,a.reference_limit),indent=2))


if __name__=='__main__':main()
