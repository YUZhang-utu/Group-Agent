"""E120 saved-pose check and time-budgeted refinement of the existing direct pool."""
import argparse
from concurrent.futures import ProcessPoolExecutor,wait,FIRST_COMPLETED
import hashlib
import json
import math
from pathlib import Path
import sqlite3
import time
from types import SimpleNamespace
import numpy as np
import kras_contact_refine as refine
import kras_soft_rank as soft

WORK_SETTINGS=None


def code_seal():
    return {p.name:hashlib.sha256(p.read_bytes().replace(b'\r\n',b'\n')).hexdigest() for p in
            (Path(__file__),Path(refine.__file__),Path(soft.__file__))}


def setup(batch,prepared,counts,settings):
    global WORK_SETTINGS
    soft.initialize(batch,prepared,counts);WORK_SETTINGS=settings


def process(task):
    index,ids=task;reader,features,states,chemistry=soft.run.WORK
    rows={t:dict(rows=[],pool_conformers=0,seeds=0,direct_poses=0,optimized_starts=0,
                 geometry_conformers=0) for t in states}
    for gid in ids:
        candidate=reader.get(gid)
        targets=[t for t in states if soft.COUNTS.execute('SELECT 1 FROM molecules WHERE target=? AND mid=?',
                                                        (t,candidate.molecule_id)).fetchone()]
        if not targets:continue
        chem=chemistry.get(gid);feature=features.get(gid)
        if any((r.molecule_id,r.conformer_id)!=(candidate.molecule_id,candidate.conformer_id) for r in (chem,feature)):
            raise ValueError('Chemical/feature identity mismatch')
        bucket=None
        for target in targets:
            row=rows[target];row['pool_conformers']+=1
            pose,stats=refine.refine_candidate(candidate,feature,chem,states[target],
                starts=WORK_SETTINGS['starts'],max_nfev=WORK_SETTINGS['max_nfev'])
            for key,value in stats.items():row[key]+=value
            if pose:
                if bucket is None:bucket=soft.diversity_bucket(chem)
                pose['bucket']=bucket;row['rows'].append(pose)
    return index,len(ids),rows


def native_checks(prepared):
    results={}
    for target in soft.run.kw.REFERENCES:
        state=soft.run.kw.load_state(prepared,target);q=state['query'];chem=SimpleNamespace(**state['record']['chemistry'])
        sets=soft.run.kw.chemical_sets(chem,q['feature_types']);initial=np.eye(4);initial[:3,3]=[.6,-.3,.2]
        matrix,evidence=refine.optimize_pose(q['shape_points'],initial,state,sets)
        if matrix is None:raise ValueError('Perturbed native constrained-refinement control failed: '+target)
        moved=q['shape_points']@matrix[:3,:3].T+matrix[:3,3]
        direct=soft.run.kw.direct_contact_masks(q['shape_points'],sets,state,matrix[None])['all'][0]
        minimum=float(state['tree'].query(moved)[0].min())
        if not direct or minimum<2:raise ValueError('Native endpoint infeasible')
        results[target]=dict(minimum_distance=minimum,all_contacts=True,evidence=evidence)
    return results


def finish_report(db,out,scheduled,scope,started,added,common,extra=None):
    result=soft.summary(db,scheduled)
    result.update(scope=scope,full_catalog=False,invocation_seconds=time.monotonic()-started,
        invocation_new_conformers=added,requested_per_target=1000000,
        recheck_source_pose_required=True,
        limitations=['Restricted to the previous direct-contact molecule pool; not a new full-library search.',
                     'Rigid optimization only; no internal torsion or receptor flexibility.',
                     'All original direct contacts and the 2 A distance floor rechecked in one pose.',
                     'Contact proximity is not hydrogen-bond angle/energy validation.',
                     'Partial coverage and quotas are explicit; no affinity or selectivity claim.'])
    if extra:result.update(extra)
    db.commit();db.close()
    result['ranking_sha256']=soft.run.kw.sha(out/'ranking.sqlite')
    soft.run.kw.save(out/'report.json',result)
    print(json.dumps(result,indent=2),flush=True)
    return result


def run_gate(args,common,batch,prepared,counts):
    kw=soft.run.kw;out=args.output
    if out.exists():raise ValueError('Use a new gate directory')
    audit=kw.read(args.examples)
    if audit.get('status')!='complete' or audit.get('scope')!='diagnostic_only_not_candidates':
        raise ValueError('Complete E119 examples required')
    if not all(r.get('matches_pilot_direct_count') and r.get('matches_pilot_geometry_count') for r in audit['targets'].values()):
        raise ValueError('E119 failed count reconciliation')
    out.mkdir(parents=True)
    kw.save(out/'protocol.json',dict(common=common,examples_sha256=kw.sha(args.examples)))
    kw.save(out/'report.json',dict(status='running',scope='saved_pose_refinement_gate'))
    started=time.monotonic();controls=native_checks(prepared);kw.save(out/'controls.json',controls)
    setup(batch,prepared,counts,refine.SETTINGS)
    reader,features,states,chemistry=soft.run.WORK
    db=sqlite3.connect(out/'ranking.sqlite');soft.create_tables(db)
    examples=[];index=0
    for target in kw.REFERENCES:
        for example in audit['targets'][target]['best_boundary_examples']:
            gid=example['gid'];candidate=reader.get(gid);chem=chemistry.get(gid);feature=features.get(gid)
            identity=(example['mid'],example['cid'])
            if any((r.molecule_id,r.conformer_id)!=identity for r in (candidate,chem,feature)):
                raise ValueError('Saved example identity mismatch')
            if not soft.COUNTS.execute('SELECT 1 FROM molecules WHERE target=? AND mid=?',(target,candidate.molecule_id)).fetchone():
                raise ValueError('Saved example not in direct pool')
            matrix=np.array(example['transform']).reshape(4,4);state=states[target]
            initial=candidate.shape_points@matrix[:3,:3].T+matrix[:3,3]
            minimum=float(state['tree'].query(initial)[0].min())
            if abs(minimum-example['minimum_distance'])>1e-6:raise ValueError('Saved example coordinates do not replay')
            sets=kw.chemical_sets(chem,candidate.feature_types)
            if not kw.direct_contact_masks(candidate.shape_points,sets,state,matrix[None])['all'][0]:
                raise ValueError('Saved example contacts do not replay')
            recovered,evidence=refine.optimize_pose(candidate.shape_points,matrix,state,sets)
            pose=refine.score_matrix(candidate,feature,chem,state,recovered) if recovered is not None else None
            if recovered is not None and pose is None:raise ValueError('Refined gate endpoint invalid')
            row=dict(rows=[],pool_conformers=1,seeds=1,direct_poses=1,optimized_starts=1,geometry_conformers=int(pose is not None))
            if pose:
                pose.update(bucket=soft.diversity_bucket(chem),refinement=evidence);row['rows'].append(pose)
            soft.commit_chunk(db,(index,1,{target:row}));index+=1
            examples.append(dict(target=target,mid=candidate.molecule_id,cid=candidate.conformer_id,
                                 rescued=pose is not None,evidence=evidence,pose=pose))
            print(json.dumps(dict(target=target,gid=gid,rescued=pose is not None,
                                  endpoint_minimum_distance=evidence['endpoint_minimum_distance'])),flush=True)
    kw.save(out/'examples.json',examples)
    result=finish_report(db,out,index,'saved_pose_refinement_gate',started,index,common,
                        dict(native_controls_passed=True,selected_examples_not_yield_estimate=True))
    soft.COUNTS.close();soft.COUNTS=None;soft.run.WORK=None
    return result


def spread_tasks(total,chunk_size):
    chunks=math.ceil(total/chunk_size);step=104729
    while math.gcd(step,chunks)!=1:step+=2
    for ordinal in range(chunks):
        index=(ordinal*step)%chunks;start=index*chunk_size
        yield index,list(range(start,min(start+chunk_size,total)))


def run_budget(args,common,batch,prepared,counts):
    kw=soft.run.kw;out=args.output
    gate=kw.read(args.gate/'report.json');gate_protocol=kw.read(args.gate/'protocol.json')
    if gate.get('status')!='complete' or gate.get('scope')!='saved_pose_refinement_gate' or gate_protocol['common']!=common:
        raise ValueError('Matching completed saved-pose gate required')
    if not gate.get('native_controls_passed') or not any(r['eligible_molecules'] for r in gate['targets'].values()):
        raise ValueError('No real saved pose rescued; stop instead of launching blind bulk optimization')
    if kw.sha(args.gate/'ranking.sqlite')!=gate['ranking_sha256']:raise ValueError('Gate ranking changed')
    signature=dict(common=common,chunk_size=args.chunk_size,gate_protocol_sha256=kw.sha(args.gate/'protocol.json'))
    out.mkdir(parents=True,exist_ok=True)
    with soft.file_lock(out/'run.lock'):
        if (out/'protocol.json').exists() and kw.read(out/'protocol.json')!=signature:
            raise ValueError('Resume mismatch; use a new output directory')
        kw.save(out/'protocol.json',signature)
        kw.save(out/'report.json',dict(status='running',scope='budgeted_direct_pool_refinement'))
        db=sqlite3.connect(out/'ranking.sqlite');soft.create_tables(db)
        # Warm-start with checked saved poses; no extra corpus coverage is claimed.
        with soft.readonly(args.gate/'ranking.sqlite') as source:
            if not db.execute('SELECT 1 FROM chunks WHERE id=-1').fetchone():
                warm={t:dict(rows=[]) for t in kw.REFERENCES}
                for target,payload in source.execute('SELECT target,payload FROM molecules'):
                    warm[target]['rows'].extend(json.loads(payload))
                soft.commit_chunk(db,(-1,0,warm))
        total=soft.run.catalog_count(kw.read(batch/'artifacts/catalog.json'))
        done={r[0] for r in db.execute('SELECT id FROM chunks')}
        tasks=(t for t in spread_tasks(total,args.chunk_size) if t[0] not in done)
        started=time.monotonic();last=started;added=0;stop_reason='plan_complete'
        def progress():
            eligible={t:db.execute('SELECT COUNT(*) FROM molecules WHERE target=?',(t,)).fetchone()[0] for t in kw.REFERENCES}
            row=dict(invocation_seconds=time.monotonic()-started,invocation_new_conformers=added,
                     eligible_molecules=eligible,requested_per_target=args.count)
            kw.save(out/'progress.json',row);print(json.dumps(row),flush=True)
            return eligible
        with ProcessPoolExecutor(max_workers=args.workers,initializer=setup,initargs=(batch,prepared,counts,refine.SETTINGS)) as pool:
            active=set();exhausted=False
            while active or not exhausted:
                while len(active)<args.workers and not exhausted:
                    if time.monotonic()-started>=args.hours*3600:
                        exhausted=True;stop_reason='time_budget';break
                    task=next(tasks,None)
                    if task is None:exhausted=True;break
                    active.add(pool.submit(process,task))
                if not active:break
                finished,active=wait(active,return_when=FIRST_COMPLETED)
                for future in finished:
                    result=future.result();soft.commit_chunk(db,result);added+=result[1]
                if time.monotonic()-last>=30:
                    eligible=progress();last=time.monotonic()
                    if all(n>=args.count for n in eligible.values()):
                        exhausted=True;stop_reason='molecule_budget'
        progress()
        return finish_report(db,out,total,'budgeted_direct_pool_refinement',started,added,common,
                             dict(stop_reason=stop_reason,requested_per_target=args.count,
                                  native_controls_passed=True,budget_finished=True))


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('mode',choices=['gate','batch'])
    p.add_argument('--direct',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--examples',type=Path);p.add_argument('--gate',type=Path)
    p.add_argument('--hours',type=float,default=4);p.add_argument('--workers',type=int,default=20)
    p.add_argument('--count',type=int,default=1000000);p.add_argument('--chunk-size',type=int,default=4)
    args=p.parse_args()
    if args.hours<=0 or not math.isfinite(args.hours) or min(args.workers,args.count,args.chunk_size)<1:p.error('Invalid budget')
    if args.mode=='gate' and not args.examples:p.error('--examples required')
    if args.mode=='batch' and not args.gate:p.error('--gate required')
    args.output=args.output.resolve();args.direct=args.direct.resolve();kw=soft.run.kw
    protocol=kw.read(args.direct/'protocol.json');report=kw.read(args.direct/'report.json')
    if report.get('status')!='complete' or not report.get('full_catalog') or report.get('scope')!='direct_contacts_only_not_final_screening_hits':
        p.error('Complete direct pool required')
    original=protocol['original'];batch=Path(original['batch']);prepared=Path(original['prepared']);counts=args.direct/'counts.sqlite'
    for protected in (args.direct,batch,prepared,args.gate):
        if protected is None:continue
        a,b=args.output,protected.resolve()
        if a.is_relative_to(b) or b.is_relative_to(a):p.error('Output overlaps protected input')
    if args.examples and args.examples.resolve().is_relative_to(args.output):p.error('Output contains input examples')
    if soft.run.seal_inputs(batch,prepared,original['chunk_size'],original['sample_size'])!=original:
        p.error('Original code/input seal mismatch')
    common=dict(original=original,settings=refine.SETTINGS,code=code_seal(),
                direct_report_sha256=kw.sha(args.direct/'report.json'),counts_sha256=kw.sha(counts))
    if args.mode=='gate':run_gate(args,common,batch,prepared,counts)
    else:run_budget(args,common,batch,prepared,counts)


if __name__=='__main__':main()
