"""Sample retained-pool molecules and measure exploratory protein-grid retention."""
import argparse
from concurrent.futures import ProcessPoolExecutor,wait,FIRST_COMPLETED
import json
import math
from pathlib import Path
import time
import numpy as np
import kras_pocket_grid as grid
import kras_soft_rank as soft
import export_kras_soft_rank as base

WORK=None


def sample_members(db,n,seed):
    rng=np.random.default_rng(seed);members={};totals={}
    for target in base.TARGETS:
        total=db.execute('SELECT COUNT(*) FROM molecules WHERE target=?',(target,)).fetchone()[0]
        totals[target]=total
        ranks=iter(sorted(rng.choice(total,size=min(n,total),replace=False).tolist()))
        wanted=next(ranks,None)
        for i,(mid,) in enumerate(db.execute('SELECT mid FROM molecules WHERE target=? ORDER BY mid',(target,))):
            if wanted is None:break
            if i==wanted:
                members.setdefault(mid,[]).append(target);wanted=next(ranks,None)
    return members,totals


def map_sample(batch,members,limit):
    """Read compact ID arrays only; coordinates are loaded only for sampled IDs."""
    if any(len(mid.encode())>16 for mid in members):raise ValueError('Unsupported artifact identity width')
    wanted=np.array(sorted(members),dtype='S16');mapped={mid:[] for mid in members}
    catalog=base.read(batch/'artifacts/catalog.json')
    for shard in sorted(catalog['shards'],key=lambda x:x['global_id_start']):
        folder=Path(shard['path'])
        if not folder.is_dir():folder=batch/'artifacts'/shard['name']
        ids=np.memmap(folder/'molecule_ids.bin',mode='r',dtype='S16')
        cids=np.memmap(folder/'conformer_ids.bin',mode='r',dtype='S16')
        if len(ids)!=shard['conformers'] or len(cids)!=len(ids):raise ValueError('ID array count mismatch')
        for start in range(0,len(ids),250000):
            chunk=ids[start:start+250000]
            for local in np.flatnonzero(np.isin(chunk,wanted)):
                i=start+int(local);mid=bytes(ids[i]).decode()
                if len(mapped[mid])<limit:
                    mapped[mid].append(dict(gid=int(shard['global_id_start'])+i,cid=bytes(cids[i]).decode()))
        del ids,cids
    return mapped


def initialize(batch,prepared,grid_root):
    global WORK
    soft.run.initialize(Path(batch),Path(prepared))
    grids={}
    for target in base.TARGETS:
        with np.load(Path(grid_root)/(target+'.npz')) as data:
            grids[target]={k:data[k].copy() for k in data.files}
    WORK=grids


def process(task):
    reader,features,states,chemistry=soft.run.WORK
    output=[]
    for mid,targets,conformers in task:
        rows={target:dict(mid=mid,target=target,conformers=0,direct_poses=0,exact_poses=0,
            outside_grid_poses=0,best=None,profile_poses={},profiles={p:False for p in grid.PROFILES}) for target in targets}
        for item in conformers:
            candidate=reader.get(item['gid']);chem=chemistry.get(item['gid'])
            if any((v.molecule_id,v.conformer_id)!=(mid,item['cid']) for v in (candidate,chem)):
                raise ValueError('Mapped artifact/chemistry identity mismatch')
            for target in targets:
                result=grid.evaluate(candidate,chem,states[target],WORK[target]);row=rows[target]
                row['conformers']+=1
                for key in ('direct_poses','exact_poses','outside_grid_poses'):row[key]+=result[key]
                for name in grid.PROFILES:row['profiles'][name] |= result['profiles'][name]
                for name,pose in result['profile_poses'].items():
                    if name not in row['profile_poses'] or pose['score']>row['profile_poses'][name]['score']:
                        row['profile_poses'][name]=pose
                if result['best'] and (row['best'] is None or result['best']['score']>row['best']['score']):
                    row['best']=result['best']
        output.extend(rows.values())
    return output


def wilson(k,n):
    if not n:return [0.,1.]
    z=1.95996398454;p=k/n;den=1+z*z/n
    center=(p+z*z/(2*n))/den;half=z*math.sqrt(p*(1-p)/n+z*z/(4*n*n))/den
    return [max(0.,center-half),min(1.,center+half)]


def run(args):
    start=time.monotonic();plan=args.plan.resolve();out=args.output.resolve()
    report=base.read(plan/'report.json');protocol=base.read(plan/'protocol.json')
    if report.get('status')!='complete' or report.get('scope')!='retained_pool_mapping':
        raise ValueError('Completed E121 plan required, not export directory')
    if base.sha(plan/'pool.sqlite')!=report['pool_sha256']:raise ValueError('Pool checksum mismatch')
    original=protocol['original'];batch=Path(original['batch']);prepared=Path(original['prepared'])
    for path in (plan,batch.resolve(),prepared.resolve()):
        if out.is_relative_to(path) or path.is_relative_to(out):raise ValueError('Output overlaps protected input')
    if out.exists():raise ValueError('Use a new output directory')
    if soft.run.seal_inputs(batch,prepared,original['chunk_size'],original['sample_size'])!=original:
        raise ValueError('Original code/input seal mismatch')
    out.mkdir(parents=True)
    base.save(out/'report.json',dict(status='running',stage='grid_setup'))
    controls={}
    for target,pdb in [('KRAS_G12D','7RPZ'),('KRAS_G12V','12FF')]:
        state=soft.run.kw.load_state(prepared,target)
        if state['record']['query_id'].split(':')[0]!=pdb:raise ValueError('Unexpected receptor')
        g=grid.make_grid(state);np.savez_compressed(out/(target+'.npz'),**g)
        controls[target]=grid.controls(state,g)
    base.save(out/'controls.json',controls)
    if not all(c['native_contacts'] and c['native_profiles']['lenient'] and c['forced_overlap_detected'] for c in controls.values()):
        base.save(out/'report.json',dict(status='blocked',reason='native_control_failed',controls=controls))
        raise ValueError('Native controls failed; no candidate run launched')
    print(json.dumps(dict(stage='controls_complete',controls=controls)),flush=True)
    with base.readonly(plan/'pool.sqlite') as db:
        members,totals=sample_members(db,args.sample,args.seed)
    print(json.dumps(dict(stage='mapping_sample_IDs',unique_sample_ids=len(members))),flush=True)
    mapped=map_sample(batch,members,grid.SETTINGS['conformers'])
    base.save(out/'sample.json',dict(members=members,conformers=mapped))
    base.save(out/'protocol.json',dict(original=original,pool_sha256=report['pool_sha256'],
        sample=args.sample,seed=args.seed,settings=grid.SETTINGS,profiles=grid.PROFILES,
        code={p.name:base.sha(p) for p in (Path(__file__),Path(grid.__file__))},
        sample_sha256=base.sha(out/'sample.json'),target_pool_molecules=totals))
    setup_seconds=time.monotonic()-start;started=time.monotonic()
    jobs=[(mid,members[mid],mapped[mid]) for mid in sorted(members)]
    tasks=iter([jobs[i:i+8] for i in range(0,len(jobs),8)])
    results=[];last=started
    with (out/'molecules.jsonl').open('w',encoding='utf-8',newline='\n') as stream:
        with ProcessPoolExecutor(max_workers=args.workers,initializer=initialize,initargs=(batch,prepared,out)) as pool:
            active=set();finished=False
            while active or not finished:
                while len(active)<args.workers*2 and not finished:
                    task=next(tasks,None)
                    if task is None:finished=True;break
                    active.add(pool.submit(process,task))
                if not active:break
                done,active=wait(active,return_when=FIRST_COMPLETED)
                for future in done:
                    rows=future.result();results.extend(rows)
                    for row in rows:stream.write(json.dumps(row,allow_nan=False)+'\n')
                if time.monotonic()-last>=15:
                    stream.flush();last=time.monotonic()
                    print(json.dumps(dict(stage='grid_scoring',completed_target_molecules=len(results),
                        compute_seconds=last-started)),flush=True)
    elapsed=time.monotonic()-started;targets={}
    for target in base.TARGETS:
        rows=[r for r in results if r['target']==target];n=len(rows)
        entry=dict(sampled_molecules=n,pool_molecules=totals[target],
            evaluated_conformers=sum(r['conformers'] for r in rows),
            missing_mapping=sum(not r['conformers'] for r in rows),
            direct_contact_molecules=sum(r['direct_poses']>0 for r in rows),
            exact_scored_molecules=sum(r['best'] is not None for r in rows),profiles={})
        for name in grid.PROFILES:
            k=sum(r['profiles'][name] for r in rows);interval=wilson(k,n)
            entry['profiles'][name]=dict(retained=k,fraction=k/n if n else 0,
                projected_pool_retained=round(totals[target]*k/n) if n else None,
                projected_95pct_interval=[round(totals[target]*v) for v in interval],
                native_control_passed=controls[target]['native_profiles'][name])
        targets[target]=entry
    full_upper=sum(totals.values())/max(1,len(results))*elapsed/3600
    result=dict(status='complete',scope='sampled_pool_grid_pilot',targets=targets,
        setup_and_ID_mapping_seconds=setup_seconds,compute_seconds=elapsed,
        full_pool_compute_hours_linear_estimate=full_upper,
        estimate_scope='Target-membership linear estimate; sharing/composition may change throughput; excludes export.',
        profiles_are_exploratory=True,production_selection_created=False,
        limitations=['Finite rigid seeds, at most three existing conformers; failure is not molecular impossibility.',
            'Exact recheck covers eight top grid poses per conformer; other poses remain unresolved.',
            'Burial is six-ray protein-occlusion proxy, not a cavity volume or ligand similarity.',
            'Profiles allow clashes and are prioritization, not validated binding poses.',
            'Wilson intervals cover sampling uncertainty only, not pose-search or model errors.',
            'Original 12FF/7RPZ pool only; no full-library recall or affinity claim.'])
    base.save(out/'report.json',result);print(json.dumps(result,indent=2),flush=True)
    return result


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--plan',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--sample',type=int,default=10000);p.add_argument('--workers',type=int,default=20)
    p.add_argument('--seed',type=int,default=122)
    args=p.parse_args()
    if min(args.sample,args.workers)<1:p.error('Positive sample and worker counts required')
    run(args)


if __name__=='__main__':main()
