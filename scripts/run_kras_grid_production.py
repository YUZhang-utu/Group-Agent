"""Resume E122 scoring on the E121 pool and export budgeted docking inputs."""
import argparse
from concurrent.futures import ProcessPoolExecutor, wait, FIRST_COMPLETED
import csv
import hashlib
import itertools
import json
from pathlib import Path
import sqlite3
import time
import numpy as np
import run_kras_grid_pilot as pilot

base, grid, soft = pilot.base, pilot.grid, pilot.soft


def schema(db):
    db.executescript('''
    CREATE TABLE IF NOT EXISTS members(target TEXT,mid TEXT,PRIMARY KEY(target,mid));
    CREATE INDEX IF NOT EXISTS members_mid ON members(mid);
    CREATE TABLE IF NOT EXISTS mapping(mid TEXT,gid INTEGER,cid TEXT,PRIMARY KEY(mid,gid));
    CREATE TABLE IF NOT EXISTS results(target TEXT,mid TEXT,balanced INTEGER,
      lenient INTEGER,strict INTEGER,score REAL,payload TEXT,PRIMARY KEY(target,mid));
    CREATE INDEX IF NOT EXISTS results_rank ON results(target,balanced,score DESC,mid);
    CREATE TABLE IF NOT EXISTS state(key TEXT PRIMARY KEY,value TEXT);
    ''')


def map_pool(db, batch):
    if db.execute("SELECT 1 FROM state WHERE key='mapping_complete'").fetchone():
        return
    db.execute('DELETE FROM mapping'); db.commit()
    counts = {mid: 0 for mid, in db.execute('SELECT DISTINCT mid FROM members')}
    for shard in sorted(base.read(batch/'artifacts/catalog.json')['shards'], key=lambda s:s['global_id_start']):
        folder = Path(shard['path'])
        if not folder.is_dir(): folder = batch/'artifacts'/shard['name']
        ids = np.memmap(folder/'molecule_ids.bin', mode='r', dtype='S16')
        cids = np.memmap(folder/'conformer_ids.bin', mode='r', dtype='S16')
        if len(ids) != shard['conformers'] or len(cids) != len(ids):
            raise ValueError('ID array count mismatch')
        rows = []
        for i, value in enumerate(ids):
            mid = bytes(value).decode()
            n = counts.get(mid)
            if n is not None and n < grid.SETTINGS['conformers']:
                rows.append((mid, int(shard['global_id_start'])+i, bytes(cids[i]).decode()))
                counts[mid] += 1
            if len(rows) >= 10000:
                db.executemany('INSERT INTO mapping VALUES (?,?,?)', rows); rows.clear()
        db.executemany('INSERT INTO mapping VALUES (?,?,?)', rows); db.commit()
        print(json.dumps(dict(stage='mapping',shard=shard['name'])), flush=True)
        del ids, cids
    db.execute("INSERT INTO state VALUES ('mapping_complete','yes')"); db.commit()


def jobs(db):
    last=''
    while True:
        mids=[r[0] for r in db.execute('SELECT DISTINCT mid FROM members WHERE mid>? ORDER BY mid LIMIT 256',(last,)).fetchall()]
        if not mids: return
        last=mids[-1]
        for mid in mids:
            targets=[r[0] for r in db.execute('SELECT target FROM members m WHERE mid=? AND NOT EXISTS '
                '(SELECT 1 FROM results r WHERE r.target=m.target AND r.mid=m.mid) ORDER BY target',(mid,)).fetchall()]
            if targets:
                conformers=[dict(gid=g,cid=c) for g,c in db.execute('SELECT gid,cid FROM mapping WHERE mid=? ORDER BY gid',(mid,)).fetchall()]
                yield mid,targets,conformers


def store_rows(db, rows):
    for row in rows:
        pose = row['profile_poses'].get('balanced')
        db.execute('INSERT INTO results VALUES (?,?,?,?,?,?,?)',
            (row['target'],row['mid'],int(row['profiles']['balanced']),
             int(row['profiles']['lenient']),int(row['profiles']['strict']),
             pose['score'] if pose else None,json.dumps(row,allow_nan=False)))
    db.commit()


def counts(db):
    result = {}
    for target in base.TARGETS:
        n = db.execute('SELECT COUNT(*) FROM members WHERE target=?',(target,)).fetchone()[0]
        done,b,l,s = db.execute('SELECT COUNT(*),COALESCE(SUM(balanced),0),'
            'COALESCE(SUM(lenient),0),COALESCE(SUM(strict),0) FROM results WHERE target=?',(target,)).fetchone()
        result[target] = dict(pool=n,processed=done,balanced=b,lenient=l,strict=s,unresolved=done-l)
    return result


def selected(db, target, cap):
    return db.execute('SELECT mid,payload FROM results WHERE target=? AND balanced=1 '
                      'ORDER BY score DESC,mid LIMIT ?', (target,cap))


def source_check(molecule, target, mid, pose):
    reader, _, states, chemistry = soft.run.WORK
    candidate, chem = reader.get(pose['gid']), chemistry.get(pose['gid'])
    if any((v.molecule_id,v.conformer_id)!=(mid,pose['cid']) for v in (candidate,chem)):
        raise ValueError('Source/artifact identity mismatch')
    heavy = [a.GetIdx() for a in molecule.GetAtoms() if a.GetAtomicNum()>1]
    points = molecule.GetConformer().GetPositions()[heavy]
    if ([molecule.GetAtomWithIdx(i).GetAtomicNum() for i in heavy] != list(chem.atomic_numbers)
        or points.shape != candidate.shape_points.shape
        or not np.isfinite(points).all()
        or np.max(np.abs(points-candidate.shape_points)) > .02):
        raise ValueError('Source/artifact geometry mismatch')
    matrix=np.array(pose['transform']).reshape(4,4); state=states[target]
    sets=soft.run.kw.chemical_sets(chem,candidate.feature_types)
    contact=bool(soft.run.kw.direct_contact_masks(points,sets,state,matrix[None])['all'][0])
    moved=points@matrix[:3,:3].T+matrix[:3,3]
    metrics=grid.metrics(state['tree'].query(moved)[0],grid.burial(moved,state['tree']))
    return contact and grid.profile_flags(metrics)['balanced']


def export_target(db, target, batch, prepared, out, cap, batch_size):
    receipt=out/(target+'_export.json')
    if receipt.exists():
        saved=base.read(receipt)
        for entry in saved['batches']:
            if base.sha(Path(entry['path']))!=entry['sha256']:
                raise ValueError('Completed export modified')
        return saved
    attempts=out/'exports'; attempts.mkdir(exist_ok=True)
    index=1
    while (attempts/(target+f'_attempt_{index:03d}')).exists(): index+=1
    root=attempts/(target+f'_attempt_{index:03d}');root.mkdir()
    db.execute('DROP TABLE IF EXISTS export_selection')
    db.execute('CREATE TABLE export_selection(mid TEXT,cid TEXT,source TEXT,idx INTEGER,digest TEXT,pose TEXT)')
    wanted=0
    with base.readonly(batch/'registry.sqlite3') as registry:
        for mid,payload in selected(db,target,cap):
            pose=json.loads(payload)['profile_poses']['balanced']
            row=registry.execute('SELECT source_path,source_record_index,content_sha256 '
                'FROM conformer WHERE molecule_id=? AND id=?',(mid,pose['cid'])).fetchone()
            if row is None: raise ValueError('Scored conformer missing from registry')
            db.execute('INSERT INTO export_selection VALUES (?,?,?,?,?,?)',
                       (mid,pose['cid'],*row,json.dumps(pose)))
            wanted+=1
    db.execute('CREATE INDEX export_source ON export_selection(source,idx)');db.commit()
    soft.run.initialize(batch,prepared)
    writer=base.BatchWriter(root/'mol2',batch_size); quarantined=0
    try:
        with (root/'molecules.csv').open('w',encoding='utf-8',newline='') as handle, \
             (root/'quarantined.csv').open('w',encoding='utf-8',newline='') as qhandle:
            csvout=csv.writer(handle);qc=csv.writer(qhandle)
            csvout.writerow(['mid','cid','batch','record','score','source','source_record','source_sha256'])
            qc.writerow(['mid','cid','reason'])
            for source, in db.execute('SELECT DISTINCT source FROM export_selection ORDER BY source'):
                cursor=db.execute('SELECT mid,cid,idx,digest,pose FROM export_selection WHERE source=? ORDER BY idx',(source,))
                pending=next(cursor,None)
                for record,raw in base.iter_mol2_blocks(Path(source)):
                    if pending is None: break
                    if record<pending[2]: continue
                    if record!=pending[2]: raise ValueError('Source record missing')
                    digest=hashlib.sha256(raw.encode('utf-8')).hexdigest()
                    while pending is not None and pending[2]==record:
                        mid,cid,_,expected,payload=pending;pose=json.loads(payload)
                        if digest!=expected: raise ValueError('Source record hash mismatch')
                        reason=None
                        try:
                            base.parse_mol2_block(Path(source),record,raw)
                            molecule,mode=base.load_rdkit_mol2(raw,cid)
                            if mode!='strict': reason='non_strict_sanitization'
                            elif not source_check(molecule,target,mid,pose): reason='original_source_profile_or_contact_failed'
                        except Exception as exc:
                            reason=type(exc).__name__+':'+str(exc)
                        if reason:
                            quarantined+=1;qc.writerow([mid,cid,reason])
                        else:
                            b,r=writer.add(base.rename_record(raw,mid+'__'+cid))
                            csvout.writerow([mid,cid,b,r,pose['score'],source,record,digest])
                        pending=next(cursor,None)
                if pending is not None: raise ValueError('Source file truncated')
                handle.flush();qhandle.flush()
                print(json.dumps(dict(stage='export',target=target,exported=writer.count,quarantined=quarantined)),flush=True)
    finally:
        writer.close()
    if writer.count+quarantined!=wanted: raise ValueError('Export accounting mismatch')
    for entry in writer.batches: entry['sha256']=base.sha(Path(entry['path']))
    result=dict(status='complete',target=target,requested_cap=cap,selected=wanted,
        exported=writer.count,quarantined=quarantined,folder=str(root),batches=writer.batches)
    base.save(receipt,result)
    return result


def run(args):
    plan=args.plan.resolve();out=args.output.resolve()
    report=base.read(plan/'report.json');original=base.read(plan/'protocol.json')['original']
    batch=Path(original['batch']);prepared=Path(original['prepared'])
    if report.get('status')!='complete' or report.get('scope')!='retained_pool_mapping':
        raise ValueError('Completed E121 plan required')
    for path in (plan,batch.resolve(),prepared.resolve()):
        if out.is_relative_to(path) or path.is_relative_to(out): raise ValueError('Output overlaps input')
    if base.sha(plan/'pool.sqlite')!=report['pool_sha256']: raise ValueError('Pool checksum mismatch')
    if soft.run.seal_inputs(batch,prepared,original['chunk_size'],original['sample_size'])!=original:
        raise ValueError('Original source/code seal mismatch')
    seal=dict(plan=str(plan),pool_sha256=report['pool_sha256'],original=original,
        registry_sha256=base.sha(batch/'registry.sqlite3'),settings=grid.SETTINGS,profiles=grid.PROFILES,
        cap=args.cap,batch_size=args.batch_size,
        code={p.name:base.sha(p) for p in (Path(__file__),Path(pilot.__file__),Path(grid.__file__),Path(base.__file__))})
    # JSON normalization preserves tuple/list equivalence across resumed processes.
    seal=json.loads(json.dumps(seal))
    if out.exists():
        if not (out/'protocol.json').exists() or base.read(out/'protocol.json')!=seal:
            raise ValueError('Resume input/code/settings mismatch; use a new output directory')
    else:
        out.mkdir(parents=True);base.save(out/'protocol.json',seal)
    db=sqlite3.connect((out/'results.sqlite').as_uri(),uri=True);schema(db)
    db.execute('ATTACH DATABASE ? AS source_pool',((plan/'pool.sqlite').as_uri()+'?mode=ro',))
    db.execute('INSERT OR IGNORE INTO members SELECT target,mid FROM source_pool.molecules');db.commit()
    db.execute('DETACH DATABASE source_pool')
    controls={}
    for target,pdb in [('KRAS_G12D','7RPZ'),('KRAS_G12V','12FF')]:
        state=soft.run.kw.load_state(prepared,target)
        if state['record']['query_id'].split(':')[0]!=pdb: raise ValueError('Unexpected receptor')
        g=grid.make_grid(state);np.savez_compressed(out/(target+'.npz'),**g)
        controls[target]=grid.controls(state,g)
    base.save(out/'controls.json',controls)
    if not all(c['native_contacts'] and all(c['native_profiles'].values()) and c['forced_overlap_detected'] for c in controls.values()):
        raise ValueError('Native controls failed')
    base.save(out/'report.json',dict(status='running',stage='mapping',targets=counts(db)))
    map_pool(db,batch)
    started=time.monotonic();last=started
    db.execute('PRAGMA journal_mode=WAL')
    # Bounded keyset reads avoid a long snapshot that would prevent WAL checkpoints.
    with base.readonly(out/'results.sqlite') as scheduler:
        iterator=iter(jobs(scheduler));finished=False;active=set()
        with ProcessPoolExecutor(max_workers=args.workers,initializer=pilot.initialize,
                                 initargs=(batch,prepared,out)) as pool:
            while active or not finished:
                while len(active)<args.workers*2 and not finished:
                    task=list(itertools.islice(iterator,8))
                    if not task: finished=True;break
                    active.add(pool.submit(pilot.process,task))
                if not active: break
                done,active=wait(active,return_when=FIRST_COMPLETED)
                for future in done: store_rows(db,future.result())
                if time.monotonic()-last>=30:
                    last=time.monotonic()
                    progress=dict(status='running',stage='grid_scoring',invocation_seconds=last-started,targets=counts(db))
                    base.save(out/'report.json',progress);print(json.dumps(progress),flush=True)
    totals=counts(db)
    if any(v['processed']!=v['pool'] for v in totals.values()): raise ValueError('Incomplete pool accounting')
    diagnostics={t:dict(no_exact_pose=0,best_pose_severe_fail=0,best_pose_overlap_fail=0,
                       best_pose_burial_fail=0) for t in base.TARGETS}
    with (out/'all_molecule_tiers.csv').open('w',encoding='utf-8',newline='') as handle:
        writer=csv.writer(handle);writer.writerow(['target','mid','tier','balanced_score','best_pose_failed_checks'])
        for t,mid,b,l,score,payload in db.execute('SELECT target,mid,balanced,lenient,score,payload FROM results ORDER BY target,mid'):
            row=json.loads(payload);best=row['best'];failed=[]
            if best is None:
                diagnostics[t]['no_exact_pose']+=1;failed=['no_exact_pose']
            elif not b:
                a,overlap,buried=grid.PROFILES['balanced']
                for name,condition in [('severe',best['severe_fraction']>a+1e-12),
                    ('overlap',best['mean_overlap']>overlap),('burial',best['buried_fraction']<buried)]:
                    if condition:
                        failed.append(name);diagnostics[t]['best_pose_'+name+'_fail']+=1
            writer.writerow([t,mid,'balanced' if b else 'lenient_reserve' if l else 'unresolved',score,';'.join(failed)])
    base.save(out/'failure_diagnostics.json',dict(targets=diagnostics,
        scope='Best exact scored pose of non-balanced molecules only; causes overlap. Not exhaustive pose or molecular rejection.'))
    base.save(out/'report.json',dict(status='running',stage='export',targets=totals))
    exports={t:export_target(db,t,batch,prepared,out,args.cap,args.batch_size) for t in base.TARGETS}
    result=dict(status='complete',scope='existing_pool_grid_docking_handoff',targets=totals,exports=exports,
        new_scientific_rejections=0,selection='balanced supporting pose score; per-target computational quota',
        limitations=['Original-source MOL2 coordinates, not receptor-aligned or docking-prepared.',
        'Profiles allow clashes; passing is not binding-pose validation.',
        'Finite seeds and three conformers; unresolved is not impossible.',
        'Non-exported IDs remain in results and original E121 exports; no diversity or affinity claim.'])
    base.save(out/'report.json',result);db.close();print(json.dumps(result,indent=2),flush=True)
    return result


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--plan',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--workers',type=int,default=20);p.add_argument('--cap',type=int,default=1000000)
    p.add_argument('--batch-size',type=int,default=100000)
    args=p.parse_args()
    if min(args.workers,args.cap,args.batch_size)<1:p.error('Positive workers, cap and batch size required')
    run(args)


if __name__=='__main__': main()
