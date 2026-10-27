"""E121: retain the entire direct pool, separating supported and unresolved IDs."""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import sqlite3
import time

import export_kras_soft_rank as base
import kras_soft_rank as soft
import run_kras_contact_refine as prior

TARGETS = base.TARGETS


def snapshot(source, destination):
    """SQLite backup preserves a consistent committed view, including WAL data."""
    with base.readonly(source) as src, sqlite3.connect(destination) as dst:
        src.backup(dst)


def protected_output(output, inputs):
    for value in inputs:
        if value is None:
            continue
        path = Path(value).resolve()
        if output.is_relative_to(path) or path.is_relative_to(output):
            raise ValueError('Output overlaps protected input')
    if output.exists():
        raise ValueError('Use a new output directory')


def map_members(db, registry, ranking):
    """Indexed pool-to-registry resolution; no library geometry scan."""
    db.executescript('''CREATE TABLE selected(target TEXT,mid TEXT,cid TEXT,tier TEXT,
        payload TEXT,source TEXT,record_index INTEGER,source_sha TEXT,
        PRIMARY KEY(target,mid));
        CREATE INDEX selected_source ON selected(source,record_index);
        CREATE TABLE mapping_issues(target TEXT,mid TEXT,reason TEXT,PRIMARY KEY(target,mid));''')
    if ranking is not None:
        for target, mid in ranking.execute('SELECT target,mid FROM molecules'):
            if not db.execute('SELECT 1 FROM molecules WHERE target=? AND mid=?', (target, mid)).fetchone():
                raise ValueError('Ranking contains molecule outside direct pool')
    count = 0
    for target, mid in db.execute('SELECT target,mid FROM molecules ORDER BY target,mid'):
        if target not in TARGETS:
            raise ValueError('Unexpected target in direct pool')
        row = ranking.execute('SELECT payload FROM molecules WHERE target=? AND mid=?',
                              (target, mid)).fetchone() if ranking is not None else None
        payload = row[0] if row else None
        if payload:
            poses = json.loads(payload)
            if not poses or any(p['mid'] != mid for p in poses):
                raise ValueError('Ranking payload identity mismatch')
            pose = poses[0]
            record = registry.execute('SELECT id,source_path,source_record_index,content_sha256 '
                'FROM conformer WHERE molecule_id=? AND id=?', (mid, pose['cid'])).fetchone()
            tier = 'priority'
            if record is None:
                raise ValueError('Saved best conformer missing from registry')
        else:
            record = registry.execute('SELECT id,source_path,source_record_index,content_sha256 '
                'FROM conformer WHERE molecule_id=? ORDER BY conformer_index,id LIMIT 1', (mid,)).fetchone()
            tier = 'unresolved'
        if record is None:
            db.execute('INSERT INTO mapping_issues VALUES (?,?,?)', (target, mid, 'missing_source_mapping'))
        else:
            cid, source, index, digest = record
            db.execute('INSERT INTO selected VALUES (?,?,?,?,?,?,?,?)',
                       (target, mid, cid, tier, payload, source, index, digest))
        count += 1
        if count % 10000 == 0:
            db.commit()
            print(json.dumps(dict(stage='mapping',mapped_pool_members=count)), flush=True)
    db.commit()


def prepare(direct, ranking, output):
    direct = Path(direct).resolve(); output = Path(output).resolve()
    ranking = Path(ranking).resolve() if ranking else None
    report = base.read(direct/'report.json'); protocol = base.read(direct/'protocol.json')
    if (report.get('status') != 'complete' or not report.get('full_catalog') or
            report.get('scope') != 'direct_contacts_only_not_final_screening_hits'):
        raise ValueError('Completed original direct pool required')
    original = protocol['original']; batch = Path(original['batch']).resolve()
    prepared = Path(original['prepared']).resolve()
    protected_output(output, [direct, ranking, batch, prepared])
    if soft.run.seal_inputs(batch, prepared, original['chunk_size'], original['sample_size']) != original:
        raise ValueError('Original source/code seal mismatch')
    if ranking:
        rp = base.read(ranking/'protocol.json'); common = rp['common']
        if (common['original'] != original or common['code'] != prior.code_seal() or
                common['settings'] != prior.refine.SETTINGS or
                common['direct_report_sha256'] != base.sha(direct/'report.json') or
                common['counts_sha256'] != base.sha(direct/'counts.sqlite')):
            raise ValueError('E120 snapshot lineage mismatch')
    output.mkdir(parents=True)
    base.save(output/'report.json', dict(status='running',stage='mapping'))
    snapshot(direct/'counts.sqlite', output/'pool.sqlite')
    if ranking:
        snapshot(ranking/'ranking.sqlite', output/'priority_snapshot.sqlite')
        base.save(output/'priority_source_protocol.json', rp)
    base.save(output/'protocol.json', dict(original=original,direct_report=report,
        direct_protocol_sha256=base.sha(direct/'protocol.json'),
        direct_counts_sha256=base.sha(direct/'counts.sqlite'),
        registry_sha256=base.sha(batch/'registry.sqlite3'),
        snapshot_sha256=base.sha(output/'priority_snapshot.sqlite') if ranking else None,
        representation='best E120 conformer or first registered conformer',
        unresolved_policy='retain all; absence of feasible pose is not rejection'))
    with sqlite3.connect(output/'pool.sqlite') as db, base.readonly(batch/'registry.sqlite3') as registry:
        ranked = base.readonly(output/'priority_snapshot.sqlite') if ranking else None
        try:
            for target in TARGETS:
                n = db.execute('SELECT COUNT(*) FROM molecules WHERE target=?',(target,)).fetchone()[0]
                if n != report['targets'][target]['distinct_molecules']:
                    raise ValueError('Direct pool/report count mismatch')
            map_members(db, registry, ranked)
            targets = {}
            for target in TARGETS:
                total = db.execute('SELECT COUNT(*) FROM molecules WHERE target=?',(target,)).fetchone()[0]
                counts = dict(db.execute('SELECT tier,COUNT(*) FROM selected WHERE target=? GROUP BY tier',(target,)))
                issues = db.execute('SELECT COUNT(*) FROM mapping_issues WHERE target=?',(target,)).fetchone()[0]
                if sum(counts.values()) + issues != total:
                    raise ValueError('Retention accounting mismatch')
                targets[target] = dict(pool_molecules=total,provisional_tiers=counts,mapping_issues=issues)
        finally:
            if ranked: ranked.close()
    result = dict(status='complete',scope='retained_pool_mapping',targets=targets,
                  pool_sha256=base.sha(output/'pool.sqlite'),new_scientific_rejections=0)
    base.save(output/'report.json',result)
    print(json.dumps(result,indent=2),flush=True)
    return result


def write_structures(plan, output, batch_size=100000, check_factory=None):
    plan=Path(plan).resolve(); output=Path(output).resolve()
    report=base.read(plan/'report.json'); protocol=base.read(plan/'protocol.json')
    if report.get('status')!='complete' or report.get('scope')!='retained_pool_mapping':
        raise ValueError('Completed retained-pool plan required')
    if base.sha(plan/'pool.sqlite')!=report['pool_sha256']:
        raise ValueError('Mapped pool checksum mismatch')
    if batch_size<1: raise ValueError('Positive batch size required')
    batch=Path(protocol['original']['batch']); prepared=Path(protocol['original']['prepared'])
    protected_output(output,[plan,batch,prepared])
    output.mkdir(parents=True)
    base.save(output/'report.json',dict(status='running',scope='retained_pool_export'))
    base.save(output/'protocol.json',dict(plan=str(plan),plan_sha256=report['pool_sha256'],
        batch_size=batch_size,selection='all direct-pool IDs; no top-N',coordinates='original_source',
        script_sha256=base.sha(Path(__file__))))
    started=time.monotonic(); check=None
    writers={}; streams={}; csvs={}; counts={t:dict(priority=0,unresolved=0,quarantined=0,downgraded=0) for t in TARGETS}
    for target in TARGETS:
        (output/target).mkdir()
        for tier in ('priority','unresolved'):
            key=(target,tier); writers[key]=base.BatchWriter(output/target/tier,batch_size)
            stream=(output/target/tier/'molecules.csv').open('w',encoding='utf-8',newline='')
            streams[key]=stream; csvs[key]=csv.writer(stream)
            csvs[key].writerow(['mid','cid','tier','reason','batch','record','source','source_record','sha256','pose_evidence'])
    db=base.readonly(plan/'pool.sqlite')
    quarantine=(output/'quarantined.csv').open('w',encoding='utf-8',newline='')
    qc=csv.writer(quarantine); qc.writerow(['target','mid','cid','reason'])
    try:
        for target,mid,reason in db.execute('SELECT target,mid,reason FROM mapping_issues'):
            counts[target]['quarantined']+=1; qc.writerow([target,mid,'',reason])
        for source, in db.execute('SELECT DISTINCT source FROM selected ORDER BY source'):
            cursor=db.execute('SELECT target,mid,cid,tier,payload,record_index,source_sha FROM selected '
                              'WHERE source=? ORDER BY record_index,target',(source,))
            pending=next(cursor,None)
            for index,raw in base.iter_mol2_blocks(Path(source)):
                if pending is None: break
                if index<pending[5]: continue
                if index!=pending[5]: raise ValueError('Selected source record missing')
                digest=hashlib.sha256(raw.encode('utf-8')).hexdigest()
                reason=None; molecule=None
                try:
                    base.parse_mol2_block(Path(source),index,raw)
                    molecule,mode=base.load_rdkit_mol2(raw,pending[2])
                    if mode!='strict': reason='non_strict_sanitization'
                    elif not all(math.isfinite(v) for p in molecule.GetConformer().GetPositions() for v in p):
                        reason='nonfinite_coordinates'
                except Exception as exc:
                    reason='source_parse_failure:'+type(exc).__name__
                while pending is not None and pending[5]==index:
                    target,mid,cid,tier,payload,_,expected=pending
                    if digest!=expected: raise ValueError('Original source record hash mismatch')
                    if reason:
                        counts[target]['quarantined']+=1; qc.writerow([target,mid,cid,reason])
                    else:
                        explanation='finite_search_unresolved_or_unprocessed'
                        if tier=='priority':
                            if check is None: check=(check_factory or base.source_pose_checker)(batch,prepared)
                            valid,minimum=check(molecule,target,json.loads(payload)[0])
                            if valid: explanation='same_pose_contacts_and_distance_floor_rechecked'
                            else:
                                tier='unresolved'; counts[target]['downgraded']+=1
                                explanation='exact_source_pose_failed_retained_for_docking'
                        b,r=writers[target,tier].add(base.rename_record(raw,mid+'__'+cid))
                        csvs[target,tier].writerow([mid,cid,tier,explanation,b,r,source,index,digest,payload or ''])
                        counts[target][tier]+=1
                    pending=next(cursor,None)
            if pending is not None: raise ValueError('Source missing selected records')
            for stream in streams.values(): stream.flush()
            quarantine.flush()
            print(json.dumps(dict(stage='export',source_complete=source,targets=counts)),flush=True)
    finally:
        for writer in writers.values(): writer.close()
        for stream in streams.values(): stream.close()
        quarantine.close(); db.close()
    for target in TARGETS:
        row=counts[target]
        if row['priority']+row['unresolved']+row['quarantined']!=report['targets'][target]['pool_molecules']:
            raise ValueError('Final pool accounting mismatch')
        row['batches']={tier:writers[target,tier].batches for tier in ('priority','unresolved')}
    result=dict(status='complete',scope='retained_pool_docking_handoff',targets=counts,
        elapsed_seconds=time.monotonic()-started,new_scientific_rejections=0,
        limitations=['Unresolved IDs have no verified feasible pose.',
                     'Parse/mapping failures are quarantined, not biological negatives.',
                     'One representative per library ID; no chemical deduplication or docking preparation.',
                     'Existing finite-seed pool does not guarantee full-library recall.'])
    base.save(output/'report.json',result)
    print(json.dumps(result,indent=2),flush=True)
    return result


def main():
    p=argparse.ArgumentParser(description=__doc__); sub=p.add_subparsers(dest='mode',required=True)
    a=sub.add_parser('prepare'); a.add_argument('--direct',type=Path,required=True)
    a.add_argument('--ranking',type=Path); a.add_argument('--output',type=Path,required=True)
    b=sub.add_parser('export'); b.add_argument('--plan',type=Path,required=True)
    b.add_argument('--output',type=Path,required=True); b.add_argument('--batch-size',type=int,default=100000)
    args=p.parse_args()
    if args.mode=='prepare': prepare(args.direct,args.ranking,args.output)
    else: write_structures(args.plan,args.output,args.batch_size)


if __name__=='__main__': main()
