"""Run two native-state KRAS screens with resumable chunks and explicit coverage."""
import argparse
import csv
from concurrent.futures import ProcessPoolExecutor, wait, FIRST_COMPLETED
import gzip
import json
import os
from pathlib import Path
import sys
import sqlite3
import time
for variable in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[variable]='1'
import numpy as np

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from aidd_agent import kras_weekend as kw
from aidd_agent.gaussian_batch import ArtifactCatalogReader
from aidd_agent.interaction_fast import InteractionFeatureReader
from aidd_agent.chemical_companion import ChemicalCompanionReader
from aidd_agent.full_library_screen import catalog_count

WORK = None


def initialize(batch, prepared):
    global WORK
    from aidd_agent.expanded_wee1 import ensure_file_descriptor_limit
    batch=Path(batch)
    if os.name=='posix':
        ensure_file_descriptor_limit(len(kw.read(batch/'artifacts/catalog.json')['shards']))
    WORK=(ArtifactCatalogReader(batch/'artifacts/catalog.json'),
          InteractionFeatureReader(batch/'artifacts/catalog.json',batch/'chemical/catalog.json'),
          {target:kw.load_state(prepared,target) for target in kw.REFERENCES},
          ChemicalCompanionReader(batch/'chemical/catalog.json'))


def process(task):
    index, ids = task
    reader, feature_reader, states, chemical_reader=WORK
    started=time.monotonic()
    outputs={target:dict(target=target,index=index,processed=len(ids),
        counts={k:0 for k in ('chemical','contact','shape','color','geometry','full')},
        seeds=0,hits=[],diagnostic_near_misses=[],direct_region_pose_counts={}) for target in states}
    for gid in ids:
        candidate=reader.get(int(gid)); features=feature_reader.get(int(gid))
        if candidate.conformer_id!=features.conformer_id or candidate.molecule_id!=features.molecule_id:
            raise ValueError('Artifact/chemical identity mismatch')
        chemistry=chemical_reader.get(int(gid))
        if chemistry.conformer_id!=candidate.conformer_id or chemistry.molecule_id!=candidate.molecule_id:
            raise ValueError('Chemical topology identity mismatch')
        sets=kw.chemical_sets(chemistry,candidate.feature_types)
        for target,state in states.items():
            result=kw.evaluate(candidate,features,state,chemistry,sets); row=outputs[target]
            row['seeds']+=result['seeds']
            for key in row['counts']: row['counts'][key]+=int(result['flags'][key])
            for key,value in result['direct_region_pose_counts'].items():
                row['direct_region_pose_counts'][key]=row['direct_region_pose_counts'].get(key,0)+value
            if result['flags']['full']: row['hits'].append(result)
            elif result['best_contact_shape'] is not None:
                row['diagnostic_near_misses'].append(result)
                row['diagnostic_near_misses']=sorted(row['diagnostic_near_misses'],
                    key=lambda r:r['best_contact_shape']['shape'],reverse=True)[:5]
    return dict(index=index,processed=len(ids),ids_sha256=kw.hashlib.sha256(
        np.asarray(ids,dtype='<i8').tobytes()).hexdigest(),targets=outputs,
        worker_seconds=time.monotonic()-started)


def tasks(total, chunk_size, sample_size):
    ids=np.unique(np.linspace(0,total-1,min(total,sample_size),dtype=np.int64)) if sample_size else None
    length=len(ids) if ids is not None else total
    for index,start in enumerate(range(0,length,chunk_size)):
        values=ids[start:start+chunk_size].tolist() if ids is not None else list(range(start,min(total,start+chunk_size)))
        yield index,values


def seal_inputs(batch, prepared, chunk_size, sample_size):
    report=kw.read(prepared/'report.json')
    for relative,digest in report['output_hashes'].items():
        path=(prepared/relative).resolve()
        if not path.is_relative_to(prepared.resolve()) or kw.sha(path)!=digest:
            raise ValueError('Reference package checksum mismatch: '+relative)
    catalogs={name:kw.sha(batch/name/'catalog.json') for name in ('artifacts','chemical')}
    companion=kw.read(batch/'chemical/catalog.json')
    if companion['artifact_v1_catalog_sha256']!=catalogs['artifacts']:
        raise ValueError('Chemical catalog is not paired to artifact catalog')
    # Content-hash shard manifests and record immutable payload metadata without
    # rereading the entire multi-terabyte library for every restart.
    shards=[]
    for name in ('artifacts','chemical'):
        for row in kw.read(batch/name/'catalog.json')['shards']:
            folder=Path(row['path'])
            if not folder.is_dir(): folder=batch/name/row['name']
            payloads={p.name:[p.stat().st_size,p.stat().st_mtime_ns] for p in sorted(folder.glob('*.bin'))}
            shards.append(dict(kind=name,path=str(folder.resolve()),manifest=kw.sha(folder/'manifest.json'),payloads=payloads))
    source=Path(kw.__file__).parent
    return dict(catalogs=catalogs,shards=shards,reference_report=kw.sha(prepared/'report.json'),
        code={p.name:kw.hashlib.sha256(p.read_bytes().replace(b'\r\n',b'\n')).hexdigest()
              for p in sorted(source.glob('*.py'))},runner_sha256=kw.sha(__file__),
        chunk_size=chunk_size,sample_size=sample_size,policy=kw.POLICY,
        numpy=np.__version__,batch=str(batch.resolve()),prepared=str(prepared.resolve()))


def execute(args):
    batch=args.batch
    if batch is None:
        batch=Path(kw.read(args.library_receipt)['batch'])
    batch=batch.resolve(); prepared=args.prepared.resolve(); output=args.output.resolve()
    if any(output.is_relative_to(p) or p.is_relative_to(output) for p in (batch,prepared)):
        raise ValueError('Output overlaps protected input directories')
    output.mkdir(parents=True,exist_ok=True)
    from aidd_agent.prompt_workflow import file_lock
    with file_lock(output/'run.lock'):
        signature=seal_inputs(batch,prepared,args.chunk_size,args.sample)
        pilot_report=getattr(args,'pilot_report',None)
        if pilot_report is not None:
            pilot=kw.read(pilot_report); pilot_protocol=kw.read(pilot_report.parent/'protocol.json')
            if pilot.get('status')!='complete' or pilot.get('scope')!='spread_pilot':
                raise ValueError('A complete pilot is required')
            for key in signature:
                if key not in ('chunk_size','sample_size') and pilot_protocol.get(key)!=signature[key]:
                    raise ValueError('Pilot/full lineage mismatch: '+key)
        protocol=output/'protocol.json'
        if protocol.exists():
            if kw.read(protocol)!=signature:
                raise ValueError('Resume input/code/settings mismatch; use a new output directory')
        else:
            kw.save(protocol,signature)
        kw.save(output/'report.json',dict(status='running',scope='spread_pilot' if args.sample else 'full_catalog'))
        controls={target:kw.self_control(kw.load_state(prepared,target)) for target in kw.REFERENCES}
        kw.save(output/'controls.json',controls)
        total=catalog_count(kw.read(batch/'artifacts/catalog.json'))
        target_total=min(total,args.sample) if args.sample else total
        chunks=output/'chunks'; chunks.mkdir(exist_ok=True)
        summary={target:dict(processed=0,total_library_conformers=total,
            scheduled_conformers=target_total,full=0,chemical=0,contact=0,shape=0,color=0,geometry=0,seeds=0)
            for target in kw.REFERENCES}
        started=time.monotonic(); deadline=started+args.hours*3600
        newly_processed=0
        database=sqlite3.connect(output/'candidates.sqlite')
        database.execute('CREATE TABLE IF NOT EXISTS candidates(target TEXT, mid TEXT, gid INTEGER, score REAL, payload TEXT, PRIMARY KEY(target,mid))')
        database.execute('DELETE FROM candidates')
        database.commit()
        def incorporate(result):
            for target,value in result['targets'].items():
                row=summary[target]
                row['processed']+=value['processed']; row['seeds']+=value['seeds']
                for key,count in value['counts'].items(): row[key]+=count
                for hit in value['hits']:
                    score=hit['passing_pose']['shape']+hit['passing_pose']['color']
                    database.execute('INSERT INTO candidates VALUES(?,?,?,?,?) ON CONFLICT(target,mid) DO UPDATE SET gid=excluded.gid,score=excluded.score,payload=excluded.payload WHERE excluded.score>candidates.score OR (excluded.score=candidates.score AND excluded.gid<candidates.gid)',
                        (target,hit['mid'],hit['gid'],score,json.dumps(hit,allow_nan=False)))
            database.commit()
        def pending_tasks():
            for task in tasks(total,args.chunk_size,args.sample):
                index,ids=task
                path=chunks/f'paired-{index:08d}.json.gz'
                receipt=path.with_suffix('.receipt.json')
                if receipt.exists():
                    check=kw.read(receipt)
                    if kw.sha(path)!=check['sha256']:
                        raise ValueError('Corrupt completed chunk: '+str(path))
                    with gzip.open(path,'rt',encoding='utf-8') as stream: result=json.load(stream)
                    if set(result['targets'])!=set(kw.REFERENCES) or result['index']!=index or result['processed']!=len(ids) or result['ids_sha256']!=kw.hashlib.sha256(np.asarray(ids,dtype='<i8').tobytes()).hexdigest():
                        raise ValueError('Chunk schedule mismatch')
                    incorporate(result)
                else:
                    yield task,path,receipt
        iterator=iter(pending_tasks()); running={}; exhausted=False
        with ProcessPoolExecutor(max_workers=args.workers,initializer=initialize,initargs=(batch,prepared)) as pool:
            while True:
                while not exhausted and len(running)<args.workers and time.monotonic()<deadline:
                    item=next(iterator,None)
                    if item is None: exhausted=True; break
                    task,path,receipt=item
                    running[pool.submit(process,task)]=(path,receipt)
                if not running: break
                done,_=wait(running,timeout=10,return_when=FIRST_COMPLETED)
                for future in done:
                    path,receipt=running.pop(future); result=future.result()
                    temporary=path.with_suffix('.partial')
                    with gzip.open(temporary,'wt',encoding='utf-8') as stream: json.dump(result,stream,allow_nan=False)
                    temporary.replace(path); kw.save(receipt,dict(sha256=kw.sha(path)))
                    incorporate(result); newly_processed+=result['processed']
                    elapsed=time.monotonic()-started
                    progress=dict(status='running',targets=summary,invocation_seconds=elapsed,
                        invocation_new_conformers=newly_processed,workers=args.workers,
                        dual_full_hours_linear_estimate=(total/newly_processed*elapsed/3600 if newly_processed else None))
                    kw.save(output/'progress.json',progress)
                    print(json.dumps(progress),flush=True)
        complete=all(row['processed']==target_total for row in summary.values())
        for target in summary:
            summary[target]['distinct_passing_molecules']=database.execute(
                'SELECT count(*) FROM candidates WHERE target=?',(target,)).fetchone()[0]
            with (output/(target+'_candidates.csv')).open('w',newline='',encoding='utf-8') as stream:
                writer=csv.writer(stream); writer.writerow(['molecule_id','global_id','conformer_id','shape','color','reference','transform'])
                for mid,gid,payload in database.execute('SELECT mid,gid,payload FROM candidates WHERE target=? ORDER BY score DESC,mid',(target,)):
                    hit=json.loads(payload); pose=hit['passing_pose']
                    writer.writerow([mid,gid,hit['cid'],pose['shape'],pose['color'],
                                     kw.REFERENCES[target][0],json.dumps(pose['transform'])])
        database.close()
        elapsed=time.monotonic()-started
        result=dict(status='complete' if complete else 'partial',
            scope='spread_pilot' if args.sample else 'full_catalog',targets=summary,
            full_library_complete=complete and not args.sample,workers=args.workers,
            invocation_seconds=elapsed,invocation_new_conformers=newly_processed,
            dual_full_hours_linear_estimate=(total/newly_processed*elapsed/3600 if newly_processed else None),
            limitations=['Finite bounded rigid poses, not all conformations/orientations.',
                'Exploratory KRAS gates, no activity enrichment or allele selectivity validation.',
                'Color and shape counters are independent after contacts; geometry includes shape; full is same-pose all gates.',
                'Each hit is retained with its best passing pose; diagnostic near misses capped at five per chunk.',
                'Pilot timing is a linear estimate, not a Monday completion guarantee.',
                'Library payload identity uses metadata plus manifest hashes, not a fresh full payload hash scan.'])
        kw.save(output/'report.json',result)
        print(json.dumps(result,indent=2))
        return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--batch',type=Path)
    parser.add_argument('--library-receipt',type=Path,default=Path.home()/'mdm2-e109/prepared/report.json')
    parser.add_argument('--prepared',type=Path,default=Path(__file__).resolve().parents[1]/'configs/targets/kras_sip_e116')
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--sample',type=int,default=0,help='0 schedules every conformer; otherwise spread sample per target')
    parser.add_argument('--workers',type=int,default=20)
    parser.add_argument('--chunk-size',type=int,default=128)
    parser.add_argument('--hours',type=float,default=60)
    parser.add_argument('--pilot-report',type=Path,help='Required for CLI full scan; verifies matching pilot lineage')
    args=parser.parse_args()
    if args.workers<1 or args.chunk_size<1 or args.sample<0 or not np.isfinite(args.hours) or args.hours<=0:
        parser.error('Invalid resource limits')
    if args.sample==0 and args.pilot_report is None:
        parser.error('Full scanning requires --pilot-report from this code, library and reference package')
    try:
        execute(args)
    except Exception as exc:
        report=args.output/'report.json'
        if report.is_file() and kw.read(report).get('status')=='running':
            kw.save(report,dict(status='failed',error_type=type(exc).__name__,error=str(exc)))
        raise


if __name__=='__main__': main()
