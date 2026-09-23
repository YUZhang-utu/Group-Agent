"""Time-bounded molecule sample through the actual consensus funnel kernels."""
import argparse
from collections import Counter
import copy
import json
import math
import multiprocessing as mp
from pathlib import Path
import random
import sqlite3
import time

import numpy as np

from . import library_acceptance as ev
from .expanded_wee1 import fingerprint
from .gaussian_batch import _atomic_json
from .screening_selection import check_hashes


def check_time(deadline):
    if time.perf_counter()>=deadline:raise TimeoutError('Pilot wall budget reached')


def reservoir(values, count, seed, deadline=float('inf')):
    rng=random.Random(seed);panel=[];total=0
    for total,value in enumerate(values,1):
        if total%4096==0:check_time(deadline)
        if total<=count:panel.append(value)
        else:
            j=rng.randrange(total)
            if j<count:panel[j]=value
    rng.shuffle(panel)
    return panel,total


def molecule_panel(batch, count, seed, deadline):
    catalog=ev.read(batch/'artifacts/catalog.json')
    from .full_library_screen import catalog_count
    catalog_count(catalog)
    with sqlite3.connect((batch/'registry.sqlite3').resolve().as_uri()+'?mode=ro',uri=True) as db:
        cursor=db.execute('SELECT id FROM molecule WHERE library_id=? ORDER BY id',(catalog['library_id'],))
        selected,total=reservoir((str(row[0]) for row in cursor),count,seed,deadline)
    if not selected:raise ValueError('No library molecules')
    keys=np.asarray(sorted(selected),dtype='S16')
    if any(len(s.encode())>16 for s in selected):raise ValueError('Unsupported molecule identity width')
    panel={mid:[] for mid in selected};sources={}
    for shard in catalog['shards']:
        check_time(deadline)
        directory=Path(shard['path'])
        if not directory.is_dir():directory=batch/'artifacts'/shard['name']
        manifest=directory/'manifest.json'
        if ev.sha(manifest)!=shard['manifest_sha256']:raise ValueError('Changed shard manifest')
        ev.verify_manifest(directory,full=False)
        sources.update(fingerprint([manifest]))
        mapped=np.memmap(directory/'molecule_ids.bin',dtype='S16',mode='r')
        if len(mapped)!=shard['conformers']:raise ValueError('Molecule-ID count mismatch')
        for start in range(0,len(mapped),262144):
            check_time(deadline);block=np.asarray(mapped[start:start+262144])
            positions=np.searchsorted(keys,block)
            selected_mask=(positions<len(keys)) & (keys[np.minimum(positions,len(keys)-1)]==block)
            for local in np.flatnonzero(selected_mask):
                mid=bytes(block[local]).decode()
                panel[mid].append(int(shard['global_id_start'])+start+int(local))
        del mapped
    if any(not ids for ids in panel.values()):raise ValueError('Sampled registry molecule has no catalog conformer')
    return panel,total,catalog,sources


def worker(task):
    from . import full_library_screen as full, preselection_full as pre
    queries,ids,target=task;root=Path(target);root.mkdir(parents=True,exist_ok=False)
    receipts=[];started=time.perf_counter();cpu_started=time.process_time()
    for index,q in enumerate(queries):
        full.initialize(q)
        receipt=pre.compute((0,len(ids),str(root/f'template-{index:03d}.npz'),ids))
        receipts.append(receipt)
        # Readers hold mmaps; release references before the next template.
        if full._FILTER_READER:full._FILTER_READER.close()
        full._STATE=None;full._FILTER_READER=None;pre._GUIDED_CACHE=None
    return dict(root=str(root),receipts=receipts,worker_wall_seconds=time.perf_counter()-started,
        cpu_seconds=time.process_time()-cpu_started)


def wilson(hits, count):
    if not count:return [0.,1.]
    z=1.959963984540054;p=hits/count;den=1+z*z/count
    center=(p+z*z/(2*count))/den
    half=z*math.sqrt(p*(1-p)/count+z*z/(4*count*count))/den
    return [max(0.,center-half),min(1.,center+half)]


def run(recommendation, batch, output, count=2000, workers=8, chunk_molecules=16, seconds=600, seed=20260923):
    from .consensus_design import adopt
    from .contact_policy import require_protein_contacts
    from .preselection_full import scaffold_key
    from .chemical_companion import ChemicalCompanionReader
    if min(count,workers,chunk_molecules,seconds)<=0:raise ValueError('Positive pilot settings required')
    batch=Path(batch).resolve();out=Path(output).resolve()
    if out.is_relative_to(batch) or batch.is_relative_to(out):raise ValueError('Pilot output overlaps library inputs')
    out.mkdir(parents=True,exist_ok=False)
    started=time.perf_counter();deadline=started+seconds
    report=dict(kind='consensus_molecule_pilot',status='running',full_library=False,
        requested_molecules=count,workers=workers,wall_budget_seconds=seconds,seed=seed,
        biological_validation='not_run',scope='Current direct-contact consensus; no new spatial definitions',
        integrity_scope='Sealed design and catalog/manifest checks, sampled artifact/chemical identity; no full-library byte checksum scan',
        timing_scope='Pilot preparation, repeated worker initialization, all templates, writes, union and scaffold grouping; not production full-run timing')
    _atomic_json(out/'report.json',report);pool=None
    try:
        design=adopt(recommendation,out/'adopted-design')
        if design['readiness']!='ready_for_consensus_funnel':raise ValueError('Crystal template controls require review')
        require_protein_contacts(design['design'],design['anchors'])
        chemical=ev.read(batch/'chemical/catalog.json');artifact=ev.read(batch/'artifacts/catalog.json')
        if chemical['library_id']!=artifact['library_id'] or chemical['artifact_v1_catalog_sha256']!=ev.sha(batch/'artifacts/catalog.json'):
            raise ValueError('Chemical/artifact catalog lineage mismatch')
        chemical_sources={}
        for shard in chemical['shards']:
            check_time(deadline);directory=Path(shard['path'])
            if not directory.is_dir():directory=batch/'chemical'/shard['name']
            if ev.sha(directory/'manifest.json')!=shard['manifest_sha256']:raise ValueError('Changed chemical shard manifest')
            ev.verify_manifest(directory,full=False)
            chemical_sources.update(fingerprint([directory/'manifest.json']))
        panel,total,catalog,manifests=molecule_panel(batch,count,seed,deadline)
        _atomic_json(out/'sample.json',dict(seed=seed,molecules=panel))
        report.update(sampled_molecules=len(panel),sampled_conformers=sum(map(len,panel.values())),
            library_molecules=total,library_conformers=catalog['conformers'],templates=len(design['templates']),
            preparation_seconds=time.perf_counter()-started,design=design['design'],
            sources={**design['sources'],**manifests,**chemical_sources,**fingerprint([batch/'artifacts/catalog.json',batch/'chemical/catalog.json',out/'sample.json'])},
            implementation_sources=fingerprint(sorted(Path(__file__).parent.glob('*.py'))))
        queries=[]
        for template in design['templates']:
            q=dict(template,anchors=copy.deepcopy(design['anchors']),consensus_npz=design['consensus_npz'],
                artifact_catalog=str(batch/'artifacts/catalog.json'),chemical_companion=str(batch/'chemical/catalog.json'),
                guided_design=design['design'],condition_policy=dict(required_anchors=design['anchor_order'],match_mode='any',
                    minimum_score=design['design']['minimum_score'],coarse_constraints=design['design']['coarse_constraints']))
            indices=sorted({a['feature_index'] for a in q['anchors']})
            for a in q['anchors']:a['score_column']=indices.index(a['feature_index'])
            queries.append(q)
        items=list(panel.items());tasks=[]
        for start in range(0,len(items),chunk_molecules):
            ids=sorted(gid for _,gids in items[start:start+chunk_molecules] for gid in gids)
            tasks.append((queries,ids,str(out/'chunks'/f'{start:08d}')))
        ctx=mp.get_context('spawn');pool=ctx.Pool(workers)
        pending=[];next_task=0;completed=[];scan_start=time.perf_counter();timings=Counter();counts=Counter()
        while next_task<len(tasks) or pending:
            check_time(deadline)
            while next_task<len(tasks) and len(pending)<workers*2:
                pending.append(pool.apply_async(worker,(tasks[next_task],)));next_task+=1
            ready=[job for job in pending if job.ready()]
            if not ready:time.sleep(.1);continue
            for job in ready:
                result=job.get();pending.remove(job);completed.append(result)
                for receipt in result['receipts']:
                    counts.update(receipt['counts']);timings.update(receipt['worker_seconds'])
                report.update(completed_chunks=len(completed),total_chunks=len(tasks),
                    counts_summed_over_templates=dict(counts),summed_worker_seconds=dict(timings),
                    elapsed_seconds=time.perf_counter()-started)
                _atomic_json(out/'report.json',report)
                print(f"Completed {len(completed)}/{len(tasks)} chunks; {report['elapsed_seconds']:.1f}s",flush=True)
        pool.close();pool.join();pool=None
        report['scan_wall_seconds']=time.perf_counter()-scan_start
        report['summed_worker_cpu_seconds']=sum(r['cpu_seconds'] for r in completed)
        merge_start=time.perf_counter();stages={mid:0 for mid in panel}
        dbpath=out/'candidates.sqlite'
        with sqlite3.connect(dbpath) as db:
            db.execute('CREATE TABLE poses(mid TEXT,template TEXT,mask INTEGER,spatial TEXT,score REAL,gid INTEGER,payload TEXT,PRIMARY KEY(mid,template,mask,spatial))')
            for result in completed:
                check_time(deadline)
                for index,receipt in enumerate(result['receipts']):
                    check_hashes(receipt['files']);path=Path(result['root'])/f'template-{index:03d}.npz'
                    with np.load(path,allow_pickle=False) as arrays:
                        for mid,level in zip(arrays['molecule_ids'],arrays['stage_levels']):
                            stages[str(mid)]=max(stages[str(mid)],int(level))
                    with path.with_suffix('.poses.jsonl').open() as stream:
                        for line in stream:
                            row=json.loads(line);row['template_id']=queries[index]['query_id']
                            spatial=json.dumps(sorted(row.get('occupied_spatial_groups',[])))
                            db.execute('INSERT INTO poses VALUES(?,?,?,?,?,?,?) ON CONFLICT(mid,template,mask,spatial) DO UPDATE SET score=excluded.score,gid=excluded.gid,payload=excluded.payload WHERE excluded.score>poses.score',
                                (row['molecule_id'],row['template_id'],row['mask'],spatial,row['composite_score'],row['global_id'],json.dumps(row)))
            db.execute('CREATE TABLE members(mid TEXT PRIMARY KEY,cluster TEXT,error TEXT)')
            reader=ChemicalCompanionReader(batch/'chemical/catalog.json')
            for mid,gid in db.execute('SELECT mid,MIN(gid) FROM poses GROUP BY mid').fetchall():
                check_time(deadline);chem=reader.get(gid)
                if chem.molecule_id!=mid:raise ValueError('Scaffold molecule identity mismatch')
                try:key,_=scaffold_key(chem);error=''
                except Exception as exc:key='unresolved-'+mid;error=type(exc).__name__
                db.execute('INSERT INTO members VALUES(?,?,?)',(mid,key,error))
            hits=db.execute('SELECT COUNT(*) FROM members').fetchone()[0]
            pose_count=db.execute('SELECT COUNT(*) FROM poses').fetchone()[0]
            clusters=db.execute('SELECT COUNT(DISTINCT cluster) FROM members').fetchone()[0]
        check_hashes(report['sources']);check_hashes(report['implementation_sources'])
        elapsed=time.perf_counter()-started;n=len(panel)
        report.update(status='complete',sample_complete=True,matching_molecules=hits,pose_records=pose_count,scaffold_groups=clusters,
            stages=[dict(stage=name,molecules=sum(v>=i for v in stages.values())) for i,name in enumerate(
                ['size','feature_count','extent','anchor_bound','pose_feasibility','gaussian','final_pose'],1)],
            merge_and_group_seconds=time.perf_counter()-merge_start,wall_seconds=elapsed,
            amortized_wall_seconds_per_molecule=elapsed/n,
            projected_matching_molecules=total*hits/n,
            approximate_95_percent_matching_interval=[total*x for x in wilson(hits,n)] if n<total else [hits,hits],
            projected_worker_cpu_hours=report['summed_worker_cpu_seconds']*total/n/3600,
            projected_processing_wall_hours=(elapsed-report['preparation_seconds'])*total/n/3600,
            projection_scope='Complete uniform molecule sample; approximate interval; same hardware/workers/templates. Repeated pilot initialization differs from production. Full-library byte integrity, active controls and nonlinear I/O/merge scaling excluded.',
            outputs=dict(candidates=str(dbpath),sample=str(out/'sample.json')),output_hashes=fingerprint([dbpath,out/'sample.json']))
    except TimeoutError as exc:
        report.update(status='budget_exhausted',sample_complete=False,error=str(exc),elapsed_seconds=time.perf_counter()-started,
            projection_scope='No population extrapolation: unfinished chunks may be systematically slower')
    except Exception as exc:
        report.update(status='failed',sample_complete=False,error=str(exc));raise
    finally:
        if pool is not None:pool.terminate();pool.join()
        _atomic_json(out/'report.json',report)
    return report


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--recommendation',required=True);p.add_argument('--output',required=True)
    location=p.add_mutually_exclusive_group(required=True)
    location.add_argument('--batch');location.add_argument('--runtime',help='Runtime JSON containing search.batch; no credentials printed')
    p.add_argument('--molecules',type=int,default=2000);p.add_argument('--workers',type=int,default=8)
    p.add_argument('--chunk-molecules',type=int,default=16);p.add_argument('--seconds',type=int,default=600);p.add_argument('--seed',type=int,default=20260923)
    a=p.parse_args()
    batch=a.batch or ev.read(a.runtime).get('search',{}).get('batch')
    if not batch:p.error('Runtime has no search.batch')
    r=run(a.recommendation,batch,a.output,a.molecules,a.workers,a.chunk_molecules,a.seconds,a.seed)
    print(json.dumps({k:r[k] for k in ('status','sample_complete','wall_seconds','matching_molecules','projected_matching_molecules','approximate_95_percent_matching_interval','projected_processing_wall_hours') if k in r},indent=2))


if __name__=='__main__':main()
