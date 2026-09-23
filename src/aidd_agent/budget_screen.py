"""Resumable molecule-budget ANN retrieval, consensus ranking and paged handoff."""
import argparse
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
import copy
import csv
import hashlib
import json
import multiprocessing as mp
import os
from pathlib import Path
import sqlite3
import time

import numpy as np

from . import library_acceptance as ev
from .expanded_wee1 import fingerprint, ensure_file_descriptor_limit
from .gaussian_batch import _atomic_json
from .screening_selection import check_hashes
from .prompt_workflow import file_lock


def readonly(path):
    return sqlite3.connect(Path(path).resolve().as_uri()+'?mode=ro',uri=True)


def unique_neighbors(ids, mids, distances, limit):
    """Best retrieved conformer per molecule, with deterministic ties."""
    best={}
    for gid,mid,d in zip(ids,mids,distances):
        mid=bytes(mid).decode() if isinstance(mid,(bytes,np.bytes_)) else str(mid)
        value=(float(d),int(gid))
        if not np.isfinite(d):raise ValueError('Nonfinite retrieval distance')
        if mid not in best or value<best[mid]:best[mid]=value
    return [(mid,*v) for mid,v in sorted(best.items(),key=lambda x:(x[1][0],x[0]))[:limit]]


def fuse(lists, limit=None, k=60):
    scores=defaultdict(float)
    for rows in lists:
        seen=set()
        for rank,mid in enumerate(rows,1):
            if mid in seen:raise ValueError('Duplicate molecule within template rank')
            seen.add(mid);scores[mid]+=1/(k+rank)
    rows=sorted(scores,key=lambda mid:(-scores[mid],mid))
    return rows[:limit] if limit is not None else rows


def query_vector(template):
    from .chemistry_prep import _ccd_molecule, enumerate_ligand_instances, standardize_parent
    from .similarity import usrcat_descriptor
    directory=Path(template['query_npz']).parent
    manifest=ev.read(directory/'native.manifest.json')
    source=manifest['source'];qid=template['query_id']
    if source['query_id']!=qid:raise ValueError('Native query identity mismatch')
    paths={key:source[key] for key in ('mmcif','ccd','query_manifest')}
    check_hashes({path:source[key+'_sha256'] for key,path in paths.items()})
    _,ccd,chain,res=qid.split(':')
    instances=[r for r in enumerate_ligand_instances(Path(paths['mmcif']),[ccd])
               if (r['chain_id'],r['residue_number'])==(chain,res)]
    if len(instances)!=1:raise ValueError('Ambiguous crystal query')
    parent,_=standardize_parent(_ccd_molecule(Path(paths['ccd']),instances[0]['atoms']))
    vector=np.asarray(usrcat_descriptor(parent),dtype=np.float32)
    if vector.shape!=(60,) or not np.isfinite(vector).all():raise ValueError('Invalid query descriptor')
    return vector


def retrieve(batch, out, design, budget, nprobe, threads):
    import faiss
    catpath=batch/'artifacts/catalog.json';catalog=ev.read(catpath)
    fm=ev.read(batch/'faiss/manifest.json')
    if fm['catalog_sha256']!=ev.sha(catpath):raise ValueError('Index/catalog lineage mismatch')
    if ev.sha(batch/'faiss/index.faiss')!=fm['index_sha256']:raise ValueError('Index checksum mismatch')
    for row in catalog['shards']:row['path']=str(ev.shard_path(row,catpath))
    corpus=ev.Corpus(catalog)
    faiss.omp_set_num_threads(threads)
    index=faiss.read_index(str(batch/'faiss/index.faiss'))
    if index.ntotal!=corpus.total or index.d!=60:raise ValueError('Wrong FAISS count/dimension')
    with np.load(batch/'faiss/transform.npz',allow_pickle=False) as z:
        mean,std=z['mean'],z['std']
    if mean.shape!=(60,) or std.shape!=(60,) or not np.isfinite(mean).all() or not np.isfinite(std).all() or np.any(std<=0):
        raise ValueError('Invalid frozen descriptor transform')
    with readonly(batch/'registry.sqlite3') as registry:
        total=registry.execute('SELECT COUNT(*) FROM molecule WHERE library_id=?',(catalog['library_id'],)).fetchone()[0]
    budget=min(budget,total);lists=[];depths=[]
    for ti,template in enumerate(design['templates']):
        q=np.ascontiguousarray((query_vector(template)-mean)/std,dtype=np.float32)
        depth=min(corpus.total,max(1024,budget*3));probe=min(nprobe,index.nlist)
        while True:
            index.nprobe=probe
            _,found=index.search(q[None],depth)
            ids=found[0][found[0]>=0]
            if len(np.unique(ids))!=len(ids):raise ValueError('Duplicate FAISS IDs')
            mids,raw=corpus.fetch(ids,vectors=True)
            distances=np.sum(((raw-mean)/std-q)**2,axis=1)
            rows=unique_neighbors(ids,mids,distances,budget)
            if len(rows)>=budget:break
            if probe<index.nlist:probe=min(index.nlist,probe*2)
            elif depth<corpus.total:depth=min(corpus.total,depth*2)
            else:break
        lists.append([r[0] for r in rows]);depths.append(dict(template=template['query_id'],molecules=len(rows),conformer_depth=depth,nprobe=probe))
        np.savez(out/f'retrieval-{ti:02d}.npz',molecule_ids=np.asarray([r[0] for r in rows]),
                 distances=np.asarray([r[1] for r in rows]),global_ids=np.asarray([r[2] for r in rows]))
        print(f"Retrieved {len(rows)} unique molecules: {template['query_id']}",flush=True)
    selected=fuse(lists,budget)
    np.save(out/'selected-molecules.npy',np.asarray(sorted(selected),dtype='S16'))
    return dict(mode='USRCAT_FAISS_ANN',library_molecules=total,library_conformers=corpus.total,
        selected_molecules=len(selected),templates=depths,recall='not_measured_on_current_target',
        scope='All indexed library entries eligible for ANN retrieval; not exhaustive pose scoring')


def collect_conformers(batch, keys):
    """Expand selected molecules to ALL stored conformers across shard boundaries."""
    if not len(keys) or not np.array_equal(keys,np.unique(keys)):raise ValueError('Sorted unique molecule IDs required')
    catalog=ev.read(batch/'artifacts/catalog.json');parts=[];seen=np.zeros(len(keys),bool)
    for row in catalog['shards']:
        directory=ev.shard_path(row,batch/'artifacts/catalog.json')
        with ev.mapped(directory/'molecule_ids.bin','S16') as mids:
            for start in range(0,len(mids),262144):
                block=mids[start:start+262144];pos=np.searchsorted(keys,block)
                keep=(pos<len(keys)) & (keys[np.minimum(pos,len(keys)-1)]==block)
                seen[pos[keep]]=True
                parts.append(np.flatnonzero(keep)+start+int(row['global_id_start']))
    ids=np.concatenate(parts).astype(np.int64)
    if not len(ids) or not seen.all():raise ValueError('Selected molecules missing stored conformers')
    return ids


def make_queries(batch, design):
    queries=[]
    for template in design['templates']:
        q=dict(template,anchors=copy.deepcopy(design['anchors']),consensus_npz=design['consensus_npz'],
            artifact_catalog=str(batch/'artifacts/catalog.json'),chemical_companion=str(batch/'chemical/catalog.json'),
            guided_design=copy.deepcopy(design['design']),selection_mode='budget',
            condition_policy=dict(required_anchors=design['anchor_order'],match_mode='any',minimum_score=.5,
                coarse_constraints=design['design']['coarse_constraints']))
        if q['guided_design'].get('exclusions'):raise ValueError('Budget mode requires an explicit exclusion-free design')
        indices=sorted({a['feature_index'] for a in q['anchors']})
        for a in q['anchors']:a['score_column']=indices.index(a['feature_index'])
        queries.append(q)
    return queries


def init_worker(queries):
    global _QUERIES, _CURRENT
    _QUERIES=queries;_CURRENT=None


def refine_chunk(task):
    global _CURRENT
    from . import full_library_screen as full, preselection_full as pre
    ti,ids,target=task
    if _CURRENT!=ti:
        if full._FILTER_READER:full._FILTER_READER.close()
        pre._GUIDED_CACHE=None;full.initialize(_QUERIES[ti]);_CURRENT=ti
    return pre.compute((0,len(ids),target,ids))


def verified_chunk(path, ids):
    receipt=path.with_suffix('.receipt.json')
    if not receipt.is_file():return False
    record=ev.read(receipt);check_hashes(record['files'])
    with np.load(path,allow_pickle=False) as z:
        if not np.array_equal(z['global_ids'],ids):raise ValueError('Changed chunk schedule')
    return True


def run_refinement(batch,out,design,ids,workers,chunk):
    queries=make_queries(batch,design);ctx=mp.get_context('spawn');done=0
    total=len(queries)*((len(ids)+chunk-1)//chunk);started=time.perf_counter()
    with ProcessPoolExecutor(workers,mp_context=ctx,initializer=init_worker,initargs=(queries,)) as pool:
        # Template-major scheduling amortizes reader initialization.
        for ti in range(len(queries)):
            pending=[];root=out/'chunks'/f'{ti:02d}';root.mkdir(parents=True,exist_ok=True)
            for start in range(0,len(ids),chunk):
                subset=ids[start:start+chunk];target=root/f'{start:010d}.npz'
                if verified_chunk(target,subset):done+=1;continue
                pending.append(pool.submit(refine_chunk,(ti,subset,str(target))))
                if len(pending)>=workers*2:
                    pending.pop(0).result();done+=1
                    print(f'Refined {done}/{total} chunks; {time.perf_counter()-started:.1f}s',flush=True)
            for job in pending:
                job.result();done+=1
                print(f'Refined {done}/{total} chunks; {time.perf_counter()-started:.1f}s',flush=True)
    return dict(chunks=total,wall_seconds=time.perf_counter()-started)


def merge_and_rank(out, design, definitions, batch, quota, k):
    from .gaussian_batch import ArtifactCatalogReader
    from .spatial_consistency import compare
    temporary=out/'ranking.partial.sqlite'
    # Only this explicitly named derived temporary database is replaced on retry.
    if temporary.exists():temporary.unlink()
    reader=ArtifactCatalogReader(batch/'artifacts/catalog.json')
    with sqlite3.connect(temporary) as db:
        db.execute('CREATE TABLE poses(mid TEXT,template TEXT,score REAL,gid INTEGER,payload TEXT,contact REAL,gaussian REAL,PRIMARY KEY(mid,template))')
        for ti,template in enumerate(design['templates']):
            for path in sorted((out/'chunks'/f'{ti:02d}').glob('*.poses.jsonl')):
                check_hashes(ev.read(path.with_name(path.name.replace('.poses.jsonl','.receipt.json')))['files'])
                with path.open() as stream:
                    for line in stream:
                        row=json.loads(line)
                        if not np.isfinite(row['composite_score']):raise ValueError('Nonfinite pose score')
                        if not np.isfinite([row['optional_score'],row['gaussian_same_pose']]).all():raise ValueError('Nonfinite ranking terms')
                        db.execute('INSERT INTO poses VALUES(?,?,?,?,?,?,?) ON CONFLICT(mid,template) DO UPDATE SET score=excluded.score,gid=excluded.gid,payload=excluded.payload,contact=excluded.contact,gaussian=excluded.gaussian WHERE (excluded.contact,excluded.gaussian,-excluded.gid)>(poses.contact,poses.gaussian,-poses.gid)',
                            (row['molecule_id'],template['query_id'],row['composite_score'],row['global_id'],line,row['optional_score'],row['gaussian_same_pose']))
            db.commit()
        db.execute('CREATE TABLE template_ranks(mid TEXT,template TEXT,rank INTEGER,PRIMARY KEY(mid,template))')
        db.execute('INSERT INTO template_ranks SELECT mid,template,ROW_NUMBER() OVER(PARTITION BY template ORDER BY contact DESC,gaussian DESC,mid) FROM poses')
        db.execute('CREATE TABLE ranking(rank INTEGER PRIMARY KEY,mid TEXT UNIQUE,rrf REAL,template_support INTEGER)')
        # Preserve the entire scored pool. A quota-priority tier precedes the
        # reserve tier, whose full RRF order remains available for later pages.
        db.execute('''INSERT INTO ranking SELECT ROW_NUMBER() OVER(ORDER BY has_contact DESC,tier DESC,rrf DESC,mid),mid,rrf,support FROM
            (SELECT t.mid,MAX(p.contact>0) has_contact,MAX(t.rank<=?) tier,SUM(1.0/(?+t.rank)) rrf,COUNT(*) support
             FROM template_ranks t JOIN poses p ON p.mid=t.mid AND p.template=t.template GROUP BY t.mid)''',(quota,k))
        db.commit()
        for mid,template,gid,payload in db.execute('SELECT mid,template,gid,payload FROM poses'):
            row=json.loads(payload);candidate=reader.get(gid)
            if candidate.molecule_id!=mid:raise ValueError('Pose identity mismatch')
            matrix=np.asarray(row['transform']).reshape(4,4)
            moved=candidate.shape_points@matrix[:3,:3].T+matrix[:3,3]
            row['region_occupancy']=(compare(moved,definitions['definitions'],definitions['ambiguity_margin'])
                                     if definitions else dict(status='unknown',reason='No reviewed region definitions for this target/frame'))
            row['template_id']=template
            db.execute('UPDATE poses SET payload=? WHERE mid=? AND template=?',(json.dumps(row),mid,template))
        db.commit()
        count=db.execute('SELECT COUNT(*) FROM ranking').fetchone()[0]
    db.close()
    temporary.replace(out/'ranking.sqlite')
    return count


def validate_regions(design, definitions):
    if definitions is None:return
    ref=design['reference'];frame=ref.get('coordinate_frame',ref.get('query_id')) if isinstance(ref,dict) else ref
    if definitions['coordinate_frame']!=frame or definitions['target']!=design['target']['accession']:
        raise ValueError('Occupancy frame/target mismatch')


def execute(args):
    from .consensus_design import adopt
    batch=args.batch.resolve();out=args.output.resolve()
    if out.is_relative_to(batch) or batch.is_relative_to(out):raise ValueError('Output overlaps library')
    out.mkdir(parents=True,exist_ok=True)
    with file_lock(out/'run.lock'):
        sourcepaths=[args.recommendation,batch/'artifacts/catalog.json',batch/'chemical/catalog.json',
                     batch/'faiss/manifest.json',batch/'faiss/transform.npz']
        if args.definitions:sourcepaths.append(args.definitions)
        protocol=dict(kind='molecule_budget_screen',ranking_policy='contact_first_v1',budget=args.retrieval_molecules,workers=args.workers,
            chunk_conformers=args.chunk_conformers,nprobe=args.nprobe,template_quota=args.template_quota,rrf_k=args.rrf_k,
            assignment_backend=os.environ.get('AIDD_ASSIGNMENT_BACKEND','python'),
            sources=fingerprint(sourcepaths),code=fingerprint(sorted(Path(__file__).parent.glob('*.py'))))
        if (out/'protocol.json').exists():
            if ev.read(out/'protocol.json')!=protocol:raise ValueError('Changed protocol; use a fresh output directory')
        else:_atomic_json(out/'protocol.json',protocol)
        try:
            _atomic_json(out/'status.json',dict(status='running'))
            if (out/'adopted-design/report.json').exists():design=ev.read(out/'adopted-design/report.json')
            else:
                supplied=ev.read(args.recommendation)
                if supplied.get('kind')=='consensus_design':
                    check_hashes(supplied['sources']);design=copy.deepcopy(supplied)
                    from .contact_policy import require_protein_contacts
                    require_protein_contacts(design['design'],design['anchors'])
                    design['sources'].update(fingerprint([args.recommendation]))
                    (out/'adopted-design').mkdir(exist_ok=True);_atomic_json(out/'adopted-design/report.json',design)
                else:design=adopt(args.recommendation,out/'adopted-design')
            check_hashes(design['sources']);definitions=ev.read(args.definitions) if args.definitions else None;validate_regions(design,definitions)
            artifact=ev.read(batch/'artifacts/catalog.json');chemical=ev.read(batch/'chemical/catalog.json')
            if chemical['artifact_v1_catalog_sha256']!=ev.sha(batch/'artifacts/catalog.json') or chemical['library_id']!=artifact['library_id']:
                raise ValueError('Chemical/artifact lineage mismatch')
            for kind,cat in [('artifacts',artifact),('chemical',chemical)]:
                for row in cat['shards']:
                    directory=ev.shard_path(row,batch/kind/'catalog.json')
                    if ev.sha(directory/'manifest.json')!=row['manifest_sha256']:raise ValueError('Changed shard manifest')
                    ev.verify_manifest(directory,full=False)
            if os.name=='posix':ensure_file_descriptor_limit(len(artifact['shards']))
            if (out/'retrieval.json').exists():
                retrieval=ev.read(out/'retrieval.json');check_hashes(retrieval['outputs'])
            else:
                retrieval=retrieve(batch,out,design,args.retrieval_molecules,args.nprobe,args.workers)
                keys=np.load(out/'selected-molecules.npy',allow_pickle=False)
                np.save(out/'selected-conformers.npy',collect_conformers(batch,keys))
                retrieval['outputs']=fingerprint([out/'selected-molecules.npy',out/'selected-conformers.npy',*out.glob('retrieval-*.npz')])
                _atomic_json(out/'retrieval.json',retrieval)
            if args.retrieve_only:
                _atomic_json(out/'status.json',dict(status='retrieval_complete',**retrieval));return
            if (out/'report.json').exists():
                report=ev.read(out/'report.json');check_hashes(report['outputs'])
                _atomic_json(out/'status.json',dict(status='complete'));return
            ids=np.load(out/'selected-conformers.npy',allow_pickle=False,mmap_mode='r')
            timing=run_refinement(batch,out,design,ids,args.workers,args.chunk_conformers)
            count=merge_and_rank(out,design,definitions,batch,args.template_quota,args.rrf_k)
            check_hashes(protocol['sources']);check_hashes(protocol['code'])
            report=dict(kind='consensus_budget_ranking',status='complete',ranked_molecules=count,retrieval=retrieval,refinement=timing,
                target=design['target'],reference=design['reference'],selected_conformers=len(ids),
                ranking='Contact-first pose and template ranks; Gaussian ties only; any-positive-contact tier before zero-contact reserve; quota tier then RRF',
                empirical_filters='annotation_only',physical_filter=design['design']['pocket'],
                biological_validation='not_run',outputs=fingerprint([out/'ranking.sqlite']))
            _atomic_json(out/'report.json',report);_atomic_json(out/'status.json',dict(status='complete'))
        except BaseException as exc:
            _atomic_json(out/'status.json',dict(status='interrupted_or_failed',error=str(exc)));raise


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--batch',type=Path,required=True);p.add_argument('--recommendation',type=Path,required=True)
    p.add_argument('--definitions',type=Path);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--retrieval-molecules',type=int,default=1000000);p.add_argument('--workers',type=int,default=24)
    p.add_argument('--chunk-conformers',type=int,default=64);p.add_argument('--nprobe',type=int,default=128)
    p.add_argument('--template-quota',type=int,default=100000);p.add_argument('--rrf-k',type=int,default=60)
    p.add_argument('--retrieve-only',action='store_true');args=p.parse_args()
    if min(args.retrieval_molecules,args.workers,args.chunk_conformers,args.nprobe,args.template_quota,args.rrf_k)<=0:p.error('Positive parameters required')
    execute(args)


if __name__=='__main__':main()
