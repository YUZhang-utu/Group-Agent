"""Uniform conformer panels across frozen partitions, with shared source exports."""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from contextlib import closing, contextmanager
import csv
import hashlib
from functools import wraps
import inspect
import json
from pathlib import Path
import random
import sqlite3
import time

from .final_work_blocks import read, sha
from .joint_spatial_profiles import save
from .boundary_pair_review import readonly as open_readonly
from .mol2 import iter_mol2_blocks, parse_mol2_block


@contextmanager
def readonly(path):
    with closing(open_readonly(path)) as db:
        yield db


def exclusive_output(function):
    """One writer per stage, including direct CLI callers outside the chat queue."""
    @wraps(function)
    def locked(*args,**kwargs):
        from .prompt_workflow import file_lock
        output=Path(inspect.signature(function).bind(*args,**kwargs).arguments['output']).resolve()
        output.parent.mkdir(parents=True,exist_ok=True)
        with file_lock(output.parent/('.'+output.name+'.e097.lock')):
            return function(*args,**kwargs)
    return locked


def validate_options(params):
    if not isinstance(params, dict) or set(params)-{'schemes', 'count', 'seed'}:
        raise ValueError('Sampling accepts schemes, count and seed only')
    schemes = params.get('schemes', ['E094', 'E095', 'E096'])
    if not isinstance(schemes, list) or not schemes or any(s not in ('E094','E095','E096') for s in schemes) or len(set(schemes)) != len(schemes):
        raise ValueError('Choose distinct E094/E095/E096 schemes')
    for name, default, low, high in [('count',100,1,10000),('seed',20261002,0,2**32-1)]:
        value=params.get(name,default)
        if type(value) is not int or not low<=value<=high: raise ValueError('Invalid '+name)
    return dict(schemes=sorted(schemes), count=params.get('count',100), seed=params.get('seed',20261002))


def verified_report(path):
    path=Path(path); r=read(path)
    if r.get('status')!='complete': raise ValueError('Incomplete report: '+str(path))
    return r


def check_outputs(directory):
    directory=Path(directory); r=verified_report(directory/'report.json')
    for name,digest in r['output_hashes'].items():
        p=(directory/name).resolve()
        if directory.resolve() not in p.parents or sha(p)!=digest: raise ValueError('Artifact changed: '+name)
    return r


def snapshot(path):
    p=Path(path).resolve(); s=p.stat()
    return dict(path=str(p), size=s.st_size, mtime_ns=s.st_mtime_ns)


@exclusive_output
def sample(profile, output, *, schemes=None, count=100, seed=20261002, resume=False):
    start=time.monotonic(); profile=Path(profile).resolve(); cfg=read(profile)
    opts=validate_options(dict(schemes=schemes or ['E094','E095','E096'],count=count,seed=seed))
    root=profile.parent
    def location(key):
        p=Path(cfg[key]); return (p if p.is_absolute() else root/p).resolve()
    profiles, descriptors, parent = (location(k) for k in ('profiles','descriptors','E094'))
    output=Path(output).resolve()
    paths={'E094':parent, **{s:location(s) for s in opts['schemes'] if s!='E094'}}
    if any(output==p or output in p.parents or p in output.parents for p in [*paths.values(),profiles,descriptors]):
        raise ValueError('Use a separate output directory')
    pr=verified_report(profiles/'report.json'); br=verified_report(parent/'report.json')
    verified_report(descriptors/'report.json')
    if pr.get('coverage')!=1 or pr['total_regular']!=br['regular_conformers']:
        raise ValueError('Complete common regular population required')
    if Path(pr['blocks']).resolve()!=parent: raise ValueError('Profile parent mismatch')
    source_seals=pr['signature']['sources']
    for p in (parent/'report.json', descriptors/'report.json'):
        if source_seals.get(str(p))!=sha(p): raise ValueError('Spatial source seal mismatch')
    receipt=location('validation'); vr=verified_report(receipt)
    if vr.get('structural_gate')!='passed' or vr.get('property_gate')!='passed' or vr.get('regular_conformers')!=br['regular_conformers']:
        raise ValueError('Prior full E094 validation receipt required')
    reports=[profiles/'report.json',parent/'report.json',descriptors/'report.json',receipt]
    databases=[profiles/'profiles.sqlite',descriptors/'descriptors.sqlite']
    populations={}
    for scheme,path in paths.items():
        r=verified_report(path/'report.json'); reports.append(path/'report.json')
        if r['regular_conformers']!=br['regular_conformers']: raise ValueError('Population mismatch')
        if scheme!='E094':
            if r.get('membership_gate')!='passed_against_complete_cached_profiles' or Path(r['parent_blocks']).resolve()!=parent:
                raise ValueError('Unvalidated property membership or different parent')
            if r['sources'].get(str(parent/'report.json'))!=sha(parent/'report.json'): raise ValueError('Parent seal mismatch')
            databases.append(path/'property_blocks.sqlite')
        csvpath=path/'blocks.csv'
        if sha(csvpath)!=r['output_hashes']['blocks.csv']: raise ValueError('Block CSV seal mismatch')
        with csvpath.open(encoding='utf-8',newline='') as f:
            rows=list(csv.DictReader(f))
        populations[scheme]={row['block_id']:int(row['conformers']) for row in rows
                             if row.get('kind','regular') not in ('special','special-exhaustive') and row['block_id']!='special-exhaustive'}
        if sum(populations[scheme].values())!=br['regular_conformers']: raise ValueError('Regular block totals disagree')
    signature=dict(options=opts, sources={str(p):sha(p) for p in reports},
                   database_snapshots=[snapshot(p) for p in databases], implementation_sha256=sha(Path(__file__)))
    if output.exists():
        if not resume or read(output/'signature.json')!=signature: raise ValueError('Use fresh output or matching --resume')
        if (output/'report.json').exists(): return check_outputs(output)
    else:
        output.mkdir(parents=True); save(output/'signature.json',signature)
    selected_path=output/'samples.sqlite'
    # Publish selection only after every full membership stream has been checked.
    if not (output/'selection.json').exists():
        tmp=output/'samples.partial.sqlite'
        if tmp.exists(): tmp.unlink()
        with closing(sqlite3.connect(tmp)) as out:
            out.executescript('''CREATE TABLE sample(scheme TEXT,block_id TEXT,cid TEXT,
                PRIMARY KEY(scheme,block_id,cid));
                CREATE TABLE population(scheme TEXT,block_id TEXT,n INTEGER,PRIMARY KEY(scheme,block_id));
                CREATE TABLE selected(cid TEXT PRIMARY KEY,alias TEXT UNIQUE,mid TEXT,file_id INTEGER,idx INTEGER,parent_block TEXT);
                CREATE TABLE files(id INTEGER PRIMARY KEY,path TEXT);''')
            for scheme in opts['schemes']:
                seen=Counter(); reservoirs={}; rngs={}
                dbpath=profiles/'profiles.sqlite' if scheme=='E094' else paths[scheme]/'property_blocks.sqlite'
                query='SELECT cid,parent FROM item' if scheme=='E094' else 'SELECT cid,block_id FROM membership'
                if scheme=='E094':
                    proxy=next((s for s in ('E095','E096') if s in paths),None)
                    if proxy:
                        # The validated child table includes unchanged E094 ownership.
                        # Avoid scanning the 317D vector payload for parent sampling.
                        dbpath=paths[proxy]/'property_blocks.sqlite'
                        query='SELECT cid,parent_block FROM membership'
                with readonly(dbpath) as db:
                    for total,(cid,bid) in enumerate(db.execute(query),1):
                        if bid not in populations[scheme]: raise ValueError('Unknown regular block')
                        if bid not in reservoirs:
                            reservoirs[bid]=[]
                            rngs[bid]=random.Random(int.from_bytes(hashlib.sha256(f'{seed}:{scheme}:{bid}'.encode()).digest(),'big'))
                        seen[bid]+=1; bucket=reservoirs[bid]
                        if len(bucket)<count: bucket.append(cid)
                        else:
                            j=rngs[bid].randrange(seen[bid])
                            if j<count: bucket[j]=cid
                        if total%1000000==0:
                            print(f'{scheme}: sampled stream {total:,}',flush=True)
                            save(output/'progress.json',dict(stage='sampling',scheme=scheme,visited=total))
                if dict(seen)!=populations[scheme]: raise ValueError('Full membership counts disagree')
                out.executemany('INSERT INTO population VALUES(?,?,?)',[(scheme,b,n) for b,n in seen.items()])
                for bid,bucket in reservoirs.items():
                    out.executemany('INSERT INTO sample VALUES(?,?,?)',[(scheme,bid,c) for c in bucket])
                out.commit()
            with readonly(profiles/'profiles.sqlite') as db:
                out.executemany('INSERT INTO files VALUES(?,?)',db.execute('SELECT id,path FROM files'))
                cids=[r[0] for r in out.execute('SELECT DISTINCT cid FROM sample ORDER BY cid')]
                for i,cid in enumerate(cids,1):
                    row=db.execute('SELECT mid,file_id,idx,parent FROM item WHERE cid=?',(cid,)).fetchone()
                    if row is None or row[0] is None: raise ValueError('Selected conformer absent from common source')
                    out.execute('INSERT INTO selected VALUES(?,?,?,?,?,?)',(cid,f'c{i:09d}',*row))
                for scheme in opts['schemes']:
                    if scheme=='E094':
                        if out.execute("SELECT 1 FROM sample s JOIN selected t USING(cid) WHERE s.scheme='E094' AND s.block_id!=t.parent_block LIMIT 1").fetchone():
                            raise ValueError('Sampled E094 parent disagrees with source profiles')
                    else:
                        with readonly(paths[scheme]/'property_blocks.sqlite') as child:
                            for cid,bid,parent_id in out.execute('SELECT s.cid,s.block_id,t.parent_block FROM sample s JOIN selected t USING(cid) WHERE s.scheme=?',(scheme,)):
                                if child.execute('SELECT block_id,parent_block FROM membership WHERE cid=?',(cid,)).fetchone()!=(bid,parent_id):
                                    raise ValueError('Sampled child ownership disagrees with common profiles')
                out.execute('CREATE UNIQUE INDEX source_locator ON selected(file_id,idx)')
                out.execute('CREATE INDEX sample_cid ON sample(cid)')
                out.commit()
        tmp.replace(selected_path); save(output/'selection.json',dict(sha256=sha(selected_path)))
    if sha(selected_path)!=read(output/'selection.json')['sha256']: raise ValueError('Selection changed')
    chunks=[]; exported=0
    (output/'ligands').mkdir(exist_ok=True); (output/'receipts').mkdir(exist_ok=True)
    with readonly(selected_path) as db, readonly(descriptors/'descriptors.sqlite') as source:
        for fid,raw in db.execute('SELECT id,path FROM files ORDER BY id'):
            wanted={i:(cid,alias,mid) for cid,alias,mid,i in db.execute('SELECT cid,alias,mid,idx FROM selected WHERE file_id=?',(fid,))}
            if not wanted: continue
            receipt_path=output/'receipts'/f'{fid}.json'
            if receipt_path.exists():
                saved=read(receipt_path)
                if saved['source']!=snapshot(raw): raise ValueError('Raw source changed after export')
                for part in saved['chunks']:
                    if sha(output/part['path'])!=part['sha256']: raise ValueError('Export changed')
                chunks.extend(saved['chunks']); exported+=sum(p['count'] for p in saved['chunks']); continue
            print(f'Export selected source records: {raw}',flush=True)
            before=snapshot(raw); parts=[]; pending=[]; metadata=[]
            def flush():
                if not pending: return
                rel=f'ligands/f{fid:06d}-{len(parts):06d}.mol2'; p=output/rel
                temp=p.with_suffix('.partial'); temp.write_text(''.join(pending),encoding='utf-8',newline='\n'); temp.replace(p)
                parts.append(dict(path=rel,sha256=sha(p),count=len(pending),records=list(metadata)))
                pending.clear(); metadata.clear()
            for index,text in iter_mol2_blocks(Path(raw)):
                item=wanted.pop(index,None)
                if item is None: continue
                cid,alias,mid=item
                row=source.execute('SELECT payload FROM descriptor WHERE cid=?',(cid,)).fetchone()
                if row is None: raise ValueError('Missing source descriptor')
                d=json.loads(row[0]); p=d['provenance']; record=parse_mol2_block(Path(raw),index,text)
                if (d['conformer_id'],d['molecule_id'])!=(cid,mid) or (record.content_sha256,record.name,record.molecule_name)!=(p['content_sha256'],p['source_record_name'],p['source_name']):
                    raise ValueError('Raw identity/hash mismatch')
                if Path(p['source_path']).resolve()!=Path(raw).resolve() or p['source_record_index']!=index: raise ValueError('Raw locator mismatch')
                lines=text.splitlines(keepends=True); lines[1]=alias+'\n'
                renamed=''.join(lines)
                if not renamed.endswith('\n'): renamed+='\n'
                # Only the title is changed. Atom types, coordinates and bonds remain verbatim.
                pending.append(renamed); metadata.append(dict(cid=cid,alias=alias,mid=mid,source_index=index,source_sha256=record.content_sha256))
                if len(pending)==100: flush()
                if not wanted: break
            if wanted: raise ValueError('Source ended before selected records')
            flush()
            if snapshot(raw)!=before: raise ValueError('Source changed during export')
            save(receipt_path,dict(source=before,chunks=parts)); chunks.extend(parts)
            exported+=sum(p['count'] for p in parts)
            save(output/'progress.json',dict(stage='export',exported=exported))
        slots=db.execute('SELECT count(*) FROM sample').fetchone()[0]
        unique=db.execute('SELECT count(*) FROM selected').fetchone()[0]
        if exported!=unique: raise ValueError('Export count mismatch')
        with (output/'block_samples.csv').open('w',encoding='utf-8',newline='') as f:
            w=csv.writer(f); w.writerow(['scheme','block_id','conformer_id','alias','molecule_id'])
            w.writerows(db.execute('SELECT s.scheme,s.block_id,s.cid,t.alias,t.mid FROM sample s JOIN selected t USING(cid) ORDER BY s.scheme,s.block_id,s.cid'))
    save(output/'exports.json',dict(chunks=chunks))
    if [snapshot(p) for p in databases]!=signature['database_snapshots']: raise ValueError('Frozen database changed')
    outputs=['samples.sqlite','block_samples.csv','exports.json','signature.json', *[p['path'] for p in chunks]]
    r=dict(status='complete',kind='block_sample',slots=slots,unique_conformers=unique,reused_slots=slots-unique,
           blocks={s:len(populations[s]) for s in opts['schemes']},options=opts,exported_conformers=exported,
           output_hashes={name:sha(output/name) for name in outputs},wall_seconds=time.monotonic()-start,
           limitations=['Uniform conformer samples, not unique-molecule samples or recall validation.',
                        'Special and deferred unsupported chemistry are excluded.',
                        'Prior full membership receipts trusted; large source databases use stat snapshots, not repeated hashes.',
                        'Raw selected record hashes verified; ligand protonation/preparation still requires review.'])
    save(output/'report.json',r); return r


def main():
    p=argparse.ArgumentParser(description=__doc__); p.add_argument('--profile',required=True); p.add_argument('--output',required=True)
    p.add_argument('--schemes',nargs='+',default=['E094','E095','E096']); p.add_argument('--count',type=int,default=100)
    p.add_argument('--seed',type=int,default=20261002); p.add_argument('--resume',action='store_true')
    a=p.parse_args(); print(json.dumps(sample(**vars(a)),indent=2))


if __name__=='__main__': main()
