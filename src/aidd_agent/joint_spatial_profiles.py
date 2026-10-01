"""E096 resumable source-verified spatial extraction, preserving E094 ownership."""
import argparse
from concurrent.futures import ProcessPoolExecutor
import json
from pathlib import Path
import sqlite3
import time

import numpy as np
from rdkit import rdBase, RDConfig

from .boundary_pair_review import readonly
from .final_work_blocks import read, sealed, sha
from .joint_spatial_descriptor import VERSION, WIDTH, compute
from .mol2 import iter_mol2_blocks

BATCH_SIZE = 256


def save(path,value):
    path=Path(path);temp=path.with_suffix(path.suffix+'.partial')
    temp.write_text(json.dumps(value,indent=2),encoding='utf-8');temp.replace(path)


def prepare(profiles, build, output, workers=8, resume=False):
    start=time.monotonic()
    profiles,build,output=[Path(p).resolve() for p in (profiles,build,output)]
    if not 1<=workers<=64:raise ValueError('Workers must be 1..64')
    for source in (profiles,build):
        if output==source or source in output.parents or output in source.parents:raise ValueError('Use independent spatial output')
    pr=sealed(profiles,['profiles.sqlite']);blocks=Path(pr['blocks']).resolve()
    br=sealed(blocks,['work_blocks.sqlite']);sr=sealed(build,['descriptors.sqlite'])
    if output==blocks or blocks in output.parents or output in blocks.parents:raise ValueError('Separate output required')
    if pr.get('width')!=69 or pr['coverage']!=1 or pr['total_regular']!=br['regular_conformers']:raise ValueError('Complete E094 69D cache required')
    if pr['signature']['blocks_report_sha256']!=sha(blocks/'report.json') or pr['signature']['build_report_sha256']!=sha(build/'report.json'):raise ValueError('E094 cache provenance mismatch')
    feature_hash=sha(Path(RDConfig.RDDataDir)/'BaseFeatures.fdef')
    if sr['schema']['rdkit']!=rdBase.rdkitVersion or sr['schema']['feature_definition_sha256']!=feature_hash:raise ValueError('Use source RDKit and feature definitions')
    codes=['joint_spatial_profiles.py','joint_spatial_descriptor.py','macrocycle_descriptors.py','macrocycle_peptide.py','macrocycle_blocks.py','mol2.py']
    signature=dict(version=VERSION,width=WIDTH,numpy=np.__version__,rdkit=rdBase.rdkitVersion,feature_definition_sha256=feature_hash,
        sources={str(p):sha(p) for p in (profiles/'report.json',blocks/'report.json',build/'report.json')},
        code_hashes={name:sha(Path(__file__).with_name(name)) for name in codes})
    if output.exists():
        if not resume or read(output/'signature.json')!=signature:raise ValueError('Matching code/inputs and --resume required')
        if (output/'report.json').exists():return sealed(output,['profiles.sqlite'])
    else:output.mkdir(parents=True);save(output/'signature.json',signature)
    with readonly(profiles/'profiles.sqlite') as old,readonly(build/'descriptors.sqlite') as descriptors,sqlite3.connect((output/'profiles.sqlite').as_uri(),uri=True) as db:
        db.executescript('''CREATE TABLE IF NOT EXISTS files(id INTEGER PRIMARY KEY,path TEXT UNIQUE);
          CREATE TABLE IF NOT EXISTS item(cid TEXT PRIMARY KEY,parent TEXT NOT NULL,file_id INTEGER,idx INTEGER,mid TEXT,vector BLOB);
          CREATE TABLE IF NOT EXISTS meta(key TEXT PRIMARY KEY,value TEXT);''')
        if not db.execute("SELECT 1 FROM meta WHERE key='planned'").fetchone():
            print('Planning complete regular spatial population from cached source locators',flush=True)
            db.execute('ATTACH DATABASE ? AS old',((profiles/'profiles.sqlite').as_uri()+'?mode=ro',))
            with db:
                db.execute('INSERT OR IGNORE INTO files SELECT id,path FROM old.files')
                db.execute('INSERT OR IGNORE INTO item(cid,parent,file_id,idx) SELECT cid,parent,file_id,idx FROM old.item')
                if db.execute('SELECT count(*) FROM item').fetchone()[0]!=br['regular_conformers']:raise ValueError('Spatial plan population mismatch')
                db.execute('CREATE INDEX IF NOT EXISTS spatial_file ON item(file_id,idx)')
                db.execute('CREATE INDEX IF NOT EXISTS spatial_parent ON item(parent)')
                db.execute("INSERT INTO meta VALUES('planned','true')")
        ready=db.execute('SELECT count(*) FROM item WHERE vector IS NOT NULL').fetchone()[0];initial_ready=ready
        print(f'Spatial profiles ready: {ready:,}/{br["regular_conformers"]:,}',flush=True)
        executor=ProcessPoolExecutor(max_workers=workers) if workers>1 else None
        try:
            for fid,path in db.execute('SELECT id,path FROM files ORDER BY id').fetchall():
                pending=iter(db.execute('SELECT idx,cid FROM item WHERE file_id=? AND vector IS NULL ORDER BY idx',(fid,)))
                wanted=next(pending,None)
                if wanted is None:continue
                print('Spatial source: '+path,flush=True);batch=[]
                def flush():
                    nonlocal ready
                    results=list(executor.map(compute,batch,chunksize=4) if executor else map(compute,batch))
                    if {r[0] for r in results}!={json.loads(t[3])['conformer_id'] for t in batch}:raise ValueError('Worker identity mismatch')
                    with db:
                        for cid,mid,vector in results:
                            if len(vector)!=WIDTH*4:raise ValueError('Spatial vector width mismatch')
                            changed=db.execute('UPDATE item SET mid=?,vector=? WHERE cid=? AND vector IS NULL',(mid,vector,cid)).rowcount
                            if changed!=1:raise ValueError('Duplicate spatial worker result')
                    ready+=len(results);elapsed=time.monotonic()-start
                    save(output/'progress.json',dict(status='running',profiles_ready=ready,total_regular=br['regular_conformers'],
                        completed_this_invocation=ready-initial_ready,wall_seconds=elapsed,last_source=path,
                        records_per_second_including_startup=(ready-initial_ready)/max(elapsed,1)))
                    print(f'Spatial profiles ready: {ready:,}/{br["regular_conformers"]:,}',flush=True);batch.clear()
                for index,text in iter_mol2_blocks(Path(path)):
                    if wanted is None:break
                    if index!=wanted[0]:continue
                    row=descriptors.execute('SELECT payload FROM descriptor WHERE cid=?',(wanted[1],)).fetchone()
                    prior=old.execute('SELECT vector FROM item WHERE cid=?',(wanted[1],)).fetchone()
                    if row is None or prior is None:raise ValueError('Missing source descriptor or broad profile')
                    batch.append((path,index,text,row[0],sr['schema']['variant'],prior[0]));wanted=next(pending,None)
                    if len(batch)>=BATCH_SIZE:flush()
                if wanted is not None:raise ValueError('Source ended before expected record')
                if batch:flush()
        finally:
            if executor:executor.shutdown(wait=True,cancel_futures=True)
        if db.execute('SELECT count(*) FROM item WHERE vector IS NULL OR mid IS NULL').fetchone()[0]:raise ValueError('Incomplete spatial profiles')
        db.execute('ATTACH DATABASE ? AS verify_old',((profiles/'profiles.sqlite').as_uri()+'?mode=ro',))
        for cid,parent,truth in db.execute('SELECT s.cid,s.parent,o.parent FROM item s LEFT JOIN verify_old.item o ON s.cid=o.cid'):
            if parent!=truth:raise ValueError('Spatial cache parent mismatch: '+cid)
    result=dict(status='complete',version=VERSION,width=WIDTH,blocks=str(blocks),total_regular=br['regular_conformers'],coverage=1,
        prior_profiles=str(profiles),signature=signature,wall_seconds_this_invocation=time.monotonic()-start,
        feature_groups={'broad_properties_and_local_sterics':[0,69],'backbone_geometry':[69,93],'typed_spatial':[93,WIDTH]},
        output_hashes={'profiles.sqlite':sha(output/'profiles.sqlite')},
        limitations=['Finite directed spatial moments can have descriptor collisions; not full atom-wise shape equivalence.',
         'Typed centers are not hydrogen-bond vectors, solvent accessibility, energies or protonation-state enumeration.',
         'Original molecule states, parent identities and source coordinates are preserved; special pool is unchanged.'])
    save(output/'report.json',result);save(output/'progress.json',dict(status='complete',profiles_ready=ready,total_regular=br['regular_conformers']))
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('profiles','build','output'):p.add_argument('--'+name,type=Path,required=True)
    p.add_argument('--workers',type=int,default=8);p.add_argument('--resume',action='store_true');a=p.parse_args()
    print(json.dumps(prepare(a.profiles,a.build,a.output,a.workers,a.resume),indent=2))
