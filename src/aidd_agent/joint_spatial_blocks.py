"""E096 bounded joint backbone, chemistry and position-sensitive spatial partitions."""
import argparse
from collections import Counter
import json
from pathlib import Path
import sqlite3
import time

import numpy as np

from .bounded_property_blocks import POLICY, CHEM, STERIC, partition
from .joint_spatial_descriptor import WIDTH, VERSION as SPATIAL_VERSION
from .macrocycle_descriptors import FAMILIES
from .boundary_pair_review import readonly
from .final_work_blocks import read, sealed, sha, identifier, write_csv

VERSION='bounded-joint-spatial-work-blocks-v1'
GROUP_SPEC=dict(broad_chemistry=.20,local_sterics=.10,backbone_geometry=.35,typed_spatial=.35)
GROUPS=[(CHEM,.20),(STERIC,.10),(np.arange(69,93),.35),(np.arange(93,WIDTH),.35)]


def scaling(db,expected):
    count=0;sums=np.zeros(WIDTH);squares=np.zeros(WIDTH)
    cursor=db.execute('SELECT vector FROM item')
    while True:
        batch=cursor.fetchmany(4096)
        if not batch:break
        if any(v is None or len(v)!=WIDTH*4 for (v,) in batch):raise ValueError('Invalid spatial cache')
        x=np.frombuffer(b''.join(v for (v,) in batch),dtype='<f4').reshape(-1,WIDTH).astype(float)
        if not np.isfinite(x).all():raise ValueError('Nonfinite spatial cache')
        count+=len(x);sums+=x.sum(0);squares+=(x*x).sum(0)
        if count//1000000!=(count-len(x))//1000000:print(f'Joint global scaling: {count:,}',flush=True)
    if count!=expected:raise ValueError('Spatial scaling population mismatch')
    mean=sums/count;std=np.sqrt(np.maximum(squares/count-mean*mean,0));weights=np.zeros(WIDTH)
    for group,budget in GROUPS:
        active=group[std[group]>=POLICY['active_std_floor']]
        if len(active):weights[active]=np.sqrt(budget/len(active))/np.maximum(std[active],POLICY['scale_floor'])
    return dict(mean=mean.tolist(),std=std.tolist(),weights=weights.tolist(),group_budgets=GROUP_SPEC,
                active_channels_by_group={name:int((weights[group]>0).sum()) for name,(group,_) in zip(GROUP_SPEC,GROUPS)},
                mean_residue_feature_counts={name:float(mean[93+i*7]*(10 if i==7 else 4)) for i,name in enumerate((*FAMILIES,'SideHeavyAtoms'))},
                scope='Global fixed group budgets; weights are engineering defaults, not calibrated screening scores')


def run(profiles, output, receipt, resume=False):
    started=time.monotonic(); profiles=Path(profiles).resolve();output=Path(output).resolve()
    pr=sealed(profiles,['profiles.sqlite']);blocks=Path(pr['blocks']).resolve()
    if pr.get('width')!=WIDTH or pr.get('version')!=SPATIAL_VERSION:raise ValueError('Completed E096 spatial cache required')
    br=sealed(blocks,['work_blocks.sqlite'])
    receipt=Path(receipt).resolve();vr=read(receipt)
    if vr.get('structural_gate')!='passed' or vr.get('property_gate')!='passed':raise ValueError('Successful E094 validation receipt required')
    if vr.get('conformers')!=br['source_conformers'] or vr.get('regular_conformers')!=br['regular_conformers']:raise ValueError('Validation population mismatch')
    if pr['coverage']!=1 or pr['total_regular']!=br['regular_conformers'] or pr['signature']['sources'][str(blocks/'report.json')]!=sha(blocks/'report.json'):raise ValueError('Complete matching profiles required')
    for source in (profiles,blocks):
        if output==source or source in output.parents or output in source.parents:raise ValueError('Separate output required')
    signature=dict(version=VERSION,policy=POLICY,implementation=sha(Path(__file__)),numpy=np.__version__,
                   partition_implementation=sha(Path(__file__).with_name('bounded_property_blocks.py')),
                   sources={str(p):sha(p) for p in (profiles/'report.json',blocks/'report.json',receipt)})
    if output.exists():
        if not resume or read(output/'signature.json')!=signature:raise ValueError('Fresh output or matching --resume required')
    else:
        output.mkdir(parents=True);(output/'signature.json').write_text(json.dumps(signature,indent=2))
    with readonly(profiles/'profiles.sqlite') as source,readonly(blocks/'work_blocks.sqlite') as work,sqlite3.connect((output/'property_blocks.sqlite').as_uri(),uri=True) as dest:
        dest.executescript('''CREATE TABLE IF NOT EXISTS membership(cid TEXT PRIMARY KEY,parent_block TEXT,block_id TEXT);
          CREATE TABLE IF NOT EXISTS block(block_id TEXT PRIMARY KEY,parent_block TEXT,n INTEGER,metadata TEXT);
          CREATE TABLE IF NOT EXISTS completed(parent TEXT PRIMARY KEY,evidence TEXT);''')
        transform_path=output/'transform.json'
        if transform_path.exists():
            transform=read(transform_path)
            if sha(transform_path)!=read(output/'transform-seal.json')['sha256']:raise ValueError('Transform changed')
        else:
            transform=scaling(source,br['regular_conformers'])
            transform_path.write_text(json.dumps(transform,indent=2))
            (output/'transform-seal.json').write_text(json.dumps(dict(sha256=sha(transform_path))))
        mean=np.array(transform['mean']);weight=np.array(transform['weights'])
        parents=work.execute("SELECT block_id,n FROM block WHERE kind!='special_exhaustive' ORDER BY block_id").fetchall()
        for number,(parent,n) in enumerate(parents,1):
            if dest.execute('SELECT 1 FROM completed WHERE parent=?',(parent,)).fetchone():continue
            print(f'Parent {number}/{len(parents)}: {parent} / {n:,}',flush=True)
            x=np.empty((n,WIDTH),dtype=np.float32);cids=[];mids=[]
            for cid,mid,blob in source.execute('SELECT cid,mid,vector FROM item WHERE parent=? ORDER BY cid',(parent,)):
                if len(cids)>=n or blob is None or len(blob)!=WIDTH*4 or not mid:raise ValueError('Invalid parent profiles')
                x[len(cids)]=(np.frombuffer(blob,dtype='<f4')-mean)*weight;cids.append(cid);mids.append(mid)
            if len(cids)!=n or not np.isfinite(x).all():raise ValueError('Incomplete parent')
            leaves,evidence=partition(x,mids)
            if len(leaves)>1 and min(map(len,leaves))<max(POLICY['minimum'],int(np.ceil(.2*n))):raise ValueError('Tiny child')
            with dest:
                for child,indices in enumerate(leaves):
                    bid=identifier([signature,parent,child]);meta=dict(block_id=bid,parent_block=parent,
                        conformers=len(indices),split=len(leaves)>1,split_evidence=evidence)
                    dest.execute('INSERT INTO block VALUES(?,?,?,?)',(bid,parent,len(indices),json.dumps(meta)))
                    dest.executemany('INSERT INTO membership VALUES(?,?,?)',((cids[int(i)],parent,bid) for i in indices))
                dest.execute('INSERT INTO completed VALUES(?,?)',(parent,json.dumps(evidence)))
            del x
        dest.execute('CREATE INDEX IF NOT EXISTS member_block ON membership(block_id)');dest.commit()
        # Full enumeration against cached source ownership; no 100-GiB frozen-model rescan.
        dest.execute('ATTACH DATABASE ? AS profiles',( (profiles/'profiles.sqlite').as_uri()+'?mode=ro',))
        counts=Counter();total=0
        for child,parent,truth,owner in dest.execute('''SELECT m.block_id,m.parent_block,p.parent,b.parent_block
            FROM membership m LEFT JOIN profiles.item p ON p.cid=m.cid LEFT JOIN block b ON b.block_id=m.block_id'''):
            if parent!=truth or parent!=owner or parent=='special-exhaustive':raise ValueError('Membership ownership mismatch')
            counts[child]+=1;total+=1
            if total%1000000==0:print(f'Validated property memberships: {total:,}',flush=True)
        expected=dict(dest.execute('SELECT block_id,n FROM block'))
        if counts!=expected or total!=br['regular_conformers']:raise ValueError('Membership conservation failed')
        rows=[json.loads(m) for (m,) in dest.execute('SELECT metadata FROM block ORDER BY block_id')]
        decisions=[dict(parent_block=p,proposals=json.loads(e)) for p,e in dest.execute('SELECT parent,evidence FROM completed ORDER BY parent')]
    write_csv(output/'blocks.csv',rows,['block_id','parent_block','conformers','split','split_evidence'])
    (output/'split-decisions.json').write_text(json.dumps(decisions,indent=2))
    reasons=Counter(e['reason'] for d in decisions for e in d['proposals'])
    report=dict(status='complete',readiness='offline_property_work_blocks_not_validated_search_dispatch',
       parent_blocks=str(blocks),profiles=str(profiles),regular_work_blocks=len(rows),parent_regular_blocks=len(parents),
       added_blocks=len(rows)-len(parents),regular_conformers=total,special_conformers=br['special_conformers'],
       special_pool='special-exhaustive in unchanged parent',minimum_regular_block_size=POLICY['minimum'],
       smallest_regular_block=min(counts.values()),largest_regular_block=max(counts.values()),
       blocks_below_minimum=sum(n<POLICY['minimum'] for n in counts.values()),logical_capacity_limit=None,
       maximum_children_per_parent=4,policy=POLICY,decision_counts=dict(reasons),
       membership_gate='passed_against_complete_cached_profiles',sources=signature['sources'],
       implementation_sha256=signature['implementation'],wall_seconds=time.monotonic()-started,
       limitations=['Joint geometry/chemistry work partitions, not discrete chemical states or atom-wise equivalence.',
        'Internal checks are grouped by source molecule ID; global scaling and adaptive reuse remain exploratory.',
        'Global scaling and repeated adaptive checks are exploratory, not unbiased evaluation.',
        'Signed residue-anchored spatial moments are position-sensitive summaries, not collision-free encodings or H-bond directionality.',
        'No old radius reuse, block rejection, docking or recall claim.',
        'Parent chemical identities and special pool remain unchanged; prior E094 validation is trusted for source ownership.'])
    report['feature_groups']=GROUP_SPEC
    report['active_channels_by_group']=transform['active_channels_by_group']
    report['mean_residue_feature_counts']=transform['mean_residue_feature_counts']
    report['check_grouping']='source_molecule_id_crc32_mod5'
    report['output_hashes']={name:sha(output/name) for name in ('property_blocks.sqlite','blocks.csv','transform.json','split-decisions.json')}
    (output/'report.json').write_text(json.dumps(report,indent=2))
    return report


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('profiles','output','receipt'):p.add_argument('--'+name,type=Path,required=True)
    p.add_argument('--resume',action='store_true');a=p.parse_args()
    print(json.dumps(run(a.profiles,a.output,a.receipt,a.resume),indent=2))
