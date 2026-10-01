"""E095: bounded backbone-parent refinement using cached broad side-chain properties."""
import argparse
from collections import Counter
import json
from pathlib import Path
import sqlite3
import time
import zlib

import numpy as np

from .final_work_blocks import identifier, read, sealed, sha, write_csv
from .boundary_pair_review import readonly

VERSION = 'bounded-property-work-blocks-v2'
CHEM = np.array([i+j for i in (0,23,46) for j in range(11)])
STERIC = np.array([i+j for i in (0,23,46) for j in range(11,23)])
POLICY = dict(minimum=5000, minimum_parent_fraction=.20, maximum_children=4,
              minimum_gain=.05, minimum_check_to_fit=.75, minimum_effect=.35,
              scale_floor=.10, active_std_floor=.02, fit_limit=8192)


def scaling(db, expected):
    total=0; sums=np.zeros(69); squares=np.zeros(69)
    cursor=db.execute('SELECT vector FROM item')
    while True:
        batch=cursor.fetchmany(8192)
        if not batch:break
        if any(v is None or len(v)!=276 for (v,) in batch):raise ValueError('Invalid profile cache')
        x=np.frombuffer(b''.join(v for (v,) in batch),dtype='<f4').reshape(-1,69).astype(float)
        if not np.isfinite(x).all():raise ValueError('Nonfinite profile')
        sums+=x.sum(0); squares+=(x*x).sum(0);total+=len(x)
        if total//1000000 != (total-len(x))//1000000:print(f'Global scaling: {total:,}',flush=True)
    if total!=expected:raise ValueError('Global profile coverage mismatch')
    mean=sums/total; std=np.sqrt(np.maximum(squares/total-mean*mean,0))
    weights=np.zeros(69)
    for group in (CHEM,STERIC):
        active=group[std[group]>=POLICY['active_std_floor']]
        if len(active):weights[active]=np.sqrt(.5/len(active))/np.maximum(std[active],POLICY['scale_floor'])
    return dict(mean=mean.tolist(),std=std.tolist(),weights=weights.tolist(),
                scope='One global transform; equal chemistry/steric group budget; fixed floor limits rare-channel amplification')


def gain(x, labels, centers, baseline):
    before=float(np.mean(np.sum((x-baseline)**2,axis=1)))
    after=float(np.mean(np.sum((x-centers[labels])**2,axis=1)))
    return (before-after)/before if before>1e-12 else 0.


def partition(x, cids, policy=None):
    """Use an internal conformer check; this is not molecule-independent evaluation."""
    p=policy or POLICY
    floor=max(p['minimum'],int(np.ceil(p['minimum_parent_fraction']*len(x))))
    check=np.array([zlib.crc32(cid.encode())%5==0 for cid in cids])
    leaves=[np.arange(len(x))]; evidence=[]

    def proposal(indices):
        note=dict(population=len(indices),minimum_child=floor)
        def reject(reason):
            note.update(accepted=False,reason=reason);evidence.append(note);return None
        if len(indices)<2*floor:return reject('population_or_parent_fraction_floor')
        train=indices[~check[indices]];test=indices[check[indices]]
        if min(len(train),len(test))<32:return reject('insufficient_internal_check')
        if len(train)>p['fit_limit']:train=train[np.linspace(0,len(train)-1,p['fit_limit'],dtype=int)]
        a=x[train].astype(float);base=a.mean(0);centered=a-base
        _,axes=np.linalg.eigh(centered.T@centered/len(a));axis=axes[:,-1]
        if axis[np.argmax(abs(axis))]<0:axis=-axis
        projection=centered@axis
        if np.ptp(projection)<1e-7:return reject('constant_projection')
        threshold=float(np.median(projection))
        fit_labels=(projection>threshold).astype(int)
        if min(np.bincount(fit_labels,minlength=2))<16:return reject('sparse_fit_child')
        centers=np.array([a[fit_labels==i].mean(0) for i in range(2)])
        labels=((x[indices]-base)@axis>threshold).astype(int)
        counts=np.bincount(labels,minlength=2)
        note.update(left_count=int(counts[0]),right_count=int(counts[1]))
        if min(counts)<floor:return reject('population_or_parent_fraction_floor')
        b=x[test].astype(float);test_labels=((b-base)@axis>threshold).astype(int)
        if min(np.bincount(test_labels,minlength=2))<16:return reject('sparse_check_child')
        fg=gain(a,fit_labels,centers,base);cg=gain(b,test_labels,centers,base)
        effect=float(np.linalg.norm(centers[0]-centers[1]))
        note.update(fit_gain=fg,check_gain=cg,standardized_effect=effect,
                    fit_conformers=len(train),check_conformers=len(test))
        if fg<p['minimum_gain'] or cg<p['minimum_gain']:return reject('insufficient_dispersion_reduction')
        if cg<p['minimum_check_to_fit']*fg:return reject('internal_check_not_stable')
        if effect<p['minimum_effect']:return reject('small_standardized_effect')
        note.update(accepted=True,reason='accepted_internal_check',axis=axis.tolist(),
                    fit_center=base.tolist(),threshold=threshold)
        evidence.append(note)
        return indices[labels==0],indices[labels==1]

    # Largest-first deterministic recursion; stop at four children, never force k.
    blocked=set()
    while len(leaves)<p['maximum_children']:
        choices=[i for i,v in enumerate(leaves) if id(v) not in blocked]
        if not choices:break
        i=max(choices,key=lambda i:len(leaves[i]));value=proposal(leaves[i])
        if value is None:blocked.add(id(leaves[i]))
        else:leaves[i:i+1]=list(value)
    return leaves,evidence


def run(profiles, output, receipt, resume=False):
    started=time.monotonic(); profiles=Path(profiles).resolve();output=Path(output).resolve()
    pr=sealed(profiles,['profiles.sqlite']);blocks=Path(pr['blocks']).resolve()
    br=sealed(blocks,['work_blocks.sqlite'])
    receipt=Path(receipt).resolve();vr=read(receipt)
    if vr.get('structural_gate')!='passed' or vr.get('property_gate')!='passed':raise ValueError('Successful E094 validation receipt required')
    if vr.get('conformers')!=br['source_conformers'] or vr.get('regular_conformers')!=br['regular_conformers']:raise ValueError('Validation population mismatch')
    if pr['coverage']!=1 or pr['total_regular']!=br['regular_conformers'] or pr['signature']['blocks_report_sha256']!=sha(blocks/'report.json'):raise ValueError('Complete matching profiles required')
    for source in (profiles,blocks):
        if output==source or source in output.parents or output in source.parents:raise ValueError('Separate output required')
    signature=dict(version=VERSION,policy=POLICY,implementation=sha(Path(__file__)),numpy=np.__version__,
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
            x=np.empty((n,69),dtype=np.float32);cids=[]
            for cid,blob in source.execute('SELECT cid,vector FROM item WHERE parent=? ORDER BY cid',(parent,)):
                if len(cids)>=n or blob is None or len(blob)!=276:raise ValueError('Invalid parent profiles')
                x[len(cids)]=(np.frombuffer(blob,dtype='<f4')-mean)*weight;cids.append(cid)
            if len(cids)!=n or not np.isfinite(x).all():raise ValueError('Incomplete parent')
            leaves,evidence=partition(x,cids)
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
       limitations=['Broad composition/steric work partitions, not discrete chemical states.',
        'Internal checks split conformers, not molecules; shared molecules can cross fit/check.',
        'Global scaling and repeated adaptive checks are exploratory, not unbiased evaluation.',
        'Position-sensitive side-chain patterns cannot be reconstructed from these summaries.',
        'No old radius reuse, block rejection, docking or recall claim.',
        'Parent chemical identities and special pool remain unchanged; prior E094 validation is trusted for source ownership.'])
    report['output_hashes']={name:sha(output/name) for name in ('property_blocks.sqlite','blocks.csv','transform.json','split-decisions.json')}
    (output/'report.json').write_text(json.dumps(report,indent=2))
    return report


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('profiles','output','receipt'):p.add_argument('--'+name,type=Path,required=True)
    p.add_argument('--resume',action='store_true');a=p.parse_args()
    print(json.dumps(run(a.profiles,a.output,a.receipt,a.resume),indent=2))
