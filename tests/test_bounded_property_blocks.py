"""E095 conservation, conservative subdivision and sealed-cache regression tests."""
import json
import sqlite3
import numpy as np
import pytest

from aidd_agent.bounded_property_blocks import POLICY, partition, run, scaling
from aidd_agent.final_work_blocks import sha


def test_constant_and_sparse_tail_do_not_fragment():
    p=dict(POLICY,minimum=50)
    x=np.zeros((1000,69),dtype=np.float32);ids=[str(i) for i in range(len(x))]
    leaves,evidence=partition(x,ids,p)
    assert len(leaves)==1 and evidence[0]['reason']=='constant_projection'
    x[:100]=10
    leaves,evidence=partition(x,ids,p)
    assert len(leaves)==1
    assert evidence[0]['reason']=='population_or_parent_fraction_floor'


def test_separated_properties_split_without_losing_parent_members():
    rng=np.random.default_rng(19)
    x=rng.normal(0,.01,(1000,69)).astype(np.float32)
    x[500:,:33]+=1
    leaves,evidence=partition(x,[str(i) for i in range(len(x))],dict(POLICY,minimum=50))
    assert 2<=len(leaves)<=4
    assert min(map(len,leaves))>=200
    assert np.array_equal(np.sort(np.concatenate(leaves)),np.arange(1000))
    assert any(e['accepted'] and e['check_gain']>=.05 for e in evidence)


def fixture(tmp_path):
    profiles=tmp_path/'profiles';blocks=tmp_path/'backbone'
    profiles.mkdir();blocks.mkdir()
    with sqlite3.connect(blocks/'work_blocks.sqlite') as db:
        db.execute('CREATE TABLE block(block_id TEXT,n INTEGER,kind TEXT)')
        db.executemany('INSERT INTO block VALUES(?,?,?)',[('parent',10000,'single_class'),('special-exhaustive',11,'special_exhaustive')])
    br=dict(status='complete',source_conformers=10011,regular_conformers=10000,special_conformers=11,
            output_hashes={'work_blocks.sqlite':sha(blocks/'work_blocks.sqlite')})
    (blocks/'report.json').write_text(json.dumps(br))
    with sqlite3.connect(profiles/'profiles.sqlite') as db:
        db.execute('CREATE TABLE item(cid TEXT PRIMARY KEY,parent TEXT,vector BLOB)')
        rows=[]
        for i in range(10000):
            x=np.zeros(69,dtype='<f4')
            if i>=5000:x[:11]=1
            rows.append((f'cid-{i:05d}','parent',x.tobytes()))
        db.executemany('INSERT INTO item VALUES(?,?,?)',rows)
        db.execute('CREATE INDEX profile_parent ON item(parent)')
    pr=dict(status='complete',coverage=1,total_regular=10000,blocks=str(blocks),
        signature={'blocks_report_sha256':sha(blocks/'report.json')},
        output_hashes={'profiles.sqlite':sha(profiles/'profiles.sqlite')})
    (profiles/'report.json').write_text(json.dumps(pr))
    receipt=tmp_path/'validation.json'
    receipt.write_text(json.dumps(dict(structural_gate='passed',property_gate='passed',conformers=10011,regular_conformers=10000)))
    return profiles,receipt


def test_complete_build_resume_and_tamper_detection(tmp_path):
    profiles,receipt=fixture(tmp_path);output=tmp_path/'v2'
    r=run(profiles,output,receipt)
    assert r['regular_work_blocks']==2 and r['smallest_regular_block']==5000
    assert r['special_conformers']==11 and r['blocks_below_minimum']==0
    assert r['membership_gate']=='passed_against_complete_cached_profiles'
    before=r['output_hashes']['property_blocks.sqlite']
    assert run(profiles,output,receipt,resume=True)['output_hashes']['property_blocks.sqlite']==before
    with pytest.raises(ValueError,match='Fresh output'):run(profiles,output,receipt)
    with sqlite3.connect(profiles/'profiles.sqlite') as db:db.execute("DELETE FROM item WHERE cid='cid-00000'")
    with pytest.raises(ValueError,match='Input hash changed'):run(profiles,output,receipt,resume=True)


def test_internal_check_failure_is_not_accepted():
    import zlib
    ids=[str(i) for i in range(1000)]
    x=np.zeros((1000,69),dtype=np.float32)
    for i,cid in enumerate(ids):
        if zlib.crc32(cid.encode())%5!=0:x[i]=(-1 if i%2 else 1)
    leaves,evidence=partition(x,ids,dict(POLICY,minimum=50))
    assert len(leaves)==1 and not evidence[0]['accepted']


def test_interrupted_parent_can_resume_and_wrong_owner_is_rejected(tmp_path,monkeypatch):
    import aidd_agent.bounded_property_blocks as module
    profiles,receipt=fixture(tmp_path);output=tmp_path/'interrupted'
    original=module.partition
    def interrupted(*args):raise RuntimeError('Simulated interruption')
    monkeypatch.setattr(module,'partition',interrupted)
    with pytest.raises(RuntimeError,match='Simulated'):run(profiles,output,receipt)
    monkeypatch.setattr(module,'partition',original)
    assert run(profiles,output,receipt,resume=True)['regular_conformers']==10000
    with sqlite3.connect(output/'property_blocks.sqlite') as db:
        db.execute("UPDATE membership SET parent_block='wrong' WHERE cid='cid-00000'")
    with pytest.raises(ValueError,match='ownership mismatch'):run(profiles,output,receipt,resume=True)
