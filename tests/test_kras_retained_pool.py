"""Retention, provenance and exact-source downgrade tests for E121."""
import hashlib
import json
from pathlib import Path
import sqlite3
import sys

import pytest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
import export_kras_retained_pool as pool

RAW='@<TRIPOS>MOLECULE\nwater\n3 2 1 0 0\nSMALL\nUSER_CHARGES\n\n@<TRIPOS>ATOM\n1 O 0 0 0 O.3 1 WAT -0.8\n2 H1 0.96 0 0 H 1 WAT 0.4\n3 H2 -0.24 0.93 0 H 1 WAT 0.4\n@<TRIPOS>BOND\n1 1 2 1\n2 2 3 1\n'
# Use valid water connectivity for strict sanitization.
RAW=RAW.replace('2 2 3 1','2 1 3 1')


def fixture(tmp_path):
    plan=tmp_path/'plan'; plan.mkdir(); source=tmp_path/'original.mol2'
    source.write_text(RAW*4,encoding='utf-8',newline='')
    db=sqlite3.connect(plan/'pool.sqlite')
    db.execute('CREATE TABLE molecules(target TEXT,mid TEXT,PRIMARY KEY(target,mid))')
    reg=sqlite3.connect(':memory:')
    reg.execute('CREATE TABLE conformer(id TEXT,molecule_id TEXT,conformer_index INTEGER,source_path TEXT,source_record_index INTEGER,content_sha256 TEXT)')
    ranked=sqlite3.connect(':memory:'); pool.soft.create_tables(ranked)
    for t in pool.TARGETS:
        db.executemany('INSERT INTO molecules VALUES (?,?)',[(t,'A'),(t,'B'),(t,'C'),(t,'MISSING')])
    digest=hashlib.sha256(RAW.encode()).hexdigest()
    for cid,mid,idx,record in [('A0','A',0,0),('A1','A',1,1),('B0','B',0,2),('C0','C',0,3)]:
        reg.execute('INSERT INTO conformer VALUES (?,?,?,?,?,?)',(cid,mid,idx,str(source),record,digest))
    for t in pool.TARGETS:
        pose=dict(mid='A',cid='A1',gid=1,score=.4)
        ranked.execute('INSERT INTO molecules VALUES (?,?,?,?,?)',(t,'A',.4,1,json.dumps([pose])))
    return plan,source,db,reg,ranked


def finish_plan(plan,db,reg,ranked):
    pool.map_members(db,reg,ranked)
    assert db.execute('SELECT COUNT(*) FROM selected').fetchone()[0]==6
    assert db.execute('SELECT COUNT(*) FROM mapping_issues').fetchone()[0]==2
    assert set(db.execute("SELECT cid FROM selected WHERE mid='A'"))=={('A1',)}
    db.close();reg.close();ranked.close()
    pool.base.save(plan/'report.json',dict(status='complete',scope='retained_pool_mapping',
        pool_sha256=pool.base.sha(plan/'pool.sqlite'),targets={t:dict(pool_molecules=4) for t in pool.TARGETS}))
    pool.base.save(plan/'protocol.json',dict(original=dict(batch=str(plan.parent/'library'),prepared=str(plan.parent/'prepared'))))


def test_partition_export_and_downgrade(tmp_path):
    plan,source,db,reg,ranked=fixture(tmp_path)
    finish_plan(plan,db,reg,ranked)
    def factory(batch,prepared):
        return lambda molecule,target,pose:(target=='KRAS_G12D',2.1)
    out=tmp_path/'export'
    result=pool.write_structures(plan,out,batch_size=2,check_factory=factory)
    d=result['targets']['KRAS_G12D'];v=result['targets']['KRAS_G12V']
    assert (d['priority'],d['unresolved'],d['quarantined'])==(1,2,1)
    assert (v['priority'],v['unresolved'],v['quarantined'],v['downgraded'])==(0,3,1,1)
    assert [b['molecules'] for b in v['batches']['unresolved']]==[2,1]
    raw=next(pool.base.iter_mol2_blocks(out/'KRAS_G12D/priority/batch_001.mol2'))[1]
    assert raw.splitlines()[2:]==RAW.splitlines()[2:]
    assert result['new_scientific_rejections']==0


def test_outside_pool_priority_rejected(tmp_path):
    plan,source,db,reg,ranked=fixture(tmp_path)
    ranked.execute('INSERT INTO molecules VALUES (?,?,?,?,?)',('KRAS_G12D','OUT',1,0,'[]'))
    with pytest.raises(ValueError,match='outside direct pool'):
        pool.map_members(db,reg,ranked)
    db.close();reg.close();ranked.close()


def test_source_tampering_fails(tmp_path):
    plan,source,db,reg,ranked=fixture(tmp_path);finish_plan(plan,db,reg,ranked)
    source.write_text(RAW.replace('-0.8','-0.7')*4,encoding='utf-8',newline='')
    with pytest.raises(ValueError,match='hash mismatch'):
        pool.write_structures(plan,tmp_path/'export')


def test_snapshot_includes_committed_wal(tmp_path):
    source=tmp_path/'live.sqlite';db=sqlite3.connect(source)
    db.execute('PRAGMA journal_mode=WAL');db.execute('CREATE TABLE x(v INTEGER)')
    db.execute('INSERT INTO x VALUES (1)');db.commit()
    pool.snapshot(source,tmp_path/'copy.sqlite')
    with sqlite3.connect(tmp_path/'copy.sqlite') as copied:
        assert copied.execute('SELECT * FROM x').fetchall()==[(1,)]
    db.close()


def test_prepare_without_ranking_preserves_all(tmp_path,monkeypatch):
    direct=tmp_path/'direct';direct.mkdir();batch=tmp_path/'library';batch.mkdir()
    original=dict(batch=str(batch),prepared=str(tmp_path/'prepared'),chunk_size=4,sample_size=0)
    monkeypatch.setattr(pool.soft.run,'seal_inputs',lambda *args:original)
    with sqlite3.connect(direct/'counts.sqlite') as db:
        db.execute('CREATE TABLE molecules(target TEXT,mid TEXT,PRIMARY KEY(target,mid))')
        db.executemany('INSERT INTO molecules VALUES (?,?)',[(t,'M') for t in pool.TARGETS])
    with sqlite3.connect(batch/'registry.sqlite3') as db:
        db.execute('CREATE TABLE conformer(id TEXT,molecule_id TEXT,conformer_index INTEGER,source_path TEXT,source_record_index INTEGER,content_sha256 TEXT)')
        db.execute('INSERT INTO conformer VALUES (?,?,?,?,?,?)',('C','M',0,'source',0,'hash'))
    pool.base.save(direct/'protocol.json',dict(original=original))
    pool.base.save(direct/'report.json',dict(status='complete',full_catalog=True,
        scope='direct_contacts_only_not_final_screening_hits',targets={t:dict(distinct_molecules=1) for t in pool.TARGETS}))
    result=pool.prepare(direct,None,tmp_path/'mapped')
    assert all(v['provisional_tiers']=={'unresolved':1} for v in result['targets'].values())
