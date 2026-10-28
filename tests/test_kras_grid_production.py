"""E123 pool membership, checkpointing, quota and original-source export tests."""
import json
import hashlib
from pathlib import Path
import sqlite3
import subprocess
import sys
import pytest

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
import run_kras_grid_production as prod
from test_kras_weekend import tiny_library,PREPARED
from test_kras_retained_pool import RAW


def test_mapping_jobs_resume_and_quota(tmp_path):
    batch=tiny_library(tmp_path/'library')
    with sqlite3.connect(':memory:') as db:
        prod.schema(db)
        db.executemany('INSERT INTO members VALUES (?,?)',[(t,m) for t in prod.base.TARGETS for m in ('MD','M2','ABSENT')])
        prod.map_pool(db,batch)
        assert db.execute('SELECT gid FROM mapping ORDER BY gid').fetchall()==[(0,),(1,),(2,)]
        assert len(list(prod.jobs(db)))==3
        for mid,score in [('MD',.2),('M2',.4)]:
            row=dict(target='KRAS_G12D',mid=mid,profiles=dict(balanced=True,lenient=True,strict=False),
                     profile_poses=dict(balanced=dict(score=score)))
            prod.store_rows(db,[row])
        assert [r[0] for r in prod.selected(db,'KRAS_G12D',1)]==['M2']
        pending={m:ts for m,ts,_ in prod.jobs(db)}
        assert pending['MD']==['KRAS_G12V'] and len(pending['ABSENT'])==2
        prod.map_pool(db,batch)
        assert db.execute('SELECT COUNT(*) FROM mapping').fetchone()[0]==3


def test_cli_export_and_resume(tmp_path):
    batch=tiny_library(tmp_path/'library');plan=tmp_path/'plan';plan.mkdir()
    source=tmp_path/'source.mol2';source.write_text(RAW,encoding='utf-8')
    with sqlite3.connect(batch/'registry.sqlite3') as reg:
        reg.execute('CREATE TABLE conformer(id TEXT,molecule_id TEXT,source_path TEXT,source_record_index INTEGER,content_sha256 TEXT)')
        reg.executemany('INSERT INTO conformer VALUES (?,?,?,?,?)',[(cid,mid,str(source),0,hashlib.sha256(RAW.encode()).hexdigest()) for cid,mid in [('C0','MD'),('C1','MD'),('C2','M2')]])
    original=prod.soft.run.seal_inputs(batch,PREPARED,4,0)
    with sqlite3.connect(plan/'pool.sqlite') as db:
        db.execute('CREATE TABLE molecules(target TEXT,mid TEXT,PRIMARY KEY(target,mid))')
        db.executemany('INSERT INTO molecules VALUES (?,?)',[(t,m) for t in prod.base.TARGETS for m in ('MD','M2')])
    prod.base.save(plan/'report.json',dict(status='complete',scope='retained_pool_mapping',pool_sha256=prod.base.sha(plan/'pool.sqlite')))
    prod.base.save(plan/'protocol.json',dict(original=original))
    out=tmp_path/'production'
    cmd=[sys.executable,str(ROOT/'scripts/run_kras_grid_production.py'),'--plan',str(plan),'--output',str(out),'--workers','2','--cap','1','--batch-size','1']
    first=subprocess.run(cmd,capture_output=True,text=True)
    assert first.returncode==0,first.stdout+first.stderr
    report=prod.base.read(out/'report.json')
    assert report['status']=='complete' and report['new_scientific_rejections']==0
    assert all(v['processed']==2 for v in report['targets'].values())
    # Water source intentionally disagrees with native artifact coordinates.
    assert sum(v['quarantined'] for v in report['exports'].values())>=1
    assert sum(v['exported'] for v in report['exports'].values())==0
    for v in report['exports'].values():
        assert v['exported']+v['quarantined']==v['selected']<=1
        for entry in v['batches']:
            assert len(list(prod.base.iter_mol2_blocks(Path(entry['path']))))==entry['molecules']==1
    second=subprocess.run(cmd,capture_output=True,text=True)
    assert second.returncode==0,second.stdout+second.stderr
    assert prod.base.read(out/'report.json')['exports']==report['exports']
    altered=subprocess.run(cmd+['--cap','2'],capture_output=True,text=True)
    assert altered.returncode!=0 and 'mismatch' in altered.stderr


def test_export_splitting_and_source_hash(tmp_path,monkeypatch):
    batch=tmp_path/'library';batch.mkdir();out=tmp_path/'output';out.mkdir()
    source=tmp_path/'source.mol2';source.write_text(RAW*2,encoding='utf-8')
    with sqlite3.connect(batch/'registry.sqlite3') as reg:
        reg.execute('CREATE TABLE conformer(id TEXT,molecule_id TEXT,source_path TEXT,source_record_index INTEGER,content_sha256 TEXT)')
        reg.executemany('INSERT INTO conformer VALUES (?,?,?,?,?)',[(str(i),str(i),str(source),i,hashlib.sha256(RAW.encode()).hexdigest()) for i in range(2)])
    monkeypatch.setattr(prod.soft.run,'initialize',lambda *a:None)
    monkeypatch.setattr(prod,'source_check',lambda *a:True)
    with sqlite3.connect(':memory:') as db:
        prod.schema(db)
        for i in range(2):
            prod.store_rows(db,[dict(target='KRAS_G12D',mid=str(i),profiles=dict(balanced=True,lenient=True,strict=True),profile_poses=dict(balanced=dict(cid=str(i),score=.4)))])
        result=prod.export_target(db,'KRAS_G12D',batch,PREPARED,out,2,1)
        assert result['exported']==2 and len(result['batches'])==2
        for b in result['batches']:
            records=list(prod.base.iter_mol2_blocks(Path(b['path'])))
            assert len(records)==1
            assert records[0][1].splitlines()[2:]==RAW.splitlines()[2:]
        source.write_text(RAW.replace('-0.8','-0.7')*2,encoding='utf-8')
        second=tmp_path/'other';second.mkdir()
        with pytest.raises(ValueError,match='hash mismatch'):
            prod.export_target(db,'KRAS_G12D',batch,PREPARED,second,2,1)


def test_source_recheck_native_and_bad_coordinates(tmp_path):
    from rdkit import Chem
    batch=tiny_library(tmp_path/'library');prod.soft.run.initialize(batch,PREPARED)
    reader,_,states,chemistry=prod.soft.run.WORK
    candidate=reader.get(0);chem=chemistry.get(0)
    mol=Chem.RWMol()
    for number in chem.atomic_numbers: mol.AddAtom(Chem.Atom(int(number)))
    conf=Chem.Conformer(len(chem.atomic_numbers))
    for i,p in enumerate(candidate.shape_points): conf.SetAtomPosition(i,tuple(float(x) for x in p))
    mol.AddConformer(conf)
    pose=dict(gid=0,cid='C0',transform=__import__('numpy').eye(4).ravel().tolist())
    assert prod.source_check(mol,'KRAS_G12D','MD',pose)
    mol.GetConformer().SetAtomPosition(0,(1000.,1000.,1000.))
    with pytest.raises(ValueError,match='geometry mismatch'):
        prod.source_check(mol,'KRAS_G12D','MD',pose)
