"""E122 native geometry, sampled-ID scheduling and CLI checks."""
from pathlib import Path
import json
import sqlite3
import subprocess
import sys
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
import kras_pocket_grid as grid
import run_kras_grid_pilot as pilot
from test_kras_weekend import tiny_library,PREPARED


def test_native_controls_and_grid_accuracy():
    for target in pilot.base.TARGETS:
        state=grid.soft.run.kw.load_state(PREPARED,target);g=grid.make_grid(state)
        row=grid.controls(state,g)
        assert row['native_contacts'] and row['forced_overlap_detected']
        assert all(row['native_profiles'].values())
        assert row['grid_distance_max_error']<np.sqrt(3)*g['step']
        p=state['query']['shape_points']+1000
        _,_,inside=grid.query_grid(p,g)
        assert not inside.any()


def test_metrics_do_not_require_all_atoms_buried():
    row=grid.metrics(np.array([3.,3.,3.,3.]),np.array([True,True,False,False]))
    assert all(grid.profile_flags(row).values())
    clash=grid.metrics(np.array([.2,.2,.2,.2]),np.ones(4))
    assert not any(grid.profile_flags(clash).values())


def test_identity_only_mapping_and_sample(tmp_path):
    batch=tiny_library(tmp_path/'library')
    members={'MD':['KRAS_G12D'],'M2':['KRAS_G12V'],'ABSENT':['KRAS_G12D']}
    mapped=pilot.map_sample(batch,members,3)
    assert [v['gid'] for v in mapped['MD']]==[0,1]
    assert [v['gid'] for v in mapped['M2']]==[2]
    assert mapped['ABSENT']==[]
    assert all(v['gid']!=3 for rows in mapped.values() for v in rows)
    with sqlite3.connect(':memory:') as db:
        db.execute('CREATE TABLE molecules(target TEXT,mid TEXT,PRIMARY KEY(target,mid))')
        db.executemany('INSERT INTO molecules VALUES (?,?)',[(t,str(i)) for t in pilot.base.TARGETS for i in range(50)])
        a,totals=pilot.sample_members(db,10,122);b,_=pilot.sample_members(db,10,122)
        assert a==b
        assert all(sum(t in ts for ts in a.values())==10 for t in pilot.base.TARGETS)
        assert set(totals.values())=={50}


def test_pilot_cli_pool_only(tmp_path):
    batch=tiny_library(tmp_path/'library');plan=tmp_path/'plan';plan.mkdir()
    original=grid.soft.run.seal_inputs(batch,PREPARED,4,0)
    with sqlite3.connect(plan/'pool.sqlite') as db:
        db.execute('CREATE TABLE molecules(target TEXT,mid TEXT,PRIMARY KEY(target,mid))')
        db.executemany('INSERT INTO molecules VALUES (?,?)',[(t,mid) for t in pilot.base.TARGETS for mid in ('MD','M2')])
    pilot.base.save(plan/'report.json',dict(status='complete',scope='retained_pool_mapping',
                                          pool_sha256=pilot.base.sha(plan/'pool.sqlite')))
    pilot.base.save(plan/'protocol.json',dict(original=original))
    out=tmp_path/'pilot'
    result=subprocess.run([sys.executable,str(ROOT/'scripts/run_kras_grid_pilot.py'),
        '--plan',str(plan),'--output',str(out),'--sample','2','--workers','2'],capture_output=True,text=True)
    assert result.returncode==0,result.stdout+result.stderr
    report=pilot.base.read(out/'report.json')
    assert report['status']=='complete' and not report['production_selection_created']
    for t,row in report['targets'].items():
        assert row['sampled_molecules']==2 and row['evaluated_conformers']==3
        assert row['profiles']['lenient']['retained']>=1
    samples=pilot.base.read(out/'sample.json')
    assert not any(v['gid']==3 for rows in samples['conformers'].values() for v in rows)
    for line in (out/'molecules.jsonl').read_text().splitlines():
        row=json.loads(line)
        for name,passed in row['profiles'].items():
            if passed:assert grid.profile_flags(row['profile_poses'][name])[name]


def test_wilson_exposes_zero_hit_uncertainty():
    assert pilot.wilson(0,100)[1]>0
    assert pilot.wilson(100,100)[0]<1
