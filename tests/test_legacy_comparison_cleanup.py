"""Cleanup is restricted to a complete legacy cohort; upstream inputs survive."""
import hashlib
import importlib.util
import itertools
import json
from pathlib import Path
import sqlite3

import pytest

spec=importlib.util.spec_from_file_location('cleanup',Path(__file__).parents[1]/'scripts/clear_legacy_block_comparison.py')
cleanup=importlib.util.module_from_spec(spec)
spec.loader.exec_module(cleanup)


def fixture(root):
    project=root/'users/u/projects/p'
    (root/'chat').mkdir()
    (project/'block-campaigns').mkdir(parents=True)
    query=project/'query.json';query.write_text('{}')
    analysis=project/'analysis.json';analysis.write_text('{}')
    with sqlite3.connect(root/'chat/conversations.sqlite3') as db:
        db.execute('CREATE TABLE jobs(id TEXT,session TEXT,request TEXT,plan TEXT,status TEXT,cancel INTEGER,error TEXT)')
        for i,(scheme,method,blocks) in enumerate(itertools.product(['E094','E095','E096'],['chemplp','equiscore'],[5,10])):
            request=dict(session='s',analysis=str(analysis),analysis_sha256='a',query=str(query),query_sha256='q',
                         conformers=100000,dock=False,search_policy='ranked_blocks_until_budget',
                         selection=dict(scheme=scheme,method=method,blocks_per_method=blocks,receptor='R'))
            payload=json.dumps(request).encode();rid=hashlib.sha256(payload).hexdigest()[:24]
            (project/'block-campaigns'/(rid+'.json')).write_bytes(payload)
            run=project/'runs'/f'PROMPT-{i:016x}';run.mkdir(parents=True)
            plan=run/'plan.json'
            plan.write_text(json.dumps(dict(plan=dict(steps=[dict(action='block_search_dock',params=dict(request_id=rid,request_sha256=hashlib.sha256(payload).hexdigest()))]))))
            db.execute('INSERT INTO jobs VALUES(?,?,?,?,?,?,?)',(str(i),'s','Search leading blocks for conformers '+rid,str(plan),'running' if i==0 else 'queued',0,None))
    return project


def test_inventory_and_apply_preserve_inputs(tmp_path,monkeypatch):
    project=fixture(tmp_path)
    report=cleanup.inventory(tmp_path,'0')
    assert len(report['tasks'])==12
    monkeypatch.setattr('sys.argv',['cleanup','--storage-root',str(tmp_path),'--anchor-task','0','--apply'])
    monkeypatch.setattr(cleanup,'require_stopped',lambda:None)
    cleanup.main()
    assert not list((project/'runs').iterdir())
    assert (project/'query.json').exists() and (project/'analysis.json').exists()
    with sqlite3.connect(tmp_path/'chat/conversations.sqlite3') as db:
        assert db.execute('SELECT count(*) FROM jobs').fetchone()[0]==0
    assert list((tmp_path/'maintenance').glob('*/inventory.json'))


def test_incomplete_cohort_is_refused(tmp_path):
    fixture(tmp_path)
    with sqlite3.connect(tmp_path/'chat/conversations.sqlite3') as db:
        db.execute("DELETE FROM jobs WHERE id='11'")
    with pytest.raises(ValueError,match='twelve-arm'):
        cleanup.inventory(tmp_path,'0')


def test_modified_request_seal_is_refused(tmp_path):
    project=fixture(tmp_path)
    next((project/'block-campaigns').glob('*.json')).write_text('{}')
    with pytest.raises(ValueError,match='seal'):
        cleanup.inventory(tmp_path,'0')
