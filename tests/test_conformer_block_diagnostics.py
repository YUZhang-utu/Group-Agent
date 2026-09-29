import json
import sqlite3

import pytest

from aidd_agent.conformer_block_diagnostics import run


def test_attributes_capacity_and_boundary_without_rebuilding(tmp_path):
    model=tmp_path/'model';model.mkdir();descriptors=tmp_path/'descriptors.sqlite'
    receipt=dict(status='complete',model_id='fixture',capacity=2,counts=dict(conformers=5,blocks=3))
    (model/'report.json').write_text(json.dumps(receipt))
    with sqlite3.connect(model/'blocks.sqlite') as db:
        db.execute('CREATE TABLE tree(node,group_id,path,n,axis)')
        db.executemany('INSERT INTO tree VALUES(?,?,?,?,?)',[
            ('a','g1','',4,0),('al','g1','0',2,None),('ar','g1','1',2,None),('b','g2','',1,None)])
        db.execute('CREATE TABLE point(cid,group_id)')
        db.executemany('INSERT INTO point VALUES(?,?)',[('one','g1'),('two','g2')])
    with sqlite3.connect(descriptors) as db:
        db.execute('CREATE TABLE descriptor(cid PRIMARY KEY,payload)')
        db.executemany('INSERT INTO descriptor VALUES(?,?)',[
            ('one',json.dumps(dict(hard_group='g1',provenance=dict(omega_states=['cis','trans'])))),
            ('two',json.dumps(dict(hard_group='g2',provenance=dict(omega_states=['boundary','trans']))))])
    before=(model/'blocks.sqlite').read_bytes(),descriptors.read_bytes()
    result=run(model,descriptors,tmp_path/'diagnostic')
    assert result['blocks_without_capacity_limit']==2
    assert result['capacity_extra_blocks']==1
    assert result['boundary_attribution']['with_boundary']['conformers']==1
    assert result['boundary_attribution']['without_boundary']['capacity_extra_blocks']==1
    assert result['maximum_uncapped_group_size']==4
    assert before==((model/'blocks.sqlite').read_bytes(),descriptors.read_bytes())
    with pytest.raises(FileExistsError):run(model,descriptors,tmp_path/'diagnostic')
    with pytest.raises(ValueError,match='separate'):run(model,descriptors,model/'diagnostic')
