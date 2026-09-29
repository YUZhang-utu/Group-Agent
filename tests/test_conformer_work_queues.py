import csv
import json
import sqlite3

import pytest

from aidd_agent.conformer_work_queues import plan, sha, table


def fixture(tmp_path):
    model=tmp_path/'model';model.mkdir();diagnostic=tmp_path/'diagnostic';diagnostic.mkdir()
    groups=[('cis_R',2,'cis;trans'),('cis_S',3,'cis;trans'),('trans_R',4,'trans;trans'),
            ('uncertain_R',1,'boundary;trans'),('uncertain_S',2,'cis;boundary'),('large',101,'trans;trans')]
    rows=[]
    with sqlite3.connect(model/'blocks.sqlite') as db:
        db.execute('CREATE TABLE tree(node,group_id,n,axis)')
        for g,n,states in groups:
            leaves=[50,51] if g=='large' else [n]
            db.executemany('INSERT INTO tree VALUES(?,?,?,NULL)',[(g+str(i),g,m) for i,m in enumerate(leaves)])
            rows.append(dict(hard_group=g,conformers=n,blocks=len(leaves),boundary_positions=states.count('boundary'),omega_pattern=states))
    (model/'report.json').write_text('{}')
    table(diagnostic/'hard_groups.csv',rows,list(rows[0]))
    receipt=dict(status='complete',model=str(model),model_id='fixture',hard_groups=6,conformers=113,blocks=7,
        provenance=dict(model_report_sha256=sha(model/'report.json')),
        output_hashes={'hard_groups.csv':sha(diagnostic/'hard_groups.csv')})
    (diagnostic/'report.json').write_text(json.dumps(receipt))
    return model,diagnostic


def test_pools_preserve_classes_and_nested_budget(tmp_path):
    model,diagnostic=fixture(tmp_path);before=(model/'blocks.sqlite').read_bytes()
    output=tmp_path/'queues';r=plan(diagnostic,output,budgets=(1,3))
    assert r['source_conformers']==113 and r['logical_classes']==6 and r['evaluation_queues']==4
    assert r['logical_capacity_limit'] is None
    assert r['queue_kinds']['boundary_review_pool']['logical_classes']==2
    with (output/'classes.csv').open() as f: classes={r['hard_group']:r for r in csv.DictReader(f)}
    assert classes['cis_R']['queue_id']==classes['cis_S']['queue_id']
    assert classes['trans_R']['queue_id']!=classes['cis_R']['queue_id']
    assert classes['cis_R']['logical_class']!=classes['cis_S']['logical_class']
    with (output/'sampling_plan.csv').open() as f: rows=list(csv.DictReader(f))
    small={r['logical_class']:int(r['planned_conformer_slots']) for r in rows if r['budget']=='1'}
    large={r['logical_class']:int(r['planned_conformer_slots']) for r in rows if r['budget']=='3'}
    assert all(small[g]<=large[g] for g in small)
    assert sum(small.values())==4 and sum(large.values())==12
    assert any(r['coverage']=='unassessed_do_not_reject' for r in rows)
    with (output/'old_block_mapping.csv').open() as f: mapping=list(csv.DictReader(f))
    assert len(mapping)==7 and sum(int(r['conformers']) for r in mapping)==113
    assert before==(model/'blocks.sqlite').read_bytes()
    with pytest.raises(FileExistsError):plan(diagnostic,output)


@pytest.mark.parametrize('changed',['csv','database','receipt'])
def test_changed_inputs_are_rejected(tmp_path,changed):
    model,diagnostic=fixture(tmp_path)
    if changed=='csv':
        with (diagnostic/'hard_groups.csv').open('a') as f:f.write('changed')
    elif changed=='database':
        with sqlite3.connect(model/'blocks.sqlite') as db:db.execute('UPDATE tree SET n=n+1')
    else:(model/'report.json').write_text('{"changed":true}')
    with pytest.raises(ValueError):plan(diagnostic,tmp_path/'output')
    assert not (tmp_path/'output').exists()
