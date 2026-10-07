"""Independent search-arm and whole-block stopping fixtures, not library benchmarks."""
import sqlite3
from pathlib import Path

import pytest

from aidd_agent import block_campaign as campaign
from aidd_agent import block_conformer_search as search
from aidd_agent.domain_tools import validate_arguments
from aidd_agent.final_work_blocks import read
from test_block_campaign import chunks, pose, imported
from test_equiscore import engine


@pytest.mark.parametrize('budget,expected,exported',[(2,['A'],2),(4,['A','B'],4),(10,['A','B'],6)])
def test_whole_blocks_stop_at_budget_and_report_shortfall(tmp_path,monkeypatch,budget,expected,exported):
    database=tmp_path/'members.sqlite'
    with sqlite3.connect(database) as db:
        db.execute('CREATE TABLE member(gid INTEGER,block_id TEXT)')
        db.executemany('INSERT INTO member VALUES(?,?)',[(i,'A' if i<3 else 'B') for i in range(6)])
    selection=dict(method='chemplp',block_ids=['B','A'],rankings={'ChemPLP':[
        dict(block_id='B',rank=2),dict(block_id='A',rank=1)]})
    request=dict(selection=selection,conformers=budget)
    calls=[]
    def refine(batch,output,design,ids,workers,chunk):
        calls.append(ids.tolist());chunks(output,[pose(int(i)) for i in ids])
    monkeypatch.setattr('aidd_agent.budget_screen.run_refinement',refine)
    rows,ranked,searched,history=search.search_ranked_blocks(request,tmp_path,tmp_path,
        dict(templates=[dict(query_id='fixture')]),database,dict(workers=1,refine_chunk=2))
    assert [h['block_id'] for h in history]==expected
    assert len(rows)==exported and ranked==searched==3*len(expected)
    assert calls[0]==[0,1,2]  # Finish the block even when only two are requested.
    assert read(tmp_path/'search-history.json')['stop_policy']=='completed_block_boundary'


def test_compare_is_twelve_independent_search_only_arms(imported,monkeypatch):
    app,sid,project,path,_=imported
    monkeypatch.setattr(campaign,'owned_report',lambda app,sid,task,project,kinds:
        (path,dict(plan=str(project/'runs/query/plan.json'))))
    calls=[]
    original_queue=campaign.queue
    def queue(app,sid,project,request,action):
        calls.append(request)
        return original_queue(app,sid,project,request,action)
    monkeypatch.setattr(campaign,'queue',queue)
    monkeypatch.setattr(campaign,'select_blocks',lambda path,scheme,receptor,method,blocks:
        dict(scheme=scheme,receptor=receptor,method=method,blocks_per_method=blocks,block_ids=['A']))
    args=dict(operation='compare',task_id='a'*16,query_task_id='b'*16,receptor='R1')
    validate_arguments('block_campaign',args)
    result=campaign.handle(app,sid,args)
    assert len(result['tasks'])==len(calls)==12
    assert len({(r['selection']['scheme'],r['selection']['method'],r['selection']['blocks_per_method']) for r in calls})==12
    assert all(r['dock'] is False and r['conformers']==100000 and r['search_policy']=='ranked_blocks_until_budget' for r in calls)
    again=campaign.handle(app,sid,args)
    assert [t['task_id'] for t in again['tasks']]==[t['task_id'] for t in result['tasks']]
    assert all(t['reused'] for t in again['tasks'])
    summary=campaign.handle(app,sid,dict(operation='comparison'))
    assert len(summary['arms'])==12 and len(app.jobs(sid))==12
    with pytest.raises(ValueError):validate_arguments('block_campaign',dict(args,dock=True))


def test_sequential_union_is_rejected():
    with pytest.raises(ValueError,match='one scoring method'):
        campaign.validate(dict(operation='start',task_id='a'*16,query_task_id='b'*16,
            scheme='E095',receptor='R1',search_policy='ranked_blocks_until_budget',method='union'))


def test_chat_compare_returns_all_receipts_without_planner_resubmission(imported,monkeypatch):
    from aidd_agent.domain_tools import DomainTools
    app,sid,_,_,_=imported
    calls=[]
    def planner(*args,**kwargs):
        calls.append(1)
        return dict(kind='tool',tool='block_campaign',arguments=dict(operation='compare',
            task_id='a'*16,query_task_id='b'*16,receptor='R1'),purpose='Compare independent arms'),dict(provider='fixture')
    app.domain_planner=planner
    monkeypatch.setattr(DomainTools,'call',lambda *args:dict(status='queued',
        tasks=[dict(task_id=f'{i:016x}',status='queued') for i in range(12)],docking_requested=False))
    answer=app.ask(sid,'Compare partitions and scorers, search only','deepseek')
    assert 'no docking' in answer and '000000000000000b' in answer and len(calls)==1
