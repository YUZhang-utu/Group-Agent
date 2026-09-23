import json
from pathlib import Path

import pytest

from aidd_agent.chat_agent import ChatAgent,validate_route,screening_feedback
from aidd_agent.budget_workflow import validate_budget,execute


def test_budget_route_and_limits():
    route=dict(intent='budget',message='',request='',task_id=None,budget={})
    assert validate_route(route)==route
    for value in ({'retrieval_molecules':True},{'export_molecules':1000001},{'path':'unsafe'}):
        with pytest.raises(ValueError):validate_route(dict(route,budget=value))
    with pytest.raises(ValueError):validate_route(dict(route,request='run'))
    assert validate_budget({'start_rank':100001},page=True)


def test_chat_budget_queues_owned_search_and_next_page(tmp_path):
    app=ChatAgent(tmp_path,start=False,router=lambda *args:dict(intent='run',request='Survey human MDM2',task_id=None,message=''))
    sid=app.new_session();app.ask(sid,'Survey human MDM2','deepseek')
    job=app.task(sid);folder=tmp_path/'PROMPT-0123456789abcdef';folder.mkdir()
    plan=folder/'plan.json';plan.write_text('{}')
    report=folder/'report.json'
    child=folder/'child.json';child.write_text(json.dumps(dict(kind='consensus_recommendation')))
    report.write_text(json.dumps({'steps':{'guided':{'action':'consensus_recommend','status':'complete','result':{'report':str(child)}}}}))
    app.update(job['id'],status='complete',plan=str(plan),report=str(report))
    app.ask(sid,'/budget','deepseek');queued=app.task(sid)
    step=json.loads(Path(queued['plan']).read_text())['plan']['steps'][0]
    assert step['action']=='consensus_budget'
    assert step['params']==dict(source_run=folder.name)
    with pytest.raises(ValueError):app.ask(app.new_session(),'/budget '+job['id'],'deepseek')
    app.router=lambda *args:dict(intent='budget',message='',request='',task_id=job['id'],budget=dict(retrieval_molecules=500000,export_molecules=100000,chunk_conformers=4048))
    app.ask(sid,'Search and deliver the top 100000 molecules','deepseek')
    natural=json.loads(Path(app.task(sid)['plan']).read_text())['plan']['steps'][0]
    assert natural['action']=='consensus_budget' and natural['params']['retrieval_molecules']==500000
    assert natural['params']['chunk_conformers']==4048
    report.write_text(json.dumps(dict(steps=dict(guided=dict(action='consensus_budget',status='complete')))))
    app.update(queued['id'],status='complete',plan=str(plan),report=str(report))
    app.ask(sid,'/budget_page','deepseek')
    step=json.loads(Path(app.task(sid)['plan']).read_text())['plan']['steps'][0]
    assert step['action']=='budget_page'
    assert 'without rescoring' in screening_feedback(dict(kind='budget_page'))
    app.close()


def test_adapter_enforces_compute_and_exports_next_page(tmp_path,monkeypatch):
    from aidd_agent import budget_screen,budget_export
    from aidd_agent.prompt_workflow import Blocked
    source=tmp_path/'source.json';run=tmp_path/'run';run.mkdir()
    (run/'report.json').write_text(json.dumps(dict(ranked_molecules=300000,ranking='contact-first')))
    from aidd_agent.expanded_wee1 import fingerprint
    source.write_text(json.dumps(dict(kind='consensus_budget',sources=fingerprint([run/'report.json']),run_directory=str(run),end_rank=100000,target={'accession':'OTHER'})))
    calls=[]
    def export(batch,path,out,start,count):
        calls.append((start,count));out.mkdir()
        result=dict(exported_molecules=count,end_rank=start+count-1,shortfall=0,review_required=0)
        (out/'report.json').write_text(json.dumps(result));return result
    monkeypatch.setattr(budget_export,'export',export)
    monkeypatch.setattr(budget_screen,'execute',lambda _:pytest.fail('Page must not rescore'))
    cfg=dict(search=dict(batch=str(tmp_path)))
    with pytest.raises(Blocked):execute('budget_page',source,tmp_path/'blocked',{},cfg,False)
    result=execute('budget_page',source,tmp_path/'out',{},cfg,True)
    assert calls==[(100001,100000)] and result['end_rank']==200000
    with pytest.raises(ValueError,match='overlaps'):execute('budget_page',source,tmp_path/'overlap',dict(start_rank=1),cfg,True)


def test_budget_search_is_target_generic(tmp_path,monkeypatch):
    from aidd_agent import budget_screen,budget_export
    from aidd_agent.expanded_wee1 import fingerprint
    evidence=tmp_path/'evidence';evidence.write_text('sealed')
    source=tmp_path/'source.json'
    source.write_text(json.dumps(dict(kind='consensus_design',sources=fingerprint([evidence]),
        target={'accession':'OTHER'},reference='OTHER:frame',anchors=[],
        design=dict(mandatory_anchors=[],alternative_groups=[],optional_weights={},evidence_ids=[],
            optional_normalization='fixed_budget',optional_budget=1))))
    calls=[]
    def screen(args):
        calls.append(args);args.output.mkdir()
        (args.output/'report.json').write_text(json.dumps(dict(ranked_molecules=1000000,ranking='contact-first')))
    def export(batch,run,out,start,count):
        out.mkdir();result=dict(exported_molecules=count,end_rank=count,shortfall=0,review_required=0)
        (out/'report.json').write_text(json.dumps(result));return result
    monkeypatch.setattr(budget_screen,'execute',screen);monkeypatch.setattr(budget_export,'export',export)
    monkeypatch.setenv('AIDD_ASSIGNMENT_BACKEND','python')
    result=execute('consensus_budget',source,tmp_path/'out',dict(chunk_conformers=4048),dict(search=dict(batch=str(tmp_path))),True)
    assert result['exported_molecules']==100000 and result['target']['accession']=='OTHER'
    assert calls[0].definitions is None and calls[0].retrieval_molecules==1000000
    assert calls[0].chunk_conformers==4048
