import pytest

from aidd_agent.chat_agent import ChatAgent
from aidd_agent.domain_tools import DomainTools
from aidd_agent.structure_chains import control,rows,tick,save
from test_structure_workflow_context import add


def start(app,sid,goal='budget'):
    return control(app,sid,'deepseek','Continue automatically',dict(operation='start',task_id='div',goal=goal,
        reference=dict(reference_query='6Q9L:HTZ:A:201',target_chain='A'),budget=dict(export_molecules=100000)))


def test_automatic_chain_waits_advances_once_and_stops_at_goal(tmp_path,monkeypatch):
    app=ChatAgent(tmp_path,start=False,allow_compute=True);sid=app.new_session()
    add(app,sid,'div','structure_diversity',dict(status='complete'))
    original=DomainTools.call;calls=[]
    def call(self,name,args):
        if name!='workflow':return original(self,name,args)
        decision=args['decision'];intent=decision['intent'];calls.append(decision)
        action={'consensus':'structure_consensus','recommend':'consensus_recommend','budget':'consensus_budget'}[intent]
        add(app,sid,intent,action,dict(status='complete',readiness='proposal_ready'),decision['task_id'])
        app.update(intent,status='queued')
        return dict(status='queued',tasks=[dict(id=intent,status='queued')])
    monkeypatch.setattr(DomainTools,'call',call)
    start(app,sid);tick(app);tick(app)
    assert len(calls)==1
    app.update('consensus',status='complete');tick(app)
    app.update('recommend',status='complete');tick(app)
    app.update('budget',status='complete');tick(app);tick(app)
    chain=rows(app,sid)[0]
    assert chain['status']=='complete' and len(calls)==3
    assert calls[-1]['budget']['export_molecules']==100000
    assert [c['task_id'] for c in calls]==['div','consensus','recommend']
    assert app.snapshot(sid)['structure_chains'][0]['current_task']=='budget'


def test_pause_block_restart_and_uncertain_dispatch(tmp_path,monkeypatch):
    app=ChatAgent(tmp_path,start=False);sid=app.new_session()
    add(app,sid,'div','structure_diversity',dict(status='complete'))
    identifier=start(app,sid,'recommendation')['chain_id']
    control(app,sid,'deepseek','pause',dict(operation='pause',chain_id=identifier));tick(app)
    assert len(app.jobs(sid))==1
    app.close();app=ChatAgent(tmp_path,start=False)
    assert rows(app,sid)[0]['status']=='paused'
    control(app,sid,'deepseek','resume',dict(operation='resume',chain_id=identifier))
    app.update('div',status='failed');tick(app)
    assert rows(app,sid)[0]['status']=='blocked'
    assert app.task(sid,'div')['status']=='failed'
    chain=rows(app,sid)[0];save(app,chain,'dispatching');app.close()
    app=ChatAgent(tmp_path,start=False);tick(app)
    assert rows(app,sid)[0]['status']=='uncertain' and len(app.jobs(sid))==1
    with pytest.raises(ValueError):control(app,sid,'deepseek','resume',dict(operation='resume',chain_id=identifier))


def test_permissions_duplicate_and_ownership(tmp_path):
    app=ChatAgent(tmp_path,start=False);sid=app.new_session()
    add(app,sid,'div','structure_diversity',dict(status='complete'))
    with pytest.raises(ValueError):start(app,sid)
    identifier=start(app,sid,'recommendation')['chain_id']
    with pytest.raises(ValueError):start(app,sid,'recommendation')
    with pytest.raises(ValueError):control(app,app.new_session(),'deepseek','pause',dict(operation='pause',chain_id=identifier))


def test_dispatch_failure_is_not_retried(tmp_path,monkeypatch):
    app=ChatAgent(tmp_path,start=False);sid=app.new_session()
    add(app,sid,'div','structure_diversity',dict(status='complete'))
    original=DomainTools.call;calls=[]
    def call(self,name,args):
        if name=='workflow':calls.append(args);raise OSError('lost receipt')
        return original(self,name,args)
    monkeypatch.setattr(DomainTools,'call',call)
    start(app,sid,'recommendation');tick(app);tick(app)
    assert len(calls)==1 and rows(app,sid)[0]['status']=='uncertain'


def test_coordinator_uses_real_validated_dispatch_and_sealed_plan(tmp_path):
    import json
    from pathlib import Path
    app=ChatAgent(tmp_path,start=False);sid=app.new_session()
    run_id='PROMPT-'+'a'*16
    add(app,sid,'div','structure_diversity',dict(status='complete'),run_id=run_id)
    start(app,sid,'recommendation');tick(app);tick(app)
    chain=rows(app,sid)[0]
    assert chain['status']=='active',chain
    assert len(app.jobs(sid))==2
    job=app.task(sid,chain['current_task'])
    plan=Path(job['plan']);envelope=json.loads(plan.read_text())
    step=envelope['plan']['steps'][0]
    assert step['action']=='structure_consensus'
    assert step['params']['source_run']==run_id
    assert step['params']['reference_query']=='6Q9L:HTZ:A:201'
    assert (plan.parent/'plan-seal.json').is_file()
