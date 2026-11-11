import io
import json
from pathlib import Path

import pytest

from aidd_agent.chat_agent import ChatAgent
from aidd_agent.research_projects import handle, memory
from aidd_agent.chat_uploads import upload, attachment
from aidd_agent.prompt_workflow import create_plan, run_plan
from aidd_agent.prompt_smoke import OfflineServices


def create(app, sid, name='CIP2A'):
    return handle(app,sid,dict(operation='create',name=name,target=name,objective='Study this target'))


def test_project_memory_shared_across_chats_and_isolated_from_other_target(tmp_path):
    app=ChatAgent(tmp_path,start=False)
    try:
        old=app.new_session(); first=create(app,old); sid=first['session']
        other_chat=app.new_session(first['project_id'])
        second=create(app,old,'MDM2')['session']
        ctx=app.context_for(sid)
        plan=create_plan(Path(ctx['db']),ctx['user_id'],ctx['project_id'],'Prepare ACDEFGHIK',local_plan=dict(
            version=1,summary='Prepare input',clarifications=[],steps=[dict(id='input',action='af3_prepare',params=dict(name='protein',protein_sequences=['ACDEFGHIK']))]))
        report,path=run_plan(Path(ctx['db']),ctx['user_id'],ctx['project_id'],plan,services=OfflineServices())
        assert report['status']=='complete' and str(plan).startswith(ctx['project_root'])
        with app.connect() as db:
            db.execute('INSERT INTO jobs(id,session,status,request,report,created) VALUES(?,?,?,?,?,?)',('known',sid,'complete','Prepared protein',str(path),1))
        handle(app,sid,dict(operation='note',kind='decision',text='Use the full supplied sequence.',task_id='known'))
        data=upload(app,sid,'input.sdf',io.BytesIO(b'example'),7)
        assert attachment(app,other_chat,data['id'])['id']==data['id']
        assert app.task(other_chat,'known')['id']=='known'
        assert memory(app,other_chat)['notes'][0]['task_id']=='known'
        assert memory(app,other_chat)['total_tasks']==1
        with pytest.raises(ValueError): app.task(second,'known')
        with pytest.raises(ValueError): attachment(app,second,data['id'])
        assert memory(app,second)['total_tasks']==0
        app.message(sid,'user','Keep the full protein sequence for the comparison.')
        assert handle(app,other_chat,dict(operation='search',query='full',limit=1))['has_more']
        assert handle(app,second,dict(operation='search',query='full'))['matches']==[]
    finally: app.close()
    reopened=ChatAgent(tmp_path,start=False)
    try:
        assert memory(reopened,other_chat)['notes'][0]['text']=='Use the full supplied sequence.'
        assert reopened.context_for(sid)['project_id']==first['project_id']
    finally: reopened.close()


def test_legacy_sessions_migrate_and_explicit_projects_do_not_switch_each_other(tmp_path):
    from concurrent.futures import ThreadPoolExecutor
    app=ChatAgent(tmp_path,start=False);sid=app.new_session();default=app.context['project_id']
    with app.connect() as db:
        db.executescript('ALTER TABLE sessions DROP COLUMN project_id;')
    app.close();app=ChatAgent(tmp_path,start=False)
    try:
        assert app.context_for(sid)['project_id']==default
        a=create(app,sid,'Alpha');b=create(app,sid,'Beta')
        def prepare(session):
            ctx=app.context_for(session)
            plan=create_plan(Path(ctx['db']),ctx['user_id'],ctx['project_id'],'Protein input ACDEFGHIK',local_plan=dict(
                version=1,summary='Prepare',clarifications=[],steps=[dict(id='protein',action='af3_prepare',params=dict(name='protein',protein_sequences=['ACDEFGHIK']))]))
            result,path=run_plan(Path(ctx['db']),ctx['user_id'],ctx['project_id'],plan,services=OfflineServices())
            assert result['status']=='complete'
            return path
        with ThreadPoolExecutor(max_workers=2) as pool:
            results=list(pool.map(prepare,[a['session'],b['session']]))
        assert results[0].is_relative_to(app.context_for(a['session'])['project_root'])
        assert results[1].is_relative_to(app.context_for(b['session'])['project_root'])
        assert app.context['project_id']==default
    finally:app.close()


def test_project_http_selection_and_prompt_memory(tmp_path):
    import threading
    from urllib.request import Request,urlopen
    from aidd_agent.chat_web import make_server
    from test_domain_agent import scripted,tool,answer
    seen=[];app=ChatAgent(tmp_path,start=False,domain_planner=scripted([
        dict(kind='tool',tool='project',arguments=dict(operation='create',name='CIP2A',target='Q8TCG1',objective='Study CIP2A'),purpose='Create the requested target project'),
        tool('project',operation='memory'),answer('The project contains its recorded constraint.',['E1'])],seen))
    sid=app.new_session();events=[]
    reply=app.ask(sid,'Create a CIP2A project to study Q8TCG1.','deepseek',on_event=events.append)
    opened=next(e for e in events if e['type']=='project_opened');assert opened['project_id'] in reply
    handle(app,opened['session'],dict(operation='note',kind='constraint',text='Preserve the complete supplied sequence.'))
    server=make_server(app,0,token='project-fixture-token');thread=threading.Thread(target=server.serve_forever,daemon=True);thread.start()
    base='http://127.0.0.1:'+str(server.server_port)
    def post(path,data):
        with urlopen(Request(base+path,data=json.dumps(data).encode(),headers={'Authorization':'Bearer project-fixture-token','Content-Type':'application/json'})) as response:return json.load(response)
    try:
        new=post('/api/session',dict(project_id=opened['project_id']))['session']
        app.ask(new,'Recall the constraints in this project.','deepseek')
        assert seen[1]['project_memory']['notes'][0]['kind']=='constraint'
        other=post('/api/project',dict(session=new,arguments=dict(operation='create',name='MDM2',target='MDM2',objective='Separate target')))
        with urlopen(Request(base+'/api/state?session='+new,headers={'Authorization':'Bearer project-fixture-token'})) as response:state=json.load(response)
        assert state['project_id']==opened['project_id'] and len(state['project_memory']['notes'])==1
        assert all(s['project_id']==opened['project_id'] for s in state['sessions'])
        assert app.snapshot(other['session'])['project_memory']['notes']==[]
        from aidd_agent.research_workflows import handle as workflow
        created=workflow(app,new,dict(operation='create',objective='Prepare ethanol',nodes=[dict(id='prep',tool='molecule_prepare',arguments=dict(smiles='CCO'),depends_on=[])]),request='Prepare CCO')
        assert post('/api/research-workflow',dict(session=new,arguments=dict(operation='pause',workflow_id=created['workflow_id'])))['status']=='paused'
        from urllib.error import HTTPError
        with pytest.raises(HTTPError) as exc:
            post('/api/research-workflow',dict(session=other['session'],arguments=dict(operation='resume',workflow_id=created['workflow_id'])))
        assert exc.value.code==400 and not app.jobs(new)
    finally:server.shutdown();server.server_close();thread.join();app.close()


def test_new_experiment_can_reuse_recorded_inputs_from_project_task(tmp_path):
    app=ChatAgent(tmp_path,start=False)
    try:
        sid=app.new_session();request='Protein sequence:\nACDEFGHIK\n\nSMILES: CCO'
        app.execute_decision(sid,request,'deepseek',dict(intent='run',message='',task_id=None,request=request))
        job=app.task(sid);app.update(job['id'],status='complete')
        other=app.new_session()
        decision=dict(intent='run',message='',task_id=job['id'],request='Prepare the same ACDEFGHIK protein and CCO ligand with seeds [47, 48].')
        app.execute_decision(other,'Reuse the inputs from task '+job['id']+' with seeds [47, 48].','deepseek',decision)
        new=app.task(other)
        assert new['id']!=job['id'] and request in new['request'].split('Original user request (preserve literal SMILES):')[1]
        foreign=create(app,sid,'Other target')['session']
        with pytest.raises(ValueError):app.execute_decision(foreign,'Reuse inputs','deepseek',decision)
    finally:app.close()
