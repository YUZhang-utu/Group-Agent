"""Reproduce the CIP2A planning failure without a live model or GPU execution."""
import copy
import io
import json
from pathlib import Path
import time
from urllib.error import HTTPError

import pytest

from aidd_agent import prompt_plan
from aidd_agent.chat_agent import ChatAgent
from aidd_agent.domain_tools import DomainTools
from aidd_agent.prompt_workflow import create_plan, initialize_context, project_root, run_plan
from aidd_agent.prompt_smoke import OfflineServices

PROMPT = (Path(__file__).parent/'fixtures/cip2a_af3_request.txt').read_text()
SEQUENCE = PROMPT.split('Protein sequence:\n')[1].split('\n')[0]
SMILES = 'c1cc4OCOc4(cc1CN2CCN(CC2)C(Nc3ccc(cc3)Cl)=S)'
PLAN = dict(version=1,summary='Prepare and predict the supplied CIP2A complex',clarifications=[],steps=[
    dict(id='input',action='af3_prepare',params=dict(name='CIP2A_ligand_seeds47_51',
         protein_sequences=[SEQUENCE],ligand_smiles=[SMILES],seeds=[47,48,49,50,51])),
    dict(id='predict',action='af3_run',params=dict(input_step='input'))])


def provider(monkeypatch, values):
    calls=[]; outputs=iter(values)
    class Opener:
        def open(self, request, timeout):
            calls.append(json.loads(request.data))
            value=next(outputs)
            if isinstance(value,Exception):raise value
            content=json.dumps(value) if not isinstance(value,str) else value
            return io.BytesIO(json.dumps(dict(id='fixture',choices=[dict(finish_reason='stop',message=dict(content=content))])).encode())
    monkeypatch.setattr(prompt_plan,'build_opener',lambda *a:Opener())
    return calls


def config(tmp_path):
    directory=tmp_path/'profiles';directory.mkdir()
    profile=directory/'deepseek.json'
    profile.write_text(json.dumps(dict(base_url='http://127.0.0.1:1/v1',model='fixture')))
    return profile


def test_cip2a_invalid_fields_repaired_before_preparation(tmp_path,monkeypatch):
    calls=provider(monkeypatch,[dict(plan=PLAN),PLAN])
    ctx=initialize_context(tmp_path/'storage');profile=config(tmp_path)
    path=create_plan(Path(ctx['db']),ctx['user_id'],ctx['project_id'],PROMPT,llm_profile=profile)
    envelope=json.loads(path.read_text())
    assert len(calls)==2 and envelope['model']['planning_attempts']==2
    assert 'missing=clarifications,steps,summary,version' in calls[1]['messages'][0]['content']
    assert calls[0]['messages'][1]['content']==calls[1]['messages'][1]['content']==PROMPT
    result,_=run_plan(Path(ctx['db']),ctx['user_id'],ctx['project_id'],path,services=OfflineServices())
    assert result['status']=='blocked' and result['steps']['input']['status']=='complete'
    payload=json.loads((path.parent/'execution/input/af3-input.json').read_text())
    assert payload['modelSeeds']==[47,48,49,50,51]
    assert payload['sequences']==[dict(protein=dict(id='A',sequence=SEQUENCE)),dict(ligand=dict(id='B',smiles=SMILES))]
    assert not (path.parent/'execution/predict/command.json').exists()


def test_bad_plan_exhaustion_creates_no_executable_run(tmp_path,monkeypatch):
    calls=provider(monkeypatch,[{},'not json',{}])
    ctx=initialize_context(tmp_path/'storage');profile=config(tmp_path)
    with pytest.raises(prompt_plan.LLMResponseError,match='after 3 validated attempts'):
        create_plan(Path(ctx['db']),ctx['user_id'],ctx['project_id'],PROMPT,llm_profile=profile)
    root=project_root(Path(ctx['db']),ctx['user_id'],ctx['project_id'])
    assert len(calls)==3 and not list(root.glob('runs/PROMPT-*/plan.json'))


def test_auth_failure_does_not_retry_or_echo_provider_body(tmp_path,monkeypatch):
    calls=provider(monkeypatch,[HTTPError('http://127.0.0.1',401,'PRIVATE',{},io.BytesIO(b'PRIVATE'))])
    ctx=initialize_context(tmp_path/'storage');profile=config(tmp_path)
    with pytest.raises(prompt_plan.LLMRequestError,match='401') as error:
        create_plan(Path(ctx['db']),ctx['user_id'],ctx['project_id'],PROMPT,llm_profile=profile)
    assert len(calls)==1 and 'PRIVATE' not in str(error.value)


def test_literal_corruption_is_repaired_without_accepting_generated_chemistry(tmp_path,monkeypatch):
    wrong=copy.deepcopy(PLAN);wrong['steps'][0]['params']['protein_sequences']=[SEQUENCE[:-1]+'A']
    calls=provider(monkeypatch,[wrong,PLAN])
    ctx=initialize_context(tmp_path/'storage');profile=config(tmp_path)
    path=create_plan(Path(ctx['db']),ctx['user_id'],ctx['project_id'],PROMPT,llm_profile=profile)
    assert len(calls)==2 and 'original user request' in calls[1]['messages'][0]['content']
    assert json.loads(path.read_text())['plan']==PLAN


@pytest.mark.parametrize('change', ['seed_list', 'shortened_sequence', 'invented_construct'])
def test_explicit_cip2a_fields_cannot_be_changed_during_planning(tmp_path,monkeypatch,change):
    wrong=copy.deepcopy(PLAN);params=wrong['steps'][0]['params']
    if change=='seed_list':params['seeds']=[1]
    elif change=='shortened_sequence':params['protein_sequences']=[SEQUENCE[:400]]
    else:params.update(start=1,end=400)
    calls=provider(monkeypatch,[wrong,PLAN])
    ctx=initialize_context(tmp_path/'storage');profile=config(tmp_path)
    path=create_plan(Path(ctx['db']),ctx['user_id'],ctx['project_id'],PROMPT,llm_profile=profile)
    assert len(calls)==2 and json.loads(path.read_text())['plan']==PLAN


def test_current_request_can_override_previous_seed_list():
    from aidd_agent.af3_inputs import literal_check
    params=copy.deepcopy(PLAN['steps'][0]['params']);params['seeds']=[9]
    literal_check(params,PROMPT+'\nCurrent user request:\nRepeat with model seeds: [9].')


def test_background_planning_failure_is_visible_without_scientific_report(tmp_path,monkeypatch):
    monkeypatch.setenv('TEST_API_KEY','private-test-value')
    def failing(*a,**k):raise ValueError('Invalid plan fields private-test-value')
    app=ChatAgent(tmp_path/'chat',planner=failing,router=lambda *a:dict(intent='run',request=PROMPT,task_id=None,message=''))
    try:
        sid=app.new_session();app.ask(sid,PROMPT,'deepseek')
        deadline=time.monotonic()+10
        while time.monotonic()<deadline and app.task(sid)['status']!='failed':time.sleep(.05)
        job=app.task(sid)
        assert job['status']=='failed' and job['plan'] is None
        assert 'planning: Invalid plan fields' in job['error'] and 'private-test-value' not in job['error']
        diagnostic=json.loads((app.root/(job['id']+'.failure.json')).read_text())
        assert diagnostic['phase']=='planning'
        assert 'Invalid plan fields' in Path(job['log']).read_text()
        state=app.snapshot(sid)
        # The worker appends the notification just after persisting the terminal state.
        while time.monotonic()<deadline and not any('failed during planning' in m['text'] for m in state['messages']):
            time.sleep(.02);state=app.snapshot(sid)
        assert any('No scientific engine was launched' in m['text'] for m in state['messages'])
        evidence=DomainTools(app,sid,'Inspect status','deepseek').call('task_report',dict(task_id=job['id']))
        assert evidence['error']==job['error'] and evidence['log']==job['log']
    finally:app.close()


def test_real_chat_queue_repairs_plan_and_retains_prepared_cip2a_input(tmp_path,monkeypatch):
    calls=provider(monkeypatch,[dict(job_name='wrong AF3 envelope'),PLAN])
    profile=config(tmp_path);decisions=[]
    def conversation(*a,**k):
        decisions.append(1)
        return dict(kind='tool',tool='workflow',purpose='Submit the supplied AF3 request',arguments=dict(
            decision=dict(intent='run',request=PROMPT,task_id=None,message=''))),{}
    app=ChatAgent(tmp_path/'chat',config_dir=profile.parent,start=False,domain_planner=conversation)
    try:
        sid=app.new_session();reply=app.ask(sid,PROMPT,'deepseek')
        original_id=app.task(sid)['id']
        assert 'not confirmed' in reply and len(decisions)==1
        app.execute(app.task(sid))
        job=app.task(sid)
        assert job['id']==original_id and len(app.jobs(sid))==1
        assert job['status']=='blocked', app.snapshot(sid)
        assert 'Compute not enabled' in job['error']
        report=json.loads(Path(job['report']).read_text())
        path=Path(report['steps']['input']['result']['input_json'])
        payload=json.loads(path.read_text())
        assert payload['modelSeeds']==[47,48,49,50,51]
        assert payload['sequences'][0]['protein']['sequence']==SEQUENCE
        assert payload['sequences'][1]['ligand']['smiles']==SMILES
        evidence=DomainTools(app,sid,'Inspect saved input','deepseek').call('task_report',dict(task_id=job['id']))
        assert any(row['path']==str(path) for row in evidence['artifacts'])
        assert str(path) in app.ask(sid,'/results '+job['id'],'deepseek')
        assert len(calls)==2
    finally:app.close()


def test_resume_failed_planning_reuses_original_task_id(tmp_path,monkeypatch):
    calls=provider(monkeypatch,[{},{},{},PLAN]);profile=config(tmp_path)
    app=ChatAgent(tmp_path/'chat',config_dir=profile.parent,router=lambda *a:dict(intent='run',request=PROMPT,task_id=None,message=''))
    def wait_status(status):
        deadline=time.monotonic()+30
        while time.monotonic()<deadline:
            if app.task(sid)['status']==status:return app.task(sid)
            time.sleep(.05)
        pytest.fail(str(app.snapshot(sid)))
    try:
        sid=app.new_session();app.ask(sid,PROMPT,'deepseek')
        failed=wait_status('failed');identifier=failed['id']
        assert failed['plan'] is None and 'after 3 validated attempts' in failed['error']
        reply=app.ask(sid,'/resume '+identifier,'deepseek')
        assert identifier in reply
        recovered=wait_status('blocked')
        assert recovered['id']==identifier and len(app.jobs(sid))==1 and len(calls)==4
        assert 'Compute not enabled' in recovered['error']
        report=json.loads(Path(recovered['report']).read_text())
        assert Path(report['steps']['input']['result']['input_json']).is_file()
    finally:app.close()
