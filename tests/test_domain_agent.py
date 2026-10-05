import json
from pathlib import Path
import sqlite3

import pytest

from aidd_agent.chat_agent import ChatAgent
from aidd_agent.domain_agent import run,validate_step
from aidd_agent.domain_tools import DomainTools,page
from aidd_agent.prompt_workflow import project_root


def tool(name,**args):return dict(kind='tool',tool=name,arguments=args,purpose='Inspect actual evidence')
def answer(text='The observed result is exploratory.',ids=None):return dict(kind='answer',answer=text,evidence_ids=ids or [])


def scripted(steps,seen=None):
    sequence=iter(steps)
    def planner(prompt,profile,**kwargs):
        if seen is not None:seen.append(json.loads(prompt))
        step=next(sequence)
        return kwargs['validator'](step),dict(source='fixture_not_live_provider')
    return planner


def completed(app,sid):
    ctx=app.context;root=project_root(Path(ctx['db']),ctx['user_id'],ctx['project_id'])
    report=root/'result.json';child=root/'child.json'
    child.write_text(json.dumps(dict(matching_molecules=226,sampled_molecules=256,validation='not_run')))
    report.write_text(json.dumps(dict(steps=dict(screen=dict(status='complete',result=dict(report=str(child)))))))
    with app.connect() as db:
        db.execute('INSERT INTO jobs(id,session,status,report,request,created) VALUES(?,?,?,?,?,?)',('known-task',sid,'complete',str(report),'MDM2 pilot',1))
    return root


def test_compound_existing_results_and_literature_without_new_jobs(tmp_path,monkeypatch):
    seen=[]
    steps=[tool('tasks'),tool('task_report',task_id='known-task',step_id='screen'),
           tool('literature_search',query='MDM2 macrocycle'),answer('226 of 256 molecules passed; validation was not run. [E2]', ['E2','E3'])]
    monkeypatch.setattr(DomainTools,'fetch_json',staticmethod(lambda url:dict(hitCount=1,resultList=dict(result=[dict(source='MED',id='123',title='A test abstract',abstractText='Binding was measured.',pubYear='2024')]))))
    app=ChatAgent(tmp_path,start=False,domain_planner=scripted(steps,seen));sid=app.new_session();completed(app,sid)
    result=app.ask(sid,'Explain my results and compare with published research.','deepseek')
    assert '226 of 256' in result and 'https://europepmc.org/article/MED/123' in result
    assert len(app.jobs(sid))==1
    assert seen[-1]['evidence'][1]['result']['content']['value']['matching_molecules']==226
    trace=json.loads(next((app.root/'agents'/sid).glob('*.json')).read_text())
    assert trace['status']=='complete' and len(trace['steps'])==4
    assert app.snapshot(sid)['agent_runs'][0]['status']=='complete'


def test_owned_report_paths_and_discovered_artifacts(tmp_path):
    app=ChatAgent(tmp_path,start=False);sid=app.new_session();root=completed(app,sid)
    tools=DomainTools(app,sid,'inspect','deepseek')
    report=tools.call('task_report',dict(task_id='known-task'))
    aid=report['artifacts'][0]['artifact_id']
    assert tools.call('read_artifact',dict(artifact_id=aid,pointer='/sampled_molecules'))['content']['value']==256
    other=app.new_session()
    with pytest.raises(ValueError):DomainTools(app,other,'inspect','deepseek').call('task_report',dict(task_id='known-task'))
    outside=tmp_path/'outside.json';outside.write_text('{}')
    (root/'result.json').write_text(json.dumps(dict(steps=dict(screen=dict(result=dict(report=str(outside)))))))
    with pytest.raises(ValueError):tools.call('task_report',dict(task_id='known-task',step_id='screen'))
    with pytest.raises(KeyError):tools.call('read_artifact',dict(artifact_id=str(outside)))


def test_library_connection_and_exact_lookup(tmp_path):
    batch=tmp_path/'library';(batch/'artifacts').mkdir(parents=True)
    (batch/'artifacts/catalog.json').write_text(json.dumps(dict(library_id='lib',conformers=3,shards=[{}])))
    with sqlite3.connect(batch/'registry.sqlite3') as db:
        db.executescript('CREATE TABLE molecule(id,library_id,source_name); CREATE TABLE conformer(id,molecule_id,conformer_index,source_record_name,source_path,source_record_index,content_sha256,atom_count,bond_count);')
        db.execute('INSERT INTO molecule VALUES(?,?,?)',('MOL-1','lib','macrocycle-cis'))
        db.execute('INSERT INTO conformer VALUES(?,?,?,?,?,?,?,?,?)',('CNF-1','MOL-1',1,'macrocycle-cis-conf1','source.mol2',17,'digest',40,42))
    runtime=tmp_path/'runtime.json';runtime.write_text(json.dumps(dict(search=dict(batch=str(batch)))))
    app=ChatAgent(tmp_path/'app',runtime=runtime,start=False);sid=app.new_session();tools=DomainTools(app,sid,'inspect','deepseek')
    assert tools.call('library_status',{})['conformers']==3
    result=tools.call('molecule_lookup',dict(source_name='macrocycle-cis'))
    assert result['molecules'][0]['conformers'][0]['source_record_index']==17
    assert not tools.call('molecule_lookup',dict(source_name="' OR 1=1 --"))['molecules']


def test_pending_job_prevents_duplicate_or_additional_mutations(tmp_path):
    seen=[]
    decision=dict(intent='run',message='',task_id=None,request='Prepare a protein input')
    steps=[tool('workflow',decision=decision),tool('workflow',decision=decision),answer('Task queued; it has not completed.', ['E1'])]
    app=ChatAgent(tmp_path,start=False,domain_planner=scripted(steps,seen));sid=app.new_session()
    result=app.ask(sid,'Prepare a protein input','deepseek')
    assert len(app.jobs(sid))==1 and 'Execution receipts:' in result
    assert seen[-1]['pending_execution'] is True
    assert seen[-1]['evidence'][-1]['status']=='failed'
    trace=json.loads(next((app.root/'agents'/sid).glob('*.json')).read_text())
    assert trace['status']=='pending'


def test_unknown_citation_repaired_and_tool_error_visible(tmp_path):
    seen=[]
    app=ChatAgent(tmp_path,start=False,domain_planner=scripted([
        tool('task_report',task_id='unknown'),answer(ids=['E99']),answer('No matching task is available.')],seen))
    sid=app.new_session();app.ask(sid,'Explain the unknown task','deepseek')
    assert seen[-1]['previous_error']=='Answer cited unknown or failed evidence'
    assert seen[-1]['evidence'][0]['status']=='failed'
    assert not app.jobs(sid)


def test_budget_exhaustion_and_persistent_trace(tmp_path):
    app=ChatAgent(tmp_path,start=False);sid=app.new_session()
    result=run(app,sid,'inspect','deepseek',planner=scripted([tool('tasks')]),max_steps=1)
    assert 'step/time limit' in result
    assert app.snapshot(sid)['agent_runs'][0]['status']=='incomplete'
    app.close();reopened=ChatAgent(tmp_path,start=False)
    assert reopened.snapshot(sid)['agent_runs'][0]['status']=='incomplete'


def test_pointer_and_schema_validation():
    assert page({'a/b':[1,2,3]},dict(pointer='/a~1b',offset=1,limit=1))['items']==[2]
    with pytest.raises(ValueError):page({},dict(offset=-1))
    with pytest.raises(ValueError):validate_step(dict(kind='tool',tool='shell',arguments={},purpose='execute'))


def test_resumed_task_is_pending_with_existing_id(tmp_path):
    app=ChatAgent(tmp_path,start=False);sid=app.new_session();completed(app,sid)
    result=DomainTools(app,sid,'Resume my task','deepseek').call('workflow',dict(
        decision=dict(intent='resume',task_id='known-task',message='',request='')))
    assert result['status']=='queued'
    assert result['tasks'][0]['id']=='known-task'
    assert len(app.jobs(sid))==1


def test_invalid_answer_schema_has_actionable_feedback():
    with pytest.raises(ValueError, match='missing=evidence_ids'):
        validate_step(dict(kind='answer',answer='No task submitted.'))
    with pytest.raises(ValueError, match='extra=decision'):
        validate_step(dict(kind='answer',answer='Done',evidence_ids=[],decision={}))
    with pytest.raises(ValueError, match='1..16000'):
        validate_step(dict(kind='answer',answer='',evidence_ids=[]))


def test_receptor_stage_rejects_sampling_options_with_repair_feedback(tmp_path):
    seen=[]
    decision=dict(intent='block_evaluation',message='',task_id=None,request='',
                  evaluation=dict(stage='adopt_receptors',count=100,seed=20261002))
    app=ChatAgent(tmp_path,start=False,domain_planner=scripted([
        tool('workflow',decision=decision),answer('The request was rejected; no task was queued.')],seen))
    sid=app.new_session()
    app.ask(sid,'Prepare the accepted receptors only. Sampling later uses 100 per block.','deepseek')
    assert 'Stage adopt_receptors does not accept count, seed' in seen[-1]['previous_error']
    assert 'Do not change the requested stage' in seen[-1]['previous_error']
    assert not app.jobs(sid)


def test_invalid_plans_do_not_bypass_time_budget(tmp_path,monkeypatch):
    from aidd_agent import domain_agent
    clock=[0.0]; calls=[]
    def planner(*args,**kwargs):
        calls.append(1);clock[0]=181.0
        raise ValueError('Malformed LLM JSON response')
    app=ChatAgent(tmp_path,start=False);sid=app.new_session()
    monkeypatch.setattr(domain_agent.time,'time',lambda:clock[0])
    result=run(app,sid,'Inspect status only','deepseek',planner=planner)
    assert len(calls)==1
    assert 'step/time limit' in result and 'do not resubmit' in result
    assert not app.jobs(sid)


def test_operation_contracts_report_exact_fields_and_repair(tmp_path):
    from aidd_agent.domain_tools import argument_contracts, validate_arguments
    contract=argument_contracts()['block_results']['operations']['attach']
    assert contract['allowed_fields']==['operation','report']
    with pytest.raises(ValueError,match='unexpected=ranking_unit,top_n'):
        validate_arguments('block_results',dict(operation='attach',report='/report.json',top_n=10,ranking_unit='molecule'))
    with pytest.raises(ValueError,match='missing=attachment_id'):
        validate_arguments('block_results',dict(operation='analyze'))
    seen=[]
    app=ChatAgent(tmp_path,start=False,domain_planner=scripted([
        tool('block_results',operation='list',top_n=10),
        tool('block_results',operation='list'),answer('No attached results yet.', ['E2'])],seen))
    sid=app.new_session();app.ask(sid,'Inspect attached results only','deepseek')
    assert 'block_results.list invalid arguments' in seen[1]['previous_error']
    assert 'allowed=operation' in seen[1]['previous_error']
    assert seen[2]['evidence'][-1]['result']['analysis_tasks']==[]
    assert not app.jobs(sid)


def test_identical_invalid_call_stops_without_exhausting_ten_steps(tmp_path):
    seen=[]
    app=ChatAgent(tmp_path,start=False,domain_planner=scripted([
        tool('tasks',limit=5),tool('tasks',limit=5)],seen))
    sid=app.new_session()
    result=app.ask(sid,'Inspect tasks','deepseek')
    assert len(seen)==2
    assert 'same tool call failed twice' in result
    assert 'unexpected=limit' in result
    assert 'step/time limit' not in result
    assert not app.jobs(sid)
