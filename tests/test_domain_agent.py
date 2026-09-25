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
