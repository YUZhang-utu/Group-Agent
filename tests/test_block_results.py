"""Completed CLI adoption: sealed reuse, session ownership and no engine dispatch."""
import json
from pathlib import Path

import pytest

from aidd_agent.block_results import handle, analyze_attachment
from aidd_agent.chat_agent import ChatAgent, screening_summary
from aidd_agent.domain_tools import DomainTools
from aidd_agent.prompt_workflow import project_root, run_plan
from aidd_agent.final_work_blocks import sha
from test_block_evaluation import fixture, plants_profile, fake_plants
from aidd_agent.block_sampling import sample
from aidd_agent.block_plants import prepare, run


@pytest.fixture
def completed(tmp_path, monkeypatch):
    app = ChatAgent(tmp_path / 'workspace', start=False)
    sid = app.new_session()
    ctx = app.context
    project = project_root(Path(ctx['db']), ctx['user_id'], ctx['project_id'])
    root = project / 'runs/cli-fixture'
    root.mkdir(parents=True)
    cfg, source = fixture(root)
    sample(cfg, root / 'sample', count=100)
    prepare(root / 'sample/report.json', plants_profile(root, source), root / 'prepared')
    with monkeypatch.context() as engine:
        engine.setattr('aidd_agent.block_plants.subprocess.run', fake_plants)
        run(root / 'prepared/report.json', root / 'docking')
    yield app, sid, project, root / 'docking/report.json'
    app.close()


def test_attach_queue_analyze_preserves_source_and_never_calls_engine(completed, monkeypatch):
    app, sid, project, source = completed
    original = sha(source)
    monkeypatch.setattr('aidd_agent.prompt_workflow.platform.platform', lambda: 'synthetic-test-platform')
    def forbidden(*args, **kwargs):
        raise AssertionError('Analysis must never execute an engine')
    monkeypatch.setattr('aidd_agent.block_plants.subprocess.run', forbidden)
    first = handle(app, sid, dict(operation='attach', report=str(source)))
    assert first == handle(app, sid, dict(operation='attach', report=str(source)))
    assert app.jobs(sid) == []
    identifier = first['attachment']['id']
    queued = handle(app, sid, dict(operation='analyze', attachment_id=identifier))
    assert handle(app, sid, dict(operation='analyze', attachment_id=identifier))['reused']
    job = app.task(sid, queued['tasks'][0]['id'])
    ctx = app.context
    result, _ = run_plan(Path(ctx['db']), ctx['user_id'], ctx['project_id'], Path(job['plan']))
    assert result['status'] == 'complete', result
    assert sha(source) == original
    child = json.loads((Path(job['plan']).parent / 'execution/analysis/blocks/report.json').read_text())
    assert child['kind'] == 'block_analyze'
    assert len(child['comparisons']) == 3
    assert child['ligand_chemistry_reviewed'] is True
    assert all(r['estimate_status'] == 'complete_panel' for r in child['comparisons'])
    assert 'block_summary.csv' in child['outputs']
    assert child['ranking']['ranking_unit'] == 'molecule'
    assert child['ranking']['top_n'] == 10
    assert screening_summary(job)['kind'] == 'block_analyze'


def test_ownership_wrong_report_and_changed_report(completed, tmp_path):
    app, sid, project, source = completed
    with pytest.raises(ValueError, match='block_plants_run'):
        handle(app, sid, dict(operation='attach', report=str(source.parent.parent / 'prepared/report.json')))
    with pytest.raises(ValueError, match='escapes'):
        handle(app, sid, dict(operation='attach', report=str(tmp_path / 'report.json')))
    row = handle(app, sid, dict(operation='attach', report=str(source)))['attachment']
    other = app.new_session()
    assert handle(app, other, dict(operation='list'))['attachments'] == []
    with pytest.raises(ValueError, match='conversation'):
        handle(app, other, dict(operation='analyze', attachment_id=row['id']))
    source.write_text(source.read_text() + ' ')
    with pytest.raises(ValueError, match='changed'):
        handle(app, sid, dict(operation='analyze', attachment_id=row['id']))


def test_tool_requires_literal_path_and_status_is_readonly(completed):
    app, sid, project, source = completed
    tools = DomainTools(app, sid, 'Show results', 'deepseek')
    with pytest.raises(ValueError, match='exact report path'):
        tools.call('block_results', dict(operation='attach', report=str(source)))
    tools = DomainTools(app, sid, 'Attach ' + str(source), 'deepseek')
    tools.call('block_results', dict(operation='attach', report=str(source)))
    assert len(tools.call('block_results', {})['attachments']) == 1
    assert not app.jobs(sid)


def test_modified_pose_is_rejected_before_statistics(completed):
    app, sid, project, source = completed
    row = handle(app, sid, dict(operation='attach', report=str(source)))['attachment']
    entry = project / 'block-results' / (row['id'] + '.json')
    pose = next(source.parent.glob('*/attempt-*/poses/docked_ligands.mol2'))
    pose.write_text(pose.read_text() + '\n')
    with pytest.raises(ValueError, match='changed'):
        analyze_attachment(project, dict(attachment_id=row['id'], attachment_sha256=sha(entry)), project / 'analysis-bad')


def test_partial_panel_remains_partial(completed):
    app, sid, project, source = completed
    report = json.loads(source.read_text())
    report.update(status='partial', failed_jobs=1)
    source.write_text(json.dumps(report))
    row = handle(app, sid, dict(operation='attach', report=str(source)))['attachment']
    assert row['summary']['status'] == 'partial'


def test_slash_commands_attach_without_provider(completed):
    app, sid, project, source = completed
    answer = json.loads(app.ask(sid, '/dock_attach ' + str(source), 'deepseek'))
    assert answer['status'] == 'attached'
    assert json.loads(app.ask(sid, '/dock_results', 'deepseek'))['attachments']
    assert not app.jobs(sid)


def test_natural_language_tool_loop_queues_one_analysis(completed):
    app, sid, project, source = completed
    calls = []
    def planner(prompt, profile, **kwargs):
        context = json.loads(prompt)
        calls.append(context)
        if len(calls) == 1:
            return dict(kind='tool', tool='block_results', arguments=dict(operation='attach', report=str(source)), purpose='Adopt the supplied completed result'), {}
        if len(calls) == 2:
            identifier = context['evidence'][-1]['result']['attachment']['id']
            return dict(kind='tool', tool='block_results', arguments=dict(operation='analyze', attachment_id=identifier), purpose='Analyze saved scores without docking'), {}
        return dict(kind='answer', answer='Analysis is queued. No docking was launched.', evidence_ids=['E1', 'E2']), {}
    app.domain_planner = planner
    answer = app.ask(sid, 'Attach and analyze the completed docking report ' + str(source), 'deepseek')
    assert 'queued' in answer
    assert len(app.jobs(sid)) == 1
    assert app.task(sid)['status'] == 'queued'


def test_active_cli_scheduler_blocks_analysis(completed):
    from aidd_agent.prompt_workflow import file_lock
    app, sid, project, source = completed
    row = handle(app, sid, dict(operation='attach', report=str(source)))['attachment']
    entry = project / 'block-results' / (row['id'] + '.json')
    with file_lock(source.parent.parent / ('.' + source.parent.name + '.e097.lock')):
        with pytest.raises(OSError):
            analyze_attachment(project, dict(attachment_id=row['id'], attachment_sha256=sha(entry)), project / 'analysis-locked')


def test_real_background_executor_runs_analysis_without_engine(completed):
    app, sid, project, source = completed
    row = handle(app, sid, dict(operation='attach', report=str(source)))['attachment']
    queued = handle(app, sid, dict(operation='analyze', attachment_id=row['id']))
    job = app.task(sid, queued['tasks'][0]['id'])
    app.execute(job)
    final = app.task(sid, job['id'])
    assert final['status'] == 'complete', app.snapshot(sid)
    assert screening_summary(final)['kind'] == 'block_analyze'


def test_owned_top_n_inspection_and_changed_options(completed):
    app, sid, project, source = completed
    row = handle(app, sid, dict(operation='attach', report=str(source)))['attachment']
    args = dict(operation='analyze', attachment_id=row['id'], top_n=2)
    queued = handle(app, sid, args)
    task = queued['tasks'][0]['id']
    app.execute(app.task(sid, task))
    tools = DomainTools(app, sid, 'Show the best blocks and their top five molecules', 'deepseek')
    groups = tools.call('block_results', dict(operation='ranks', task_id=task))
    assert groups['status'] == 'select_group'
    group = groups['groups'][0]
    ranks = tools.call('block_results', dict(operation='ranks', task_id=task, **group))
    assert ranks['blocks'] and ranks['top_n'] == 2
    top_args = dict(operation='top', task_id=task, **group, block_id=ranks['blocks'][0]['block_id'])
    top = tools.call('block_results', top_args)
    assert len(top['candidates']) == len({r['molecule_id'] for r in top['candidates']})
    assert top['candidates'][0]['score'] <= top['candidates'][-1]['score']
    answer = json.loads(app.ask(sid, f'/block_ranks {task}', 'deepseek'))
    assert answer['groups'] == groups['groups']
    with pytest.raises(ValueError):
        handle(app, app.new_session(), top_args)
    assert handle(app, sid, args)['reused']
    changed = handle(app, sid, dict(args, top_n=5))
    assert changed['tasks'][0]['id'] != task


def test_repair_attach_arguments_then_return_queue_receipt_without_more_planning(completed):
    app, sid, project, source = completed
    calls=[]
    def planner(prompt, profile, **kwargs):
        context=json.loads(prompt);calls.append(context)
        assert kwargs['capabilities']['argument_contracts']['version']=='top-n-arguments-v1'
        if len(calls)==1:
            arguments=dict(operation='attach', report=str(source), top_n=10)
        elif len(calls)==2:
            assert 'unexpected=top_n' in context['previous_error']
            arguments=dict(operation='attach', report=str(source))
        elif len(calls)==3:
            identifier=context['evidence'][-1]['result']['attachment']['id']
            arguments=dict(operation='analyze', attachment_id=identifier, top_n=10, ranking_unit='molecule')
        else:
            raise AssertionError('A queued analysis must return its receipt without another planner call')
        return dict(kind='tool', tool='block_results', arguments=arguments, purpose='Adopt and rank saved results'), {}
    app.domain_planner=planner
    reply=app.ask(sid,'Attach and rank the existing report '+str(source),'deepseek')
    assert len(calls)==3 and len(app.jobs(sid))==1
    job=app.jobs(sid)[0]
    assert job['id'] in reply and 'queued' in reply
    listing=json.loads(app.ask(sid,'/dock_results','deepseek'))
    assert listing['analysis_tasks'][0]['id']==job['id']
    assert listing['argument_contract_version']=='top-n-arguments-v1'
    assert handle(app,app.new_session(),dict(operation='list'))['analysis_tasks']==[]
