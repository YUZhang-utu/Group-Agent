import json
from pathlib import Path

import pytest

from aidd_agent.chat_agent import ChatAgent
from aidd_agent.domain_agent import run
from aidd_agent.domain_tools import DomainTools
from aidd_agent.pdb_review import coordinate, note_requested, request_task, render
from aidd_agent.prompt_workflow import create_plan, run_plan
from aidd_agent.research_projects import handle as project
from test_pdb_inspection import Services


TASK = 'fe88df6f8f664528'
REQUEST = f'''Read the saved inspection report for completed task {TASK}.
Do not repeat the inspection or launch any calculation.
Summarize:
1. Protein identity and deposited evidence for G12V.
2. Deposited sequence length and coordinate coverage for each protein chain.
3. Asymmetric-unit chains versus annotated biological assemblies.
4. Ligands and ions, including chain and residue identifiers.
5. Database mapping and unresolved isoform assignment.
6. Saved mmCIF and inspection-report paths.
Distinguish deposited annotations from independently verified conclusions.
Record a concise project note linked to this task and report whether it was saved.'''


@pytest.fixture
def saved(tmp_path):
    app = ChatAgent(tmp_path / 'chat', start=False)
    sid = app.new_session()
    ctx = app.context_for(sid)
    plan = dict(version=1, summary='Inspect', clarifications=[], steps=[
        dict(id='inspect', action='pdb_inspect', params=dict(pdb_id='1ABC'))])
    path = create_plan(Path(ctx['db']), ctx['user_id'], ctx['project_id'], 'Inspect 1ABC', local_plan=plan)
    report, _ = run_plan(Path(ctx['db']), ctx['user_id'], ctx['project_id'], path, services=Services())
    with app.connect() as db:
        db.execute('INSERT INTO jobs(id,session,status,plan,report,request,created) VALUES(?,?,?,?,?,?,?)',
            (TASK, sid, 'complete', str(path), str(path.parent / 'execution/report.json'), 'Inspect 1ABC', 1))
    yield app, sid, report['steps']['inspect']['result']
    app.close()


def test_native_route_without_provider_or_planner_and_idempotent_note(saved, monkeypatch):
    app, sid, result = saved
    from aidd_agent import domain_agent, llm_profiles
    monkeypatch.setattr(domain_agent, 'chat_plan', lambda *a, **k: pytest.fail('No model call'))
    monkeypatch.setattr(llm_profiles, 'select_llm_profile', lambda *a, **k: pytest.fail('No provider required'))
    first = run(app, sid, REQUEST, 'deepseek')
    assert 'Project note: saved' in first and 'G12V' not in first
    assert 'GLY' in first and 'VAL' in first and 'monomeric' in first
    assert 'not_independently_verified' in first
    assert result['report'] in first and result['path'] in first
    other = app.new_session()
    second = run(app, other, REQUEST, 'deepseek')
    assert 'Project note: already_saved' in second
    with app.connect() as db:
        assert db.execute('SELECT count(*) FROM research_notes').fetchone()[0] == 1
        assert db.execute('SELECT count(*) FROM jobs').fetchone()[0] == 1
    trace = json.loads(next((app.root / 'agents' / other).glob('*.json')).read_text())
    assert trace['status'] == 'complete'
    assert trace['answer_mode'] == 'local_pdb_report_review_no_model_call'
    assert trace['steps'][0]['result']['structure_sha256'] == result['sha256']


def test_note_is_opt_in_and_tool_cannot_invent_authorization(saved):
    app, sid, _ = saved
    request = f'Read the saved inspection report for task {TASK}. Do not save a project note.'
    assert not note_requested(request)
    output = coordinate(app, sid, request)
    assert 'Project note: not_requested' in output
    with pytest.raises(ValueError, match='explicit current user'):
        DomainTools(app, sid, request, 'deepseek').call('pdb_review', dict(task_id=TASK, save_note=True))
    with app.connect() as db:
        assert db.execute('SELECT count(*) FROM research_notes').fetchone()[0] == 0


def test_chat_entry_stores_answer_and_emits_receipts_without_provider(saved, monkeypatch):
    app, sid, _ = saved
    from aidd_agent import llm_profiles
    monkeypatch.setattr(llm_profiles, 'select_llm_profile', lambda *a, **k: pytest.fail('No provider needed'))
    events = []
    answer = app.ask(sid, REQUEST, 'deepseek', on_event=events.append)
    assert 'Project note: saved' in answer
    assert [e['type'] for e in events] == ['tool_start', 'tool_end']
    with app.connect() as db:
        stored = db.execute("SELECT text FROM messages WHERE session=? AND role='assistant' ORDER BY id DESC", (sid,)).fetchone()[0]
    assert stored == answer


def test_missing_structure_and_report_escape_stop_reading(saved, tmp_path):
    app, sid, result = saved
    structure = Path(result['path'])
    structure.rename(structure.with_suffix('.retained'))
    answer = coordinate(app, sid, REQUEST)
    assert 'could not be completed' in answer and 'Project note: saved' not in answer
    job = app.task(sid, TASK)
    task_report = Path(job['report'])
    data = json.loads(task_report.read_text())
    outside = tmp_path / 'foreign.json'
    outside.write_text('{"secret":"DO_NOT_DISCLOSE"}')
    data['steps']['inspect']['result']['report'] = str(outside)
    task_report.write_text(json.dumps(data))
    answer = coordinate(app, sid, REQUEST)
    assert 'could not be completed' in answer and 'DO_NOT_DISCLOSE' not in answer


@pytest.mark.parametrize('field', ['report', 'path'])
def test_changed_report_or_structure_never_summarized_or_noted(saved, field):
    app, sid, result = saved
    path = Path(result[field])
    if field == 'report':
        data = json.loads(path.read_text()); data['pdb_id'] = '9XYZ'
        path.write_text(json.dumps(data))
    else:
        path.write_text('changed structure')
    output = coordinate(app, sid, REQUEST)
    assert 'hash mismatch' in output and 'No scientific task was submitted' in output
    with app.connect() as db:
        assert db.execute('SELECT count(*) FROM research_notes').fetchone()[0] == 0


def test_foreign_project_and_incomplete_task(saved):
    app, sid, _ = saved
    other = project(app, sid, dict(operation='create', name='Other', target='Other', objective='Other'))['session']
    output = coordinate(app, other, REQUEST)
    assert 'No matching task in this project' in output and 'GLY' not in output
    with app.connect() as db:
        db.execute("UPDATE jobs SET status='running' WHERE id=?", (TASK,))
    output = coordinate(app, sid, REQUEST)
    assert 'running' in output and 'No completed inspection' in output


@pytest.mark.parametrize('prompt', [
    f'Read report for task {TASK} and run AF3.',
    f'Review report for task {TASK} and search literature.',
    f'Compare task {TASK} to a new structure.',
    f'Read reports for tasks {TASK} and aaaaaaaaaaaaaaaa.',
    f'Resume task {TASK}.',
])
def test_mixed_or_ambiguous_requests_keep_planner(prompt):
    assert request_task(prompt) is None


def test_exact_request_and_chinese_read_request_match():
    assert request_task(REQUEST) == TASK
    assert request_task('\u8bfb\u53d6\u68c0\u67e5\u62a5\u544a ' + TASK) == TASK


def test_display_coverage_and_missing_annotations():
    report = dict(pdb_id='1ABC', structure='/saved/1ABC.cif',
        polymers=[dict(entity_id=str(i), description='x' * 300) for i in range(25)])
    answer, coverage = render(report, TASK, '/saved/report.json')
    assert coverage['polymers'] == dict(total=25, displayed=20, shortened_fields=20)
    assert '(20/25 rows)' in answer and '[field shortened]' in answer
    assert 'not evidence of absence' in answer


def test_other_artifact_types_are_not_claimed_as_pdb_reviews(saved):
    app, sid, _ = saved
    job = app.task(sid, TASK)
    path = Path(job['report']); data = json.loads(path.read_text())
    data['steps']['inspect']['action'] = 'af3_run'; path.write_text(json.dumps(data))
    assert coordinate(app, sid, REQUEST) is None
