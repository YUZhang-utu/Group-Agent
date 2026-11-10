"""Confidence label diagnostics preserve raw scientific outputs."""
import copy
import hashlib
import json
from pathlib import Path

import pytest

from aidd_agent.af3_analysis import assess_chain_mapping, compact_confidence, explain_confidence
from aidd_agent.chat_agent import ChatAgent
from aidd_agent.domain_tools import DomainTools
from aidd_agent.prediction import inspect_alphafold3_output
from aidd_agent.prompt_workflow import project_root
from aidd_agent.structure_review import confidence as review_confidence


def repeated():
    return dict(chain_ids=['A'] * 905 + ['B'] * 26, chain_iptm=[.84, .84],
                chain_ptm=[.62, .3], chain_pair_iptm=[[.62, .84], [.84, .3]],
                chain_pair_pae_min=[[.76, 4.18], [6.6, .77]], iptm=.84, ptm=.63)


def test_repeated_labels_are_compact_without_inventing_mapping():
    raw = repeated(); before = copy.deepcopy(raw)
    analysis = explain_confidence(raw, ['A', 'B'])
    mapping = analysis['chain_mapping']
    assert mapping['status'] == 'unverified' and mapping['matrix_chain_ids'] is None
    assert mapping['label_counts'] == {'A': 905, 'B': 26}
    assert mapping['metric_dimensions']['chain_pair_iptm'] == 2
    assert len(mapping['warnings']) == 2
    assert analysis['metrics']['chain_ids']['value']['total_entries'] == 931
    assert analysis['metrics']['chain_pair_iptm']['value'] == raw['chain_pair_iptm']
    assert raw == before
    assert len(json.dumps(compact_confidence(raw))) < 600


def test_unique_output_order_is_preserved_over_input_order():
    raw = repeated(); raw['chain_ids'] = ['B', 'A']
    analysis = explain_confidence(raw, ['A', 'B'])
    assert analysis['chain_mapping']['matrix_chain_ids'] == ['B', 'A']
    assert analysis['chain_mapping']['warnings'] == []


@pytest.mark.parametrize('change', [
    {'chain_ids': None}, {'chain_ids': []}, {'chain_ids': ['A', None]},
    {'chain_ptm': [.1]}, {'chain_pair_iptm': [[.1, .2]]},
    {'chain_ptm': [[.1], [.2]]}, {'chain_pair_pae_min': []},
])
def test_invalid_mapping_never_assigns_chain_labels(change):
    raw = repeated(); raw['chain_ids'] = ['A', 'B']; raw.update(change)
    assessment = assess_chain_mapping(raw)
    assert assessment['warnings'] and assessment['matrix_chain_ids'] is None


def test_full_array_labels_do_not_override_summary_mapping():
    raw = dict(chain_ids=['A', 'B'], token_chain_ids=['A'] * 10 + ['B'] * 2,
               atom_chain_ids=['A'] * 90 + ['B'] * 2, chain_ptm=[.5, .3])
    assert assess_chain_mapping(raw)['matrix_chain_ids'] == ['A', 'B']
    display = compact_confidence(raw)
    assert display['token_chain_ids']['counts'] == {'A': 10, 'B': 2}
    assert display['atom_chain_ids']['total_entries'] == 92


def test_parser_preserves_raw_summary_and_hash(tmp_path):
    (tmp_path/'job_model.cif').write_text('data_job\n#\n')
    path = tmp_path/'job_summary_confidences.json'
    path.write_text(json.dumps(repeated()))
    before = path.read_bytes()
    manifest = inspect_alphafold3_output(tmp_path, construct_name='job', model_version='fixture', chain_ids=['A', 'B'])
    assert manifest['confidence'] == repeated()
    assert manifest['provenance']['confidence_sha256'] == hashlib.sha256(before).hexdigest()
    assert path.read_bytes() == before


def test_existing_report_checked_in_chat_and_tools_without_reexecution(tmp_path):
    app = ChatAgent(tmp_path/'chat', start=False)
    try:
        sid = app.new_session(); ctx = app.context
        root = project_root(Path(ctx['db']), ctx['user_id'], ctx['project_id'])
        path = root/'existing.json'
        path.write_text(json.dumps(dict(status='complete', steps=dict(predict=dict(
            action='af3_run', status='complete', result=dict(confidence=repeated()))))))
        before = path.read_bytes()
        with app.connect() as db:
            db.execute('INSERT INTO jobs(id,session,status,report,plan,request,created) VALUES(?,?,?,?,?,?,?)',
                       ('existing-af3', sid, 'complete', str(path), str(root/'plan.json'), 'Existing AF3 output', 1))
        reply = app.ask(sid, '/results existing-af3', 'deepseek')
        assert 'repeated labels' in reply and "'A': 905" in reply
        assert "'A', 'A', 'A'" not in reply and len(reply) < 2000
        tools = DomainTools(app, sid, 'Read existing results', 'deepseek')
        for extra in ({}, {'step_id': 'predict'}):
            result = tools.call('task_report', dict(task_id='existing-af3', **extra))
            assert result['confidence_chain_mapping']['predict']['matrix_chain_ids'] is None
        analysis = review_confidence(app.task(sid, 'existing-af3'), root)
        assert analysis['chain_mapping']['warnings']
        assert path.read_bytes() == before and len(app.jobs(sid)) == 1
        assert app.task(sid, 'existing-af3')['status'] == 'complete'
    finally:
        app.close()
