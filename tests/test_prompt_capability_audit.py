"""Real Chat workers with explicitly synthetic planner and AF3 runner."""
import json
from pathlib import Path
import sys

import pytest

from aidd_agent.chat_agent import ChatAgent
from aidd_agent.domain_tools import DomainTools
from aidd_agent.prompt_workflow import create_plan
from test_domain_agent import answer, scripted, tool
from test_general_prompt_workflows import molecular_files, params


@pytest.mark.parametrize('heterodimer', [False, True])
@pytest.mark.parametrize('with_ligand', [False, True])
def test_prompt_dimer_reaches_real_worker_and_external_runner(tmp_path, heterodimer, with_ligand):
    a, b = 'ACDEFGHIKLMNPQRSTVWY', 'MNPQRSTVWYACDEFGHIKL'
    prompt = (f'Run AF3 for a heterodimer with one copy each. Sequence A: {a}. Sequence B: {b}.'
              if heterodimer else f'Run AF3 for a homodimer: exactly two copies of {a}.')
    prompt += '\nModel seeds: [47, 48, 49, 50, 51].'
    if with_ligand: prompt += '\nInclude exactly one ligand with SMILES CCO.'
    else: prompt += '\nDo not include any ligand.'
    fields = dict(name='synthetic_dimer', protein_sequences=[a, b] if heterodimer else [a],
                  seeds=[47, 48, 49, 50, 51])
    if not heterodimer: fields['protein_copies'] = 2
    if with_ligand: fields['ligand_smiles'] = ['CCO']
    plan = dict(version=1, summary='Synthetic dimer dispatch audit', clarifications=[], steps=[
        dict(id='input', action='af3_prepare', params=fields),
        dict(id='predict', action='af3_run', params=dict(input_step='input'))])
    def planner(db, user, project, request, **kwargs):
        return create_plan(db, user, project, request, local_plan=plan)
    runner = tmp_path/'synthetic_runner.py'
    runner.write_text("""import json, sys
from pathlib import Path
args = dict(arg[2:].split('=', 1) for arg in sys.argv[1:] if arg.startswith('--') and '=' in arg)
payload = json.loads(Path(args['json_path']).read_text())
output = Path(args['output_dir']); output.mkdir(parents=True, exist_ok=True)
(output/'observed-input.json').write_text(json.dumps(payload))
(output/'synthetic_dimer_model.cif').write_text('data_synthetic_fixture_not_prediction\\n#\\n')
(output/'synthetic_dimer_summary_confidences.json').write_text(json.dumps({'fixture_not_scientific_prediction': True}))
""")
    models = tmp_path/'models'; models.mkdir()
    databases = tmp_path/'databases'; databases.mkdir()
    profile = tmp_path/'af3.json'
    profile.write_text(json.dumps(dict(backend='alphafold3', python=sys.executable, runner=str(runner),
        model_parameters=str(models), databases=str(databases), license_acknowledged=True)))
    runtime = tmp_path/'runtime.json'; runtime.write_text(json.dumps(dict(af3_profile=str(profile))))
    domain = scripted([tool('workflow', decision=dict(intent='run', request=prompt, task_id=None, message=''))])
    app = ChatAgent(tmp_path/'chat', runtime=runtime, allow_compute=True, start=False,
                    planner=planner, domain_planner=domain)
    try:
        sid = app.new_session(); app.ask(sid, prompt, 'deepseek')
        assert len(app.jobs(sid)) == 1
        app.execute(app.task(sid))
        job = app.task(sid); assert job['status'] == 'complete', job.get('error')
        report = json.loads(Path(job['report']).read_text())
        payload = json.loads(Path(report['steps']['input']['result']['input_json']).read_text())
        assert payload['modelSeeds'] == fields['seeds']
        entities = payload['sequences']
        assert [r['protein']['sequence'] for r in entities[:2]] == [a, b if heterodimer else a]
        assert [next(iter(r.values()))['id'] for r in entities] == (['A', 'B', 'C'] if with_ligand else ['A', 'B'])
        if with_ligand: assert entities[2]['ligand']['smiles'] == 'CCO'
        output = Path(report['steps']['predict']['result']['output_directory'])
        assert json.loads((output/'observed-input.json').read_text()) == payload
        assert len(app.jobs(sid)) == 1
    finally:
        app.close()


def test_full_library_prompt_dispatches_generic_tool(tmp_path):
    receptor, ligand = molecular_files(tmp_path/'inputs')
    request = dict(operation='start', **params(receptor, ligand, 'PROMPT_AUDIT'))
    prompt = f'Search the entire configured local library using receptor {receptor}, chain B, and posed reference ligand {ligand}. Name it PROMPT_AUDIT.'
    app = ChatAgent(tmp_path/'chat', start=False, result_import_roots=[receptor.parent],
                    domain_planner=scripted([dict(kind='tool', tool='target_search', arguments=request,
                                                 purpose='Dispatch the full-library request')]))
    try:
        sid = app.new_session(); reply = app.ask(sid, prompt, 'deepseek')
        assert len(app.jobs(sid)) == 1
        assert app.task(sid)['id'] in reply and 'not confirmed' in reply
        plan = json.loads(Path(app.task(sid)['plan']).read_text())['plan']
        assert [s['action'] for s in plan['steps']] == ['target_query_prepare', 'target_library_search']
        assert app.task(sid)['status'] == 'queued'
    finally:
        app.close()


def test_target_search_cannot_dispatch_after_pending_tool_receipt(tmp_path, monkeypatch):
    calls = []
    def pending_call(self, name, args):
        calls.append(name)
        return dict(status='running', task_id='already-active')
    monkeypatch.setattr(DomainTools, 'call', pending_call)
    request = dict(operation='prepare', name='DO_NOT_SUBMIT', pdb_id='1ABC',
                   receptor_chains=['A'], ligand_instance=dict(ccd='ABC', chain='B', residue='1'))
    app = ChatAgent(tmp_path/'chat', start=False, domain_planner=scripted([
        tool('tasks'), dict(kind='tool', tool='target_search', arguments=request, purpose='Attempt another dispatch'),
        answer('The existing task is still active.', ['E1'])]))
    try:
        sid = app.new_session(); app.ask(sid, 'Inspect the active task before searching 1ABC.', 'deepseek')
        assert calls == ['tasks'] and not app.jobs(sid)
        trace = json.loads(next((app.root/'agents'/sid).glob('*.json')).read_text())
        assert 'pending' in trace['steps'][1]['error']
    finally:
        app.close()
