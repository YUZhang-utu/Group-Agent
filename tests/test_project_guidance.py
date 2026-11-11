"""Mission checks and learning are local records, not evidence of model fine-tuning."""
import json
from pathlib import Path

import pytest

from aidd_agent.chat_agent import ChatAgent
from aidd_agent import project_guidance as g, project_execution as execution
from aidd_agent.molecule_preparation import handle as prepare
from aidd_agent.research_projects import handle as project


@pytest.fixture
def app(tmp_path):
    a=ChatAgent(tmp_path,start=False)
    yield a
    a.close()


def update_mission(app,sid,**changes):
    original=g.mission(app,sid)
    value={k:original[k] for k in ('objective','target','requirements','rules','excluded_work','constraints')}
    return g.revise(app,sid,{**value,**changes},original['revision'],'Test direct user action')


def test_mission_baseline_proposals_revisions_and_manual_edits(app):
    sid=app.new_session();original=g.mission(app,sid)
    changed={k:original[k] for k in ('objective','target','requirements','rules','excluded_work','constraints')}
    changed['objective']='Deliver candidate MOL2 files'
    p=g.handle(app,sid,dict(operation='propose',mission=changed))
    assert g.mission(app,sid)['objective']==original['objective']
    with pytest.raises(ValueError):g.handle(app,sid,dict(operation='apply',mission=changed))
    accepted=g.user_control(app,sid,dict(operation='adopt',proposal_id=p['proposal_id']))['mission']
    assert accepted['revision']==2 and accepted['baseline']['objective']==original['objective']
    with pytest.raises(ValueError,match='changed'):g.revise(app,sid,changed,1,'stale edit')
    path=Path(g.render(app,sid)['path']);path.write_text('User manual edit')
    assert g.render(app,sid)['status']=='manual_edit_conflict' and path.read_text()=='User manual edit'


def test_pinned_feedback_survives_history_and_isolates_targets(app):
    sid=app.new_session();quote='Use complete sequences and keep independent ligand runs separate.'
    app.message(sid,'user',quote)
    entry=g.handle(app,sid,dict(operation='feedback',kind='correction',quote=quote),quote)
    for i in range(120):app.message(sid,'user','Unrelated later progress message '+str(i))
    other=app.new_session();brief=g.brief(app,other)
    assert brief['pinned_feedback'][0]['text']==quote
    assert len(g.encoded(brief))<22000
    with pytest.raises(ValueError,match='actual user'):g.handle(app,sid,dict(operation='feedback',kind='constraint',quote='fabricated'),request='fabricated')
    foreign=project(app,sid,dict(operation='create',name='Other',objective='Other work',target='Other'))['session']
    assert not g.brief(app,foreign)['pinned_feedback']
    with pytest.raises(ValueError):g.user_control(app,foreign,dict(operation='resolve_feedback',entry_id=entry['entry_id']))
    g.user_control(app,other,dict(operation='resolve_feedback',entry_id=entry['entry_id']))
    assert not g.brief(app,sid)['pinned_feedback']


def test_real_preparation_reuses_content_and_keeps_evidence(app):
    sid=app.new_session();first=prepare(app,sid,dict(smiles='CCO'));job=app.task(sid,first['task_id']);app.execute(job)
    assert app.task(sid,job['id'])['status']=='complete'
    assert any(e['kind']=='evidence' and e.get('value',{}).get('work',{}).get('prepared_molecules')==1 for e in g.entries(app,sid))
    other=app.new_session();second=prepare(app,other,dict(smiles='CCO'));new=app.task(other,second['task_id']);app.execute(new)
    result=app.task(other,new['id'])
    assert result['status']=='complete' and result['report']==app.task(sid,job['id'])['report']
    data=g.handle(app,other,dict(operation='read'))
    card=next(c for c in data['task_cards'] if c['task_id']==new['id'])
    assert card['reuse_task_id']==job['id'] and not (Path(new['plan']).parent/'execution').exists()
    assert any(e.get('reused_from_task')==job['id'] for e in data['learning'])
    # Evidence is sourced by field; a changed report cannot remain fresh.
    e=g.handle(app,sid,dict(operation='evidence',task_id=job['id'],pointer='/steps/work/result/prepared_molecules'))
    assert e['value']==1
    report=Path(job['report']);report.write_text(report.read_text()+' ')
    assert any(v.get('freshness')=='stale_or_unavailable' for v in g.brief(app,sid,job['id'])['evidence'])
    third=prepare(app,sid,dict(smiles='CCO'));app.execute(app.task(sid,third['task_id']))
    assert app.task(sid,third['task_id'])['status']=='blocked'


def test_forbidden_plan_is_blocked_without_calculation_and_can_be_revised(app):
    sid=app.new_session();update_mission(app,sid,rules=dict(forbidden_actions=['research_task']))
    receipt=prepare(app,sid,dict(smiles='CCO'));job=app.task(sid,receipt['task_id']);app.execute(job)
    assert app.task(sid,job['id'])['status']=='blocked'
    assert 'forbids' in app.task(sid,job['id'])['error'] and not Path(job['report']).exists()
    update_mission(app,sid,rules={});app.execute(app.task(sid,job['id']))
    assert app.task(sid,job['id'])['status']=='complete'


def test_actual_af3_parameters_and_coverage_are_checked(app):
    sid=app.new_session();m=update_mission(app,sid,rules=dict(protein_sequences=['ACDEFGHIK'],protein_copies=2,seeds=[47,48],ligand_smiles=['CCO'],full_library=True))
    root=Path(app.context_for(sid)['project_root'])
    p=dict(name='test',protein_sequences=['ACDEFGHIK'],protein_copies=2,seeds=[47,48],ligand_smiles=['CCO'])
    plan=dict(steps=[dict(id='input',action='af3_prepare',params=p)])
    assert execution.validate(root,plan,m,{})['status']=='typed_rules_passed'
    for key,value in [('protein_copies',1),('seeds',[1]),('protein_sequences',['ACDEF']),('ligand_smiles',['CCC'])]:
        with pytest.raises(ValueError):execution.validate(root,dict(steps=[dict(id='i',action='af3_prepare',params={**p,key:value})]),m,{})
    with pytest.raises(ValueError,match='full-catalog'):execution.validate(root,dict(steps=[dict(id='s',action='consensus_budget',params={})]),m,{})


def test_requirements_constrain_actions_and_replicates_are_explicit(app):
    sid=app.new_session();update_mission(app,sid,requirements=[dict(id='R1',text='Prepare protein',actions=['af3_prepare'],acceptance='JSON input')])
    receipt=prepare(app,sid,dict(smiles='CCO'));job=app.task(sid,receipt['task_id']);app.execute(job)
    assert app.task(sid,job['id'])['status']=='blocked'
    update_mission(app,sid,requirements=[dict(id='R1',text='Prepare molecule',actions=['research_task'],acceptance='Prepared 3D record')])
    app.execute(app.task(sid,job['id']));assert app.task(sid,job['id'])['status']=='complete'
    receipt=prepare(app,sid,dict(smiles='CCO'));new=app.task(sid,receipt['task_id'])
    g.user_control(app,sid,dict(operation='replicate',task_id=new['id'],reason='Independent preparation requested'))
    app.execute(new);assert app.task(sid,new['id'])['report']!=app.task(sid,job['id'])['report']
    assert Path(new['report']).exists()


def test_context_never_silently_drops_essential_constraints(app):
    sid=app.new_session();update_mission(app,sid,constraints=['Always retain this constraint'])
    with pytest.raises(ValueError,match='not silently dropped'):g.brief(app,sid,budget=20)


def test_budget_is_enforced_after_optional_evidence_and_task_metadata(app):
    sid=app.new_session()
    for i in range(40):g.add_entry(app,sid,'observation',dict(text='Observation '+str(i)+' '+'x'*200))
    for budget in (1800,2500,4000,7000):
        value=g.brief(app,sid,budget=budget)
        assert value['context_characters']==len(g.encoded(value))<=budget
        assert value['mission']['objective']==g.mission(app,sid)['objective']


def test_http_user_controls_restart_and_revision_bound_acceptance(tmp_path):
    import threading
    from urllib.request import Request,urlopen
    from urllib.error import HTTPError
    from aidd_agent.chat_web import make_server
    app=ChatAgent(tmp_path,start=False);sid=app.new_session()
    server=make_server(app,0,token='guidance-fixture-token');thread=threading.Thread(target=server.serve_forever,daemon=True);thread.start()
    def call(session,arguments):
        with urlopen(Request('http://127.0.0.1:'+str(server.server_port)+'/api/project-guidance',
            data=json.dumps(dict(session=session,arguments=arguments)).encode(),headers={
                'Authorization':'Bearer guidance-fixture-token','Content-Type':'application/json'})) as response:return json.load(response)
    try:
        initial=call(sid,dict(operation='read'))['mission']
        value={k:initial[k] for k in ('objective','target','requirements','rules','excluded_work','constraints')}
        value['constraints']=['Keep artifacts local.']
        call(sid,dict(operation='apply',base_revision=1,mission=value))
        with pytest.raises(HTTPError):call(sid,dict(operation='apply',base_revision=1,mission=value))
        receipt=prepare(app,sid,dict(smiles='CCO'));app.execute(app.task(sid,receipt['task_id']))
        evidence=next(e for e in call(sid,dict(operation='read'))['learning'] if e['kind']=='evidence')
        args=dict(operation='accept_requirement',base_revision=2,requirement_id='R1',evidence_ids=[evidence['id']],reason='Reviewed the prepared molecule report.')
        other=project(app,sid,dict(operation='create',name='Other',objective='Other',target='Other'))['session']
        with pytest.raises(HTTPError):call(other,args)
        assert call(sid,args)['status']=='user_accepted'
        assert g.brief(app,sid)['accepted_requirements'][0]['requirement_id']=='R1'
        call(sid,dict(operation='apply',base_revision=2,mission=value))
        assert not g.brief(app,sid)['accepted_requirements']
        with pytest.raises(HTTPError):call(sid,args)
    finally:server.shutdown();server.server_close();thread.join();app.close()
    reopened=ChatAgent(tmp_path,start=False)
    try:
        assert g.mission(reopened,sid)['revision']==3
        assert g.mission(reopened,sid)['baseline']==initial['baseline']
        assert any(e['kind']=='acceptance' for e in g.entries(reopened,sid))
        again=prepare(reopened,sid,dict(smiles='CCO'));reopened.execute(reopened.task(sid,again['task_id']))
        assert reopened.task(sid,again['task_id'])['report']==reopened.task(sid,receipt['task_id'])['report']
    finally:reopened.close()


def test_scripted_chat_records_progress_then_uses_it_in_a_later_conversation(tmp_path):
    from test_domain_agent import scripted,tool,answer
    seen=[];quote='Exclude unstable compounds before the next export.'
    app=ChatAgent(tmp_path,start=False,domain_planner=scripted([
        tool('project_guidance',operation='feedback',kind='correction',quote=quote),answer('Recorded your correction.',['E1']),
        answer('The next export must account for the recorded exclusion.',[])],seen))
    try:
        sid=app.new_session();app.ask(sid,quote,'deepseek')
        for i in range(100):app.message(sid,'user','Later progress '+str(i))
        other=app.new_session();app.ask(other,'What should the next task account for?','deepseek')
        assert seen[-1]['project_memory']['pinned_feedback'][0]['text']==quote
        assert g.mission(app,sid)['revision']==1
    finally:app.close()


def test_changed_export_count_reuses_ranking_lane_and_failed_duplicate_blocks(tmp_path):
    import io
    from test_ligand_library import tiny_library
    from aidd_agent import ligand_library,candidate_sets
    from aidd_agent.chat_uploads import upload
    batch,raw=tiny_library(tmp_path);runtime=tmp_path/'runtime.json'
    runtime.write_text(json.dumps(dict(search=dict(batch=str(batch),workers=1,refine_chunk=1))))
    app=ChatAgent(tmp_path/'chat',start=False,runtime=runtime,allow_compute=True)
    try:
        sid=app.new_session();file=upload(app,sid,'query.mol2',io.BytesIO(raw.encode()),len(raw.encode()))
        args=dict(operation='start',attachment_id=file['id'],export_molecules=1)
        first=ligand_library.handle(app,sid,args);job=app.task(sid,first['task_id']);app.execute(job)
        assert app.task(sid,job['id'])['status']=='complete'
        second=ligand_library.handle(app,sid,{**args,'export_molecules':2});app.execute(app.task(sid,second['task_id']))
        blocked=app.task(sid,second['task_id'])
        assert blocked['status']=='blocked' and 'Only the export count changed' in blocked['error']
        assert not Path(blocked['report']).exists()
        imported=candidate_sets.handle(app,sid,dict(operation='import',task_id=job['id'],name='Full saved ranking'))
        app.execute(app.task(sid,imported['task_id']));source=app.task(sid,imported['task_id'])
        assert source['status']=='complete',source['error']
        set_id=json.loads(Path(source['report']).read_text())['steps']['work']['result']['set_id']
        exported=candidate_sets.handle(app,sid,dict(operation='export',set_id=set_id,count=2))
        app.execute(app.task(sid,exported['task_id']));done=app.task(sid,exported['task_id'])
        assert done['status']=='complete',done['error']
        # A failed equivalent preparation is diagnosed rather than retried as a new task.
        one=prepare(app,sid,dict(smiles='CCO'));a=app.task(sid,one['task_id']);execution.check(app,a)
        app.update(a['id'],status='failed',error='Recorded interruption fixture')
        two=prepare(app,sid,dict(smiles='CCO'));app.execute(app.task(sid,two['task_id']))
        assert a['id'] in app.task(sid,two['task_id'])['error']
    finally:app.close()


def test_changed_mission_blocks_started_task_and_preserves_historical_link(app):
    sid=app.new_session();receipt=prepare(app,sid,dict(smiles='CCO'));job=app.task(sid,receipt['task_id'])
    execution.check(app,job);execution.started(app,job)
    with pytest.raises(ValueError,match='historical'):g.handle(app,sid,dict(operation='link',task_id=job['id'],requirement_ids=['R1'],contribution='Changed after launch'))
    update_mission(app,sid,constraints=['New scope after this task started'])
    app.execute(job);assert 'Mission changed since' in app.task(sid,job['id'])['error']
    assert not Path(job['report']).exists()


def test_missing_artifact_prevents_reuse_even_when_report_is_unchanged(app):
    sid=app.new_session();receipt=prepare(app,sid,dict(smiles='CCO'));job=app.task(sid,receipt['task_id']);app.execute(job)
    card=next(c for c in g.handle(app,sid,dict(operation='read'))['task_cards'] if c['task_id']==job['id'])
    artifact=next(Path(p) for p in card['output_hashes'] if p.endswith('.sdf'));artifact.unlink()
    second=prepare(app,sid,dict(smiles='CCO'));app.execute(app.task(sid,second['task_id']))
    assert app.task(sid,second['task_id'])['status']=='blocked'


def test_empty_allowlist_denies_execution_and_live_af3_is_not_implicitly_cached(app):
    sid=app.new_session();m=update_mission(app,sid,rules=dict(allowed_actions=[]));root=Path(app.context_for(sid)['project_root'])
    plan=dict(steps=[dict(id='input',action='af3_prepare',params=dict(name='test',protein_sequences=['ACDEFGHIK']))])
    with pytest.raises(ValueError,match='allowlist'):execution.validate(root,plan,m,{})
    assert execution.reusable(root,plan)
    plan['steps'].append(dict(id='predict',action='af3_run',params=dict(input_step='input')))
    assert not execution.reusable(root,plan)


def test_existing_project_notes_are_in_both_the_dossier_and_pinned_context(app):
    sid=app.new_session()
    project(app,sid,dict(operation='note',kind='constraint',text='Retain every source compound identifier.'))
    dossier=Path(g.render(app,sid)['path']).read_text(encoding='utf-8')
    assert 'Retain every source compound identifier.' in dossier
    assert 'not independently verified' in dossier
    assert g.brief(app,sid)['notes'][0]['text']=='Retain every source compound identifier.'
