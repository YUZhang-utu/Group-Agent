"""Fixture evidence only: no claim of real MDM2 retrieval or docking quality."""
import json
import os
from pathlib import Path
import sqlite3

import numpy as np
import pytest

from aidd_agent import block_campaign as campaign
from aidd_agent import block_conformer_search as search
from aidd_agent import equiscore as es
from aidd_agent.block_sampling import sample
from aidd_agent.chat_agent import ChatAgent
from aidd_agent.domain_tools import DomainTools, validate_arguments, argument_contracts
from aidd_agent.expanded_wee1 import fingerprint
from aidd_agent.final_work_blocks import read, sha
from aidd_agent.joint_spatial_profiles import save
from aidd_agent.prompt_workflow import project_root, run_plan
from test_equiscore import engine, source_panel
from test_block_evaluation import fixture, raw, plants_profile, fake_plants


@pytest.fixture
def imported(engine, tmp_path, monkeypatch):
    monkeypatch.setenv('PYTHONPATH',str(Path(__file__).resolve().parents[1]/'src'))
    _, profile, calls = engine
    app = ChatAgent(tmp_path/'chat', start=False)
    sid = app.new_session()
    c = app.context
    project = project_root(Path(c['db']),c['user_id'],c['project_id'])
    root = project/'runs/cli';root.mkdir(parents=True)
    source = source_panel(root)
    es.run(source,profile,root/'pilot')
    es.run(source,profile,root/'full',scope='full',pilot_report=root/'pilot/report.json')
    es.analyze(root/'full/report.json',root/'analysis')
    path = root/'analysis/report.json'
    yield app,sid,project,path,calls
    app.close()


def test_import_runs_in_real_queue_without_inference_and_inspects(imported):
    app,sid,project,path,calls = imported
    result = campaign.handle(app,sid,dict(operation='import',report=str(path)))
    assert result['status']=='queued'
    assert campaign.handle(app,sid,dict(operation='import',report=str(path)))['task_id']==result['task_id']
    app.execute(app.task(sid,result['task_id']))
    job = app.task(sid,result['task_id'])
    assert job['status']=='complete', app.snapshot(sid)
    assert len(calls)==2
    ranks = campaign.handle(app,sid,dict(operation='ranks',task_id=job['id'],scheme='E095',receptor='R1'))
    assert ranks['blocks'][0]['block_id']=='B'
    old = campaign.handle(app,sid,dict(operation='ranks',task_id=job['id'],scheme='E095',receptor='R1',score='ChemPLP'))
    assert old['blocks'][0]['block_id']=='A'
    preview = campaign.handle(app,sid,dict(operation='preview',task_id=job['id'],scheme='E095',receptor='R1'))
    assert preview['candidate_unit']=='conformer' and preview['block_population_conformers']==2000
    assert len(campaign.handle(app,sid,{})['tasks'])==1
    assert campaign.select_blocks(path,'E095','R1','intersection',1)['block_ids']==[]
    assert campaign.select_blocks(path,'E095','R1','union',1)['block_ids']==['A','B']
    other=app.new_session()
    with pytest.raises(ValueError):
        campaign.handle(app,other,dict(operation='ranks',task_id=job['id']))


def test_import_rejects_changed_seals_and_ungrounded_paths(imported):
    app,sid,project,path,calls=imported
    tools=DomainTools(app,sid,'Import saved scores','deepseek')
    with pytest.raises(ValueError,match='exact analysis'):
        tools.call('block_campaign',dict(operation='import',report=str(path)))
    request=dict(report=str(path),report_sha256=sha(path))
    (path.parent/'block_ranking.sqlite').write_bytes(b'changed')
    output=project/'failed-import';output.mkdir()
    with pytest.raises(ValueError,match='changed'):
        campaign.adopt(request,output,project)
    assert len(calls)==2


@pytest.mark.parametrize('arguments',[
    {'operation':'start','task_id':'a'*16,'query_task_id':'b'*16,'scheme':'E095','receptor':'R1','conformers':True},
    {'operation':'start','task_id':'a'*16,'query_task_id':'b'*16,'scheme':'E095','receptor':'R1','molecules':100000},
    {'operation':'preview','task_id':'a'*16,'scheme':'E095','receptor':'R1','blocks':6},
    {'operation':'import','report':'x','dock':True},
])
def test_strict_tool_arguments(arguments):
    with pytest.raises(ValueError):validate_arguments('block_campaign',arguments)
    assert 'start' in argument_contracts()['block_campaign']


def member_fixture(root):
    profile, source = fixture(root)
    sample(profile,root/'sample',count=100)
    selection=dict(scheme='E095',block_ids=['b0'])
    database=root/'members.sqlite'
    locations=search.members(profile,selection,root/'sample/report.json',database)
    batch=root/'batch';shard=batch/'artifacts/shard';shard.mkdir(parents=True)
    (batch/'chemical').mkdir()
    np.array([f'c{i}' for i in range(12)],dtype='S16').tofile(shard/'conformer_ids.bin')
    save(shard/'manifest.json',dict(conformers=12,global_id_start=0,files={
        'conformer_ids.bin':dict(bytes=192,sha256=sha(shard/'conformer_ids.bin'))}))
    catalog=dict(library_id='fixture',shards=[dict(path=str(shard),name='shard',manifest_sha256=sha(shard/'manifest.json'))])
    save(batch/'artifacts/catalog.json',catalog)
    save(batch/'chemical/catalog.json',dict(library_id='fixture',artifact_v1_catalog_sha256=sha(batch/'artifacts/catalog.json'),shards=[]))
    with sqlite3.connect(batch/'registry.sqlite3') as db:
        db.execute('CREATE TABLE conformer(id TEXT,source_record_index INTEGER,source_path TEXT,content_sha256 TEXT)')
        import hashlib
        db.executemany('INSERT INTO conformer VALUES(?,?,?,?)',
            ((f'c{i}',i,str(source),hashlib.sha256(raw(i).encode()).hexdigest()) for i in range(12)))
    return profile,source,database,locations,batch


def pose(gid,contact=0.8):
    return dict(global_id=gid,conformer_id=f'c{gid}',molecule_id=f'artifact-m{gid//2}',
                composite_score=contact,optional_score=contact,gaussian_same_pose=0.5,transform=np.eye(4).tolist())


def chunks(output,rows):
    root=output/'chunks/00';root.mkdir(parents=True,exist_ok=True)
    path=root/'0000000000.poses.jsonl'
    path.write_text(''.join(json.dumps(row)+'\n' for row in rows))
    save(root/'0000000000.receipt.json',dict(files=fingerprint([path])))


def test_membership_mapping_preserves_conformers_and_docking_handoff(tmp_path,monkeypatch):
    profile,source,database,locations,batch=member_fixture(tmp_path)
    ids=search.map_artifacts(batch,database,locations['descriptors'])
    assert ids.tolist()==list(range(6))  # No conformer from b1 is admitted.
    chunks(tmp_path,[pose(i) for i in ids.tolist()]+[pose(0,0.9)])
    rows,total=search.rank_conformers(tmp_path,dict(templates=[dict(query_id='fixture')]),database,4)
    assert total==6 and [r['gid'] for r in rows]==[0,1,2,3]
    result=search.export_selected(tmp_path/'selected',database,rows,'E095')
    assert result['unique_conformers']==4 and result['unique_molecules']==2
    # Ranking again must allow deterministic resume without rewriting source chemistry.
    rows2,_=search.rank_conformers(tmp_path,dict(templates=[dict(query_id='fixture')]),database,4)
    assert search.export_selected(tmp_path/'selected',database,rows2,'E095')==result
    from aidd_agent.block_plants import prepare,run
    plants=plants_profile(tmp_path,source)
    prepare(tmp_path/'selected/report.json',plants,tmp_path/'prepare')
    monkeypatch.setattr('aidd_agent.block_plants.subprocess.run',fake_plants)
    report=run(tmp_path/'prepare/report.json',tmp_path/'dock')
    assert report['scored']==4 and report['failed_jobs']==0
    top=search.docking_candidates(tmp_path,tmp_path/'selected',tmp_path/'dock')
    assert len(top)==4 and {r['block_id'] for r in top}=={'b0'}
    assert top[0]['score']<top[-1]['score']


def test_missing_mapping_does_not_expand_outside_block(tmp_path):
    _,_,database,locations,batch=member_fixture(tmp_path)
    with sqlite3.connect(batch/'registry.sqlite3') as db:db.execute("DELETE FROM conformer WHERE id='c0'")
    with pytest.raises(ValueError,match='lacks 1/6'):
        search.map_artifacts(batch,database,locations['descriptors'])


def test_foreign_pose_is_rejected(tmp_path):
    _,_,database,locations,batch=member_fixture(tmp_path)
    search.map_artifacts(batch,database,locations['descriptors'])
    chunks(tmp_path,[pose(7)])
    with pytest.raises(ValueError,match='outside'):
        search.rank_conformers(tmp_path,dict(templates=[dict(query_id='fixture')]),database,100000)


def test_changed_source_hash_is_rejected(tmp_path):
    _,_,database,locations,batch=member_fixture(tmp_path)
    with sqlite3.connect(batch/'registry.sqlite3') as db:db.execute("UPDATE conformer SET content_sha256='wrong' WHERE id='c0'")
    with pytest.raises(ValueError,match='hash mismatch'):
        search.map_artifacts(batch,database,locations['descriptors'])


@pytest.mark.parametrize('use_tools',[False,True])
@pytest.mark.parametrize('do_dock',[False,True])
@pytest.mark.parametrize('threshold',[False,True])
def test_search_to_docking_pipeline_reports_shortfall_and_reuses_receipt(tmp_path,monkeypatch,use_tools,do_dock,threshold):
    import aidd_agent.consensus_design  # Load SciPy before replacing the docking subprocess.
    profile,source,_,locations,batch=member_fixture(tmp_path)
    plants=plants_profile(tmp_path,source)
    from aidd_agent.block_plants import prepare, run as plants_run
    prepare(tmp_path/'sample/report.json',plants,tmp_path/'original-prepared')
    monkeypatch.setattr('aidd_agent.block_plants.subprocess.run',fake_plants)
    plants_run(tmp_path/'original-prepared/report.json',tmp_path/'original-docking')
    original=tmp_path/'original-docking/report.json'
    full=tmp_path/'full.json';save(full,dict(sample_report=str(tmp_path/'sample/report.json'),source_report=str(original),identity=dict(source_sha256=sha(original),profile_inputs=fingerprint([tmp_path/'sample/report.json']))))
    analysis=tmp_path/'analysis.json';save(analysis,dict(source_report=str(full),ranking=dict(source_report_sha256=sha(full))))
    evidence=tmp_path/'evidence.json';save(evidence,dict(fixture=True))
    query=tmp_path/'query.json'
    save(query,dict(kind='consensus_design',status='complete',target=dict(accession='Q00987'),
        templates=[dict(query_id='fixture')],anchors=[],sources=fingerprint([evidence]),
        design=dict(mandatory_anchors=[],alternative_groups=[],optional_weights={})))
    request=dict(analysis=str(analysis),query=str(query),selection=dict(scheme='E095',receptor='MDM2_fixture',block_ids=['b0']),conformers=100000,dock=True)
    if threshold:
        request['search_engine']='consensus-threshold-v1'
        query_data=read(query)
        query_data['readiness']='ready_for_consensus_funnel'
        query_data['design'].update(minimum_score=.73,minimum_pose_score=.61,coarse_constraints={})
        save(query,query_data)
    cfg=dict(search=dict(batch=str(batch),workers=1,refine_chunk=64),block_evaluation=dict(sampling_profile=str(profile),plants_profile=str(plants)))
    if use_tools:
        tools=tmp_path/'receptor-tools.json'
        settings=read(plants)
        save(tools,dict(plants_executable=settings['executable'],workers=2,timeout_seconds=600))
        cfg['block_evaluation'].pop('plants_profile')
        cfg['block_evaluation']['receptor_tools_profile']=str(tools)
    calls=[]
    if not do_dock:
        request['dock']=False
        if use_tools:
            request['search_policy']='ranked_blocks_until_budget'
            request['selection'].update(method='equiscore',rankings={'EquiScore':[dict(block_id='b0',rank=1)]})
        cfg['block_evaluation'].pop('plants_profile',None)
        cfg['block_evaluation'].pop('receptor_tools_profile',None)
        def forbidden(*args,**kwargs):
            raise AssertionError('Search-only must not inspect or run PLANTS')
        monkeypatch.setattr(search,'inherited_plants_profile',forbidden)
        monkeypatch.setattr('aidd_agent.block_plants.prepare',forbidden)
    def refine(batch_,output,design,ids,workers,chunk):
        if threshold:
            assert design['block_search_mode']=='threshold'
            assert design['design']['minimum_score']==.73
        calls.append(ids.tolist());chunks(output,[pose(int(i)) for i in ids])
    monkeypatch.setattr('aidd_agent.budget_screen.run_refinement',refine)
    monkeypatch.setattr('aidd_agent.block_plants.subprocess.run',fake_plants)
    output=tmp_path/'continuation';output.mkdir()
    result=search.run(request,output,cfg)
    assert result['status']=='complete'
    if do_dock:
        assert result['scored']==6 and len(result['top_candidates'])==6
    else:
        assert result['planned_pairs']==0 and result['receptors']==0 and 'scored' not in result
        assert not (output/'docking').exists()
        assert list(Path(result['mol2_directory']).glob('*.mol2'))
    assert result['unique_molecules']==3 and result['exported_conformers']==6
    assert result['shortfall']==99994 and result['planned_pairs']==(6 if do_dock else 0)
    assert search.run(request,output,cfg)==result and len(calls)==1
    if use_tools and do_dock:
        inherited,_=search.inherited_plants_profile(read(full),tools,'MDM2_fixture',receptor_tools=True)
        assert inherited['workers']==2 and inherited['timeout_seconds']==600
        assert inherited['receptors'][0]['center']==settings['receptors'][0]['center']
        changed=tmp_path/'different-engine';changed.write_text('different executable')
        save(tools,dict(plants_executable=str(changed)))
        with pytest.raises(ValueError,match='engine differs'):
            search.inherited_plants_profile(read(full),tools,'MDM2_fixture',receptor_tools=True)


def test_config_retry_preserves_receipt_and_is_idempotent(tmp_path):
    app=ChatAgent(tmp_path/'chat',start=False)
    try:
        sid=app.new_session();ctx=app.context
        project=project_root(Path(ctx['db']),ctx['user_id'],ctx['project_id'])
        original=campaign.queue(app,sid,project,dict(session=sid,dock=True),'block_search_dock')
        job=app.task(sid,original['task_id'])
        report=Path(job['report']);report.parent.mkdir()
        save(report,dict(status='blocked',steps=dict(campaign=dict(status='blocked',
            error='Configure search.batch, block_evaluation.sampling_profile and plants_profile'))))
        with app.connect() as db: db.execute("UPDATE jobs SET status='blocked' WHERE id=?",(job['id'],))
        # Real execute_step creates this directory even when runtime checks block.
        (report.parent/'campaign/blocks').mkdir(parents=True)
        args=dict(operation='retry_config',task_id=job['id'])
        validate_arguments('block_campaign',args)
        retry=campaign.handle(app,sid,args)
        assert retry['task_id']!=job['id'] and retry['status']=='queued'
        assert campaign.handle(app,sid,args)['task_id']==retry['task_id']
        assert app.task(sid,job['id'])['status']=='blocked'
        with pytest.raises(ValueError): campaign.handle(app,app.new_session(),args)
        with pytest.raises(ValueError): campaign.handle(app,sid,dict(operation='retry_config',task_id=retry['task_id']))
        (report.parent/'campaign/blocks/protocol.json').write_text('{}')
        with pytest.raises(ValueError,match='artifacts'): campaign.handle(app,sid,args)
        (report.parent/'campaign/blocks/protocol.json').unlink()
        (report.parent/'campaign.stage.json').write_text('{}')
        with pytest.raises(ValueError,match='artifacts'): campaign.handle(app,sid,args)
        (report.parent/'campaign.stage.json').unlink()
        (report.parent/'campaign/protocol.json').write_text('{}')
        with pytest.raises(ValueError,match='artifacts'): campaign.handle(app,sid,args)
        (report.parent/'campaign/protocol.json').unlink()
        save(report,dict(status='blocked',steps=dict(campaign=dict(status='blocked',error='Review query chemistry'))))
        with pytest.raises(ValueError,match='Only missing'): campaign.handle(app,sid,args)
    finally:
        app.close()


def test_queued_import_keeps_failure_receipt_instead_of_resubmitting(imported):
    app,sid,project,path,_=imported
    args=dict(operation='import',report=str(path))
    first=campaign.handle(app,sid,args)
    with app.connect() as db:db.execute("UPDATE jobs SET status='failed' WHERE id=?",(first['task_id'],))
    second=campaign.handle(app,sid,args)
    assert second['status']=='failed' and second['reused'] and first['task_id']==second['task_id']
    assert len(app.jobs(sid))==1


def test_prompt_import_stops_after_receipt_and_persists_task(imported):
    app,sid,project,path,_=imported
    calls=[]
    def planner(*args,**kwargs):
        calls.append(1)
        return dict(kind='tool',tool='block_campaign',arguments=dict(operation='import',report=str(path)),purpose='Record completed analysis'),dict(provider='fixture')
    app.domain_planner=planner
    answer=app.ask(sid,'Import the completed EquiScore analysis '+str(path),'deepseek')
    assert 'queued' in answer and 'without inference or docking' in answer
    assert len(calls)==1 and len(app.jobs(sid))==1
    result=app.ask(sid,'/block_campaign {"operation":"list"}','deepseek')
    assert json.loads(result)['tasks'][0]['id']==app.jobs(sid)[0]['id']


def test_start_binds_query_and_conformer_budget_without_bypassing_compute_gate(imported):
    app,sid,project,path,_=imported
    adopted=campaign.handle(app,sid,dict(operation='import',report=str(path)))
    app.execute(app.task(sid,adopted['task_id']))
    query_root=project/'runs/PROMPT-0123456789abcdef';query_root.mkdir()
    execution=query_root/'execution';execution.mkdir()
    save(query_root/'plan.json',dict(fixture=True))
    save(execution/'query.json',dict(kind='consensus_design',status='complete',target=dict(accession='Q00987')))
    save(execution/'report.json',dict(steps=dict(query=dict(action='consensus_design',status='complete',result=dict(report=str(execution/'query.json'))))))
    query_id='0123456789abcdef'
    with app.connect() as db:
        db.execute('INSERT INTO jobs(id,session,request,provider,profile,status,plan,report,log,created) VALUES(?,?,?,?,?,?,?,?,?,?)',
                   (query_id,sid,'Fixture confirmed query','local','','complete',str(query_root/'plan.json'),str(execution/'report.json'),'',1))
    discovered=campaign.handle(app,sid,{})
    assert discovered['query_tasks'][0]['task_id']==query_id
    args=dict(operation='start',task_id=adopted['task_id'],query_task_id=query_id,scheme='E095',receptor='R1',conformers=100000)
    queued=campaign.handle(app,sid,args)
    assert campaign.handle(app,sid,args)['task_id']==queued['task_id']
    job=app.task(sid,queued['task_id'])
    plan=json.loads(Path(job['plan']).read_text())['plan']['steps'][0]
    assert plan['action']=='block_search_dock'
    request=json.loads((project/'block-campaigns'/(plan['params']['request_id']+'.json')).read_text())
    assert request['conformers']==100000 and request['query_run']==query_root.name
    assert request['selection']['candidate_unit']=='conformer' and request['dock'] is True
    app.execute(job)
    final=app.task(sid,job['id'])
    assert final['status']=='blocked', Path(final['log']).read_text()
