"""Recovery of original CLI ligand-query references, with no new science run."""
import json
from pathlib import Path

import pytest

from aidd_agent.block_campaign import handle
from aidd_agent.block_saved_inputs import adopt_query
from aidd_agent.chat_agent import ChatAgent
from aidd_agent.expanded_wee1 import fingerprint
from aidd_agent.final_work_blocks import sha
from aidd_agent.gaussian_batch import write_gaussian_query
from aidd_agent.joint_spatial_profiles import save
from aidd_agent.prompt_workflow import project_root
from aidd_agent.domain_tools import DomainTools


def query_package(root):
    root.mkdir(parents=True)
    source={}
    for key in ('mmcif','ccd','query_manifest'):
        file=root/(key+'.json');file.write_text('{}')
        source[key]=str(file);source[key+'_sha256']=sha(file)
    source['query_id']='7NA2:1I0:A:201'
    values=dict(shape_points=[[0,0,0],[1,0,0]],feature_points=[[0,0,0]],
                feature_types=[1],anchored_weights=[1.],anchor_feature_indices=[0])
    write_gaussian_query(root/'native.npz',**values,source=source)
    write_gaussian_query(root/'aligned.npz',**values,source=dict(query_id=source['query_id'],protein_alignment=dict(transform=[[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]])))
    write_gaussian_query(root/'consensus.npz',**values,source=dict(fixture=True))
    report=dict(kind='consensus_design',status='complete',target=dict(accession='Q00987'),
        readiness='ready_for_consensus_funnel',templates=[dict(query_id=source['query_id'],query_npz=str(root/'aligned.npz'))],
        consensus_npz=str(root/'consensus.npz'),sources=fingerprint([Path(source['mmcif'])]),
        anchors=[],anchor_order=[],design=dict(mandatory_anchors=[],alternative_groups=[],optional_weights={}))
    save(root/'report.json',report)
    return root/'report.json'


def test_discover_old_cli_query_and_import_as_owned_task(tmp_path,monkeypatch):
    monkeypatch.setenv('PYTHONPATH',str(Path(__file__).resolve().parents[1]/'src'))
    library_parent=tmp_path/'old-workstation'
    source=query_package(library_parent/'e059-budget/adopted-design')
    runtime=tmp_path/'runtime.json';save(runtime,dict(search=dict(batch=str(library_parent/'library'),workers=1,refine_chunk=64)))
    app=ChatAgent(tmp_path/'chat',runtime=runtime,start=False)
    try:
        sid=app.new_session()
        assert handle(app,sid,{})['query_tasks']==[]
        found=handle(app,sid,dict(operation='discover'))
        assert len(found['candidates'])==1 and found['computation_launched'] is False
        row=found['candidates'][0]
        assert row['details']['query_ids']==['7NA2:1I0:A:201']
        args=dict(operation='import_query',candidate_id=row['candidate_id'])
        # No literal external path is needed after trusted discovery.
        queued=DomainTools(app,sid,'Reuse the original co-crystal ligands','deepseek').call('block_campaign',args)
        assert handle(app,sid,args)['task_id']==queued['task_id']
        app.execute(app.task(sid,queued['task_id']))
        job=app.task(sid,queued['task_id'])
        assert job['status']=='complete',Path(job['log']).read_text()
        queries=handle(app,sid,{})['query_tasks']
        assert len(queries)==1 and queries[0]['task_id']==queued['task_id']
        result=json.loads(Path(queries[0]['report']).read_text())
        assert result['templates']==json.loads(source.read_text())['templates']
        refs=json.loads(Path(result['outputs']['reference_ligands']).read_text())['references']
        assert refs[0]['ligand_ccd']=='1I0' and refs[0]['ligand_chain']=='A'
        assert refs[0]['native_query']==str(source.parent/'native.npz')
        with pytest.raises((ValueError,FileNotFoundError)):
            handle(app,app.new_session(),args)
    finally:
        app.close()


def test_import_rejects_modified_cocrystal_coordinates(tmp_path):
    source=query_package(tmp_path/'old/adopted-design')
    request=dict(saved_query=str(source),saved_query_sha256=sha(source))
    (source.parent/'native.npz').write_bytes(b'changed')
    output=tmp_path/'out';output.mkdir()
    with pytest.raises(ValueError,match='checksum'):
        adopt_query(request,output)
    assert not (output/'report.json').exists()


def test_discover_score_path_without_requiring_user_to_repeat_it(tmp_path):
    app=ChatAgent(tmp_path/'chat',start=False)
    try:
        sid=app.new_session();ctx=app.context
        project=project_root(Path(ctx['db']),ctx['user_id'],ctx['project_id'])
        root=project/'runs/PROMPT-0123456789abcdef/execution/equiscore-v3/equiscore-analysis';root.mkdir(parents=True)
        save(root/'report.json',dict(kind='block_analyze',status='complete',source_report=str(root.parent/'full/report.json'),
            ranking=dict(score='EquiScore',better='higher',ranking_unit='molecule',top_n=10)))
        candidates=handle(app,sid,dict(operation='discover'))['candidates']
        assert len(candidates)==1 and candidates[0]['category']=='scores'
        (root/'report.json').write_text('{}')
        with pytest.raises(ValueError,match='changed'):
            handle(app,sid,dict(operation='import',candidate_id=candidates[0]['candidate_id']))
        assert app.jobs(sid)==[]
    finally:
        app.close()


def test_literal_query_path_must_come_from_user(tmp_path):
    app=ChatAgent(tmp_path/'chat',start=False)
    try:
        sid=app.new_session();source=query_package(tmp_path/'external/adopted-design')
        args=dict(operation='import_query',report=str(source))
        with pytest.raises(ValueError,match='exact analysis'):
            DomainTools(app,sid,'Import an old query','deepseek').call('block_campaign',args)
        result=DomainTools(app,sid,'Reuse '+str(source),'deepseek').call('block_campaign',args)
        assert result['status']=='queued'
    finally:
        app.close()
