"""Local CPU geometry, durable orchestration and saved-library reuse checks."""
import io
import json
from pathlib import Path

import numpy as np
import pytest

from aidd_agent.chat_agent import ChatAgent
from aidd_agent.chat_uploads import upload
from aidd_agent import candidate_sets, molecule_preparation, research_workflows as flows
from aidd_agent.research_projects import handle as project, memory
from aidd_agent.ligand_library import handle as search, queries
from test_ligand_library import tiny_library


def execute(app,sid,receipt):
    job=app.task(sid,receipt['task_id']);app.execute(job)
    assert app.task(sid,job['id'])['status']=='complete',app.task(sid,job['id'])['error']
    report=json.loads(Path(job['report']).read_text())
    return next(iter(report['steps'].values()))['result']


def test_preparation_records_failures_and_preserves_3d(tmp_path):
    from rdkit import Chem
    data=tmp_path/'input.csv';data.write_text('name,smiles\nacetone,CC(=O)C\nbad,not_smiles\nisomer,CC(O)C(=O)O\nsalt,C[NH3+].[Cl-]\n')
    result=molecule_preparation.execute(dict(input_path=str(data)),tmp_path/'prepared')
    assert result['prepared_molecules']==1 and result['failed_molecules']==3
    assert len(result['prepared_queries'])==1
    details=json.loads(Path(result['details']).read_text())
    sdf=Path(details[0]['conformers'][0]['sdf'])
    again=molecule_preparation.execute(dict(input_path=str(sdf),conformers=3),tmp_path/'preserved')
    second=json.loads(Path(again['details']).read_text())[0]
    assert second['coordinates']=='original_3d_preserved' and len(second['conformers'])==1
    a=Chem.SDMolSupplier(str(sdf),removeHs=False)[0]
    b=Chem.SDMolSupplier(second['conformers'][0]['sdf'],removeHs=False)[0]
    np.testing.assert_array_equal(a.GetConformer().GetPositions(),b.GetConformer().GetPositions())
    repeated=molecule_preparation.execute(dict(smiles='CC(=O)C'),tmp_path/'repeat')
    c=Chem.SDMolSupplier(repeated['prepared_queries'][0]['sdf'],removeHs=False)[0]
    np.testing.assert_array_equal(a.GetConformer().GetPositions(),c.GetConformer().GetPositions())


def test_project_workflow_real_worker_and_candidate_reuse(tmp_path,monkeypatch):
    batch,raw=tiny_library(tmp_path)
    runtime=tmp_path/'runtime.json';runtime.write_text(json.dumps(dict(search=dict(batch=str(batch),workers=1,refine_chunk=1))))
    app=ChatAgent(tmp_path/'workspace',start=False,runtime=runtime,allow_compute=True)
    try:
        sid=app.new_session();sid=project(app,sid,dict(operation='create',name='CIP2A',target='Q8TCG1',objective='Find candidates'))['session']
        attached=upload(app,sid,'ligand.smi',io.BytesIO(b'CC(=O)C acetone'),14)
        nodes=[
            dict(id='prepare',tool='molecule_prepare',arguments=dict(attachment_id=attached['id']),depends_on=[]),
            dict(id='search',tool='ligand_library',arguments=dict(operation='start',task_id='{{task:prepare}}',query_id='{{query:prepare}}',export_molecules=1),depends_on=['prepare']),
            dict(id='collect',tool='candidate_sets',arguments=dict(operation='import',task_id='{{task:search}}',name='All candidates'),depends_on=['search']),
            dict(id='export',tool='candidate_sets',arguments=dict(operation='export',set_id='{{set:collect}}',count=2),depends_on=['collect'])]
        receipt=flows.handle(app,sid,dict(operation='create',objective='Prepare, search, retain and export',nodes=nodes))
        flows.tick(app);row=flows.rows(app,sid)[0];assert row['nodes'][0]['status']=='submitted'
        execute(app,sid,dict(task_id=row['nodes'][0]['task_id']))
        # New chat and process retain the goal, target and prepared query.
        new_sid=app.new_session(app.context_for(sid)['project_id']);app.close()
        app=ChatAgent(tmp_path/'workspace',start=False,runtime=runtime,allow_compute=True)
        assert len(queries(app,new_sid,row['nodes'][0]['task_id']))==1
        for _ in range(3):
            flows.tick(app);row=flows.rows(app,new_sid)[0]
            assert row['status']=='active',row
            pending=next(n for n in row['nodes'] if n['status']=='submitted')
            execute(app,new_sid,dict(task_id=pending['task_id']))
        flows.tick(app);assert flows.rows(app,new_sid)[0]['status']=='complete'
        assert len(app.project_jobs(new_sid))==4
        flows.tick(app);assert len(app.project_jobs(new_sid))==4
        item=candidate_sets.handle(app,new_sid,dict(operation='list'))['sets'][0];assert item['count']==2
        def change(**args): return execute(app,new_sid,candidate_sets.handle(app,new_sid,args))
        subset=change(operation='filter',set_id=item['set_id'],name='First',limit=1)
        assert subset['count']==1
        assert change(operation='union',set_id=subset['set_id'],other_set_id=item['set_id'],name='Union')['count']==2
        assert change(operation='intersection',set_id=subset['set_id'],other_set_id=item['set_id'],name='Intersection')['count']==1
        difference=change(operation='difference',set_id=item['set_id'],other_set_id=subset['set_id'],name='Remaining')
        assert difference['count']==1
        assert change(operation='diverse',set_id=item['set_id'],name='Scaffolds')['count']==1
        member=candidate_sets.handle(app,new_sid,dict(operation='show',set_id=subset['set_id']))['members'][0]['molecule_id']
        candidate_sets.handle(app,new_sid,dict(operation='label',set_id=item['set_id'],molecule_ids=[member],label='tested'))
        assert change(operation='filter',set_id=item['set_id'],name='Untested',exclude_labels=['tested'])['count']==1
        # Export reuses stored ranking/poses; forbid another search in this process.
        import aidd_agent.ligand_library as library
        monkeypatch.setattr(library,'execute',lambda *a,**k:pytest.fail('Unexpected rescan'))
        result=change(operation='export',set_id=difference['set_id'],count=1)
        assert result['exported_molecules']==1
        assert len(list(Path(result['original_mol2_directory']).glob('*.mol2')))==1
        assert memory(app,new_sid)['candidate_sets']
        foreign=project(app,new_sid,dict(operation='create',name='MDM2',target='MDM2',objective='Other target'))['session']
        with pytest.raises(ValueError):candidate_sets.handle(app,foreign,dict(operation='show',set_id=item['set_id']))
        with pytest.raises(ValueError):flows.handle(app,foreign,dict(operation='pause',workflow_id=receipt['workflow_id']))
    finally:app.close()


def test_flow_pause_restart_and_uncertain_dispatch(tmp_path,monkeypatch):
    app=ChatAgent(tmp_path,start=False)
    try:
        sid=app.new_session();args=dict(operation='create',objective='Prepare a molecule',nodes=[dict(id='prep',tool='molecule_prepare',arguments=dict(smiles='CCO'),depends_on=[])])
        receipt=flows.handle(app,sid,args,request='Prepare CCO')
        control=dict(workflow_id=receipt['workflow_id'])
        flows.handle(app,sid,dict(operation='pause',**control));flows.tick(app);assert not app.jobs(sid)
        flows.handle(app,sid,dict(operation='resume',**control));flows.tick(app);assert len(app.jobs(sid))==1
        row=flows.rows(app,sid)[0];node=row['nodes'][0];jid=node.pop('task_id');node['status']='dispatching';flows.save(app,row)
        app.close();app=ChatAgent(tmp_path,start=False)
        flows.tick(app);row=flows.rows(app,sid)[0];assert row['status']=='blocked' and 'Uncertain' in row['error']
        assert len(app.jobs(sid))==1
        with pytest.raises(ValueError):flows.handle(app,sid,dict(operation='resume',**control))
        flows.handle(app,sid,dict(operation='reconcile',node_id='prep',task_id=jid,**control))
        flows.tick(app);assert len(app.jobs(sid))==1
        flows.handle(app,sid,dict(operation='cancel',**control));flows.tick(app);assert len(app.jobs(sid))==1
    finally:app.close()


def test_multiple_prepared_queries_block_until_explicit_resolution(tmp_path):
    app=ChatAgent(tmp_path,start=False)
    try:
        sid=app.new_session();payload=b'CCO ethanol\nCC(=O)C acetone\n'
        attached=upload(app,sid,'molecules.smi',io.BytesIO(payload),len(payload))
        receipt=flows.handle(app,sid,dict(operation='create',objective='Prepare and search selected molecule',nodes=[
            dict(id='prep',tool='molecule_prepare',arguments=dict(attachment_id=attached['id']),depends_on=[]),
            dict(id='search',tool='ligand_library',arguments=dict(operation='start',task_id='{{task:prep}}',query_id='{{query:prep}}',export_molecules=1),depends_on=['prep'])]))
        flows.tick(app);node=flows.rows(app,sid)[0]['nodes'][0];execute(app,sid,node)
        flows.tick(app);row=flows.rows(app,sid)[0]
        assert row['status']=='blocked' and 'explicit query_id' in row['error'] and len(app.jobs(sid))==1
        query=queries(app,sid,node['task_id'])[1]['query_id']
        flows.handle(app,sid,dict(operation='resolve',workflow_id=receipt['workflow_id'],node_id='search',
            arguments=dict(operation='start',task_id='{{task:prep}}',query_id=query,export_molecules=1)))
        flows.tick(app);assert len(app.jobs(sid))==2 and flows.rows(app,sid)[0]['nodes'][1]['status']=='submitted'
    finally:app.close()
