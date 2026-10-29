import copy
import json
from pathlib import Path

import numpy as np
import pytest

from aidd_agent.chat_agent import ChatAgent
from aidd_agent.prompt_workflow import create_plan, run_plan, initialize_context
from aidd_agent.prompt_plan import validate_plan
from aidd_agent.prompt_smoke import OfflineServices
from aidd_agent.target_query import prepare, design
from aidd_agent.target_search_workflow import handle


def af3_plan(**params):
    return dict(version=1,summary='Prepare user-supplied complex',clarifications=[],steps=[
        dict(id='input',action='af3_prepare',params=dict(name='custom_complex',**params))])


def execute_plan(tmp_path, plan, prompt):
    ctx=initialize_context(tmp_path/'storage')
    path=create_plan(Path(ctx['db']),ctx['user_id'],ctx['project_id'],prompt,local_plan=plan)
    result,_=run_plan(Path(ctx['db']),ctx['user_id'],ctx['project_id'],path,services=OfflineServices())
    return path,result


def test_literal_heteromer_seeds_ligand_no_identity_lookup(tmp_path):
    sequences=['ACDEFGHIK','MNPQRSTVWY','ACDEFGHIK']
    smiles='C[C@H](O)C(=O)O'
    plan=af3_plan(protein_sequences=sequences,seed_count=4,seed_start=7,ligand_smiles=[smiles])
    path,result=execute_plan(tmp_path,plan,'Predict '+ ' '.join(sequences)+' with SMILES '+smiles)
    assert result['status']=='complete'
    payload=json.loads((path.parent/'execution/input/af3-input.json').read_text())
    assert payload['modelSeeds']==[7,8,9,10]
    assert [row['protein']['sequence'] for row in payload['sequences'][:3]]==sequences
    assert [next(iter(row.values()))['id'] for row in payload['sequences']]==list('ABCD')
    evidence=json.loads((path.parent/'execution/input/input-provenance.json').read_text())
    assert all(row['accession'] is None and row['provenance']=='literal_user_sequence_not_database_verified' for row in evidence['proteins'])


def test_sequence_cannot_be_invented_by_router(tmp_path):
    plan=af3_plan(protein_sequences=['ACDEFGHIK'])
    with pytest.raises(ValueError,match='original user'):
        execute_plan(tmp_path,plan,'ACDEFGHIK\nOriginal user request (preserve literal SMILES):\nPredict an unknown protein')


@pytest.mark.parametrize('params',[
    dict(protein_sequences=['ACDX']),dict(protein_sequences=[]),
    dict(protein_sequences=['ACD'],protein_steps=['unknown']),dict(protein_steps=['unknown']),
    dict(protein_sequences=['ACD','EFG'],protein_copies=2),
    dict(protein_sequences=['ACD'],seeds=[1],seed_count=2),
    dict(protein_sequences=['ACD'],seed_count=True),dict(protein_sequences=['ACD'],seed_start=3),
    dict(protein_sequences=['ACD','EFG'],start=1,end=3),
])
def test_invalid_generic_af3_inputs(params):
    with pytest.raises(ValueError):validate_plan(af3_plan(**params))


def molecular_files(tmp_path):
    from rdkit import Chem
    from rdkit.Chem import AllChem
    tmp_path.mkdir(parents=True,exist_ok=True)
    mol=Chem.AddHs(Chem.MolFromSmiles('CCO'));assert AllChem.EmbedMolecule(mol,randomSeed=17)==0
    ligand=tmp_path/'reference.sdf'
    with Chem.SDWriter(str(ligand)) as writer:writer.write(mol)
    o=np.asarray(mol.GetConformer().GetAtomPosition(2))
    atom=o+np.array([0.,0.,3.])
    receptor=tmp_path/'receptor.pdb'
    receptor.write_text(f'ATOM      1  O   GLY B  47    {atom[0]:8.3f}{atom[1]:8.3f}{atom[2]:8.3f}  1.00 20.00           O  \nTER\nEND\n')
    return receptor,ligand


def params(receptor,ligand,name='OTHER_TARGET'):
    return dict(name=name,receptor=str(receptor),reference_ligand=str(ligand),receptor_chains=['B'])


def test_arbitrary_target_query_and_same_pose_contacts(tmp_path):
    receptor,ligand=molecular_files(tmp_path/'inputs')
    report=prepare(params(receptor,ligand,'UNRELATED_ENZYME'),tmp_path/'query')
    assert report['target']=='UNRELATED_ENZYME'
    assert any(a['feature_class']=='HBD' for a in report['query']['anchors'])
    result=design(tmp_path/'query/report.json',tmp_path/'design',contacts=[dict(chain='B',residue='47',feature_class='HBD',distance=4.)])
    assert result['reference_control']==dict(pocket_passed=True,anchor_rule_passed=True)
    from aidd_agent.gaussian_batch import _load_query
    from aidd_agent.guided_filters import GuidedFilter
    _,query=_load_query(Path(report['query']['query_npz']))
    expanded=dict(query,anchor_feature_indices=np.asarray(sorted({a['feature_index'] for a in report['query']['anchors']})))
    gate=GuidedFilter(report['query'],expanded,result['design'])
    transforms=np.asarray([np.eye(4),np.eye(4)]);transforms[1,:3,3]=[50.,0.,0.]
    assert gate.pocket_mask(query['shape_points'],transforms).tolist()==[True,False]
    with pytest.raises(ValueError,match='no reference-supported'):
        design(tmp_path/'query/report.json',tmp_path/'bad',contacts=[dict(chain='B',residue='999')])


def test_reference_must_be_in_receptor_coordinates(tmp_path):
    from rdkit import Chem
    receptor,ligand=molecular_files(tmp_path/'inputs')
    mol=next(iter(Chem.SDMolSupplier(str(ligand),removeHs=False)))
    conf=mol.GetConformer()
    for i,p in enumerate(conf.GetPositions()):conf.SetAtomPosition(i,p+np.array([100.,0.,0.]))
    with Chem.SDWriter(str(ligand)) as writer:writer.write(mol)
    with pytest.raises(ValueError,match='outside'):
        prepare(params(receptor,ligand),tmp_path/'query')


def test_required_residues_cannot_be_combined_across_poses(tmp_path):
    from aidd_agent.guided_filters import GuidedFilter
    receptor=tmp_path/'receptor.npz';np.savez(receptor,points=np.array([[-100.,0,0]]))
    settings=dict(mandatory_anchors=[],alternative_groups=[],receptor_npz=str(receptor),
        pocket=dict(minimum_heavy_atom_distance=1.2,maximum_clashing_fraction=0.),
        direct_contact_groups=[dict(points=[[0.,0,2]],distance=2.5),dict(points=[[10.,0,2]],distance=2.5)])
    gate=GuidedFilter(dict(anchors=[]),{},settings)
    transforms=np.asarray([np.eye(4),np.eye(4)]);transforms[1,0,3]=10.
    assert gate.pocket_mask(np.array([[0.,0,0]]),transforms).tolist()==[False,False]


def test_explicit_cofactor_and_residue_identity(tmp_path):
    from aidd_agent.target_query import load_receptor
    receptor,ligand=molecular_files(tmp_path/'inputs')
    text=receptor.read_text().replace('TER\nEND\n','')
    text+='HETATM    2 MG    MG B 201      30.000  30.000  30.000  1.00 20.00          MG  \nEND\n'
    receptor.write_text(text)
    assert len(load_receptor(receptor,['B'])[1])==1
    assert len(load_receptor(receptor,['B'],[dict(chain='B',residue='201',name='MG')])[1])==2
    with pytest.raises(ValueError,match='not found'):
        load_receptor(receptor,['B'],[dict(chain='B',residue='202',name='GDP')])
    prepare(params(receptor,ligand),tmp_path/'query')
    with pytest.raises(ValueError,match='no reference-supported'):
        design(tmp_path/'query/report.json',tmp_path/'design',contacts=[dict(chain='B',residue='47',residue_name='HIS')])


def test_literal_complex_uses_existing_af3_execution_adapter(tmp_path):
    import sys
    ctx=initialize_context(tmp_path/'storage')
    payload=af3_plan(protein_sequences=['ACDEFGHIK','MNPQRSTVWY'],seed_count=3,ligand_smiles=['CCO'])
    payload['steps'].append(dict(id='predict',action='af3_run',params=dict(input_step='input')))
    path=create_plan(Path(ctx['db']),ctx['user_id'],ctx['project_id'],'Predict ACDEFGHIK and MNPQRSTVWY with SMILES CCO',local_plan=payload)
    runner=tmp_path/'runner.py';runner.write_text('# offline fixture')
    models=tmp_path/'models';models.mkdir();database=tmp_path/'databases';database.mkdir()
    profile=tmp_path/'af3.json';profile.write_text(json.dumps(dict(backend='alphafold3',python=sys.executable,
        runner=str(runner),model_parameters=str(models),databases=str(database),license_acknowledged=True)))
    runtime=tmp_path/'runtime.json';runtime.write_text(json.dumps(dict(af3_profile=str(profile))))
    called=[]
    class FakeAF3(OfflineServices):
        @staticmethod
        def run_command(argv,log):
            called.append(argv)
            input_path=Path(next(s.split('=',1)[1] for s in argv if s.startswith('--json_path=')))
            document=json.loads(input_path.read_text())
            assert document['modelSeeds']==[1,2,3] and len(document['sequences'])==3
            output=Path(next(s.split('=',1)[1] for s in argv if s.startswith('--output_dir=')));output.mkdir(parents=True)
            (output/'fixture_model.cif').write_text('data_fixture\n')
            (output/'fixture_summary_confidences.json').write_text(json.dumps(dict(iptm=.7)))
            log.write_text('Offline simulated inference, not a GPU run')
    result,_=run_plan(Path(ctx['db']),ctx['user_id'],ctx['project_id'],path,runtime=runtime,services=FakeAF3(),allow_compute=True)
    assert result['status']=='complete' and len(called)==1


def test_pdb_ligand_instance_path_uses_explicit_identity(tmp_path,monkeypatch):
    from rdkit import Chem
    from aidd_agent import chemistry_prep
    receptor,ligand=molecular_files(tmp_path/'inputs')
    mol=next(iter(Chem.SDMolSupplier(str(ligand),removeHs=False)))
    app=ChatAgent(tmp_path/'chat',start=False)
    try:
        sid=app.new_session()
        receipt=handle(app,sid,dict(operation='prepare',name='SYNTHETIC_PDB',pdb_id='1ABC',receptor_chains=['B'],
            ligand_instance=dict(ccd='ABC',chain='C',residue='501')))
        class FakePDB(OfflineServices):
            @staticmethod
            def pdb_fetch(code,output):
                assert code=='1ABC'
                return dict(path=str(receptor))
        monkeypatch.setattr(chemistry_prep,'enumerate_ligand_instances',lambda *_:[dict(chain_id='C',residue_number='501',altloc='',insertion_code='',model='1',atoms=[])])
        monkeypatch.setattr(chemistry_prep,'fetch_ccd',lambda *_:tmp_path/'fixture.cif')
        monkeypatch.setattr(chemistry_prep,'_ccd_molecule',lambda *_:mol)
        ctx=app.context
        result,_=run_plan(Path(ctx['db']),ctx['user_id'],ctx['project_id'],Path(receipt['plan']),services=FakePDB())
        assert result['status']=='complete'
    finally:app.close()


def test_general_search_queue_compute_gate_and_existing_engine(tmp_path,monkeypatch):
    receptor,ligand=molecular_files(tmp_path/'inputs')
    app=ChatAgent(tmp_path/'chat',start=False,result_import_roots=[receptor.parent])
    try:
        sid=app.new_session()
        receipt=handle(app,sid,dict(operation='start',**params(receptor,ligand),required_contacts=[dict(chain='B',residue='47')]))
        task=app.task(sid,receipt['task_id']);ctx=app.context
        result,_=run_plan(Path(ctx['db']),ctx['user_id'],ctx['project_id'],Path(task['plan']),services=OfflineServices())
        assert result['steps']['query']['status']=='complete' and result['status']=='blocked'
        # A fresh configuration gets a new sealed plan; never bypass resume identity.
        receipt=handle(app,sid,dict(operation='start',**params(receptor,ligand,'CUSTOM_MUTANT')))
        task=app.task(sid,receipt['task_id'])
        runtime=tmp_path/'runtime.json';runtime.write_text(json.dumps(dict(search=dict(batch=str(tmp_path/'library'),workers=2,refine_chunk=64))))
        calls=[]
        from aidd_agent import guided_workflow
        def scan(action,source,output,options,cfg,compute):
            saved=json.loads(Path(source).read_text());calls.append(saved)
            assert action=='guided_funnel' and compute and cfg['search']['workers']==2
            assert saved['query']['query_id']=='CUSTOM_MUTANT'
            assert saved['reference_control']['pocket_passed']
            output.mkdir(parents=True)
            value=dict(kind='guided_funnel',status='complete',full_coverage=True,outputs={},counts=dict(input_conformers=3))
            (output/'report.json').write_text(json.dumps(value));return value
        monkeypatch.setattr(guided_workflow,'execute',scan)
        result,_=run_plan(Path(ctx['db']),ctx['user_id'],ctx['project_id'],Path(task['plan']),services=OfflineServices(),runtime=runtime,allow_compute=True)
        assert result['status']=='complete' and len(calls)==1
        assert result['steps']['screen']['result']['full_coverage'] is True
    finally:app.close()


def test_foreign_inputs_and_changed_imports_are_rejected(tmp_path):
    receptor,ligand=molecular_files(tmp_path/'inputs')
    app=ChatAgent(tmp_path/'chat',start=False)
    try:
        sid=app.new_session()
        with pytest.raises(ValueError,match='outside'):
            handle(app,sid,dict(operation='prepare',**params(receptor,ligand)))
        app.result_import_roots=(receptor.parent,)
        receipt=handle(app,sid,dict(operation='prepare',**params(receptor,ligand)))
        plan=Path(receipt['plan']);envelope=json.loads(plan.read_text())
        identifier=envelope['plan']['steps'][0]['params']['request_id']
        root=plan.parent.parent.parent/'target-search/requests'/identifier
        (root/'reference_ligand.sdf').write_text('changed')
        ctx=app.context
        result,_=run_plan(Path(ctx['db']),ctx['user_id'],ctx['project_id'],plan,services=OfflineServices())
        assert result['status']=='failed'
    finally:app.close()
