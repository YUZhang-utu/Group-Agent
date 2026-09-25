import copy
import hashlib
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from aidd_agent.af3_analysis import ligand_evidence, explain_confidence
from aidd_agent.chat_agent import ChatAgent, validate_route
from aidd_agent.prompt_plan import validate_plan
from aidd_agent.prompt_workflow import create_plan, run_plan, initialize_context, load_runtime
from aidd_agent.prompt_smoke import OfflineServices
from aidd_agent.pymol_bridge import execute_command, selection, validate_view
from aidd_agent.structure_review import structures, contacts


def dump(path,value):
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(value),encoding='utf-8')
    return path


def smiles_plan(smiles):
    return dict(version=1,summary='Prepare supplied ligand complex',clarifications=[],steps=[
        dict(id='protein',action='protein_fetch',params=dict(accession='P30291')),
        dict(id='fold',action='af3_prepare',params=dict(protein_step='protein',name='complex',ligand_smiles=[smiles]))])


def test_smiles_prepare_end_to_end_preserves_chirality_and_charge(tmp_path):
    smiles='C[C@H]([NH3+])C(=O)[O-]'
    ctx=initialize_context(tmp_path/'home')
    plan=create_plan(Path(ctx['db']),ctx['user_id'],ctx['project_id'],'Prepare SMILES '+smiles,
                     local_plan=smiles_plan(smiles))
    report,_=run_plan(Path(ctx['db']),ctx['user_id'],ctx['project_id'],plan,services=OfflineServices())
    assert report['status']=='complete'
    payload=json.loads((plan.parent/'execution/fold/af3-input.json').read_text())
    assert payload['sequences'][1]=={'ligand':{'id':'B','smiles':smiles}}
    evidence=json.loads((plan.parent/'execution/fold/ligand-chemistry.json').read_text())
    assert evidence[0]['unassigned_tetrahedral_centers']==[]
    with pytest.raises(ValueError,match='exactly'):
        create_plan(Path(ctx['db']),ctx['user_id'],ctx['project_id'],'Prepare aspirin',local_plan=smiles_plan(smiles))
    with pytest.raises(ValueError,match='exactly'):
        create_plan(Path(ctx['db']),ctx['user_id'],ctx['project_id'],
                    'Router supplied '+smiles+'\nOriginal user request (preserve literal SMILES):\nPrepare aspirin',
                    local_plan=smiles_plan(smiles))


@pytest.mark.parametrize('smiles',['bad_smiles','C.C',''])
def test_invalid_chemistry(smiles):
    with pytest.raises(ValueError):ligand_evidence(smiles)


def test_confidence_does_not_fabricate_metrics():
    result=explain_confidence(dict(iptm=.72,chain_ids=['A','B'],has_clash=True))
    assert result['metrics']['iptm']['value']==.72
    assert 'ptm' not in result['metrics']
    assert 'affinity' in result['metrics']['iptm']['meaning']


@pytest.mark.parametrize('view',[
    dict(operation='run',code='print(1)'),dict(operation='color',color='red; quit'),
    dict(operation='sticks',chain='A or all'),dict(operation='zoom',objects=['all']),
    dict(operation='contacts',cutoff=float('nan')),dict(operation='rotate',angle=float('inf'))])
def test_view_rejects_unstructured_code(view):
    with pytest.raises(ValueError):validate_view(view)


class FakeCmd:
    def __init__(self):self.loaded=[];self.calls=[];self.xyz=np.array([[1.,2.,3.]])
    def delete(self,obj):self.calls.append(('delete',obj))
    def load(self,path,obj):self.loaded.append(obj)
    def get_names(self,kind):return self.loaded
    def get_coords(self,*a,**k):return self.xyz
    def load_coords(self,xyz,*a,**k):self.xyz=xyz
    def get_distance(self,*a,**k):return 3.2
    def count_atoms(self,*a):return 1
    def find_pairs(self,*a,**k):return [(('v001',1),('v001',2))]
    def get_model(self,*a,**k):return SimpleNamespace(atom=[SimpleNamespace(name='N',chain='A',resi='23',resn='TRP')])
    def __getattr__(self,name):return lambda *a,**k:self.calls.append((name,a,k))


def test_pymol_uses_recorded_transform_and_exports_pairs(tmp_path):
    path=tmp_path/'model.cif';path.write_text('fixture')
    matrix=np.eye(4);matrix[0,3]=10
    catalog=[dict(id='v001',path=str(path),sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                  label='test',ligand=dict(chain='B'),transform=matrix.tolist())]
    cmd=FakeCmd()
    execute_command(cmd,dict(id='1-abcdefabcdef',view=dict(operation='open')),catalog,tmp_path)
    assert np.array_equal(cmd.xyz,[[11.,2.,3.]])
    result=execute_command(cmd,dict(id='2-abcdefabcdef',view=dict(operation='polar_contacts')),catalog,tmp_path)
    assert result['contact_count']==1
    assert Path(result['artifacts']['contacts_csv']).is_file()
    assert result['examples'][0]['kind']=='pymol_polar_candidate'
    path.write_text('changed')
    with pytest.raises(ValueError,match='changed'):
        execute_command(cmd,dict(id='3-abcdefabcdef',view=dict(operation='open')),catalog,tmp_path)


def test_catalog_ownership_and_confidence_chat(tmp_path,monkeypatch):
    app=ChatAgent(tmp_path/'home',router=lambda *a:dict(intent='run',request='Prepare protein',task_id=None,message=''),start=False)
    sid=app.new_session();app.ask(sid,'Prepare protein','deepseek');job=app.task(sid)
    from aidd_agent.prompt_workflow import project_root
    ctx=app.context;project=project_root(Path(ctx['db']),ctx['user_id'],ctx['project_id'])
    root=project/'runs/PROMPT-aaaaaaaaaaaaaaaa';root.mkdir(parents=True)
    model=root/'prediction.cif';model.write_text('fixture')
    report=dump(root/'execution/report.json',dict(steps=dict(fold=dict(action='af3_run',status='complete',result=dict(model=str(model),confidence=dict(iptm=.8))))))
    app.update(job['id'],status='complete',plan=str(root/'plan.json'),report=str(report))
    assert '.8' in app.ask(sid,'/confidence '+job['id'],'deepseek')
    with pytest.raises(ValueError):app.ask(app.new_session(),'/confidence '+job['id'],'deepseek')
    from aidd_agent import pymol_bridge
    requests=[]
    monkeypatch.setattr(pymol_bridge,'submit',lambda root,view,catalog,executable:requests.append(view) or dict(status='queued'))
    app.ask(sid,'/pymol '+job['id'],'deepseek')
    app.ask(sid,'/view {"operation":"surface","target":"protein"}','deepseek')
    assert requests[-1]['target']=='protein'
    from aidd_agent import pymol_agent
    planned=[]
    monkeypatch.setattr(pymol_agent,'run_agent',lambda request,*args:planned.append(request) or dict(status='complete'))
    assert 'complete' in app.ask(sid,'/pymol_agent Show the ligand pocket with mixed colors','deepseek')
    assert planned==['Show the ligand pocket with mixed colors']
    with pytest.raises(ValueError,match='Open a completed task'):
        app.ask(app.new_session(),'/pymol_agent Show all ligands','deepseek')
    external=tmp_path/'outside.cif';external.write_text('fixture')
    dump(report,dict(steps=dict(fold=dict(action='af3_run',status='complete',result=dict(model=str(external))))))
    with pytest.raises(ValueError):structures(app.task(sid,job['id']),project)
    app.close()


def test_runtime_and_experimental_seed_route(tmp_path):
    runtime=dump(tmp_path/'runtime.json',dict(pymol=dict(executable='/opt/pymol')))
    assert load_runtime(runtime)[0]['pymol']['executable']=='/opt/pymol'
    route=dict(intent='budget',request='',message='',task_id=None,
               budget=dict(seed_search=dict(backend='batched',max_pair_seeds=512)))
    assert validate_route(route)==route
    route['budget']['seed_search']['backend']='unknown'
    with pytest.raises(ValueError):validate_route(route)


def test_pocket_elements_chain_removal_and_contact_overview(tmp_path):
    cmd=FakeCmd();cmd.loaded=['v001']
    catalog=[dict(id='v001',label='fixture',ligand=dict(chain='B'))]
    def run(view):
        return execute_command(cmd,dict(id='7-abcdefabcdef',view=view),catalog,tmp_path)
    result=run(dict(operation='pocket_view',radius=5))
    assert result['residues']['v001'][0]['residue']=='23'
    assert Path(result['artifacts']['residues_json']).exists()
    assert any(c[0]=='show' and 'within 5' in str(c) for c in cmd.calls)
    assert any(c[0]=='color' and 'elem O' in str(c) and 'red' in str(c) for c in cmd.calls)
    assert run(dict(operation='chains'))['chains']['v001']==[dict(chain='A',near_ligand=True)]
    result=run(dict(operation='remove_chain',objects=['v001'],chain='A'))
    assert result['removed_atoms']=={'v001':1}
    removed=[c for c in cmd.calls if c[0]=='remove']
    assert 'polymer.protein' in str(removed) and 'chain A' in str(removed)


@pytest.mark.parametrize('view',[
    dict(operation='remove_chain',chain='A'),dict(operation='remove_chain',objects=['v001']),
    dict(operation='color',scheme='exec'),dict(operation='pocket_view',radius=float('nan'))])
def test_new_view_parameters_reject_ambiguity(view):
    with pytest.raises(ValueError):validate_view(view)


def test_consensus_catalog_and_filtered_ledger(tmp_path):
    structure=tmp_path/'structures/model.cif'
    structure.parent.mkdir();structure.write_text('fixture')
    query=tmp_path/'queries/Q/aligned.npz'
    digest=hashlib.sha256(structure.read_bytes()).hexdigest()
    dump(query.parent/'native.manifest.json',dict(source=dict(mmcif=str(structure),mmcif_sha256=digest)))
    qid='3JZK:YIN:A:1'; transform=np.eye(4).tolist()
    ledger=dump(tmp_path/'contacts.json',dict(complexes=[dict(query_id=qid,pairs=[
        dict(target_residue=96,protein_author_residue='96',classes=['heavy_atom_proximity']),
        dict(target_residue=23,protein_author_residue='23',classes=['heavy_atom_proximity'])])]))
    child=dump(tmp_path/'consensus/report.json',dict(kind='structure_consensus',
        proposed_template_ids=[qid],prepared_complexes=[dict(query_id=qid,query_npz=str(query))],
        cohort=dict(admitted=[dict(query_id=qid,admission=dict(transform=transform))]),
        contact_evidence=dict(path=str(ledger)),sources={str(ledger):hashlib.sha256(ledger.read_bytes()).hexdigest()}))
    report=dump(tmp_path/'execution/report.json',dict(steps=dict(s=dict(action='structure_consensus',status='complete',result=dict(report=str(child))))))
    job=dict(status='complete',report=str(report),plan=str(tmp_path/'plan.json'))
    catalog=structures(job,tmp_path)
    assert catalog[0]['transform']==transform and catalog[0]['ligand']['resn']=='YIN'
    result=contacts(job,tmp_path,residue='96')
    assert result['total_matching_pairs']==1 and result['pairs'][0]['target_residue']==96
    ledger.write_text('{}')
    with pytest.raises(ValueError,match='hash'):contacts(job,tmp_path)
