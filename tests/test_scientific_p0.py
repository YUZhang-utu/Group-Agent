"""E136 software checks with real CPU chemistry and explicitly synthetic AF3 output."""
import io
import json
from pathlib import Path

import numpy as np
import pytest

from aidd_agent import af3_batch, molecule_preparation, pose_quality, screening_evaluation
from aidd_agent.workflow_errors import Blocked


def ligand(smiles='CCO'):
    from rdkit import Chem
    from rdkit.Chem import AllChem
    mol=Chem.AddHs(Chem.MolFromSmiles(smiles));assert AllChem.EmbedMolecule(mol,randomSeed=47)==0
    AllChem.MMFFOptimizeMolecule(mol)
    return Chem.RemoveHs(mol)


def test_second_chain_clash_and_same_pose_contacts():
    mol=ligand();xyz=mol.GetConformer().GetPositions()
    atoms=[dict(chain='A',residue='1',name='ALA',element='C',xyz=[100,0,0]),
           dict(chain='B',residue='2',name='ALA',element='C',xyz=xyz[0].tolist())]
    result=pose_quality.assess(mol,atoms,'CCO',[dict(chain='A',residue='1',distance=4)])
    assert result['qc_status']=='failed' and result['assembly_clashes'][0]['chain']=='B'
    assert set(result['failures'])=={'assembly_clashes','required_contact'}
    assert result['experimental_status']=='not_assessed'


def test_pose_invalid_identity_and_missing_receptor_not_passed():
    result=pose_quality.assess(ligand(),expected_smiles='CCN')
    assert result['checks']['assembly_clashes'] is None
    assert result['checks']['chemical_identity'] is False
    assert result['qc_status']=='failed'


def test_pose_bond_and_chiral_geometry():
    mol=ligand('C[C@H](O)C(=O)O')
    xyz=mol.GetConformer().GetPositions();xyz[:,0]*=-1
    for i,p in enumerate(xyz):mol.GetConformer().SetAtomPosition(i,p)
    result=pose_quality.assess(mol)
    assert result['wrong_stereo_atoms']
    mol.GetConformer().SetAtomPosition(0,[100,100,100])
    assert pose_quality.assess(mol)['bad_bonds']


def test_benchmark_labels_and_false_negatives():
    result=screening_evaluation.evaluate(['a','b','c','d'],[('coarse',['a','c']),('refine',['c'])],
                                         dict(a=1,b=1,c=0,d=0),2,.9)
    a,b=result['stages']
    assert a['positive_recall']==.5 and a['enrichment_factor']==1
    assert a['exact_reference_top_k_recovery']==.5 and a['lost_positive_ids']==['b']
    assert b['lost_positive_ids']==['a'] and b['recovery_gate']=='failed'
    partial=screening_evaluation.evaluate(['a','b'],[('x',['a'])],{'a':1},1)
    assert partial['stages'][0]['enrichment_factor'] is None


@pytest.mark.parametrize('stages,labels',[([('x',['a','a'])],{}),([('x',['unknown'])],{}),([('x',['a'])],{'unknown':1}),([('x',['a'])],{'a':True})])
def test_benchmark_rejects_invalid_membership(stages,labels):
    with pytest.raises(ValueError):screening_evaluation.evaluate(['a','b'],stages,labels)


def test_explicit_fragment_and_state_lineage(tmp_path):
    result=molecule_preparation.execute(dict(smiles='CC[NH3+].[Cl-]',fragment_policy='largest',
        state_smiles=['CCN'],max_states=2),tmp_path/'states')
    rows=json.loads(Path(result['details']).read_text())
    assert result['input_records']==1 and len(rows)==2
    assert {r['formal_charge'] for r in rows}=={0,1}
    assert all(r['lineage']['parent_smiles']=='CC[NH3+].[Cl-]' for r in rows)
    assert rows[0]['lineage']['removed_fragments']==['[Cl-]']
    assert rows[0]['lineage']['state_id']!=rows[1]['lineage']['state_id']
    with pytest.raises(ValueError,match='connectivity'):
        molecule_preparation.execute(dict(smiles='CCO',state_smiles=['CCN']),tmp_path/'wrong')


def test_tautomer_and_adaptive_macrocycle_preparation(tmp_path):
    report=molecule_preparation.execute(dict(smiles='CC(=O)C',tautomer_policy='enumerate',max_states=8),tmp_path/'tautomers')
    assert report['prepared_molecules']>=2
    report=molecule_preparation.execute(dict(smiles='C1CCCCCCCCCCC1',conformer_policy='adaptive'),tmp_path/'macrocycle')
    row=json.loads(Path(report['details']).read_text())[0]
    assert row['macrocycle'] and row['requested_conformers']==32 and row['actual_conformers']>0
    assert row['conformer_coverage']=='finite_sample_not_exhaustive'


def batch_args():
    return dict(name='test',protein_sequences=['ACDEFG'],protein_copies=2,seeds=[47,48],
                ligands=[dict(name='first',smiles='CCO'),dict(name='second',smiles='CCN')])


class FakeAF3:
    mode='fixture_not_live'
    def __init__(self,fail=None):self.calls=[];self.fail=fail
    def run_command(self,argv,log):
        self.calls.append(argv)
        input_path=Path(next(s.split('=',1)[1] for s in argv if s.startswith('--json_path=')))
        output=Path(next(s.split('=',1)[1] for s in argv if s.startswith('--output_dir=')))
        data=json.loads(input_path.read_text());output.mkdir(parents=True,exist_ok=True)
        if '--run_inference=false' in argv:
            for item in data['sequences']:
                if 'protein' in item:item['protein'].update(unpairedMsa='',pairedMsa='',templates=[])
            (output/'test_data.json').write_text(json.dumps(data));return
        if data['name']==self.fail:raise RuntimeError('Synthetic child failure')
        (output/(data['name']+'_model.cif')).write_text('data_fixture\n#\n')
        (output/(data['name']+'_summary_confidences.json')).write_text(json.dumps(dict(iptm=.5)))
        (output/'ranking_scores.csv').write_text('seed,sample,ranking_score\n47,0,0.5\n48,0,0.6\n')
        for seed in (47,48):
            sample=output/f'seed-{seed}_sample-0';sample.mkdir()
            (sample/'model.cif').write_text('data_fixture\n#\n')
            (sample/'summary_confidences.json').write_text(json.dumps(dict(iptm=.5,ranking_score=.5)))


def configured(tmp_path,monkeypatch):
    path=tmp_path/'profile.json'
    path.write_text(json.dumps(dict(backend='alphafold3',python='python',runner='fake.py',model_parameters='weights',databases='db',license_acknowledged=True)))
    monkeypatch.setattr(af3_batch,'validate_af3_installation',lambda p:None)
    return dict(af3_profile=str(path))


def test_batch_prepares_before_configuration_block(tmp_path):
    with pytest.raises(Blocked,match='af3_profile'):
        af3_batch.execute(batch_args(),tmp_path/'batch',{},True,FakeAF3())
    state=json.loads((tmp_path/'batch/batch-state.json').read_text())
    assert len(state['jobs'])==2
    data=json.loads(Path(state['jobs'][0]['input']).read_text())
    assert data['modelSeeds']==[47,48]
    assert [r['protein']['sequence'] for r in data['sequences'] if 'protein' in r]==['ACDEFG']*2
    assert data['sequences'][-1]['ligand']['smiles']=='CCO'


def test_batch_shared_features_and_partial_resume(tmp_path,monkeypatch):
    cfg=configured(tmp_path,monkeypatch);service=FakeAF3(fail='ligand-00002');root=tmp_path/'batch'
    with pytest.raises(RuntimeError,match='failed children'):af3_batch.execute(batch_args(),root,cfg,True,service)
    assert len(service.calls)==3
    service.fail=None
    result=af3_batch.execute(batch_args(),root,cfg,True,service)
    assert len(service.calls)==4 and result['status']=='complete'
    assert result['jobs'][0]['attempt']==1 and result['jobs'][1]['attempt']==2
    assert result['jobs'][0]['sample_summary']['best_recorded_sample']['seed']=='48'
    assert result['jobs'][0]['sample_summary']['sample_count']==2
    assert result['jobs'][0]['pose_quality']['qc_status']=='not_assessed'
    af3_batch.execute(batch_args(),root,cfg,True,service);assert len(service.calls)==4
    Path(result['jobs'][0]['model']).write_text('changed')
    with pytest.raises(ValueError):af3_batch.execute(batch_args(),root,cfg,True,service)


def test_batch_uncertain_dispatch_blocks_and_reconciliation(tmp_path,monkeypatch):
    cfg=configured(tmp_path,monkeypatch);args=dict(batch_args(),shared_msa=False);root=tmp_path/'batch'
    service=FakeAF3();af3_batch.execute(args,root,cfg,True,service)
    path=root/'batch-state.json';state=json.loads(path.read_text());state['jobs'][1]['status']='running';path.write_text(json.dumps(state))
    with pytest.raises(Blocked,match='Uncertain'):af3_batch.execute(args,root,cfg,True,service)
    assert len(service.calls)==2
    with pytest.raises(ValueError):af3_batch.reconcile(path,'ligand-00002')
    af3_batch.reconcile(path,'ligand-00002',True)
    af3_batch.execute(args,root,cfg,True,service);assert len(service.calls)==3


def test_shared_protein_identity_mismatch():
    from aidd_agent.af3_shared import transfer
    source=dict(sequences=[dict(protein=dict(id='A',sequence='ACD',unpairedMsa='',pairedMsa='',templates=[]))])
    with pytest.raises(ValueError):transfer(dict(sequences=[dict(protein=dict(id='A',sequence='ACE'))]),source)


def test_missing_requested_seed_is_not_prediction_success():
    summary=dict(ranking_rows=[dict(seed='47',sample='0',ranking_score='.5')],samples=[])
    with pytest.raises(ValueError,match='seed coverage'):
        af3_batch.validate_sample_coverage(summary,dict(modelSeeds=[47,48]))


def test_live_acceptance_rejects_fixture_and_wrong_expectations():
    from aidd_agent.live_acceptance import check_report
    case=dict(expected_actions=['af3_run'],expect={'/steps/fold/result/status':'prediction_completed'})
    report=dict(status='complete',service_mode='offline',steps=dict(fold=dict(status='complete',action='af3_run',result=dict(status='prediction_completed'))))
    with pytest.raises(ValueError,match='Fixture'):check_report(report,case)
    report['service_mode']='live';assert check_report(report,case)['status']=='passed'
    report['steps']['fold']['result']['status']='prepared_not_predicted'
    with pytest.raises(ValueError,match='mismatch'):check_report(report,case)


def test_prompt_batch_queues_and_preserves_project_isolation(tmp_path):
    from aidd_agent.chat_agent import ChatAgent
    from aidd_agent.domain_tools import DomainTools
    from test_domain_agent import scripted,tool
    args=batch_args();app=ChatAgent(tmp_path,start=False,domain_planner=scripted([dict(kind='tool',tool='af3_batch',arguments=args,purpose='Run independent ligand predictions')]))
    try:
        sid=app.new_session();answer=app.ask(sid,'Run ACDEFG homodimer independently with CCO and CCN. seeds: [47,48]','deepseek')
        assert 'queued' in answer and len(app.jobs(sid))==1
        job=app.jobs(sid)[0]
        from aidd_agent.prompt_workflow import run_plan
        ctx=app.context_for(sid)
        report,_=run_plan(Path(ctx['db']),ctx['user_id'],ctx['project_id'],Path(job['plan']),allow_compute=True)
        assert report['status']=='blocked' and 'af3_profile' in report['steps']['work']['error']
        assert list(Path(job['plan']).parent.glob('execution/work/inputs/*.json'))
    finally:app.close()


def test_saved_candidates_benchmark_and_qc_through_real_worker(tmp_path):
    from aidd_agent.chat_agent import ChatAgent
    from aidd_agent.chat_uploads import upload
    from aidd_agent import candidate_sets,ligand_library,scientific_tools
    from test_ligand_library import tiny_library
    from test_research_p0 import execute
    batch,raw=tiny_library(tmp_path)
    runtime=tmp_path/'runtime.json';runtime.write_text(json.dumps(dict(search=dict(batch=str(batch),workers=1,refine_chunk=1))))
    app=ChatAgent(tmp_path/'chat',runtime=runtime,allow_compute=True,start=False)
    try:
        sid=app.new_session();query=upload(app,sid,'query.mol2',io.BytesIO(raw.encode()),len(raw.encode()))
        receipt=ligand_library.handle(app,sid,dict(operation='start',attachment_id=query['id'],export_molecules=1))
        result=execute(app,sid,receipt)
        assert result['timing']['fresh_latency_eligible'] and result['timing']['wall_seconds']>0
        reference=execute(app,sid,candidate_sets.handle(app,sid,dict(operation='import',task_id=receipt['task_id'],name='Reference')))
        stage=execute(app,sid,candidate_sets.handle(app,sid,dict(operation='filter',set_id=reference['set_id'],limit=1,name='Top one')))
        labels=b'molecule_id,label\nM0,1\nM1,0\n';uploaded=upload(app,sid,'labels.csv',io.BytesIO(labels),len(labels))
        report=execute(app,sid,scientific_tools.handle(app,sid,'screening_benchmark',dict(reference_set_id=reference['set_id'],
            stage_set_ids=[stage['set_id']],labels_attachment_id=uploaded['id'],top_k=1,minimum_recovery=1)))
        assert report['stages'][0]['enrichment_factor']==2
        assert report['stages'][0]['recovery_gate']=='passed'
        qc=execute(app,sid,scientific_tools.handle(app,sid,'pose_qc',dict(set_id=stage['set_id'],count=1)))
        assert qc['checked']==1 and qc['poses'][0]['checks']['chemical_identity']
        assert qc['poses'][0]['checks']['assembly_clashes'] is None
        root=Path(app.context_for(sid)['project_root'])
        resolved=scientific_tools.resolve_batch(dict(set_id=stage['set_id'],count=1),tmp_path/'resolved',root)
        assert len(resolved['ligands'])==1 and resolved['ligands'][0]['molecule_id']=='M0'
        from aidd_agent.research_projects import handle
        other=handle(app,sid,dict(operation='create',name='Other',target='OTHER',objective='Isolation'))['session']
        with pytest.raises(ValueError):scientific_tools.handle(app,other,'pose_qc',dict(set_id=stage['set_id'],count=1))
    finally:app.close()


def test_receptor_reader_keeps_partner_chain_metal_and_excludes_native(tmp_path):
    pdb=tmp_path/'assembly.pdb'
    pdb.write_text('ATOM      1  CA  ALA A   1       1.000   0.000   0.000  1.00 20.00           C  \n'
                   'ATOM      2  CA  ALA B   2       2.000   0.000   0.000  1.00 20.00           C  \n'
                   'HETATM    3 ZN    ZN B   3       3.000   0.000   0.000  1.00 20.00          ZN  \n'
                   'HETATM    4  C1  LIG A   4       4.000   0.000   0.000  1.00 20.00           C  \nEND\n')
    atoms=pose_quality.receptor_atoms(pdb,exclude=[('A','4','LIG')])
    assert len(atoms)==3 and {a['chain'] for a in atoms}=={'A','B'}
    assert 'ZN' in {a['element'] for a in atoms}
    with pytest.raises(ValueError,match='missing'):pose_quality.receptor_atoms(pdb,['Z'])


def test_pose_consistency_alignment_and_separate_binding_modes(monkeypatch):
    from aidd_agent import pose_consistency as pc
    a={i:np.array(x,dtype=float) for i,x in enumerate([[0,0,0],[1,0,0],[0,1,0],[0,0,1]])}
    b={0:np.array([3.,3.,3.]),1:np.array([4.,3.,3.])}
    def coords(path,*unused):
        shift=np.array([10,20,30]);extra=np.array([8,0,0]) if path.name=='different' else 0
        return ({k:v+shift for k,v in a.items()},{k:v+shift+extra for k,v in b.items()})
    monkeypatch.setattr(pc,'coordinates',coords)
    data=dict(sequences=[dict(protein=dict(id='A')),dict(ligand=dict(id='B'))])
    result=pc.summarize([dict(sample=s,model=s) for s in ['one','same','different']],data)
    assert result['clusters']==2 and result['samples'][1]['ligand_rmsd_to_first']<1e-6


def test_batch_shared_failure_can_retry_fresh_attempt(tmp_path,monkeypatch):
    cfg=configured(tmp_path,monkeypatch);service=FakeAF3();original=service.run_command
    def fail(argv,log):raise RuntimeError('shared pipeline failed')
    service.run_command=fail;root=tmp_path/'batch'
    with pytest.raises(RuntimeError):af3_batch.execute(batch_args(),root,cfg,True,service)
    service.run_command=original;report=af3_batch.execute(batch_args(),root,cfg,True,service)
    assert report['shared_dispatch']['attempt']==2 and report['status']=='complete'


def test_batch_respects_project_sequence_and_action_constraints(tmp_path):
    from aidd_agent.chat_agent import ChatAgent
    from aidd_agent.scientific_tools import handle
    from aidd_agent.project_execution import validate
    app=ChatAgent(tmp_path,start=False)
    try:
        sid=app.new_session();receipt=handle(app,sid,'af3_batch',batch_args(),'ACDEFG CCO CCN seeds: [47,48]')
        job=app.task(sid,receipt['task_id']);envelope=json.loads(Path(job['plan']).read_text())
        root=Path(app.context_for(sid)['project_root'])
        plan=envelope['plan']
        with pytest.raises(ValueError,match='forbids'):validate(root,plan,dict(rules=dict(forbidden_actions=['af3_run'])),{})
        with pytest.raises(ValueError,match='sequence'):validate(root,plan,dict(rules=dict(protein_sequences=['AAAA'])),{})
    finally:app.close()


def test_live_suite_records_blocked_preflight_without_dispatch(tmp_path,monkeypatch):
    from aidd_agent import live_acceptance as live
    suite=tmp_path/'suite.json';suite.write_text(json.dumps(dict(cases=[dict(id='test',prompt='Prepare CCO',
        expected_actions=['research_task'],result_expect={'research_task':{'kind':'molecule_prepare'}})])))
    monkeypatch.setattr(live,'preflight',lambda *a:dict(provider=dict(status='blocked',error='No credential')))
    result=live.run(suite,tmp_path/'output',tmp_path/'storage',execute=True)
    assert result['status']=='blocked' and not result['cases']
    assert not (tmp_path/'storage').exists()


def test_live_acceptance_checks_actual_input_stoichiometry(tmp_path):
    from aidd_agent.live_acceptance import check_report
    result=af3_batch.prepare(batch_args(),tmp_path)
    report=dict(status='complete',service_mode='live',steps=dict(work=dict(status='complete',action='research_task',result=result)))
    case=dict(expected_actions=['research_task'],af3_inputs=[dict(protein_sequences=['ACDEFG']*2,seeds=[47,48],ligand_count=1)])
    check_report(report,case)
    case['af3_inputs'][0]['protein_sequences']=['ACDEFG']
    with pytest.raises(ValueError,match='stoichiometry'):check_report(report,case)


def test_batch_reuses_exact_owned_protein_input_and_rejects_mutation(tmp_path):
    from aidd_agent.chat_agent import ChatAgent
    from aidd_agent.scientific_tools import handle
    from aidd_agent.library_acceptance import sha
    app=ChatAgent(tmp_path,start=False)
    try:
        sid=app.new_session();root=Path(app.context_for(sid)['project_root'])
        source=af3_batch.prepare(batch_args(),root/'source')
        input_path=Path(source['jobs'][0]['input']);report_path=root/'source-report.json'
        report_path.write_text(json.dumps(dict(steps=dict(prep=dict(status='complete',result=dict(input_json=str(input_path),input_sha256=sha(input_path)))))))
        with app.connect() as db:
            db.execute('INSERT INTO jobs(id,session,status,report,request,created) VALUES(?,?,?,?,?,?)',
                       ('prepared-protein',sid,'complete',str(report_path),'Prepared input fixture only',1))
        args=dict(name='reuse',protein_task_id='prepared-protein',seeds=[47,48],ligands=[dict(smiles='CCO')])
        receipt=handle(app,sid,'af3_batch',args,'Use existing assembly with CCO seeds: [47,48]')
        job=app.task(sid,receipt['task_id']);plan=json.loads(Path(job['plan']).read_text())['plan']
        request=root/'research-requests'/plan['steps'][0]['params']['request_id']/'request.json'
        assert json.loads(request.read_text())['arguments']['protein_sequences']==['ACDEFG']*2
        input_path.write_text('{}')
        with pytest.raises(ValueError,match='changed'):handle(app,sid,'af3_batch',args,'CCO seeds: [47,48]')
    finally:app.close()
