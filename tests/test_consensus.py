import copy
import json
from types import SimpleNamespace

import numpy as np
import pytest

from aidd_agent.consensus_admission import fit_protein, admit
from aidd_agent.consensus_model import consensus, transform_query
from aidd_agent.consensus_design import validate, adaptive_rules
from aidd_agent.guided_filters import adaptive_stage, same_pose_gaussian, pose_rank, GuidedFilter
from test_preselection_full import geometry


def test_protein_fit_recovers_frame_without_ligand_alignment():
    fixed=np.array([[0.,0,0],[1,0,0],[0,2,0],[0,0,3]])
    rotation=np.array([[0,-1,0],[1,0,0],[0,0,1.]])
    moved=fixed@rotation.T+[10,20,30]
    matrix,rmsd,_=fit_protein(moved,fixed)
    assert rmsd<1e-12
    assert np.allclose(moved@matrix[:3,:3].T+matrix[:3,3],fixed)
    assert admit({}, {}, '', {'chains':{}}, {'query_id':'wrong'}, 'A')['status']=='pending'
    with pytest.raises(ValueError):fit_protein([[0,0,0],[1,0,0],[2,0,0]],[[0,0,0],[1,0,0],[2,0,0]])


def test_authoritative_ccd_uppercase_halogen_symbols(tmp_path):
    pytest.importorskip('gemmi');pytest.importorskip('rdkit')
    from aidd_agent.chemistry_prep import _ccd_molecule
    ccd=tmp_path/'LIG.cif'
    ccd.write_text('data_LIG\nloop_\n_chem_comp_atom.atom_id\n_chem_comp_atom.type_symbol\n_chem_comp_atom.charge\nC1 C 0\nCL1 CL 0\nloop_\n_chem_comp_bond.atom_id_1\n_chem_comp_bond.atom_id_2\n_chem_comp_bond.value_order\n_chem_comp_bond.pdbx_aromatic_flag\nC1 CL1 SING N\n')
    molecule=_ccd_molecule(ccd,[dict(atom_id='C1',x=0,y=0,z=0),dict(atom_id='CL1',x=1.8,y=0,z=0)])
    assert [a.GetAtomicNum() for a in molecule.GetAtoms()]==[6,17]


def test_consensus_deduplicates_crystal_copies_and_splits_spatial_modes():
    observations=[]
    for qid,pdb,chem,x in [('a','P1','C1',0),('b','P1','C1',.1),('c','P2','C2',.2),('d','P3','C3',4)]:
        observations.append(dict(key=[20,'N','hydrogen_bond',1,0],query_id=qid,pdb_id=pdb,chemotype=chem,
            source_anchor=qid,point=[x,0,0],direction=[0,0,0],kind=0))
    prepared=[dict(pdb_id=f'P{i}',observed_partner_atoms=['20:N']) for i in range(1,13)]
    anchors,_=consensus(observations,prepared)
    assert len(anchors)==2
    assert sorted(a['distinct_structures'] for a in anchors)==[1,2]
    assert all(a['eligible_structures']==12 and not a['mandatory_proposal_eligible'] for a in anchors)


def survey_design():
    survey=dict(anchors=[dict(anchor_id='A',evidence_id='eA',mandatory_proposal_eligible=False),
                        dict(anchor_id='B',evidence_id='eB',mandatory_proposal_eligible=True)],templates=[dict(query_id='T1'),dict(query_id='T2')])
    design=dict(mandatory_anchors=[],alternative_groups=[],optional_weights={'A':.5,'B':1.},template_ids=['T1'],evidence_ids=[],
        rationale='Exploratory optional contacts.',exclusions=[],permissiveness=1.,gaussian_weight=.7,optional_weight=.3,minimum_pose_score=0.)
    return survey,design


def test_rare_mandatory_blocked_and_templates_independent():
    s,d=survey_design();assert validate(d,s)==d
    other=dict(d,template_ids=['T2']);assert validate(other,s)['optional_weights']==d['optional_weights']
    with pytest.raises(ValueError,match='support'):
        validate(dict(d,mandatory_anchors=['A'],optional_weights={'B':1},evidence_ids=['eA']),s)
    with pytest.raises(ValueError):validate(dict(d,template_ids=['unknown']),s)
    with pytest.raises(ValueError):validate(dict(d,optional_weights={'A':float('nan')}),s)


def test_adaptive_rules_retain_reference_range_and_reject_size_outlier():
    q,_,_=geometry();other=copy.deepcopy(q);other['shape_points']=np.vstack([q['shape_points'],[1,1,1]])
    rule=adaptive_rules([q,other],1)
    assert adaptive_stage(SimpleNamespace(**q),rule)==3
    assert adaptive_stage(SimpleNamespace(**other),rule)==3
    outlier=dict(q,shape_points=np.zeros((100,3)))
    assert adaptive_stage(SimpleNamespace(**outlier),rule)==0


def test_same_pose_gaussian_and_exclusions(tmp_path):
    q,_,seeds=geometry();candidate=SimpleNamespace(**q)
    values=same_pose_gaussian(q,candidate,seeds)
    assert values[0]==pytest.approx(1.) and values[1]<.1
    s,d=survey_design();d.update(receptor_npz=str(tmp_path/'r.npz'),minimum_score=.5,
        pocket=dict(minimum_heavy_atom_distance=1,maximum_clashing_fraction=0),
        exclusions=[dict(mode='hard',center=[0,0,0],radius=.5,weight=1,evidence='user-region',rationale='Reviewed region')])
    np.savez(d['receptor_npz'],points=np.array([[100,100,100.]]))
    gate=GuidedFilter(dict(anchors=[dict(anchor_id='A',feature_index=0),dict(anchor_id='B',feature_index=1)]),q,d)
    assert not gate.pocket_mask(q['shape_points'],np.eye(4)[None])[0]
    d['exclusions'][0]['mode']='soft'
    assert gate.pocket_mask(q['shape_points'],np.eye(4)[None])[0]
    rank=pose_rank(np.array([1.,.2]),np.array([0,1]),['A','B'],.8,q['shape_points'],np.eye(4),d)
    assert rank['exclusion_penalty']==pytest.approx(1/3)
    assert rank['optional_score']==pytest.approx(.7/1.5)
    with pytest.raises(ValueError,match='invent'):validate({k:d[k] for k in survey_design()[1]},s,llm=True)


def test_transformed_shape_does_not_change_consensus(tmp_path):
    q,_,_=geometry();before=copy.deepcopy(q);m=np.eye(4);m[0,3]=5
    out=transform_query(q,m,tmp_path/'aligned.npz',{})
    assert np.allclose(out['shape_points'],q['shape_points']+[5,0,0])
    assert all(np.array_equal(q[k],before[k]) for k in q)


def test_chat_consensus_route_and_action_contract():
    from aidd_agent.chat_agent import validate_route
    from aidd_agent.prompt_plan import validate_plan
    decision=dict(intent='consensus',message='Build consensus.',task_id=None,request='',reference=dict(reference_query='5C5A:NUT:A:201',target_chain='A'))
    assert validate_route(decision)==decision
    plan=dict(version=1,summary='Reference pocket',clarifications=[],steps=[dict(id='c',action='structure_consensus',params=dict(source_run='PROMPT-0123456789abcdef',**decision['reference']))])
    assert validate_plan(plan)
    plan['steps'][0]['params']['reference_query']='../../other'
    with pytest.raises(ValueError):validate_plan(plan)


def test_active_records_require_exact_target_assay_and_units():
    from aidd_agent.active_controls import binding_records
    protein=dict(accession='P1',organism=dict(taxonId=9606))
    def fetch(url):
        if '/target.json?' in url:return dict(targets=[dict(target_type='SINGLE PROTEIN',tax_id=9606,target_components=[dict(accession='P1')],target_chembl_id='T')])
        if '/activity.json?' in url:return dict(activities=[dict(standard_relation='=',standard_units='nM',standard_type='Ki',standard_value=10,
            assay_chembl_id='A',molecule_chembl_id='M',activity_id=1,canonical_smiles='CCO')])
        if '/assay.json?' in url:return dict(assays=[dict(assay_chembl_id='A',confidence_score=9,target_chembl_id='T')])
        return dict(molecule_structures=dict(canonical_smiles='CCO'))
    rows,_=binding_records(protein,fetch);assert len(rows)==1 and rows[0]['type']=='Ki'
    def wrong(url):
        result=fetch(url)
        if '/assay.json?' in url:result['assays'][0]['target_chembl_id']='OTHER'
        return result
    assert binding_records(protein,wrong)[0]==[]


def test_same_pocket_rejects_displaced_ligand_even_with_chain_identity(monkeypatch):
    from aidd_agent import consensus_admission as module
    monkeypatch.setattr(module,'assembly_compatible',lambda *a:True)
    rng=np.random.default_rng(7)
    atoms=[dict(auth_asym_id='A',auth_seq_id=str(i),auth_comp_id='ALA',group_PDB='ATOM',
                auth_atom_id='CA',canonical_residue=i,xyz=p,label_alt_id='',occupancy='1',pdbx_PDB_ins_code='')
           for i,p in enumerate(rng.normal(size=(20,3)))]
    lig=[dict(auth_asym_id='L',auth_seq_id='1',auth_comp_id='LIG',group_PDB='HETATM',xyz=p) for p in [[0.,0,0],[1,0,0],[0,1,0]]]
    ref=dict(atoms=atoms+lig,chains={'A':atoms});row=dict(chain='L',residue='1',ccd_id='LIG',query_id='P:LIG:L:1')
    assert admit(ref,row,'A',ref,row,'A')['status']=='admitted'
    moved=copy.deepcopy(ref)
    for a in moved['atoms'][-3:]:a['xyz']=np.array(a['xyz'])+[30,0,0]
    assert admit(ref,row,'A',moved,row,'A')['status']=='rejected'
    alternate=copy.deepcopy(ref);alternate['chains']['A'][0]['label_alt_id']='B'
    assert admit(ref,row,'A',alternate,row,'A')['status']=='pending'


def test_multi_template_union_preserves_provenance_and_same_pose_selection(tmp_path,monkeypatch):
    import sqlite3
    from aidd_agent import consensus_funnel as module, preselection_full, active_controls
    from aidd_agent.expanded_wee1 import fingerprint
    from aidd_agent.guided_workflow import select_candidates
    batch=tmp_path/'library'
    for name in ['artifacts','chemical']:
        (batch/name).mkdir(parents=True);(batch/name/'catalog.json').write_text('{}')
    source=tmp_path/'design.json';_,design=survey_design()
    design.update(minimum_score=.5,coarse_constraints={})
    source.write_text(json.dumps(dict(readiness='ready_for_consensus_funnel',sources=fingerprint([batch/'artifacts/catalog.json']),design=design,anchor_order=['A','B'],
        anchors=[],consensus_npz='pocket',templates=[dict(query_id='T1'),dict(query_id='T2')],independent_active_validation='not_run',
        reference=dict(coordinate_frame='ref'))))
    monkeypatch.setattr(module.ev,'accept_library',lambda *a,**k:(None,dict(library_conformers=3,library_molecules=2)))
    monkeypatch.setattr(active_controls,'run',lambda *a:dict(status='unavailable',acceptance='not_validated'))
    calls=[]
    def scan(selection,output,**kwargs):
        q=json.loads(selection.read_text())['query'];calls.append(q['query_id']);output.mkdir(parents=True)
        database=output/'candidates.sqlite';mask=1 if q['query_id']=='T1' else 2
        with sqlite3.connect(database) as db:
            preselection_full.schema(db)
            preselection_full.merge_pose(db,dict(molecule_id='M',mask=mask,min_matched_score=.8,global_id=1,
                scores=[.8,0] if mask==1 else [0,.8],assignments=[0,1],transform=np.eye(4).reshape(-1).tolist()))
            db.execute('INSERT INTO members VALUES (?,?,?,?,?,?)',('M','scaffold','C','',mask,1))
        report=dict(full_coverage=True,counts=dict(input_conformers=3),outputs={'candidates.sqlite':str(database)},output_hashes=fingerprint([database]))
        (output/'report.json').write_text(json.dumps(report));return report
    monkeypatch.setattr(preselection_full,'run',scan)
    result=module.run(source,tmp_path/'output',dict(batch=str(batch),workers=1,refine_chunk=2))
    assert calls==['T1','T2'] and result['matching_molecules']==1 and result['pose_records']==2
    with sqlite3.connect(result['outputs']['candidates.sqlite']) as db:
        assert {r[0] for r in db.execute('SELECT template FROM poses')}=={'T1','T2'}
    selected=select_candidates(tmp_path/'output/report.json',tmp_path/'selected',dict(required_anchors=['A','B'],minimum_score=.5,match_mode='all'))
    assert selected['matching_molecules']==0


def test_consensus_compute_scores_each_retained_pose(tmp_path,monkeypatch):
    from aidd_agent import preselection_full as pre
    from aidd_agent.necessary_conditions import NecessaryConditions
    query,features,seeds=geometry();candidate=SimpleNamespace(**vars(features),shape_points=query['shape_points'])
    _,design=survey_design();design.update(minimum_score=.5,adaptive_coarse=adaptive_rules([query],1))
    policy=dict(required_anchors=['A','B'],minimum_score=.5,coarse_constraints={})
    q=dict(consensus_npz='separate-pocket',condition_policy=policy,guided_design=design,
           anchors=[dict(anchor_id='A',score_column=0),dict(anchor_id='B',score_column=1)])
    monkeypatch.setattr(pre.full,'_STATE',(SimpleNamespace(get=lambda _:candidate),query,query,q))
    monkeypatch.setattr(pre.full,'_FILTER_READER',SimpleNamespace(get=lambda _:features))
    monkeypatch.setattr(pre.full,'_BOUND',NecessaryConditions(query,[0,1],'any',.5))
    monkeypatch.setattr(pre,'prepare_seeds',lambda *a,**k:(seeds,2))
    gate=SimpleNamespace(design=design,geometry_possible=lambda _:True,geometry_seeds=lambda f,s,p:p,pocket_mask=lambda p,t:np.ones(len(t),bool))
    monkeypatch.setattr(pre,'_GUIDED_CACHE',(json.dumps(design,sort_keys=True),gate))
    def forbidden(*a,**k):raise AssertionError('Conformer winner cannot replace same-pose scoring')
    monkeypatch.setattr(pre,'_score_ids',forbidden)
    receipt=pre.compute((0,1,str(tmp_path/'chunk.npz')))
    rows=[json.loads(line) for line in (tmp_path/'chunk.poses.jsonl').read_text().splitlines()]
    assert receipt['counts']['gaussian_evaluated_seeds']==2 and len(rows)==2
    assert rows[0]['gaussian_same_pose']!=rows[1]['gaussian_same_pose']
    assert all('composite_score' in r for r in rows)
