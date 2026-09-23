import copy
import json

import numpy as np
import pytest

from aidd_agent.contact_evidence import pair_ledger, quality
from aidd_agent.contact_groups import resolve_regions, grouped_terms, selected_anchors
from aidd_agent.consensus_design import validate
from aidd_agent.expanded_wee1 import fingerprint
from aidd_agent.guided_filters import pose_rank, rule_passes
from test_consensus import survey_design


def atom(name, point, element='C', residue=301, chain='Z', **kw):
    return dict(atom=name, point=point, element=element, quality_issues=[], chain=chain,
                author_residue=str(residue), canonical_residue=residue, **kw)


@pytest.mark.parametrize('residue,chain',[(301,'Z'),(905,'B')])
def test_atom_ledger_preserves_aromatic_halogen_and_all_polar_partners(residue,chain):
    ligand=[atom('CL1',[0,0,0],'CL'),atom('N1',[0,1,0],'N',chemistry={'donor':True})]
    receptor=[atom('O',[2,0,0],'O',residue,chain,acceptor=True),
              atom('OD1',[2,1,0],'O',residue+1,chain,acceptor=True),
              atom('CG',[0,0,3],'C',residue+2,chain)]
    rows=pair_ledger(ligand,receptor)
    assert len(rows)==6
    assert len([r for r in rows if 'donor_acceptor_proximity' in r['classes']])==2
    assert len([r for r in rows if 'halogen_proximity_not_halogen_bond' in r['classes']])==3
    assert {r['protein_chain'] for r in rows}=={chain}
    assert all(r['status']=='geometric_hypothesis' and r['angle_status']=='not_evaluated' for r in rows)


def test_uncertain_contact_is_not_promoted_and_input_is_preserved():
    ligand=[atom('C1',[0,0,0])];receptor=[atom('CG',[1,0,0])]
    receptor[0]['quality_issues']=['alternate_location']
    before=copy.deepcopy(receptor)
    assert pair_ledger(ligand,receptor)==[] and receptor==before


def test_optional_family_counts_max_once_and_is_not_a_hard_or():
    survey,design=survey_design()
    design.update(optional_weights={},optional_groups=[dict(id='ring',anchor_ids=['A','B'],weight=.6)])
    validate(design,survey)
    result=pose_rank(np.array([.8,.7]),np.array([0,1]),['A','B'],.5,np.zeros((1,3)),np.eye(4),design)
    assert result['optional_score']==pytest.approx(.8)
    assert result['composite_score']==pytest.approx(.7*.5+.3*.8)
    assert rule_passes(0,['A','B'],design)
    assert selected_anchors(design)==['A','B']
    missing=pose_rank(np.array([.8,.7]),np.array([-1,-1]),['A','B'],.5,np.zeros((1,3)),np.eye(4),design)
    assert missing['optional_score']==0


@pytest.mark.parametrize('bad', ['duplicate','overlap','unknown','nan','required'])
def test_invalid_optional_family_rejected(bad):
    survey,design=survey_design()
    design.update(optional_weights={},optional_groups=[dict(id='ring',anchor_ids=['A','B'],weight=.6)])
    if bad=='duplicate':design['optional_groups'][0]['anchor_ids']=['A','A']
    if bad=='overlap':design['optional_weights']={'A':.4}
    if bad=='unknown':design['optional_groups'][0]['anchor_ids']=['missing']
    if bad=='nan':design['optional_groups'][0]['weight']=float('nan')
    if bad=='required':design['mandatory_anchors']=['B']
    with pytest.raises(ValueError):validate(design,survey)


def spatial_fixture(tmp_path):
    survey,design=survey_design()
    source=tmp_path/'source';source.write_text('verified structure fixture')
    ledger=dict(target='P_OTHER',coordinate_frame='F',sources=fingerprint([source]),complexes=[
        dict(query_id='TEST:LIG:Z:501',reference_state_check=dict(status='no_severe_overlap'),
             ligand_atoms=[atom('C1',[0,0,0]),atom('C2',[6,0,0]),atom('C3',[0,6,0])])])
    path=tmp_path/'contacts.json';path.write_text(json.dumps(ledger))
    survey.update(target=dict(accession='P_OTHER'),cohort=dict(reference=dict(coordinate_frame='F')),
                  contact_evidence=dict(path=str(path)),sources=fingerprint([path]))
    groups=[dict(id='region'+str(i),reference_query='TEST:LIG:Z:501',ligand_atoms=[name],radius=2.,minimum_atoms=1)
            for i,name in enumerate(['C1','C2','C3'])]
    design.update(spatial_groups=groups,occupancy_rewards=[0,.1,.3,1.],occupancy_weight=.5,spatial_ambiguity=.5)
    return survey,design,path,ledger


def test_spatial_regions_use_sealed_coordinates_and_nonlinear_same_pose_score(tmp_path):
    survey,design,_,_=spatial_fixture(tmp_path)
    validate(design,survey)
    design['resolved_spatial_groups']=resolve_regions(design,survey)
    points=np.array([[0.,0,0],[6,0,0],[0,6,0]])
    def score(p):return pose_rank(np.ones(2),np.arange(2),['A','B'],.5,p,np.eye(4),design)
    one,two,three=[score(points[:n]) for n in [1,2,3]]
    assert [x['occupancy_score'] for x in [one,two,three]]==[.1,.3,1.]
    assert three['composite_score']-two['composite_score'] > two['composite_score']-one['composite_score']
    assert score(points[[1]])['occupied_group_count']==1  # A separate pose cannot accumulate the previous pose's count.
    assert design['occupancy_rewards']==[0,.1,.3,1.]


def test_overlap_does_not_double_count_an_atom(tmp_path):
    survey,design,_,_=spatial_fixture(tmp_path)
    regions=resolve_regions(design,survey)
    regions[1]['points']=[[2,0,0]]
    design['resolved_spatial_groups']=regions
    result=grouped_terms([1,1],[0,1],['A','B'],np.array([[1.,0,0]]),design)
    assert result['occupied_group_count']==0 and result['ambiguous_spatial_atoms']==1


@pytest.mark.parametrize('bad',['source_tamper','ledger_tamper','target','frame','unknown_atom','uncertain_atom','duplicate_atom','rewards','invented_points','unknown_query','weight','receptor_state'])
def test_spatial_validation_rejects_untrusted_or_invalid_regions(tmp_path,bad):
    survey,design,path,ledger=spatial_fixture(tmp_path)
    if bad=='source_tamper':(tmp_path/'source').write_text('changed')
    if bad=='ledger_tamper':path.write_text('{}')
    if bad=='target':survey['target']['accession']='WRONG'
    if bad=='frame':survey['cohort']['reference']['coordinate_frame']='OTHER'
    if bad=='unknown_atom':design['spatial_groups'][0]['ligand_atoms']=['C99']
    if bad=='uncertain_atom':
        ledger['complexes'][0]['ligand_atoms'][0]['quality_issues']=['alternate_location']
        path.write_text(json.dumps(ledger));survey['sources']=fingerprint([path])
    if bad=='duplicate_atom':design['spatial_groups'][1]['ligand_atoms']=['C1']
    if bad=='rewards':design['occupancy_rewards']=[0,.8,.2,1.]
    if bad=='invented_points':design['spatial_groups'][0]['points']=[[9,9,9]]
    if bad=='unknown_query':design['spatial_groups'][0]['reference_query']='unknown'
    if bad=='weight':design['occupancy_weight']=float('nan')
    if bad=='receptor_state':
        ledger['complexes'][0]['reference_state_check']['status']='incompatible'
        path.write_text(json.dumps(ledger));survey['sources']=fingerprint([path])
    with pytest.raises(ValueError):validate(design,survey)


def test_legacy_scores_unchanged():
    _,design=survey_design()
    result=pose_rank(np.array([.4,.9]),np.array([0,1]),['A','B'],.6,np.zeros((1,3)),np.eye(4),design)
    expected_optional=(.5*.4+.9)/1.5
    assert result['composite_score']==pytest.approx(.7*.6+.3*expected_optional)
    assert 'occupancy_score' not in result


def test_reference_state_clash_is_not_reclassified_as_ligand_inactivity():
    from aidd_agent.consensus_admission import reference_state_check
    receptor=[dict(xyz=[0.,0.,0.],canonical_residue=301,auth_atom_id='NZ',label_alt_id='',occupancy='1')]
    failed=reference_state_check(np.array([[.6,0,0]]),receptor)
    assert failed['status']=='incompatible' and failed['clashes'][0]['reference_residue']==301
    assert reference_state_check(np.array([[3.,0,0]]),receptor)['status']=='no_severe_overlap'
    receptor[0]['label_alt_id']='A'
    assert reference_state_check(np.array([[.6,0,0]]),receptor)['status']=='unknown'


def test_spatial_score_common_rotation_and_translation_preserves_group_counts(tmp_path):
    survey,design,_,_=spatial_fixture(tmp_path)
    design['resolved_spatial_groups']=resolve_regions(design,survey)
    points=np.array([[0.,0,0],[6,0,0],[0,6,0]])
    before=grouped_terms([1.,1.],[0,1],['A','B'],points,design)
    r=np.array([[0,-1,0],[1,0,0],[0,0,1]])
    t=np.array([10.,-3.,4.]);other=copy.deepcopy(design)
    for group in other['resolved_spatial_groups']:group['points']=(np.array(group['points'])@r.T+t).tolist()
    after=grouped_terms([1.,1.],[0,1],['A','B'],points@r.T+t,other)
    assert before==after


def test_plan_and_chat_accept_new_fields_without_accepting_coordinates():
    from aidd_agent.chat_agent import validate_route
    from aidd_agent.prompt_plan import validate_plan
    edit=dict(optional_groups=[dict(id='ring',anchor_ids=['A','B'],weight=.6)])
    route=dict(intent='design',message='Review edited design',task_id=None,request='',design=edit)
    assert validate_route(route)==route
    plan=dict(version=1,summary='Review groups',clarifications=[],steps=[dict(id='d',action='consensus_design',
        params=dict(source_run='PROMPT-0123456789abcdef',design=edit))])
    assert validate_plan(plan)==plan


def test_sql_retains_distinct_spatial_occupancy_without_unioning_poses():
    import sqlite3
    from aidd_agent.preselection_full import schema,merge_pose
    with sqlite3.connect(':memory:') as db:
        schema(db)
        row=dict(molecule_id='m',mask=3,min_matched_score=.8,global_id=1,occupied_spatial_groups=['siteA'])
        merge_pose(db,row)
        merge_pose(db,dict(row,global_id=2,occupied_spatial_groups=['siteB'],min_matched_score=.7))
        merge_pose(db,dict(row,global_id=3,min_matched_score=.9))
        records=[json.loads(r[0]) for r in db.execute('SELECT payload FROM poses')]
        assert len(records)==2
        assert {r['global_id'] for r in records}=={2,3}
        assert all(len(r['occupied_spatial_groups'])==1 for r in records)


def test_seed_representatives_preserve_same_anchor_mask_with_different_regions(tmp_path,monkeypatch):
    from aidd_agent import preselection_full as pre
    from test_preselection_full import geometry
    survey,design,_,_=spatial_fixture(tmp_path)
    design['resolved_spatial_groups']=resolve_regions(design,survey)
    design['minimum_pose_score']=0.
    query,features,seeds=geometry()
    seeds[1].transform_matrix=np.eye(4).reshape(-1)
    seeds[1].transform_matrix[3]=6.
    def match(query,*args):
        return None,np.array([[0,1],[0,1]]),np.array([[.8,.8],[.8,.8]])
    monkeypatch.setattr(pre,'match_batch',match)
    rows,count=pre.pose_representatives(features,seeds,np.ones(2,bool),query,[0,1],.5,
        design,['A','B'],np.array([.6,.6]),np.array([[0.,0.,0.]]))
    assert count==2 and len(rows)==2 and {r['mask'] for r in rows}=={3}
    assert {tuple(r['occupied_spatial_groups']) for r in rows}=={('region0',),('region1',)}
