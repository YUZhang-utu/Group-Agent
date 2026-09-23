import copy

import numpy as np
import pytest

from aidd_agent.contact_groups import contact_semantics, grouped_terms, resolve_regions
from aidd_agent.consensus_design import validate
from aidd_agent.guided_filters import pose_rank
from test_consensus import survey_design
from test_contact_groups import spatial_fixture


def rank(design,values=(1.,0.)):
    return pose_rank(np.array(values),np.array([0,1]),['A','B'],.5,np.zeros((1,3)),np.eye(4),design)


def test_fixed_budget_deleting_unmatched_contact_does_not_inflate_other_score():
    survey,design=survey_design()
    design.update(optional_normalization='fixed_budget',optional_budget=1.5)
    validate(design,survey)
    before=rank(design)
    del design['optional_weights']['B']
    validate(design,survey)
    assert rank(design)['composite_score']==before['composite_score']
    design['optional_weights']['A']=1.
    assert rank(design)['optional_score']==pytest.approx(2*before['optional_score'])
    assert rank(design)['optional_denominator']==1.5


def test_family_weight_uses_frozen_budget_not_member_count():
    survey,design=survey_design()
    design.update(optional_weights={},optional_groups=[dict(id='family',anchor_ids=['A','B'],weight=.6)],
        optional_normalization='fixed_budget',optional_budget=2.)
    validate(design,survey)
    assert rank(design,(.8,.9))['optional_score']==pytest.approx(.6*.9/2)
    design['optional_groups'][0]['anchor_ids']=['B']
    assert rank(design,(.8,.9))['optional_score']==pytest.approx(.6*.9/2)


@pytest.mark.parametrize('update',[
    dict(optional_budget=0),dict(optional_budget=float('nan')),dict(optional_budget=1.),
    dict(optional_normalization='unknown'),dict(gaussian_weight=.8),dict(optional_budget=True)])
def test_invalid_budget_or_dimension_contract_rejected(update):
    survey,design=survey_design()
    design.update(optional_normalization='fixed_budget',optional_budget=1.5)
    design.update(update)
    with pytest.raises(ValueError):validate(design,survey)


def test_backbone_and_sidechain_are_distinct_without_changing_ids():
    assert contact_semantics('O','HBD')==dict(protein_part='backbone',ligand_role='donor',interaction_type='HBD')
    assert contact_semantics('CG1','hydrophobic')['protein_part']=='sidechain'
    assert contact_semantics('CG/ND1/CE1/NE2/CD2','pi_stacking')['protein_part']=='sidechain'
    assert contact_semantics('N','HBA')['ligand_role']=='acceptor'


def test_region_specific_margin_considers_all_competitors_and_order(tmp_path):
    survey,design,_,_=spatial_fixture(tmp_path)
    for group in design['spatial_groups']:group['ambiguity_margin']=.1
    validate(design,survey)
    regions=resolve_regions(design,survey)
    for group,x in zip(regions,[0.,.3,.6]):group['points']=[[x,0,0]]
    design['resolved_spatial_groups']=regions
    def result(d):return grouped_terms([1,1],[0,1],['A','B'],np.array([[0.,0,0]]),d)
    assert result(design)['occupied_spatial_groups']==['region0']
    regions[2]['ambiguity_margin']=.7
    assert result(design)['ambiguous_spatial_atoms']==1
    other=copy.deepcopy(design);other['resolved_spatial_groups'].reverse()
    assert result(other)['ambiguous_spatial_atoms']==1
    assert result(other)['occupied_group_count']==0


def test_plan_accepts_fixed_budget_contract():
    from aidd_agent.chat_agent import validate_route
    route=dict(intent='design',message='Freeze score budget',task_id=None,request='',
        design=dict(optional_normalization='fixed_budget',optional_budget=9.75))
    assert validate_route(route)==route
