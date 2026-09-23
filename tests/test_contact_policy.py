import copy

import pytest

from aidd_agent.contact_policy import without_mediators, require_protein_contacts
from aidd_agent.consensus_design import validate
from test_consensus import survey_design


def test_water_removed_from_groups_and_budget_preserved():
    survey,design=survey_design()
    survey['anchors'][0]['feature_class']='water_bridge'
    survey['anchors'][1]['feature_class']='hydrophobic'
    original=copy.deepcopy(design)
    updated,removed=without_mediators(design,survey['anchors'])
    assert removed==['A'] and design==original
    assert updated['optional_weights']=={'B':1.}
    assert updated['optional_budget']==1.5
    validate(updated,survey)
    with pytest.raises(ValueError,match='evidence-only'):validate(design,survey,llm=True)
    validate(design,survey)  # Historical reports remain auditable.
    with pytest.raises(ValueError,match='evidence-only'):require_protein_contacts(design,survey['anchors'])


def test_partial_family_retains_weight_and_no_water_can_satisfy_rule():
    survey,design=survey_design()
    survey['anchors'][0]['feature_class']='water_bridge'
    design.update(optional_weights={},optional_groups=[dict(id='mixed',anchor_ids=['A','B'],weight=.6)],
        optional_normalization='fixed_budget',optional_budget=2.)
    updated,_=without_mediators(design,survey['anchors'])
    assert updated['optional_groups']==[dict(id='mixed',anchor_ids=['B'],weight=.6)]
    assert updated['optional_budget']==2.


def test_water_only_design_cannot_be_migrated_into_empty_screen():
    survey,design=survey_design()
    for a in survey['anchors']:a['feature_class']='water_bridge'
    updated,_=without_mediators(design,survey['anchors'])
    with pytest.raises(ValueError,match='1 to 20'):validate(updated,survey)
