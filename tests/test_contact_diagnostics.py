import itertools

import numpy as np
import pytest

from aidd_agent.contact_diagnostics import (duplicate_relations, optional_group_assignment,
    match_details, distributions, contributions, source_index)


def test_identity_and_controlled_angular_response():
    points=np.array([[0.,0,0],[4,0,0],[8,0,0]])
    directions=np.array([[0.,0,0],[1,0,0],[1,0,0]])
    query=[points,np.array([5,1,6]),directions,np.array([0,1,2])]
    candidate=dict(zip(('feature_points','feature_types','feature_directions','feature_direction_kinds'),query))
    labels=[dict(layer='identity') for _ in points]
    rows,_,_=match_details(query,candidate,labels)
    assert all(r['score']==pytest.approx(1) for r in rows)
    candidate=dict(candidate,feature_directions=np.array([[0.,0,0],[-1,0,0],[-1,0,0]]))
    rows,_,_=match_details(query,candidate,labels)
    assert [r['score'] for r in rows]==[1.,0.,1.]
    summary=distributions(rows)
    assert sum(r['unmatched'] for r in summary)==1
    assert next(r for r in summary if r['direction_kind']==1)['score']['quantiles']==[0.]*5


def evidence(name,x,query='Q',indices=None):
    return dict(query_id=query,source_anchor=name,anchor_id=name,point=[x,0,0],
        ligand_atom_indices=[1] if indices is None else indices,key=[54,'O','HBD'])


def test_relations_are_source_scoped_pairwise_not_transitive():
    rows=[evidence('A',0),evidence('B',.4),evidence('C',.8),evidence('D',0,query='OTHER'),evidence('E',0,indices=[])]
    result=duplicate_relations(rows,.5)
    assert {(r['left'],r['right']) for r in result['pairs']}=={('A','B'),('B','C')}
    assert len(result['unresolved'])==1
    rows[1]['key']=[94,'NZ','salt_bridge']
    assert not duplicate_relations(rows,.5)['pairs'][0]['same_interaction_type']


def test_source_mapping_requires_instance_namespace():
    manifest=dict(source=dict(anchor_mapping={'A0:HBD:3':2}))
    assert source_index(dict(query_id='Q',source_anchor='Q/A0:HBD:3'),manifest)==2
    assert source_index(dict(query_id='Q',source_anchor='Q/F7:salt_bridge'),manifest)==7
    assert source_index(dict(query_id='Q',source_anchor='OTHER/F7:salt_bridge'),manifest) is None


def test_group_assignment_avoids_wasted_candidate_and_matches_exhaustive_optimum():
    matrix=np.array([[1.,0.],[0.,.9],[0.,.8]])
    design=dict(mandatory_anchors=[],alternative_groups=[],optional_weights={'C':1.},
        optional_groups=[dict(id='family',anchor_ids=['A','B'],weight=1.)])
    r=optional_group_assignment(matrix,['A','B','C'],design)
    assert r['optional_numerator']==pytest.approx(1.8)
    assert len({s['candidate_feature'] for s in r['selections']})==2
    # Exhaust all group/member/candidate selections, including unmatched choices.
    possibilities=[]
    for a,b in itertools.product([-1,0,1],repeat=2):
        if a>=0 and a==b:continue
        possibilities.append((max(matrix[0,a],matrix[1,a]) if a>=0 else 0)+(matrix[2,b] if b>=0 else 0))
    assert r['optional_numerator']==max(possibilities)
    design['mandatory_anchors']=['C']
    assert optional_group_assignment(matrix,['A','B','C'],design)['status']=='not_applicable_hard_rules'


@pytest.mark.parametrize('seed',range(8))
def test_group_assignment_random_small_exhaustive(seed):
    matrix=np.random.default_rng(seed).random((3,3))
    design=dict(mandatory_anchors=[],alternative_groups=[],optional_weights={'C':.3},
        optional_groups=[dict(id='family',anchor_ids=['A','B'],weight=.7)])
    r=optional_group_assignment(matrix,['A','B','C'],design)
    brute=max((.7*max(matrix[0,a],matrix[1,a]) if a>=0 else 0)+(.3*matrix[2,b] if b>=0 else 0)
        for a,b in itertools.product([-1,0,1,2],repeat=2) if a<0 or a!=b)
    assert r['optional_numerator']==pytest.approx(brute)


def test_type_contributions_split_group_ties_and_include_unmatched():
    design=dict(optional_weights={'C':.5},optional_groups=[dict(id='g',anchor_ids=['A','B'],weight=1.)],
        optional_normalization='fixed_budget',optional_budget=2.,gaussian_weight=.7,optional_weight=.3)
    anchors={a:dict(feature_class=k) for a,k in zip('ABC',['hydrophobic','HBD','water_bridge'])}
    rows=contributions([.8,.8,.9],[0,1,-1],list('ABC'),design,anchors)
    assert sum(r['contact_contribution'] for r in rows)==pytest.approx(.4)
    assert sum(r['composite_contribution'] for r in rows)==pytest.approx(.12)
    assert next(r for r in rows if r['anchor_id']=='C')['numerator']==0
