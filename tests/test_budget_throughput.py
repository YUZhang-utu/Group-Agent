from types import SimpleNamespace

import numpy as np
import pytest

from aidd_agent import preselection_full as pre
from aidd_agent.guided_filters import same_pose_gaussian
from test_preselection_full import geometry


@pytest.mark.parametrize('mode',['unique','ties','zero','masked'])
def test_lazy_gaussian_preserves_exact_representative(mode):
    query,features,seeds=geometry()
    candidate=SimpleNamespace(**vars(features),shape_points=query['shape_points'])
    if mode=='zero':features.feature_points+=100
    design=dict(optional_weights={'A':1.,'B':1. if mode in ('ties','zero') else .3},
        optional_normalization='fixed_budget',optional_budget=3.,gaussian_weight=.7,
        optional_weight=.3,minimum_pose_score=.2,optional_groups=[dict(id='g',anchor_ids=['A'],weight=.1)])
    if mode=='ties':design['optional_groups']=[]
    seeds.append(seeds[0])
    possible=np.ones(len(seeds),bool)
    if mode=='masked':possible[0]=False
    eager=pre.pose_representatives(features,seeds,possible,query,[0,1],.5,design,['A','B'],
        same_pose_gaussian(query,candidate,seeds),candidate.shape_points,budget_mode=True)
    lazy=pre.contact_first_representatives(features,seeds,possible,query,[0,1],.5,design,['A','B'],candidate,query)
    assert lazy[0]==eager[0] and lazy[1]==eager[1]
    assert 1<=lazy[2]<=possible.sum()
    if mode=='unique':assert lazy[2]<possible.sum()


def test_compute_uses_lazy_gaussian_and_labels_missing_all_pose_max(tmp_path,monkeypatch):
    query,features,seeds=geometry();candidate=SimpleNamespace(**vars(features),shape_points=query['shape_points'])
    design=dict(optional_weights={'A':1.,'B':.3},gaussian_weight=.7,optional_weight=.3,minimum_pose_score=.2)
    q=dict(selection_mode='budget',consensus_npz='fixture',guided_design=design,
        condition_policy=dict(required_anchors=['A','B'],minimum_score=.5),
        anchors=[dict(anchor_id='A',score_column=0),dict(anchor_id='B',score_column=1)])
    from aidd_agent import guided_filters
    monkeypatch.setattr(guided_filters,'GuidedFilter',lambda *a:SimpleNamespace(pocket_mask=lambda p,m:np.ones(len(m),bool)))
    monkeypatch.setattr(pre,'_GUIDED_CACHE',None)
    monkeypatch.setattr(pre.full,'_STATE',(SimpleNamespace(get=lambda _:candidate),query,query,q))
    monkeypatch.setattr(pre.full,'_FILTER_READER',SimpleNamespace(get=lambda _:features))
    monkeypatch.setattr(pre,'prepare_seeds',lambda *a,**kw:(seeds,2))
    path=tmp_path/'chunk.npz';pre.compute((0,1,str(path),[7]))
    import json
    receipt=json.loads(path.with_suffix('.receipt.json').read_text())
    pose=json.loads(path.with_suffix('.poses.jsonl').read_text())
    assert receipt['counts']['gaussian_evaluated_seeds']==1
    assert receipt['counts']['exact_seed_annotations']==2
    assert pose['conformer_gaussian_best'] is None and 'Not evaluated' in pose['conformer_gaussian_scope']
