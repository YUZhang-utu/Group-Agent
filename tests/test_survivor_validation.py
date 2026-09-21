import numpy as np
import pytest

from aidd_agent.survivor_validation import anchor_diagnostics, self_control, validate_survivors


def test_crystal_feature_self_control_matches_identity_and_gaussian():
    from test_pose_feasibility import fixture
    q,_,_=fixture()
    policy=dict(required_anchors=['a','b'],match_mode='all',minimum_score=.5,
        coarse_constraints=dict(heavy_atom_ratio=[.7,1.3],maximum_extent_distance=.35,minimum_feature_coverage=.7))
    definition=dict(condition_policy=policy,anchors=[dict(anchor_id='a',feature_index=0,score_column=0),
                                                    dict(anchor_id='b',feature_index=1,score_column=1)])
    result=self_control(q,{**q,'anchor_feature_indices':np.array([0,1])},definition)
    for k in ('joint_eligible','anchor_bound_passed','seed_feasibility_passed','identity_passed','gaussian_selected_pose_passed'):
        assert result[k]
    assert result['identity_scores']==[1.,1.]


def test_diagnostics_do_not_merge_different_poses():
    rows=anchor_diagnostics(np.array([[.9,.1],[.1,.9]]),np.array([[0,1],[0,1]]),
                            [dict(anchor_id='a'),dict(anchor_id='b')],.5)
    assert [r['individual_passes'] for r in rows]==[1,1]
    assert rows[-1]['cumulative_all_passes']==0


@pytest.mark.parametrize('lose_positive',[False,True])
def test_survivor_validation_detects_lost_library_positive(tmp_path,monkeypatch,lose_positive):
    from aidd_agent import full_library_screen as full, funnel_benchmark as bench, survivor_validation as mod
    definition=dict(condition_policy=dict(required_anchors=['a'],minimum_score=.5),anchors=[dict(anchor_id='a',score_column=0)])
    monkeypatch.setattr(full,'_STATE',(None,None,None,definition))
    monkeypatch.setattr(mod,'self_control',lambda *args:dict.fromkeys(
        ['joint_eligible','anchor_bound_passed','seed_feasibility_passed','identity_passed','gaussian_selected_pose_passed'],True))
    def work(ids):
        optimized=definition['pose_feasibility']
        selected=ids[ids==0] if optimized else ids
        if optimized and lose_positive:selected=selected[:0]
        return dict(global_ids=selected,checked_global_ids=ids,condition_passed=selected==0,
            **{full.OBJECTIVE+'__anchor_scores':(selected==0).astype(float)[:,None],
               full.OBJECTIVE+'__anchor_assignments':np.zeros((len(selected),1),dtype=int)}),{}
    monkeypatch.setattr(bench,'work',work)
    result=validate_survivors(definition,[0],[1,2],tmp_path)
    assert result['reference_positive_ids']==[0]
    assert result['status']==('failed' if lose_positive else 'passed')
    assert result['checks'][0]['false_rejections']==int(lose_positive)
    assert result['anchor_diagnostics'][0]['individual_passes']==1
