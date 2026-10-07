"""Exact matching, inherited funnel semantics and sealed shared-block reuse."""
import copy
import json
from types import SimpleNamespace

import numpy as np
import pytest

from aidd_agent import budget_screen, preselection_full as pre
from aidd_agent.block_refinement_cache import refine_cached
from aidd_agent.expanded_wee1 import fingerprint
from aidd_agent.final_work_blocks import read
from aidd_agent.joint_spatial_profiles import save
from aidd_agent.interaction_matching import _assignment_reference, _assignment_vectorized
from aidd_agent.necessary_conditions import NecessaryConditions
from test_preselection_full import geometry


def test_numpy_assignment_preserves_exact_indices_including_ties():
    rng = np.random.default_rng(106)
    for rows in (0, 1, 4, 8, 20):
        for cols in (0, 3, 12, 40):
            for _ in range(20):
                weights = rng.integers(0, 4, rows).astype(float)
                for scores in (rng.integers(0, 4, (rows, cols))/3,
                               rng.random((rows, cols)), np.zeros((rows, cols))):
                    np.testing.assert_array_equal(_assignment_vectorized(scores, weights),
                                                  _assignment_reference(scores, weights))


def test_queries_inherit_thresholds_without_changing_legacy_budget(tmp_path):
    design = dict(templates=[dict(query_id='fixture')], anchors=[dict(feature_index=0)],
                  consensus_npz='query.npz', anchor_order=['A'],
                  design=dict(minimum_score=.73, minimum_pose_score=.61,
                              coarse_constraints={'fixture': True}))
    legacy = budget_screen.make_queries(tmp_path, design)[0]
    assert legacy['selection_mode'] == 'budget'
    assert legacy['condition_policy']['minimum_score'] == .5
    q = budget_screen.make_queries(tmp_path, dict(design, block_search_mode='threshold'))[0]
    assert q['selection_mode'] == 'threshold'
    assert q['condition_policy']['minimum_score'] == .73
    assert q['guided_design'] == design['design']
    assert q['condition_policy']['coarse_constraints'] == {'fixture': True}


@pytest.mark.parametrize('minimum_pose_score,large,expected',[(.1,False,1),(.99,False,0),(.1,True,0)])
def test_block_threshold_equals_existing_full_funnel(tmp_path, monkeypatch, minimum_pose_score, large, expected):
    from aidd_agent import guided_filters
    query, features, seeds = geometry()
    candidate = SimpleNamespace(**vars(features), shape_points=query['shape_points'])
    if large:
        candidate.shape_points = np.tile(candidate.shape_points, (5, 1))
    receptor = tmp_path/'receptor.npz'
    np.savez(receptor, points=np.array([[100.,100.,100.]]))
    rule = dict(heavy_atom_ratio=[.7,1.3], minimum_feature_coverage=.5, maximum_extent_distance=.35)
    design = dict(mandatory_anchors=['A'], alternative_groups=[], optional_weights={'A':1.},
                  minimum_score=.73, minimum_pose_score=minimum_pose_score, coarse_constraints=rule,
                  receptor_npz=str(receptor), pocket=dict(minimum_heavy_atom_distance=1.,maximum_clashing_fraction=0.),
                  gaussian_weight=1., optional_weight=1., exclusions=[])
    q = dict(condition_policy=dict(required_anchors=['A','B'], minimum_score=.73, match_mode='any', coarse_constraints=rule),
             anchors=[dict(anchor_id='A',score_column=0,feature_index=0),
                      dict(anchor_id='B',score_column=1,feature_index=1)],
             guided_design=design, consensus_npz='fixture')
    monkeypatch.setattr(pre.full,'_FILTER_READER',SimpleNamespace(get=lambda gid:features))
    monkeypatch.setattr(pre.full,'_BOUND',NecessaryConditions(query,[0,1],'any',.73))
    monkeypatch.setattr(pre.full,'_JOINT',None)
    monkeypatch.setattr(pre,'prepare_seeds',lambda *a,**kw:(seeds,2))
    monkeypatch.setattr(guided_filters,'same_pose_gaussian',lambda original,candidate,seeds:np.full(len(seeds),.8))
    results = []
    for name, backend, mode in [('full','python',{}),('block','numpy',{'selection_mode':'threshold'})]:
        monkeypatch.setenv('AIDD_ASSIGNMENT_BACKEND',backend)
        monkeypatch.setattr(pre,'_GUIDED_CACHE',None)
        monkeypatch.setattr(pre.full,'_STATE',(SimpleNamespace(get=lambda gid:candidate),query,query,dict(q,**mode)))
        target = tmp_path/(name+'.npz')
        receipt = pre.compute((0,1,str(target),np.array([0])))
        results.append((receipt['counts'], target.with_suffix('.poses.jsonl').read_text()))
    assert results[0] == results[1]
    assert results[1][0]['matching_conformers'] == expected
    if expected:
        assert json.loads(results[1][1])['matched_anchors'] == ['A']
    if large:
        assert results[1][0]['original_seeds'] == 0


def test_complete_block_cache_reuses_only_identical_inputs_and_checks_copies(tmp_path):
    calls = []
    def refine(batch, output, design, ids, workers, chunk):
        calls.append(ids.tolist())
        for start in range(0,len(ids),chunk):
            p=output/'chunks/00'/f'{start:010d}.npz'
            p.parent.mkdir(parents=True,exist_ok=True)
            np.savez(p,global_ids=ids[start:start+chunk])
            rows=p.with_suffix('.poses.jsonl');rows.write_text('')
            save(p.with_suffix('.receipt.json'),dict(files=fingerprint([p,rows])))
    design=dict(templates=[dict(query_id='A')],block_search_mode='threshold',
                _block_cache=dict(root=str(tmp_path/'cache'),identity={'code':'v1'}))
    settings=dict(workers=1,refine_chunk=2)
    for i, expected in enumerate([False,True]):
        assert refine_cached(tmp_path,tmp_path/str(i),design,np.arange(3),settings,refine) == expected
    assert len(calls)==1
    local=next((tmp_path/'1').glob('chunks/*/*.receipt.json'))
    assert all(str(tmp_path/'1') in p for p in read(local)['files'])
    changed=copy.deepcopy(design);changed['_block_cache']['identity']['code']='v2'
    assert not refine_cached(tmp_path,tmp_path/'2',changed,np.arange(3),settings,refine)
    assert len(calls)==2
    seal=next((tmp_path/'cache').glob('*/complete.json'))
    artifact=next(iter(read(seal)['files']))
    from pathlib import Path
    Path(artifact).write_bytes(b'corrupt')
    original = design if read(seal)['identity']['inputs']['code']=='v1' else changed
    with pytest.raises(ValueError):
        refine_cached(tmp_path,tmp_path/'3',original,np.arange(3),settings,refine)
