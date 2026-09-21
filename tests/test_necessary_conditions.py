from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from aidd_agent.necessary_conditions import NecessaryConditions, distinct_assignment
from aidd_agent.interaction_matching import interaction_match
from aidd_agent import full_library_screen as full
from aidd_agent.expanded_wee1 import fingerprint
from aidd_agent.gaussian_batch import _atomic_json, _atomic_savez

O = full.OBJECTIVE


def query(points, types=None, kinds=None):
    n=len(points)
    return dict(feature_points=np.array(points,dtype=float),feature_types=np.array(types or [1]*n),
                feature_direction_kinds=np.array(kinds or [0]*n))


def candidate(points, types=None, kinds=None):
    q=query(points,types,kinds)
    return SimpleNamespace(feature_points=q['feature_points'].reshape(-1,3),
        feature_types=q['feature_types'],feature_kinds=q['feature_direction_kinds'])


def test_type_cardinality_and_pair_distance_rejections_are_not_rank_budgets():
    q=query([[0,0,0],[10,0,0]])
    bound=NecessaryConditions(q,[0,1],'all',.5)
    assert bound.check(candidate([[0,0,0]]))[1]=='distinct_assignment'
    assert bound.check(candidate([[0,0,0],[30,0,0]]))[1]=='pair_distance'
    assert bound.check(candidate([[200,0,0],[210,0,0]]))[0]
    assert NecessaryConditions(q,[0,1],'any',.5).check(candidate([[0,0,0]]))[0]
    assert NecessaryConditions(q,[0,0],'all',.5).check(candidate([[0,0,0]]))[0]
    directional=NecessaryConditions(query([[0,0,0]],[1],[1]),[0],'all',.5)
    assert not directional.check(candidate([[0,0,0]],[1],[0]))[0]
    assert not distinct_assignment(np.array([[1,0,0],[1,0,0],[0,1,1]],dtype=bool))


@pytest.mark.parametrize('threshold',[.25,.5,.75,1.])
def test_passing_rigid_poses_survive_random_rotations_and_threshold_boundary(threshold):
    rng=np.random.default_rng(631)
    radius=np.sqrt(-2*np.log(threshold))
    for n in range(2,7):
        points=rng.normal(size=(n,3))*8
        types=np.arange(n)+1
        q=query(points,types.tolist())
        bound=NecessaryConditions(q,range(n),'all',threshold)
        for _ in range(15):
            noise=rng.normal(size=(n,3))
            noise*=radius*.8/np.maximum(np.linalg.norm(noise,axis=1)[:,None],1e-12)
            moved=points+noise
            exact=interaction_match(points,types,np.zeros((n,3)),np.zeros(n),np.ones(n),
                moved,types,np.zeros((n,3)),np.zeros(n),sigma=1.,cutoff=4.5,angular_power=2.)
            assert np.all(exact['anchor_scores']>=threshold-1e-12)
            rotation,_=np.linalg.qr(rng.normal(size=(3,3)))
            raw=(moved-rng.normal(size=3)*100)@rotation
            assert bound.check(candidate(raw,types.tolist()))[0]
    q=query([[0,0,0],[10,0,0]],[1,2])
    assert NecessaryConditions(q,[0,1],'all',threshold).check(
        candidate([[-radius,0,0],[10+radius,0,0]],[1,2]))[0]


def test_funnel_skips_impossible_rows_and_preserves_uncapped_pass_set(tmp_path,monkeypatch):
    from test_full_library_screen import fixture, InlinePool
    source=fixture(tmp_path)
    policy=dict(required_anchors=['TEST:LIG:A:1/F0','TEST:LIG:A:1/F1'],match_mode='all',minimum_score=.5)
    selection=tmp_path/'selection.json'
    _atomic_json(selection,dict(kind='selection_preview',policy=policy,evidence_report=str(source),
                               sources=fingerprint([source])))
    records=[]
    for gid, points in enumerate(([[0,0,0]],[[0,0,0],[30,0,0]],[[0,0,0],[10,0,0]],[[100,0,0],[110,0,0]])):
        c=candidate(points,[1] if gid==0 else [1,2]);c.molecule_id=['m0','m0','m1','m2'][gid]
        c.conformer_id='c'+str(gid);records.append(c)
    original=query([[0,0,0],[10,0,0]],[1,2])
    original['anchor_feature_indices']=np.array([0,1])
    def initialize(q):
        q['pose_feasibility']=False  # Fixture substitutes scoring; real seed path is tested separately.
        full._STATE=(None,original,original,q)
        full._FILTER_READER=SimpleNamespace(get=lambda gid:records[gid])
        full._BOUND=NecessaryConditions(original,[0,1],'all',.5)
    scored=[]
    def pose(reader,q,ids,**kwargs):
        scored.extend(ids.tolist())
        return dict(global_ids=ids,molecule_ids=np.array([records[i].molecule_id for i in ids]),
            conformer_ids=np.array([records[i].conformer_id for i in ids]),
            **{O+'__objective':np.ones(len(ids)),O+'__transform':np.tile(np.eye(4).reshape(1,16),(len(ids),1))})
    def score(*args,**kwargs):
        ids=args[3]
        values=np.array([[1.,1.] if i==2 else [0.,0.] for i in ids])
        return {O:values.mean(axis=1)},{O:np.tile([0,1],(len(ids),1))},{O:values}
    class Pool(InlinePool):
        def __init__(self,**kwargs):kwargs['initializer'](*kwargs['initargs'])
    monkeypatch.setattr(full,'initialize',initialize)
    monkeypatch.setattr(full,'ProcessPoolExecutor',Pool)
    monkeypatch.setattr(full,'_score_ids',pose)
    monkeypatch.setattr(full,'score_batched',score)
    result=full.run_funnel(selection,tmp_path/'funnel',workers=1,chunk_size=2)
    assert scored==[2,3]
    counts=result['queries'][0]['counts']
    assert counts['evaluated_conformers']==4 and counts['prefilter_rejected']==2
    assert counts['pose_evaluated_conformers']==2 and counts['matching_molecules']==1
    from aidd_agent.screening_selection import preview
    selected=preview(tmp_path/'funnel/report.json',tmp_path/'preview',policy)
    assert selected['counts']['matching_conformers']==selected['counts']['matching_molecules']==1
    assert '"global_id": 2' in Path(selected['representatives_jsonl']).read_text()
    with pytest.raises(ValueError,match='Changed rule'):
        preview(tmp_path/'funnel/report.json',tmp_path/'changed',dict(policy,match_mode='any'))
    # Independent unfiltered fixture reference passes only ID 2; filtered output agrees.
    all_values=score(None,None,None,np.arange(4))[2][O]
    assert np.flatnonzero((all_values>=.5).all(axis=1)).tolist()==[2]

    # Reusing a fully computed legacy block avoids Gaussian work but audits the bound.
    legacy=tmp_path/'legacy';legacy.mkdir()
    full._STATE[3]['reuse_directory']=str(legacy);full._STATE[3]['reuse_chunk_size']=2
    arrays=dict(global_ids=np.arange(2),molecule_ids=np.array(['m0','m0']),conformer_ids=np.array(['c0','c1']),
        query_anchor_feature_indices=np.array([0,1]),**{O+'__objective':np.ones(2),
        O+'__transform':np.tile(np.eye(4).reshape(1,16),(2,1)),
        O+'__anchor_scores':np.zeros((2,2)),O+'__anchor_assignments':np.tile([0,1],(2,1))})
    old=legacy/'chunk-00000000.npz';_atomic_savez(old,arrays)
    _atomic_json(old.with_suffix('.receipt.json'),dict(start=0,stop=2,sha256=full.ev.sha(old)))
    receipt=full.compute_chunk((0,2,str(tmp_path/'reused.npz')))
    assert receipt['legacy_chunk_reused'] and scored==[2,3]
    arrays[O+'__anchor_scores'][1]=1.
    _atomic_savez(old,arrays);_atomic_json(old.with_suffix('.receipt.json'),dict(start=0,stop=2,sha256=full.ev.sha(old)))
    with pytest.raises(ValueError,match='rejected an existing passing pose'):
        full.compute_chunk((0,2,str(tmp_path/'invalid.npz')))


@pytest.mark.parametrize('kinds',[[1,1],[2,2],[1,2]])
@pytest.mark.parametrize('threshold',[.25,.5,.75,1.])
def test_direction_bounds_preserve_true_rigid_matches(kinds,threshold):
    rng=np.random.default_rng(462)
    for _ in range(30):
        points=rng.normal(size=(2,3))
        directions=rng.normal(size=(2,3));directions/=np.linalg.norm(directions,axis=1)[:,None]
        q=query(points,[1,2],kinds);q['feature_directions']=directions
        rot,_=np.linalg.qr(rng.normal(size=(3,3)))
        c=candidate(points@rot,[1,2],kinds);c.feature_directions=directions@rot
        assert NecessaryConditions(q,[0,1],'all',threshold).check(c)[0]
        # Perturb both directions within the score's angular cone.
        perturbed=directions+rng.normal(size=(2,3))*.1
        perturbed/=np.linalg.norm(perturbed,axis=1)[:,None]
        c.feature_directions=perturbed@rot
        scores=interaction_match(points,[1,2],directions,kinds,[1,1],points,[1,2],perturbed,kinds)['anchor_scores']
        if np.all(np.array(scores)>=threshold):
            assert NecessaryConditions(q,[0,1],'all',threshold).check(c)[0]


def test_direction_pair_can_reject_geometry_that_distance_alone_retains():
    q=query([[0,0,0],[2,0,0]],[1,2],[1,1]);q['feature_directions']=np.array([[1.,0,0],[1.,0,0]])
    c=candidate([[0,0,0],[2,0,0]],[1,2],[1,1]);c.feature_directions=np.array([[1.,0,0],[-1.,0,0]])
    assert not NecessaryConditions(q,[0,1],'all',.75).check(c)[0]
    assert NecessaryConditions(q,[0,1],'any',.75).check(c)[0]
