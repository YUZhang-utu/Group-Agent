from types import SimpleNamespace
import numpy as np
import pytest
from aidd_agent.gaussian_overlay import RigidSeed
from aidd_agent.gaussian_batch import prepare_seeds,_score_ids,OBJECTIVES
from aidd_agent.necessary_conditions import NecessaryConditions
from aidd_agent.pose_feasibility import possible_seed_mask,prepare_survivors
from aidd_agent.interaction_fast import match_batch
from aidd_agent.full_library_screen import POSE_PARAMETERS
from aidd_agent.interaction_review import OBJECTIVE


def seed(translation=(0,0,0)):
    m=np.eye(4);m[:3,3]=translation
    return RigidSeed('test',(0,1),(0,1),0,tuple(m.ravel()))


def feature_record(points,types,kinds=None,directions=None):
    n=len(points)
    return SimpleNamespace(feature_points=np.array(points,dtype=float).reshape(-1,3),
        feature_types=np.array(types),feature_kinds=np.array(kinds if kinds is not None else [0]*n),
        feature_directions=np.array(directions,dtype=float) if directions is not None else np.zeros((n,3)))


def bound(points,types,mode='all',threshold=.5,kinds=None,directions=None):
    c=feature_record(points,types,kinds,directions)
    q=dict(feature_points=c.feature_points,feature_types=c.feature_types,
        feature_direction_kinds=c.feature_kinds,feature_directions=c.feature_directions)
    return NecessaryConditions(q,range(len(types)),mode,threshold)


def test_all_must_be_possible_in_same_pose_and_boundary_survives():
    b=bound([[0,0,0],[10,0,0]],[1,2])
    c=feature_record([[0,0,0],[30,0,0]],[1,2])
    assert not possible_seed_mask(b,c,[seed(),seed((-20,0,0))]).any()
    b.mode='any'
    assert possible_seed_mask(b,c,[seed(),seed((-20,0,0))]).all()
    b=bound([[0,0,0]],[1]);r=np.sqrt(-2*np.log(.5))
    assert possible_seed_mask(b,feature_record([[r,0,0]],[1]),[seed()])[0]


@pytest.mark.parametrize('mode',['all','any'])
@pytest.mark.parametrize('threshold',[.25,.5,.75,1.])
def test_actual_assignment_passes_survive_optimistic_seed_check(mode,threshold):
    rng=np.random.default_rng(470)
    for _ in range(25):
        qp=rng.normal(size=(4,3));types=np.array([1,2,1,6]);kinds=np.array([1,1,0,2])
        qd=rng.normal(size=(4,3));qd/=np.linalg.norm(qd,axis=1)[:,None]
        b=bound(qp,types,mode,threshold,kinds,qd)
        points=qp+rng.normal(size=(4,3))*.2
        c=feature_record(points,types,kinds,qd)
        seeds=[seed(),seed((1,0,0)),seed((-1,0,0))]
        matrices=np.array([s.transform_matrix for s in seeds])
        _,assignments,scores=match_batch((qp,types,qd,kinds,np.ones(4)),
            np.repeat(points[None],3,axis=0),np.repeat(types[None],3,axis=0),
            np.repeat(qd[None],3,axis=0),np.repeat(kinds[None],3,axis=0),matrices)
        hits=(scores>=threshold)&(assignments>=0)
        actual=hits.all(axis=1) if mode=='all' else hits.any(axis=1)
        assert not np.any(actual & ~possible_seed_mask(b,c,seeds))


def fixture():
    rng=np.random.default_rng(471)
    shape=rng.normal(size=(12,3))*2
    f=shape[:6].copy();types=np.array([1,2,1,2,5,6])
    q=dict(shape_points=shape,feature_points=f,feature_types=types,anchored_weights=np.ones(6),
        anchor_feature_indices=np.array([0,1]),feature_direction_kinds=np.zeros(6,dtype=int),
        feature_directions=np.zeros((6,3)))
    records=[];artifacts=[]
    for i in range(3):
        points=f.copy()
        if i:points[4:]+=i*100
        c=feature_record(points,types);c.molecule_id='m'+str(i);c.conformer_id='c'+str(i)
        records.append(c)
        artifacts.append(SimpleNamespace(shape_points=shape,feature_points=points,feature_types=types,
            molecule_id=c.molecule_id,conformer_id=c.conformer_id))
    return q,SimpleNamespace(get=lambda i:artifacts[i]),records


def test_reuse_all_seeds_and_avoid_gaussian_work_on_impossible_conformers():
    q,reader,records=fixture();b=NecessaryConditions(q,[4,5],'all',.75)
    ids=np.arange(3);invariant=np.array([b.check(c)[0] for c in records])
    assert invariant.all()
    keep,prepared,stats=prepare_survivors(reader,q,b,records,ids,invariant,POSE_PARAMETERS)
    assert keep.tolist()==[True,False,False]
    assert stats['pose_feasibility_rejected']==2
    ref=_score_ids(reader,q,ids,**POSE_PARAMETERS)
    fast=_score_ids(reader,q,ids[keep],prepared=prepared,backend='numpy',**POSE_PARAMETERS)
    for key in fast:
        expected=ref[key] if key=='objective_names' else ref[key][keep]
        np.testing.assert_array_equal(expected,fast[key],err_msg=key)
    # Every original seed, including nonpassing seeds, still competes on survivors.
    seeds,count=prepare_seeds(reader.get(0),q,**POSE_PARAMETERS)
    assert prepared[0][1]==seeds and prepared[0][2]==count
    # Reference actual passes: identity match is positive, shifted selected features fail.
    passed=[]
    for i,c in enumerate(records):
        _,a,s=match_batch((q['feature_points'][[4,5]],q['feature_types'][[4,5]],
            q['feature_directions'][[4,5]],q['feature_direction_kinds'][[4,5]],np.ones(2)),
            c.feature_points[None],c.feature_types[None],c.feature_directions[None],c.feature_kinds[None],
            ref[OBJECTIVE+'__transform'][i:i+1])
        passed.append(bool(((s>=.75)&(a>=0)).all()))
    assert passed==[True,False,False]


def test_filtered_comparison_detects_missing_positive_and_keeps_checked_ids():
    from aidd_agent.funnel_benchmark import compare_filtered
    ref=dict(global_ids=np.array([0,1]),checked_global_ids=np.array([0,1]),condition_passed=np.array([True,False]),
             score=np.array([.9,.1]))
    actual=dict(global_ids=np.array([0]),checked_global_ids=np.array([0,1]),condition_passed=np.array([True]),score=np.array([.9]))
    assert compare_filtered(ref,actual)['passed']
    bad=dict(global_ids=np.array([1]),checked_global_ids=np.array([0,1]),condition_passed=np.array([False]),score=np.array([.1]))
    check=compare_filtered(ref,bad)
    assert not check['passed'] and check['false_rejections']==1



def test_full_worker_keeps_coverage_but_skips_proven_negative_gaussian(tmp_path,monkeypatch):
    from aidd_agent import full_library_screen as full
    q,reader,records=fixture()
    expanded={**q,'anchor_feature_indices':np.array([4,5])}
    definition=dict(condition_policy=dict(required_anchors=['f4','f5'],match_mode='all',minimum_score=.75),
        anchors=[dict(anchor_id='f4',score_column=0),dict(anchor_id='f5',score_column=1)],
        pose_feasibility=True,pose_backend='numpy',artifact_catalog='fixture',chemical_companion='fixture')
    monkeypatch.setattr(full,'_STATE',(reader,q,expanded,definition))
    monkeypatch.setattr(full,'_FILTER_READER',SimpleNamespace(get=lambda i:records[i]))
    monkeypatch.setattr(full,'_BOUND',NecessaryConditions(q,[4,5],'all',.75))
    def annotate(artifact,chemical,query,ids,mids,cids,objectives,transforms,**kwargs):
        assert ids.tolist()==[0]
        indices=query['anchor_feature_indices'];c=records[0]
        values=tuple(query[k][indices] for k in ('feature_points','feature_types','feature_directions','feature_direction_kinds','anchored_weights'))
        score,assign,anchor=match_batch(values,c.feature_points[None],c.feature_types[None],
            c.feature_directions[None],c.feature_kinds[None],transforms[OBJECTIVE])
        return {OBJECTIVE:score},{OBJECTIVE:assign},{OBJECTIVE:anchor}
    monkeypatch.setattr(full,'score_batched',annotate)
    path=tmp_path/'chunk.npz';receipt=full.compute_chunk((0,3,str(path)))
    arrays=full.load_chunk(path,0,3,[4,5])
    assert arrays['global_ids'].tolist()==[0,1,2]
    assert arrays['prefilter_passed'].all()
    assert arrays['pose_evaluated'].tolist()==[True,False,False]
    assert arrays['condition_passed'].tolist()==[True,False,False]
    assert receipt['pose_candidates']==1 and receipt['pose_feasibility_rejected']==2
    counts=full.Counts(tmp_path/'counts.sqlite',2)
    try:
        counts.add(arrays)
        assert counts.processed==3 and counts.pose_candidates==1 and counts.matching_conformers==1
    finally:counts.db.close()


def test_evidence_panel_includes_saved_passes_and_boundary(tmp_path):
    from aidd_agent.funnel_benchmark import evidence_panel
    from aidd_agent.gaussian_batch import _atomic_savez
    path=tmp_path/'side.npz';rigid=tmp_path/'rigid.npz'
    _atomic_savez(rigid,dict(global_ids=np.array([11,12,13])))
    _atomic_savez(path,dict(global_ids=np.array([11,12,13]),**{
        OBJECTIVE+'__anchor_scores':np.array([[.9],[.499],[.1]]),
        OBJECTIVE+'__anchor_assignments':np.array([[0],[0],[0]])}))
    q=dict(sidecar=str(path),rigid=str(rigid),anchors=[dict(anchor_id='a',score_column=0)])
    ids,panel=evidence_panel(q,dict(required_anchors=['a'],match_mode='all',minimum_score=.5),1000,4)
    assert panel['saved_positive_ids']==[11]
    assert {11,12,13}<=set(ids) and 0 in ids and 999 in ids
