from types import SimpleNamespace
import numpy as np
import pytest
from aidd_agent.gaussian_batch import _score_ids
from aidd_agent.full_library_screen import POSE_PARAMETERS
from aidd_agent.pose_acceleration import hardware_options, cross_overlaps, winner_indices
from aidd_agent.gaussian_overlay import RigidSeed, gaussian_overlap, apply_transform


def fixture():
    rng=np.random.default_rng(460)
    points=rng.normal(size=(12,3))*3
    features=points[:6]
    query=dict(shape_points=points, feature_points=features,feature_types=np.array([1,2,1,2,5,6]),
        anchored_weights=np.array([2.,3.,1.,1.,1.,1.]),anchor_feature_indices=np.array([0,1,2]))
    records=[]
    for i in range(5):
        records.append(SimpleNamespace(molecule_id='m'+str(i),conformer_id='c'+str(i),
            shape_points=points+rng.normal(size=points.shape)*.3,
            feature_points=features+rng.normal(size=features.shape)*.3,
            feature_types=query['feature_types']))
    return query,SimpleNamespace(get=lambda i:records[i])


def test_complete_arrays_and_winning_poses_match_reference():
    query,reader=fixture()
    reference=_score_ids(reader,query,np.arange(5),**POSE_PARAMETERS)
    batched=_score_ids(reader,query,np.arange(5),backend='numpy',**POSE_PARAMETERS)
    for key in reference:
        np.testing.assert_array_equal(reference[key],batched[key],err_msg=key)


@pytest.mark.parametrize('cutoff',[None,4.5])
def test_bounded_batches_reuse_color_without_changing_overlap(cutoff):
    query,reader=fixture();candidate=reader.get(0)
    seeds=[]
    for i in range(7):
        m=np.eye(4);m[:3,3]=[i*.1,0,0]
        seeds.append(RigidSeed(str(i),(0,1),(0,1),0,tuple(m.ravel())))
    values=cross_overlaps(query,candidate,seeds,sigma=1.,cutoff=cutoff,batch_size=2)
    for i,s in enumerate(seeds):
        moved=apply_transform(candidate.feature_points,s.transform_matrix)
        expected=gaussian_overlap(query['feature_points'],moved,sigma=1.,cutoff=cutoff,
            types_a=query['feature_types'],types_b=candidate.feature_types,weights_a=query['anchored_weights'])
        assert values[i,2] == expected


def test_hardware_override_and_gpu_single_owner(monkeypatch):
    monkeypatch.setenv('AIDD_POSE_BACKEND','numpy')
    monkeypatch.setenv('AIDD_SCREEN_WORKERS','16')
    assert hardware_options()==('numpy',16)
    with pytest.raises(ValueError,match='one GPU owner'):hardware_options('cupy',16)
    assert hardware_options('cupy',1)==('cupy',1)
    with pytest.raises(ValueError):hardware_options('bad',1)


def test_near_ties_are_reference_rechecked_in_original_order():
    q=dict(shape=1.,color_unweighted=1.,color_anchored=1.)
    c=dict(shape=1.,color=1.)
    cross=np.array([[1.,1.,1.],[1.+1e-12,1.,1.],[.1,.1,.1]])
    assert winner_indices(cross,q,c,contenders=True).tolist()==[0,1]


def test_cuda_optional_full_array_equivalence():
    cp=pytest.importorskip('cupy')
    if not cp.cuda.runtime.getDeviceCount():pytest.skip('No CUDA device')
    query,reader=fixture()
    ref=_score_ids(reader,query,np.arange(5),**POSE_PARAMETERS)
    gpu=_score_ids(reader,query,np.arange(5),backend='cupy',**POSE_PARAMETERS)
    for key in ref:np.testing.assert_array_equal(ref[key],gpu[key],err_msg=key)




def test_pilot_sampling_and_equivalence_detect_membership_and_pose_change():
    from aidd_agent.funnel_benchmark import sample_ids,compare
    assert sample_ids(1000000,256)[[0,-1]].tolist()==[0,999999]
    assert len(sample_ids(3,256))==3
    ref=dict(condition_passed=np.array([False,True]),transform=np.eye(4))
    assert compare(ref,ref)['passed']
    bad={**ref,'condition_passed':np.array([True,True])}
    assert not compare(ref,bad)['passed']


def test_benchmark_is_trusted_chat_action(tmp_path):
    import json
    from pathlib import Path
    from aidd_agent.chat_agent import ChatAgent,validate_route
    from test_screening_selection import attach_fixture
    from aidd_agent.gaussian_batch import _atomic_json
    app=ChatAgent(tmp_path,start=False)
    try:
        sid=app.new_session();jid=attach_fixture(app,sid);job=app.task(sid,jid)
        report=json.loads(Path(job['report']).read_text())
        for step in report['steps'].values():step['action']='select_screening'
        _atomic_json(Path(job['report']),report)
        reply=app.ask(sid,'/benchmark '+jid,'deepseek')
        assert 'bounded hardware' in reply
        reply=app.ask(sid,'/benchmark','deepseek')
        assert jid in reply
        plan=json.loads(Path(app.task(sid)['plan']).read_text())['plan']
        assert plan['steps'][0]['action']=='benchmark_funnel'
        validate_route(dict(intent='benchmark',message='Test hardware',request='',task_id=jid))
    finally:app.close()



def test_benchmark_report_uses_spread_rows_and_rejects_disagreement(tmp_path,monkeypatch):
    import json
    from aidd_agent import funnel_benchmark as bench
    from test_full_library_screen import fixture,data
    from aidd_agent.gaussian_batch import _atomic_json
    from aidd_agent.expanded_wee1 import fingerprint
    source=fixture(tmp_path)
    doc=json.loads(source.read_text());qpath=source.parent/'query.npz';qpath.write_bytes(b'fixture')
    doc['queries'][0]['query_npz']=str(qpath);_atomic_json(source,doc)
    selection=tmp_path/'selection';selection.mkdir()
    selected=selection/'report.json'
    policy=dict(required_anchors=['TEST:LIG:A:1/F0'],match_mode='all',minimum_score=.5)
    _atomic_json(selected,dict(kind='selection_preview',policy=policy,evidence_report=str(source),sources=fingerprint([source])))
    state={}
    def initialize(q):state.update(q)
    class Pool:
        def __init__(self,**kw):kw['initializer'](*kw['initargs'])
        def __enter__(self):return self
        def __exit__(self,*args):pass
        def map(self,fn,batches):return map(fn,batches)
    def work(ids):
        values=data(0,4);values['condition_passed']=np.array([True,False,False,True])
        return values,dict(prefilter_seconds=.01,pose_seconds=.1,annotation_seconds=.01,prefilter_passed=4,rejection_reasons={})
    monkeypatch.setattr(bench.full,'initialize',initialize)
    monkeypatch.setattr(bench,'ProcessPoolExecutor',Pool)
    monkeypatch.setattr(bench,'work',work)
    monkeypatch.setattr(bench,'hardware_options',lambda *a:('numpy',2))
    result=bench.run(selected,tmp_path/'bench',count=4,repeats=1,include_gpu=False)
    assert result['status']=='complete' and result['sample_count']==4
    assert result['sample_ids']==[0,1,2,3]
    assert all(r['checks'][0]['passed'] for r in result['scenarios'])
    original=work
    def wrong(ids):
        arrays,stats=original(ids)
        if state['pose_backend']=='numpy':arrays['condition_passed'][0]=False
        return arrays,stats
    monkeypatch.setattr(bench,'work',wrong)
    result=bench.run(selected,tmp_path/'bad-bench',count=4,repeats=1,include_gpu=False)
    assert result['status']=='failed'
    assert result['recommendation']['backend']=='reference'
