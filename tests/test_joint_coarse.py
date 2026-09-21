from types import SimpleNamespace
import copy
import numpy as np
import pytest

from aidd_agent.joint_coarse import JointCoarse, validate_constraints
from aidd_agent.full_library_screen import normalize_policy
from aidd_agent.necessary_conditions import pose_mask
from aidd_agent.screening_selection import select_rows


RULE = dict(heavy_atom_ratio=[.7, 1.3], maximum_extent_distance=.2, minimum_feature_coverage=.75)


def fixture():
    points = np.array([[0., 0, 0], [2, 0, 0], [0, 3, 0], [0, 0, 4]])
    query = dict(shape_points=points, feature_types=np.array([1, 1, 2, 3]))
    return query, SimpleNamespace(**query)


def test_joint_predicates_rigid_invariance_and_independent_reasons():
    q, candidate = fixture()
    gate = JointCoarse(q, RULE)
    assert gate.check(candidate)[0]
    candidate.shape_points = candidate.shape_points[:, [1, 2, 0]] + 100
    assert gate.check(candidate)[0]
    candidate.feature_types = np.array([1, 2])
    assert gate.check(candidate) == (False, 'joint_feature_coverage')
    candidate.feature_types = np.array([1, 1, 2])
    assert gate.check(candidate)[0]  # Inclusive count-coverage threshold.
    candidate.shape_points *= 2
    assert gate.check(candidate) == (False, 'joint_shape_extent')
    candidate.shape_points = candidate.shape_points[:2]
    assert gate.check(candidate) == (False, 'joint_heavy_atom_ratio')


@pytest.mark.parametrize('field,value', [('heavy_atom_ratio',[2,1]), ('heavy_atom_ratio',[True,2]),
    ('maximum_extent_distance',float('nan')), ('minimum_feature_coverage',0), ('unknown',1)])
def test_bad_joint_rules_fail_closed(field, value):
    rule = copy.deepcopy(RULE); rule[field] = value
    with pytest.raises(ValueError): validate_constraints(rule)


def test_policy_round_trip_and_same_membership_in_preview_and_pose_gate():
    policy = dict(required_anchors=['a'], match_mode='all', minimum_score=.5,
                  coarse_constraints=copy.deepcopy(RULE))
    assert normalize_policy(policy) == policy
    r = dict(global_ids=np.array([0,1]), molecule_ids=np.array(['a','b']),
        conformer_ids=np.array(['c','d']), joint_eligible=np.array([True,False]),
        x__anchor_scores=np.array([[.9],[.9]]), x__anchor_assignments=np.array([[0],[0]]))
    assert pose_mask(r,[0],policy,'x').tolist() == [True,False]
    from aidd_agent.interaction_review import OBJECTIVE
    r.update({OBJECTIVE+'__objective':np.array([.5,.6]),
              OBJECTIVE+'__anchor_scores':r['x__anchor_scores'],
              OBJECTIVE+'__anchor_assignments':r['x__anchor_assignments']})
    rows, count = select_rows(r,r,[0],policy,'query')
    assert count == 1 and rows[0]['global_id'] == 0
    del r['joint_eligible']
    with pytest.raises(ValueError,match='Missing joint'): pose_mask(r,[0],policy,'x')


def test_preseed_joint_rejection_bypasses_anchor_and_pose_work(monkeypatch):
    from aidd_agent import full_library_screen as full
    q,c = fixture(); c.molecule_id='m'; c.conformer_id='c'
    rejected = copy.deepcopy(c); rejected.shape_points *= 3
    reader = SimpleNamespace(get=lambda gid: [c,rejected][gid])
    class AnchorGate:
        def check(self, record):
            assert record is c
            return True, 'possible'
    monkeypatch.setattr(full,'_STATE',(reader,q,q,{}))
    monkeypatch.setattr(full,'_JOINT',JointCoarse(q,RULE))
    monkeypatch.setattr(full,'_BOUND',AnchorGate())
    decisions,joint=full.coarse_decisions(np.array([0,1]),[c,rejected])
    assert joint.tolist() == [True,False]
    assert decisions == [(True,'possible'),(False,'joint_shape_extent')]


def test_joint_funnel_writes_all_ids_and_skips_seed_generation_for_rejects(tmp_path,monkeypatch):
    from test_pose_feasibility import fixture as pose_fixture
    from aidd_agent import full_library_screen as full
    from aidd_agent.necessary_conditions import NecessaryConditions
    from aidd_agent.interaction_review import OBJECTIVE
    q,reader,records=pose_fixture()
    reader.get(1).shape_points = reader.get(1).shape_points * 3
    reader.get(2).shape_points = reader.get(2).shape_points * 4
    policy=dict(required_anchors=['a'],match_mode='all',minimum_score=.5,coarse_constraints=RULE)
    expanded={**q,'anchor_feature_indices':np.array([0])}
    definition=dict(condition_policy=policy,anchors=[dict(anchor_id='a',score_column=0)],
                    pose_feasibility=False,pose_backend='numpy',artifact_catalog='x',chemical_companion='y')
    monkeypatch.setattr(full,'_STATE',(reader,q,expanded,definition))
    monkeypatch.setattr(full,'_FILTER_READER',SimpleNamespace(get=lambda gid:records[gid]))
    monkeypatch.setattr(full,'_JOINT',JointCoarse(q,RULE))
    monkeypatch.setattr(full,'_BOUND',NecessaryConditions(q,[0],'all',.5))
    original=full._score_ids
    def score(reader, query, ids, **kwargs):
        assert ids.tolist()==[0]
        return original(reader,query,ids,**kwargs)
    monkeypatch.setattr(full,'_score_ids',score)
    def annotate(*args,**kwargs):
        return {},{OBJECTIVE:np.array([[0]])},{OBJECTIVE:np.array([[1.]])}
    monkeypatch.setattr(full,'score_batched',annotate)
    path=tmp_path/'chunk.npz'; receipt=full.compute_chunk((0,3,str(path)))
    data=full.load_chunk(path,0,3,[0])
    assert data['global_ids'].tolist()==[0,1,2]
    assert data['joint_eligible'].tolist()==[True,False,False]
    assert data['condition_passed'].tolist()==[True,False,False]
    assert receipt['rejection_reasons']=={'joint_shape_extent':2}


def test_coarse_audit_never_calls_scoring(tmp_path,monkeypatch):
    import json
    from test_full_library_screen import fixture as evidence_fixture
    from aidd_agent import funnel_benchmark as bench
    from aidd_agent.gaussian_batch import _atomic_json
    from aidd_agent.expanded_wee1 import fingerprint
    source=evidence_fixture(tmp_path)
    doc=json.loads(source.read_text()); query_path=source.parent/'query.npz'; query_path.write_bytes(b'fixture')
    doc['queries'][0]['query_npz']=str(query_path); _atomic_json(source,doc)
    selected=tmp_path/'selection'/'report.json'; selected.parent.mkdir()
    policy=dict(required_anchors=['TEST:LIG:A:1/F0'],match_mode='all',minimum_score=.5,coarse_constraints=RULE)
    _atomic_json(selected,dict(kind='selection_preview',policy=policy,evidence_report=str(source),sources=fingerprint([source])))
    monkeypatch.setattr(bench.full,'initialize',lambda q:None)
    monkeypatch.setattr(bench.full,'_FILTER_READER',SimpleNamespace(get=lambda gid:gid))
    monkeypatch.setattr(bench.full,'coarse_decisions',lambda ids,records:
        ([(i==0,'possible' if i==0 else 'joint_shape_extent') for i in ids],ids==0))
    def forbidden(*args,**kwargs): raise AssertionError('Coarse audit must not score poses')
    monkeypatch.setattr(bench,'work',forbidden)
    monkeypatch.setattr(bench.full,'_score_ids',forbidden)
    result=bench.run(selected,tmp_path/'audit',count=4,coarse_only=True)
    assert result['kind']=='coarse_audit' and result['status']=='complete'
    assert result['retained_conformers']==1 and result['rejection_fraction']==.75
    assert result['seed_generation_count']==result['gaussian_evaluated']==0
    assert not result['sample_rejection_target_met']


def test_coarse_chat_uses_latest_completed_selection(tmp_path):
    import json
    from pathlib import Path
    from aidd_agent.chat_agent import ChatAgent
    from test_screening_selection import attach_fixture
    from aidd_agent.gaussian_batch import _atomic_json
    app=ChatAgent(tmp_path,start=False)
    try:
        sid=app.new_session(); jid=attach_fixture(app,sid); job=app.task(sid,jid)
        report=json.loads(Path(job['report']).read_text())
        for step in report['steps'].values(): step['action']='select_screening'
        _atomic_json(Path(job['report']),report)
        assert jid in app.ask(sid,'/coarse','deepseek')
        step=json.loads(Path(app.task(sid)['plan']).read_text())['plan']['steps'][0]
        assert step['action']=='benchmark_funnel' and step['params']['coarse_only'] is True
    finally: app.close()
