import json
from concurrent.futures import Future
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from aidd_agent import full_library_screen as full
from aidd_agent import screening_selection as selection
from aidd_agent.expanded_wee1 import fingerprint
from aidd_agent.gaussian_batch import _atomic_json, _atomic_savez
from aidd_agent.library_acceptance import sha

O = full.OBJECTIVE


def data(start, stop):
    # m0 meets two criteria in DIFFERENT poses; only the last conformer meets both.
    scores = np.array([[.9,.1],[.1,.9],[.1,.1],[.8,.8]])[start:stop]
    return dict(global_ids=np.arange(start,stop),
        molecule_ids=np.array(['m0','m0','m1','m2'])[start:stop],
        conformer_ids=np.array(['c0','c1','c2','c3'])[start:stop],
        query_anchor_feature_indices=np.array([0,1]),
        **{O+'__anchor_scores':scores, O+'__anchor_assignments':np.tile([0,1],(stop-start,1)),
           O+'__objective':np.ones(stop-start), O+'__transform':np.tile(np.eye(4),(stop-start,1,1))})


def save_chunk(path, start, stop):
    _atomic_savez(path,data(start,stop))
    result=dict(start=start,stop=stop,sha256=sha(path),wall_seconds=.1)
    _atomic_json(path.with_suffix('.receipt.json'),result)
    return result


def test_molecule_counts_merge_across_chunks_without_double_counting(tmp_path):
    counts=full.Counts(tmp_path/'counts.sqlite',2)
    try:
        counts.add(data(0,1));counts.add(data(1,4))
        assert counts.processed==4 and counts.total_molecules()==3
        assert counts.molecules.tolist()==[[2,2],[2,2],[2,2]]
        assert counts.conformers.tolist()==[[2,2],[2,2],[2,2]]
    finally:counts.db.close()


def fixture(tmp_path):
    source=tmp_path/'evidence';source.mkdir()
    shard=source/'shard';shard.mkdir()
    _atomic_json(shard/'manifest.json',dict(files={}))
    catalog=source/'catalog.json'
    _atomic_json(catalog,dict(conformers=4,shards=[dict(name='shard',path=str(shard),
        global_id_start=0,conformers=4,manifest_sha256=sha(shard/'manifest.json'))]))
    anchors=[dict(anchor_id=f'TEST:LIG:A:1/F{i}',feature_index=i,score_column=i,
        feature_class='HBD',interaction_class='hydrogen_bond',ligand_atom_indices=[i],
        evidence={},evidence_class='fixture',interpretation='Fixture only',diagnostic_counts=[]) for i in range(2)]
    report=dict(kind='screening_evidence',status='complete',classification='fixture',
        sources=fingerprint([catalog]),library=dict(library_conformers=4,library_molecules=3),
        classes={'hydrogen_bond':[a['anchor_id'] for a in anchors]},unsupported_classes=['halogen_bond'],
        queries=[dict(query_id='TEST:LIG:A:1',anchors=anchors,artifact_catalog=str(catalog),
                      chemical_companion=str(catalog),query_npz='unused-fixture')])
    _atomic_json(source/'report.json',report)
    return source/'report.json'


class InlinePool:
    def __init__(self, **kwargs):pass
    def __enter__(self):return self
    def __exit__(self,*args):pass
    def submit(self, fn, task):
        future=Future()
        try:future.set_result(fn(task))
        except Exception as exc:future.set_exception(exc)
        return future


def test_full_run_partial_resume_last_chunk_selection_and_export(tmp_path,monkeypatch):
    source=fixture(tmp_path);output=tmp_path/'full'
    calls=[]
    def compute(task):
        start,stop,path=task;calls.append((start,stop))
        return save_chunk(Path(path),start,stop)
    monkeypatch.setattr(full,'ProcessPoolExecutor',InlinePool)
    monkeypatch.setattr(full,'compute_chunk',compute)
    partial=full.run(source,output,workers=1,chunk_size=1,max_chunks=2)
    assert partial['status']=='partial' and not partial['queries'][0]['full_coverage']
    policy=dict(required_anchors=['TEST:LIG:A:1/F0','TEST:LIG:A:1/F1'],match_mode='all',minimum_score=.5)
    with pytest.raises(ValueError,match='completed full-library'):
        selection.preview(output/'report.json',tmp_path/'bad',policy)
    report=full.run(source,output,workers=2,chunk_size=1)
    assert calls==[(0,1),(1,2),(2,3),(3,4)]
    assert report['status']=='complete'
    assert report['queries'][0]['counts']['evaluated_conformers']==4
    assert all(a['diagnostic_counts'][1]['molecules']==2 for a in report['queries'][0]['anchors'])
    chosen=selection.preview(output/'report.json',tmp_path/'selected',policy)
    rows=[json.loads(line) for line in Path(chosen['representatives_jsonl']).read_text().splitlines()]
    assert chosen['counts']['matching_molecules']==1 and rows[0]['global_id']==3
    any_result=selection.preview(output/'report.json',tmp_path/'any',dict(policy,match_mode='any'))
    assert any_result['counts']['matching_molecules']==2
    zero=selection.preview(output/'report.json',tmp_path/'zero',dict(policy,minimum_score=1.))
    exported=selection.export(tmp_path/'zero/report.json',tmp_path/'export')
    assert exported['counts']['selected_molecules']==0
    assert json.loads(Path(exported['ids']).read_text())==[]
    from aidd_agent import chemical_companion
    def record(gid):
        return SimpleNamespace(global_id=gid,molecule_id='m2',conformer_id='c3',
            shape_points=np.array([[1.,2.,3.]]),atomic_numbers=[6],formal_charges=[0],bonds=[])
    reader=lambda *args:SimpleNamespace(get=record)
    monkeypatch.setattr(full,'ArtifactCatalogReader',reader)
    monkeypatch.setattr(chemical_companion,'ChemicalCompanionReader',reader)
    exported=selection.export(tmp_path/'selected/report.json',tmp_path/'nonempty-export')
    assert json.loads(Path(exported['ids']).read_text())[0]['global_id']==3
    assert 'AIDD_GLOBAL_ID>\n3' in Path(exported['sdf']).read_text()
    assert 'AIDD_MOLECULE_ID>\nm2' in Path(exported['sdf']).read_text()
    chunk=output/'query-0/chunk-00000003.npz';chunk.write_bytes(b'corrupt')
    with pytest.raises(ValueError,match='receipt'):
        full.run(source,output,workers=1,chunk_size=1)


def test_catalog_gaps_and_chunk_feature_changes_rejected(tmp_path):
    with pytest.raises(ValueError,match='gaps'):
        full.catalog_count(dict(conformers=2,shards=[dict(global_id_start=1,conformers=2)]))
    path=tmp_path/'chunk.npz';save_chunk(path,0,2)
    with pytest.raises(ValueError,match='IDs/features'):
        full.load_chunk(path,0,2,[1,0])
    source=fixture(tmp_path)
    report=json.loads(source.read_text());report['queries']=[]
    _atomic_json(source,report)
    with pytest.raises(ValueError,match='nonempty query'):
        full.run(source,tmp_path/'empty')


def test_compute_evaluates_every_id_and_preserves_gaussian_pose(tmp_path,monkeypatch):
    from test_gaussian_batch import _multi_artifact
    from aidd_agent.gaussian_batch import write_gaussian_query, _score_ids, ArtifactCatalogReader, _load_query
    catalog=_multi_artifact(tmp_path,3)
    query_path=tmp_path/'query.npz'
    write_gaussian_query(query_path,shape_points=[[0,0,0],[1,0,0],[0,2,0]],
        feature_points=[[0,0,0],[1,0,0]],feature_types=[1,2],anchored_weights=[1,1],anchor_feature_indices=[0,1])
    q=dict(artifact_catalog=str(catalog),chemical_companion='mocked',query_npz=str(query_path),
           anchors=[dict(feature_index=0),dict(feature_index=1)])
    _, query=_load_query(query_path)
    expected=_score_ids(ArtifactCatalogReader(catalog),query,np.arange(3),**full.POSE_PARAMETERS)
    def score(*args,**kwargs):
        ids=args[3];assert np.array_equal(ids,np.arange(3))
        assert np.array_equal(args[7][O],expected[O+'__transform'])
        return {O:np.ones(3)},{O:np.tile([0,1],(3,1))},{O:np.ones((3,2))}
    monkeypatch.setattr(full,'score_batched',score)
    full.initialize(q);target=tmp_path/'computed.npz'
    full.compute_chunk((0,3,str(target)))
    result=full.load_chunk(target,0,3,[0,1])
    assert np.array_equal(result[O+'__transform'],expected[O+'__transform'])
    assert np.array_equal(result['global_ids'],np.arange(3))


@pytest.mark.parametrize('intent,source_action,action',[('full_count','classify_screening','full_library_screen'),
                                                      ('funnel','select_screening','condition_funnel')])
def test_full_count_is_local_chat_action_without_threshold_question(tmp_path,intent,source_action,action):
    from aidd_agent.chat_agent import ChatAgent, validate_route
    from test_screening_selection import attach_fixture
    app=ChatAgent(tmp_path,start=False)
    try:
        sid=app.new_session();jid=attach_fixture(app,sid);job=app.task(sid,jid)
        report=json.loads(Path(job['report']).read_text())
        for step in report['steps'].values():step['action']=source_action
        _atomic_json(Path(job['report']),report)
        reply=app.ask(sid,'/'+intent+' '+jid,'deepseek')
        assert 'Top-K' in reply and 'No new search' not in reply
        plan=json.loads(Path(app.task(sid)['plan']).read_text())['plan']
        assert plan['steps'][0]['action']==action
        validate_route(dict(intent=intent,message='Evaluate every conformer',request='',task_id=jid))
    finally:app.close()
