import json
import sqlite3
from types import SimpleNamespace

import numpy as np
import pytest

from aidd_agent import preselection_full as pre
from aidd_agent.necessary_conditions import NecessaryConditions


def geometry():
    query = dict(shape_points=np.array([[0.,0,0],[1,0,0],[0,1,0]]),
                 feature_points=np.array([[0.,0,0],[1,0,0]]),
                 feature_types=np.array([1,2]), feature_directions=np.zeros((2,3)),
                 feature_direction_kinds=np.zeros(2,dtype=int),
                 anchored_weights=np.ones(2),anchor_feature_indices=np.array([0,1]))
    features = SimpleNamespace(feature_points=np.array([[0.,0,0],[10,0,0]]),
                               feature_types=np.array([1,2]),feature_directions=np.zeros((2,3)),
                               feature_kinds=np.zeros(2,dtype=int),molecule_id='m0',conformer_id='c0')
    matrices = [np.eye(4),np.eye(4)]
    matrices[1][0,3] = -9
    seeds = [SimpleNamespace(transform_matrix=m.reshape(-1),seed_id=str(i)) for i,m in enumerate(matrices)]
    return query,features,seeds


def test_alternative_pose_combinations_preserved_without_false_union():
    query, features, seeds = geometry()
    bound = NecessaryConditions(query,[0,1],'any',.5)
    possible = pre.possible_seed_mask(bound,features,seeds)
    rows, matched = pre.pose_representatives(features,seeds,possible,query,[0,1],.5)
    assert matched == 2
    assert [r['mask'] for r in rows] == [1,2]
    assert all(r['matched_anchor_count'] == 1 for r in rows)
    # Compare the accelerated seed subset with unpruned exact annotation.
    reference, _ = pre.pose_representatives(features,seeds,np.ones(2,dtype=bool),query,[0,1],.5)
    assert rows == reference


def test_best_representative_keeps_real_pose_and_honors_required_column_order():
    query, features, seeds = geometry()
    seeds.append(seeds[0])
    rows, matched = pre.pose_representatives(features,seeds,np.ones(3,dtype=bool),query,[1,0],.5)
    assert matched == 3 and len(rows) == 2
    by_mask = {r['mask']:r for r in rows}
    assert by_mask[2]['seed_index'] == 0
    assert by_mask[1]['seed_index'] == 1


def test_molecule_merge_preserves_separate_masks_and_best_conformer():
    db = sqlite3.connect(':memory:'); pre.schema(db)
    for mask, quality, gid in [(1,.7,3),(2,.8,2),(1,.9,4),(1,.8,1)]:
        pre.merge_pose(db,dict(molecule_id='m0',mask=mask,min_matched_score=quality,global_id=gid))
    assert db.execute('SELECT mask,gid FROM poses ORDER BY mask').fetchall() == [(1,4),(2,2)]
    db.close()


def test_coarse_any_does_not_require_two_feature_types():
    query, features, _ = geometry()
    candidate = SimpleNamespace(shape_points=query['shape_points'],feature_types=np.array([1]))
    features.feature_points=features.feature_points[:1]
    features.feature_types=np.array([1]); features.feature_directions=np.zeros((1,3)); features.feature_kinds=np.zeros(1,dtype=int)
    rule=dict(heavy_atom_ratio=[.7,1.3],minimum_feature_coverage=.5,maximum_extent_distance=.35)
    assert pre.coarse_stages(candidate,features,query,rule,NecessaryConditions(query,[0,1],'any',.5)) == 4
    assert pre.coarse_stages(candidate,features,query,rule,NecessaryConditions(query,[0,1],'all',.5)) == 3
    assert pre.coarse_stages(candidate,features,query,dict(rule,minimum_feature_coverage=.7),NecessaryConditions(query,[0,1],'any',.5)) == 1


def test_compute_seals_chunk_and_preserves_alternative_poses(tmp_path,monkeypatch):
    query, features, seeds = geometry()
    candidate = SimpleNamespace(**vars(features),shape_points=query['shape_points'])
    policy=dict(required_anchors=['A','B'],minimum_score=.5,match_mode='any',coarse_constraints={})
    q=dict(condition_policy=policy,anchors=[dict(anchor_id='A',score_column=0),dict(anchor_id='B',score_column=1)])
    monkeypatch.setattr(pre.full,'_STATE',(SimpleNamespace(get=lambda gid:candidate),query,query,q))
    monkeypatch.setattr(pre.full,'_FILTER_READER',SimpleNamespace(get=lambda gid:features))
    monkeypatch.setattr(pre.full,'_BOUND',NecessaryConditions(query,[0,1],'any',.5))
    monkeypatch.setattr(pre,'coarse_stages',lambda *args,**kwargs:4)
    monkeypatch.setattr(pre,'prepare_seeds',lambda *args,**kwargs:(seeds,2))
    monkeypatch.setattr(pre,'_score_ids',lambda *args,**kwargs:{pre.full.OBJECTIVE+'__objective':np.array([.9])})
    path=tmp_path/'chunk.npz'
    receipt=pre.compute((0,1,str(path)))
    assert receipt['counts']['matching_conformers']==1
    assert receipt['counts']['pose_combination_records']==2
    assert pre.verify_chunk(path,0,1)==receipt
    rows=[json.loads(line) for line in path.with_suffix('.poses.jsonl').read_text().splitlines()]
    assert [r['matched_anchors'] for r in rows]==[['A'],['B']]
    path.with_suffix('.poses.jsonl').write_text('corrupt')
    with pytest.raises(ValueError):pre.verify_chunk(path,0,1)


def test_scaffold_grouping_does_not_collapse_acyclic_structures():
    pytest.importorskip('rdkit')
    dtype=[('begin','i4'),('end','i4'),('order','i4')]
    methane=SimpleNamespace(atomic_numbers=[6],formal_charges=[0],bonds=np.array([],dtype=dtype))
    ethane=SimpleNamespace(atomic_numbers=[6,6],formal_charges=[0,0],bonds=np.array([(0,1,1)],dtype=dtype))
    assert pre.scaffold_key(methane)[0] != pre.scaffold_key(ethane)[0]


def test_aggregate_retains_members_and_pose_combinations(tmp_path,monkeypatch):
    from aidd_agent import chemical_companion
    root=tmp_path; (root/'chunks').mkdir()
    path=root/'chunks/chunk-00000000.npz'
    pre.full._atomic_savez(path,dict(global_ids=np.arange(3),molecule_ids=np.array(['m0','m0','m1']),stage_levels=np.array([7,7,7])))
    rows=[dict(molecule_id='m0',global_id=0,mask=1,min_matched_score=.8),
          dict(molecule_id='m0',global_id=1,mask=2,min_matched_score=.9),
          dict(molecule_id='m1',global_id=2,mask=1,min_matched_score=.7)]
    path.with_suffix('.poses.jsonl').write_text(''.join(json.dumps(r)+'\n' for r in rows))
    pre.record_json(path.with_suffix('.receipt.json'),dict(start=0,stop=3,
        files=pre.full.fingerprint([path,path.with_suffix('.poses.jsonl')])))
    monkeypatch.setattr(chemical_companion,'ChemicalCompanionReader',lambda _:SimpleNamespace(
        get=lambda gid:SimpleNamespace(molecule_id='m0' if gid<2 else 'm1')))
    monkeypatch.setattr(pre,'scaffold_key',lambda _:('shared','C'))
    result=pre.aggregate(root,dict(chemical_companion='unused'),3,3,2)
    assert result['matching_molecules']==2 and result['clusters']==1
    assert result['molecules_by_stage']==[2]*8
    assert result['pose_combination_records']==3
    assert result['combinations']==[dict(mask=1,molecules=2),dict(mask=2,molecules=1)]
    with sqlite3.connect(root/'candidates.sqlite') as db:
        assert db.execute("SELECT union_mask,pose_count FROM members WHERE mid='m0'").fetchone()==(3,2)
    with pytest.raises(ValueError,match='coverage mismatch'):
        pre.aggregate(root,dict(chemical_companion='unused'),3,3,4)


def test_full_orchestration_resume_any_policy_and_no_pilot(tmp_path,monkeypatch):
    import sys
    from concurrent.futures import Future
    from aidd_agent import chemical_companion
    monkeypatch.setitem(sys.modules,'rdkit',SimpleNamespace(rdBase=SimpleNamespace(rdkitVersion='fixture')))
    source=tmp_path/'source';source.mkdir()
    library=tmp_path/'library';library.mkdir()
    shard=library/'shard';shard.mkdir()
    pre.record_json(shard/'manifest.json',{})
    catalog=library/'catalog.json'
    pre.record_json(catalog,dict(conformers=3,shards=[dict(path=str(shard),name='shard',global_id_start=0,
        conformers=3,manifest_sha256=pre.full.ev.sha(shard/'manifest.json'))]))
    query=source/'query.npz';query.write_bytes(b'fixture')
    q=dict(query_id='TEST:LIG:A:1',query_npz=str(query),artifact_catalog=str(catalog),chemical_companion=str(catalog),
           anchors=[dict(anchor_id='A',feature_index=0),dict(anchor_id='B',feature_index=1)])
    evidence=source/'report.json'
    pre.record_json(evidence,dict(status='complete',kind='screening_evidence',queries=[q],sources={},
                                library=dict(library_conformers=3,library_molecules=2)))
    selection=source/'selection.json'
    policy=dict(required_anchors=['A','B'],match_mode='all',minimum_score=.5,
                coarse_constraints=dict(heavy_atom_ratio=[.7,1.3],minimum_feature_coverage=.7,maximum_extent_distance=.35))
    pre.record_json(selection,dict(kind='selection_preview',policy=policy,evidence_report=str(evidence),
                                  sources=pre.full.fingerprint([evidence])))
    monkeypatch.setattr(pre.full.ev,'verify_manifest',lambda *a,**k:{})
    initialized=[]
    class Pool:
        def __init__(self,**kwargs):initialized.append(kwargs['initargs'][0])
        def __enter__(self):return self
        def __exit__(self,*args):pass
        def submit(self,fn,task):
            f=Future()
            try:f.set_result(fn(task))
            except Exception as exc:f.set_exception(exc)
            return f
    monkeypatch.setattr(pre,'ProcessPoolExecutor',Pool)
    calls=[]
    def compute(task):
        start,stop,target=task;target=pre.Path(target);calls.append((start,stop))
        mids=np.array(['m0','m0','m1'])[start:stop]
        pre.full._atomic_savez(target,dict(global_ids=np.arange(start,stop),molecule_ids=mids,stage_levels=np.full(stop-start,7)))
        target.with_suffix('.poses.jsonl').write_text(''.join(json.dumps(dict(molecule_id=str(mid),global_id=gid,
            mask=1 if gid!=1 else 2,min_matched_score=.8))+'\n' for gid,mid in zip(range(start,stop),mids)))
        counts={key:stop-start for key in ('input_conformers','size_passed','coverage_passed','extent_passed',
            'anchor_bound_passed','pose_feasibility_passed','gaussian_evaluated_conformers','matching_conformers')}
        receipt=dict(start=start,stop=stop,counts=counts,worker_seconds={'gaussian':.1},
                     files=pre.full.fingerprint([target,target.with_suffix('.poses.jsonl')]))
        pre.record_json(target.with_suffix('.receipt.json'),receipt)
        return receipt
    monkeypatch.setattr(pre,'compute',compute)
    monkeypatch.setattr(chemical_companion,'ChemicalCompanionReader',lambda _:SimpleNamespace(
        get=lambda gid:SimpleNamespace(molecule_id='m0' if gid<2 else 'm1')))
    monkeypatch.setattr(pre,'scaffold_key',lambda _:('shared','C'))
    output=tmp_path/'output'
    result=pre.run(selection,output,workers=2,chunk_size=2)
    assert result['status']=='complete' and result['full_coverage']
    assert result['counts']['input_conformers']==3
    assert result['aggregation']['matching_molecules']==2
    assert result['stage_counts'][-1]['output_molecules']==2
    assert initialized[0]['condition_policy']['match_mode']=='any'
    assert pre.full.ev.read(selection)['policy']['match_mode']=='all'
    assert calls==[(0,2),(2,3)]
    resumed=pre.run(selection,output,workers=2,chunk_size=2)
    assert calls==[(0,2),(2,3)] and resumed['reused_conformers']==3
    with pytest.raises(ValueError,match='Changed protocol'):
        pre.run(selection,output,workers=2,chunk_size=1)
