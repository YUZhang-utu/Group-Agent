from types import SimpleNamespace

import numpy as np
import pytest

from aidd_agent.gaussian_overlay import pair_alignment_seeds
from aidd_agent.gaussian_batch import prepare_seeds
from aidd_agent.seed_search import prepare, policy
from aidd_agent.full_library_screen import POSE_PARAMETERS


@pytest.mark.parametrize('case',range(8))
def test_batched_ordered_prefix(case):
    rng=np.random.default_rng(case)
    c,q=rng.normal(size=(12,3)),rng.normal(size=(5,3))
    if case==0:c[:]=0
    if case==1:c[1:]=c[0]+np.arange(11)[:,None]
    if case==2:q=-c[:5]
    if case==3:c[1]=c[0];q[1]=q[0]
    ct,qt=np.ones(len(c)),np.ones(len(q))
    if case==4:ct[:]=9
    for cap in (0,1,17,128,512):
        assert pair_alignment_seeds(c,ct,q,qt,max_seeds=cap,backend='batched') == pair_alignment_seeds(c,ct,q,qt,max_seeds=cap)


def fixture():
    rng=np.random.default_rng(193)
    c=SimpleNamespace(shape_points=rng.normal(size=(30,3)),feature_points=rng.normal(size=(12,3)),feature_types=np.ones(12))
    q=dict(shape_points=rng.normal(size=(25,3)),feature_points=rng.normal(size=(6,3)),feature_types=np.ones(6),anchor_feature_indices=np.arange(6))
    return c,q


@pytest.mark.parametrize('backend',['reference','batched'])
def test_generated_budget_preserves_default_seeds_and_mask(backend):
    c,q=fixture()
    guided=SimpleNamespace(pocket_mask=lambda p,m:m[:,0,3]>0)
    original,_=prepare_seeds(c,q,**POSE_PARAMETERS)
    seeds,mask,stats,_=prepare(c,q,guided,dict(backend=backend,seed_batch=17),POSE_PARAMETERS)
    assert seeds==original
    np.testing.assert_array_equal(mask,guided.pocket_mask(c.shape_points,np.asarray([s.transform_matrix for s in original]).reshape(-1,4,4)))
    assert stats['pair_seeds']==512 and stats['stop_reason']=='generated_cap'


@pytest.mark.parametrize('target',[1,3,100,200])
def test_target_returns_first_surviving_prefix(target):
    c,q=fixture();guided=SimpleNamespace(pocket_mask=lambda p,m:m[:,0,3]>0)
    all_seeds,all_mask,_,_=prepare(c,q,guided,dict(max_pair_seeds=1024),POSE_PARAMETERS)
    seeds,mask,stats,_=prepare(c,q,guided,dict(max_pair_seeds=1024,survivor_target=target,seed_batch=31),POSE_PARAMETERS)
    assert seeds==all_seeds[:len(seeds)]
    np.testing.assert_array_equal(mask,all_mask[:len(seeds)])
    assert mask.sum()==target and mask[-1] and stats['stop_reason']=='survivor_target'
    assert len(seeds)==np.flatnonzero(all_mask)[target-1]+1


def test_zero_survivors_terminates_at_cap_and_bounds_query_batches():
    c,q=fixture();widths=[]
    def reject(p,m):widths.append(len(m));return np.zeros(len(m),bool)
    seeds,mask,stats,_=prepare(c,q,SimpleNamespace(pocket_mask=reject),
        dict(max_pair_seeds=128,survivor_target=200,seed_batch=17),POSE_PARAMETERS)
    assert not mask.any() and stats['pair_seeds']==128 and stats['stop_reason']=='generated_cap'
    assert max(widths)<=17 and sum(widths)==len(seeds)


@pytest.mark.parametrize('bad',[dict(max_pair_seeds=-1),dict(survivor_target=True),dict(seed_batch=0),dict(backend='other'),dict(unknown=1)])
def test_reject_invalid_policy(bad):
    with pytest.raises(ValueError):policy(bad)


def test_worker_experimental_default_preserves_pose_payload(tmp_path,monkeypatch):
    import json
    from aidd_agent import preselection_full as pre, guided_filters
    from test_preselection_full import geometry
    query,features,_=geometry()
    features.feature_points=query['feature_points'].copy()
    candidate=SimpleNamespace(**vars(features),shape_points=query['shape_points'])
    design=dict(optional_weights={'A':1.,'B':.3},gaussian_weight=.7,optional_weight=.3,minimum_pose_score=.2)
    q=dict(selection_mode='budget',consensus_npz='fixture',guided_design=design,
        condition_policy=dict(required_anchors=['A','B'],minimum_score=.5),
        anchors=[dict(anchor_id='A',score_column=0),dict(anchor_id='B',score_column=1)])
    monkeypatch.setattr(guided_filters,'GuidedFilter',lambda *a:SimpleNamespace(pocket_mask=lambda p,m:m[:,0,3]>=0))
    monkeypatch.setattr(pre,'_GUIDED_CACHE',None)
    monkeypatch.setattr(pre.full,'_STATE',(SimpleNamespace(get=lambda _:candidate),query,query,q))
    monkeypatch.setattr(pre.full,'_FILTER_READER',SimpleNamespace(get=lambda _:features))
    old=tmp_path/'old.npz';pre.compute((0,1,str(old),[7]))
    q['seed_search']=dict(backend='batched')
    new=tmp_path/'new.npz';pre.compute((0,1,str(new),[7]))
    assert old.with_suffix('.poses.jsonl').read_text()==new.with_suffix('.poses.jsonl').read_text()
    a=json.loads(old.with_suffix('.receipt.json').read_text())
    b=json.loads(new.with_suffix('.receipt.json').read_text())
    assert a['counts']==b['counts']
    assert b['seed_search']['conformers'][0]['surviving']==b['counts']['pocket_tested_seeds']-b['counts']['pocket_rejected_seeds']


def test_tail_statistics_deduplicate_conformers_and_expose_sample_sizes():
    from aidd_agent.block_tail_audit import summarize
    rows=[dict(block_id=b,molecule_id=m,score=s) for b,m,s in
          [('a','x',1),('a','x',3),('a','y',2),('a','z',0),('b','q',8)]]
    report=summarize(rows,2,2)
    a,b=report['blocks']
    assert a['sampled_molecules']==3 and a['maximum']==3 and a['top_k_mean']==2.5
    assert a['high_score_fraction']==2/3 and b['top_k_used']==1
    assert not report['comparable_sample_sizes']
    with pytest.raises(ValueError):summarize([dict(block_id='a',molecule_id='x',score=float('nan'))])


def test_rank_comparison_reports_missing_and_improved_candidates():
    from aidd_agent.seed_budget_audit import compare
    a=(['a','b'],{('a','t'):(.5,.2),('b','t'):(.3,.1)},{('a','t'):'a'})
    b=(['b','c'],{('b','t'):(.4,.1),('c','t'):(.1,.2)},{('b','t'):'b'})
    row=compare(a,b,2)
    assert row['reference_top_recall']==.5 and row['improved_contact_pairs']==1
    assert row['missing_molecule_templates']==1 and row['new_molecule_templates']==1
    assert not row['exact_pose_payloads']


def test_audit_freezes_one_molecule_panel_and_never_modifies_source(tmp_path,monkeypatch):
    import json
    import sqlite3
    from aidd_agent import seed_budget_audit as audit
    from aidd_agent.expanded_wee1 import fingerprint
    from aidd_agent.gaussian_batch import _atomic_json
    batch=tmp_path/'batch';(batch/'artifacts').mkdir(parents=True)
    catalog=batch/'artifacts/catalog.json';_atomic_json(catalog,{})
    source=tmp_path/'source';(source/'adopted-design').mkdir(parents=True)
    np.save(source/'selected-molecules.npy',np.asarray(['m1','m2','m3'],dtype='S16'))
    _atomic_json(source/'protocol.json',dict(sources=fingerprint([catalog]),template_quota=20,rrf_k=60))
    _atomic_json(source/'retrieval.json',dict(outputs=fingerprint([source/'selected-molecules.npy'])))
    _atomic_json(source/'adopted-design/report.json',dict(sources=fingerprint([catalog]),templates=[]))
    before=fingerprint(list(source.rglob('*.*')));panels=[];calls=[]
    def collect(b,keys):
        panels.append(keys.copy());return np.array([3,7,9])
    def refine(b,out,design,ids,workers,chunk,settings):
        calls.append((ids.copy(),dict(settings)))
        root=out/'chunks/00';root.mkdir(parents=True)
        _atomic_json(root/'0.receipt.json',dict(files=fingerprint([catalog]),worker_seconds={'seed_generation':1},counts={},
            seed_search=dict(conformers=[dict(global_id=3,molecule_id='m1',generated=5,pair_seeds=0,surviving=2,stop_reason='exhausted_pairs')])) )
        return dict(wall_seconds=1)
    def merge(out,*args):
        with sqlite3.connect(out/'ranking.sqlite') as db:
            db.executescript('CREATE TABLE ranking(rank INTEGER,mid TEXT); INSERT INTO ranking VALUES(1,"m1"); CREATE TABLE poses(mid TEXT,template TEXT,contact REAL,gaussian REAL,payload TEXT); INSERT INTO poses VALUES("m1","t",0.2,0.3,"same");')
        return 1
    monkeypatch.setattr(audit,'collect_conformers',collect)
    monkeypatch.setattr(audit,'run_refinement',refine);monkeypatch.setattr(audit,'merge_and_rank',merge)
    args=SimpleNamespace(run=source,output=tmp_path/'audit',molecules=2,seed=5,workers=2,chunk_conformers=32,top=1)
    report=audit.run(args)
    assert report['equivalence_passed'] and len(calls)==8 and len(panels)==1
    assert all(np.array_equal(ids,[3,7,9]) for ids,_ in calls)
    assert calls[-1][1]['survivor_target']==200
    assert before==fingerprint(list(source.rglob('*.*')))
    assert json.loads((args.output/'status.json').read_text())['status']=='complete'
    with pytest.raises(ValueError,match='fresh'):audit.run(args)
    _atomic_json(catalog,dict(changed=True))
    args.output=tmp_path/'changed-audit'
    with pytest.raises(ValueError,match='Upstream evidence changed'):audit.run(args)
    assert not args.output.exists()
