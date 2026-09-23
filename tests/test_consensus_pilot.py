import json
from types import SimpleNamespace

import numpy as np
import pytest

from aidd_agent.consensus_pilot import reservoir, molecule_panel, wilson
from aidd_agent import preselection_full as pre
from aidd_agent.necessary_conditions import NecessaryConditions
from test_preselection_full import geometry
from test_library_acceptance import library


def test_reservoir_reproducible_without_replacement_and_not_prefix():
    one,total=reservoir(map(str,range(10000)),100,42)
    assert one==reservoir(map(str,range(10000)),100,42)[0]
    assert len(set(one))==100 and total==10000 and max(map(int,one))>9000
    assert set(reservoir(['a','b'],100,42)[0])=={'a','b'}
    with pytest.raises(TimeoutError):reservoir(range(10000),10,42,deadline=0)


def test_sample_includes_all_conformers_across_shards(library):
    batch,_,molecules=library
    panel,total,_,_=molecule_panel(batch,3,42,float('inf'))
    assert len(panel)==3 and total==len(set(molecules))
    for mid,ids in panel.items():assert ids==np.flatnonzero(molecules==mid.encode()).tolist()


def test_zero_hit_interval_is_not_zero_population_claim():
    lo,hi=wilson(0,100)
    assert lo==pytest.approx(0) and .03<hi<.04
    assert wilson(100,100)[1]==pytest.approx(1)


def test_explicit_pilot_ids_reuse_production_compute(tmp_path,monkeypatch):
    query,features,seeds=geometry()
    candidate=SimpleNamespace(**vars(features),shape_points=query['shape_points'])
    policy=dict(required_anchors=['A','B'],minimum_score=.5,match_mode='any',coarse_constraints={})
    q=dict(condition_policy=policy,anchors=[dict(anchor_id='A',score_column=0),dict(anchor_id='B',score_column=1)])
    monkeypatch.setattr(pre.full,'_STATE',(SimpleNamespace(get=lambda gid:candidate),query,query,q))
    monkeypatch.setattr(pre.full,'_FILTER_READER',SimpleNamespace(get=lambda gid:features))
    monkeypatch.setattr(pre.full,'_BOUND',NecessaryConditions(query,[0,1],'any',.5))
    monkeypatch.setattr(pre,'_GUIDED_CACHE',None)
    monkeypatch.setattr(pre,'coarse_stages',lambda *args,**kwargs:4)
    monkeypatch.setattr(pre,'prepare_seeds',lambda *args,**kwargs:(seeds,2))
    monkeypatch.setattr(pre,'_score_ids',lambda *args,**kwargs:{pre.full.OBJECTIVE+'__objective':np.array([.9])})
    path=tmp_path/'pilot.npz'
    result=pre.compute((0,2,str(path),[7,100]))
    assert result['counts']['input_conformers']==2 and result['counts']['matching_conformers']==2
    with np.load(path) as a:assert a['global_ids'].tolist()==[7,100]
    rows=[json.loads(line) for line in path.with_suffix('.poses.jsonl').read_text().splitlines()]
    assert {r['global_id'] for r in rows}=={7,100}


@pytest.mark.parametrize('timeout',[False,True])
def test_orchestrator_all_templates_union_and_budget_no_projection(library,tmp_path,monkeypatch,timeout):
    from aidd_agent import consensus_pilot as pilot, consensus_design, chemical_companion
    from aidd_agent.expanded_wee1 import fingerprint
    batch,_,_=library
    design=dict(minimum_score=.5,coarse_constraints={},mandatory_anchors=[],alternative_groups=[],optional_weights={'A':1})
    adopted=dict(readiness='ready_for_consensus_funnel',design=design,sources={},
        anchors=[dict(anchor_id='A',feature_index=0)],anchor_order=['A'],consensus_npz='unused',
        templates=[dict(query_id='T1'),dict(query_id='T2')])
    monkeypatch.setattr(consensus_design,'adopt',lambda *args:adopted)
    submitted=[];terminated=[]
    class Job:
        def __init__(self,task):self.task=task
        def ready(self):
            if timeout:raise TimeoutError('fixture budget')
            return True
        def get(self):
            queries,ids,target=self.task;root=pilot.Path(target);root.mkdir(parents=True)
            receipts=[]
            for i,q in enumerate(queries):
                path=root/f'template-{i:03d}.npz'
                from aidd_agent.gaussian_batch import _atomic_savez
                mids=['MOL-000000000000']*len(ids)
                _atomic_savez(path,dict(global_ids=np.array(ids),molecule_ids=np.array(mids),stage_levels=np.zeros(len(ids))))
                path.with_suffix('.poses.jsonl').write_text('')
                receipts.append(dict(counts=dict(input_conformers=len(ids)),worker_seconds={},files=fingerprint([path,path.with_suffix('.poses.jsonl')])))
            return dict(root=str(root),receipts=receipts,cpu_seconds=.01)
    class Pool:
        def apply_async(self,fn,args):submitted.append(args[0]);return Job(args[0])
        def close(self):pass
        def join(self):pass
        def terminate(self):terminated.append(True)
    monkeypatch.setattr(pilot.mp,'get_context',lambda _:SimpleNamespace(Pool=lambda _:Pool()))
    monkeypatch.setattr(chemical_companion,'ChemicalCompanionReader',lambda _:None)
    result=pilot.run('fixture',batch,tmp_path/'pilot',count=100,workers=2,chunk_molecules=2)
    assert all(len(t[0])==2 for t in submitted)
    if timeout:
        assert result['status']=='budget_exhausted' and terminated
        assert 'projected_matching_molecules' not in result
    else:
        assert result['status']=='complete' and result['sample_complete']
        assert result['matching_molecules']==0 and result['approximate_95_percent_matching_interval']==[0,0]
