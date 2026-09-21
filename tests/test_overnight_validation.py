import json
from pathlib import Path

import pytest

from aidd_agent import overnight_validation as suite
from aidd_agent.gaussian_batch import _atomic_json


@pytest.mark.parametrize('coverage',[True,False])
def test_suite_requires_full_coverage_and_never_caps_scan(tmp_path,monkeypatch,coverage):
    source=tmp_path/'source'; source.mkdir()
    selection=source/'report.json'
    evidence=source/'evidence.json'; evidence.write_text('{}')
    policy=dict(required_anchors=['a'],match_mode='all',minimum_score=.5,
        coarse_constraints=dict(heavy_atom_ratio=[.7,1.3],maximum_extent_distance=.35,minimum_feature_coverage=.7))
    _atomic_json(selection,dict(kind='selection_preview',policy=policy,sources=suite.full.fingerprint([evidence]),evidence_report=str(evidence),
        query=dict(artifact_catalog=str(source/'catalog.json'),chemical_companion=str(source/'chemical.json'))))
    monkeypatch.setattr(suite,'pilot',lambda *a,**k:dict(status='complete'))
    def run(*args,**kwargs):
        assert kwargs==dict(workers=22,chunk_size=2048,backend='numpy',pose_feasibility=True)
        return dict(status='complete' if coverage else 'partial',library={'library_conformers':100},
            queries=[dict(query_id='q',full_coverage=coverage,counts=dict(total_conformers=100,matching_molecules=0))])
    monkeypatch.setattr(suite.full,'run_funnel',run)
    def preview(source,output,rule):
        assert 'max_molecules' not in rule
        output.mkdir(); path=output/'representatives.jsonl'; path.write_text('')
        return dict(representatives_jsonl=str(path))
    monkeypatch.setattr(suite.full,'preview',preview)
    output=tmp_path/'out'
    if coverage:
        result=suite.run(selection,output)
        assert result['status']=='complete' and result['candidate_top_k'] is None
        assert result['final_hit_validation']['status']=='not_covered_no_final_hits'
    else:
        with pytest.raises(ValueError,match='coverage'): suite.run(selection,output)
        result=json.loads((output/'suite-report.json').read_text())
        assert result['status']=='failed' and 'all_matches' not in result['stages']


def test_reservoir_samples_across_result_list_reproducibly(tmp_path):
    path=tmp_path/'rows.jsonl'
    path.write_text(''.join(json.dumps(dict(global_id=i))+'\n' for i in range(1000)))
    a=suite.sample_representatives(path)
    assert len(a)==128 and len(set(a))==128 and max(a)>900
    assert a==suite.sample_representatives(path)
