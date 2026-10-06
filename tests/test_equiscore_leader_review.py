"""Saved-result review validation, using sealed synthetic ranking tables only."""
import importlib.util
import json
from pathlib import Path
import pytest

spec=importlib.util.spec_from_file_location('leader_review',Path(__file__).resolve().parents[1]/'scripts/review_equiscore_leaders.py')
module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)


def fixture(tmp_path):
    analysis=tmp_path/'analysis';analysis.mkdir()
    source=tmp_path/'full.json'
    source.write_text(json.dumps(dict(scope='full',status='complete',failed_pairs=0,pairs=200,scored=200,worker_exit_code=0)))
    rankings=[];candidates=[];comparisons=[]
    for receptor in ('A','B'):
        for rank in range(1,13):
            block='block-'+str(rank);score=1-rank/100
            for n in (5,10):
                rankings.append(dict(scheme='E095',receptor=receptor,block_id=block,top_n=n,rank=rank,top_n_mean=score,status='eligible',missing=0,population=1000))
            for i in range(1,11):
                candidates.append(dict(scheme='E095',receptor=receptor,block_id=block,candidate_rank=i,molecule_id='m'+str(i),cid='c'+str(i),score=score,pose_file='/saved/poses.mol2',pose_name='fixture',pose_index=i))
            comparisons.append(dict(scheme='E095',receptor=receptor,block_id=block,equiscore_rank=rank,chemplp_rank=13-rank,rank_improvement=13-2*rank,equiscore_top10_mean=score,chemplp_top10_mean=-80,shared_top10_molecules=3))
    for name,data in [('block_rankings.csv',rankings),('block_top_candidates.csv',candidates),('rank_comparison.csv',comparisons)]:module.write_csv(analysis/name,data)
    report=dict(status='complete',source_report=str(source),ranking=dict(score='EquiScore',better='higher',ranking_unit='molecule',top_n=10,ranked_blocks=24,source_report_sha256=module.sha(source)),output_hashes={p.name:module.sha(p) for p in analysis.iterdir()})
    (analysis/'report.json').write_text(json.dumps(report))
    return analysis


@pytest.mark.parametrize('limit',[5,10])
def test_separate_groups_bounded_candidates_and_no_reexecution(tmp_path,limit):
    analysis=fixture(tmp_path);before={p.name:module.sha(p) for p in analysis.iterdir()}
    output=tmp_path/'review'
    result=module.review(analysis,output,limit)
    assert result['selected_blocks']==2*limit and result['selected_candidate_rows']==20*limit
    assert not result['model_executed'] and not result['pose_geometry_inspected']
    assert all(g['top5_block_overlap_between_top10_and_top5_means']==5 for g in result['groups'])
    assert (tmp_path/'review.zip').exists()
    assert before=={p.name:module.sha(p) for p in analysis.iterdir()}
    with pytest.raises(ValueError,match='new review'):module.review(analysis,output,limit)


def test_tamper_stops_before_output(tmp_path):
    analysis=fixture(tmp_path)
    with (analysis/'rank_comparison.csv').open('a') as f:f.write('tampered')
    with pytest.raises(ValueError,match='artifact'):module.review(analysis,tmp_path/'review')
    assert not (tmp_path/'review').exists()


def test_candidate_identity_must_remain_distinct(tmp_path):
    analysis=fixture(tmp_path);path=analysis/'block_top_candidates.csv'
    data=module.rows(path);data[1]['molecule_id']=data[0]['molecule_id'];module.write_csv(path,data)
    report=json.loads((analysis/'report.json').read_text());report['output_hashes'][path.name]=module.sha(path)
    (analysis/'report.json').write_text(json.dumps(report))
    with pytest.raises(ValueError,match='distinct'):module.review(analysis,tmp_path/'review')


def test_incomplete_full_receipt_rejected_even_with_analysis_complete(tmp_path):
    analysis=fixture(tmp_path);source=tmp_path/'full.json'
    run=json.loads(source.read_text());run['failed_pairs']=1;source.write_text(json.dumps(run))
    path=analysis/'report.json';report=json.loads(path.read_text());report['ranking']['source_report_sha256']=module.sha(source);path.write_text(json.dumps(report))
    with pytest.raises(ValueError,match='successful full'):module.review(analysis,tmp_path/'review')
