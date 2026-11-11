"""Original-output retrieval is read-only and supports official filename variants."""
import json
from pathlib import Path

import pytest

from aidd_agent.af3_batch import sample_summary,validate_sample_coverage
from aidd_agent.af3_analysis import assess_chain_mapping
from aidd_agent.library_acceptance import sha


def files(root,prefixed=True):
    root.mkdir(parents=True,exist_ok=True)
    for seed in (47,48):
        name=f'seed-{seed}_sample-0';folder=root/name;folder.mkdir()
        prefix='job_'+name+'_' if prefixed else ''
        (folder/(prefix+'summary_confidences.json')).write_text(json.dumps(dict(iptm=.84)))
        (folder/(prefix+'model.cif')).write_text('data_fixture\n#\n')
    (root/'job_ranking_scores.csv').write_text('seed,sample,ranking_score\n47,0,0.8\n48,0,0.9\n')


@pytest.mark.parametrize('prefixed',[True,False])
def test_sample_filename_variants_and_coverage(tmp_path,prefixed):
    files(tmp_path/'prediction',prefixed)
    result=sample_summary(tmp_path/'prediction')
    assert result['sample_count']==2 and all(s['model'] for s in result['samples'])
    assert result['sample_model_count']==2 and result['missing_sample_models']==[]
    validate_sample_coverage(result,dict(modelSeeds=[47,48]))
    assert result['seed_coverage']['status']=='verified_ranked_samples'
    assert result['best_recorded_sample']['seed']=='48'


def test_ambiguous_sample_summary_is_rejected(tmp_path):
    files(tmp_path)
    (tmp_path/'seed-47_sample-0/summary_confidences.json').write_text('{}')
    with pytest.raises(ValueError,match='Ambiguous'):sample_summary(tmp_path)


def test_consistent_metrics_do_not_claim_dimensions_disagree():
    value=dict(chain_ids=['A']*905+['B']*26,chain_ptm=[.62,.3],chain_iptm=[.84,.84],
               chain_pair_iptm=[[.62,.84],[.84,.3]])
    result=assess_chain_mapping(value)
    assert result['metric_dimensions_consistent'] and result['status']=='unverified'
    assert not any('dimensions disagree with' in w for w in result['warnings'])
    assert result['matrix_chain_ids'] is None


def test_original_output_review_no_mutation_no_task_and_isolation(tmp_path):
    from aidd_agent.chat_agent import ChatAgent
    from aidd_agent.domain_tools import DomainTools
    from aidd_agent.research_projects import handle
    app=ChatAgent(tmp_path/'workspace',start=False)
    try:
        sid=app.new_session();root=Path(app.context_for(sid)['project_root'])
        output=root/'prediction';files(output)
        final=output/'job_model.cif';final.write_text('data_fixture\n#\n')
        inp=root/'af3-input.json';inp.write_text(json.dumps(dict(modelSeeds=[47,48])))
        report=root/'old-report.json'
        report.write_text(json.dumps(dict(status='complete',steps=dict(
            prepare=dict(status='complete',action='af3_prepare',result=dict(input_json=str(inp),input_sha256=sha(inp))),
            predict=dict(status='complete',action='af3_run',result=dict(model=str(final),structure_sha256=sha(final),output_directory=str(output)))))))
        with app.connect() as db:
            db.execute('INSERT INTO jobs(id,session,status,report,request,created) VALUES(?,?,?,?,?,?)',
                       ('old-task',sid,'complete',str(report),'Old prediction fixture',1))
        before={str(p):sha(p) for p in root.rglob('*') if p.is_file()}
        tools=DomainTools(app,sid,'Review old-task without running predictions','deepseek')
        result=tools.call('af3_review',dict(task_id='old-task',limit=1))
        assert result['sample_count']==2 and result['best_recorded_sample']['seed']=='48'
        assert result['best_recorded_sample_model'].endswith('_model.cif')
        assert result['final_model_check']=='hash_verified'
        assert result['seed_coverage']['status']=='verified_ranked_samples'
        assert result['samples']['total']==2
        assert len(app.jobs(sid))==1
        assert before=={str(p):sha(p) for p in root.rglob('*') if p.is_file()}
        other=handle(app,sid,dict(operation='create',name='Other',target='OTHER',objective='Isolation'))['session']
        with pytest.raises(ValueError):DomainTools(app,other,'Read old-task','deepseek').call('af3_review',dict(task_id='old-task'))
        final.write_text('changed')
        with pytest.raises(ValueError,match='CIF changed'):tools.call('af3_review',dict(task_id='old-task'))
    finally:app.close()
