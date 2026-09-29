import copy
import json

import numpy as np
import pytest

from aidd_agent import receptor_advice as ra
from aidd_agent import pocket_states as ps


def rows(n, entries=None):
    return [dict(id=f'{i}:A', pdb_id=str(entries[i] if entries else i),
                 target_chain='A', transform=np.eye(4).tolist(), resolution=2.,
                 queries=[dict(query_id=f'{i}:LIG:A:1', smiles='c1ccccc1' if i%2 else 'C1CCC1')])
            for i in range(n)]


def test_sparse_single_does_not_claim_rigidity():
    result = ra.recommend(rows(2), np.array([[0,.8],[.8,0]]))
    assert result['mode'] == 'single_insufficient_evidence'
    assert result['recommended_k'] == 1
    assert result['coverage']['uncovered']
    assert not result['coverage']['target_met']
    assert 'not proof of rigidity' in result['limitations'][0]


def test_near_repeats_can_support_unresolved_variation():
    result = ra.recommend(rows(8,[0,0,1,1,2,2,3,3]), np.full((8,8),.04)-np.eye(8)*.04)
    assert result['mode'] == 'single_unresolved_variation'
    assert result['repeat_baseline']['entries'] == 4
    assert result['coverage']['target_met']


def test_supported_partition_survives_entry_removal():
    d = np.full((8,8),.8)
    d[:4,:4] = d[4:,4:] = .03
    np.fill_diagonal(d,0)
    result = ra.recommend(rows(8),d)
    assert result['mode'] == 'partition_representatives'
    assert result['recommended_k'] == 2
    assert any(x.get('leave_entry_out_stable') for x in result['plateau_diagnostics'])
    limited = ra.recommend(rows(8),d,dict(maximum_representatives=1))
    assert limited['recommended_k'] == 1
    assert 'exceeds' in limited['recommendation']


def test_continuous_coverage_budget_and_rare_retention():
    x = np.arange(10)/10
    d = abs(x[:,None]-x)
    result = ra.recommend(rows(10),d,dict(maximum_representatives=2))
    assert result['mode'] == 'coverage_representatives'
    assert result['recommended_k'] == 2
    assert not result['coverage']['target_met']
    assert result['coverage']['uncovered']
    assert all(p['role'] == 'local_coverage_reference' for p in result['representatives'])
    assert 'not an adequate-coverage claim' in ra.summary(result)


def test_chain_multiplicity_and_unclassified_ligands():
    r = rows(20, [0]*20)
    result = ra.recommend(r,np.zeros((20,20)))
    assert result['evidence']['independent_pdbs'] == 1
    assert result['repeat_baseline']['entries'] == 1
    assert not result['repeat_baseline']['sufficient']
    for row in r:
        row['queries'][0]['smiles'] = ''
    result = ra.recommend(r,np.zeros((20,20)))
    assert result['evidence']['independent_series_proxy'] == 0
    assert result['evidence']['unclassified_ligand_instances'] == 20


def test_fusion_preserves_zero_and_records_scales():
    overlap = np.array([[1,.7],[.7,1]])
    local = np.array([[0,.4],[.4,0]])
    d, calibration = ra.fuse(overlap,local,np.zeros((2,2)))
    assert d[0,0] == 0 and d[0,1] == d[1,0]
    assert calibration['shape_weights'] == [.5,.5]
    assert np.allclose(calibration['scales'],[.3,.4])


def test_advice_default_and_individual_adoption(tmp_path):
    r = rows(8)
    d = np.abs(np.arange(8)[:,None]-np.arange(8))/10
    original = dict(kind='pocket_states',descriptor_version=2,readiness='needs_user_adoption',
        sources={},policy=ps.policy(),structures=r,clusters=[],sensitivity=[],limitations=[],
        diversity_report='diversity.json',pairwise=dict(distance=d.tolist(),overlap=(1-d).tolist(),
                    local_difference=d.tolist(),chemical_distance=np.zeros((8,8)).tolist()))
    original_copy = copy.deepcopy(original)
    report = ra.attach(original)
    evidence = tmp_path/'source.txt'; evidence.write_text('immutable evidence')
    report['sources'] = ps.fingerprint([evidence])
    assert report['pairwise']['legacy_distance'] == original_copy['pairwise']['distance']
    assert len(report['receptor_options']) == 8
    path=tmp_path/'report.json';path.write_text(json.dumps(report))
    adopted = ps.adopt(path,tmp_path/'default')
    assert adopted['selected_cluster_ids'] == [c['id'] for c in report['clusters']]
    choice = report['receptor_options'][-1]
    ps.adopt(path,tmp_path/'individual',[choice['id']])
    _,cohort,_ = ps.state_cohort(tmp_path/'individual/report.json')
    assert [q['query_id'] for q in cohort['admitted']] == ['7:LIG:A:1']
    assert cohort['receptor_selection_role'] == 'user_selected_single_structure'
    with pytest.raises(ValueError):
        ps.adopt(path,tmp_path/'invalid',['pocket-ffffffffffff'])


def test_chat_shows_advice_instead_of_json():
    from aidd_agent.chat_agent import screening_feedback
    a = ra.recommend(rows(2), np.array([[0,.4],[.4,0]]))
    message = screening_feedback(dict(kind='pocket_states',receptor_advice=a,clusters=[{'sentinel':True}]))
    assert 'Evidence: 2 qualified chains' in message
    assert 'sentinel' not in message
    assert 'Literature has not yet been reviewed' in message


def test_receptor_budget_is_validated_in_agent_plan():
    from aidd_agent.prompt_plan import validate_plan
    plan = dict(version=1,summary='Review receptors',clarifications=[],steps=[dict(
        id='pockets',action='pocket_states',params=dict(source_run='PROMPT-'+'a'*16,
        maximum_representatives=8,coverage_fraction=.9))])
    assert validate_plan(plan)
    plan['steps'][0]['params']['coverage_fraction'] = 0
    with pytest.raises(ValueError):
        validate_plan(plan)


def test_viewer_follows_individual_adoption(tmp_path):
    from aidd_agent.structure_review import structures
    files = [tmp_path/f'{i}.cif' for i in range(2)]
    for file in files:
        file.write_text('structure fixture')
    r = rows(2)
    for row,file in zip(r,files):
        row['structure_path'] = str(file)
    pocket = tmp_path/'pockets.json'
    pocket.write_text(json.dumps(dict(kind='pocket_states',structures=r,
        sources=ps.fingerprint(files),clusters=[dict(id='pocket-000000000001',representative='0:A')],
        receptor_options=[dict(id='pocket-000000000002',representative='1:A')])) )
    adoption = tmp_path/'adoption.json'
    adoption.write_text(json.dumps(dict(kind='pocket_adoption',pocket_report=str(pocket),
        sources=ps.fingerprint([pocket,*files]), selected_cluster_ids=['pocket-000000000002'])))
    execution=tmp_path/'execution.json'
    execution.write_text(json.dumps(dict(steps=dict(adopt=dict(action='pocket_adopt',status='complete',result=dict(report=str(adoption)))))))
    catalog=structures(dict(status='complete',report=str(execution),plan=str(tmp_path/'plan.json')),tmp_path,collection='diverse')
    assert len(catalog)==1
    assert str(files[1]) in json.dumps(catalog).replace('\\\\','\\')


@pytest.mark.parametrize('d',[np.array([[0,-1],[-1,0]]),np.array([[0,.2],[.3,0]]),np.array([[0,np.nan],[np.nan,0]])])
def test_bad_matrices_rejected(d):
    with pytest.raises(ValueError):
        ra.recommend(rows(2),d)
