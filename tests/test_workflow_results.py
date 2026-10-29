import json
import threading
from urllib.error import HTTPError
from urllib.request import Request, urlopen

import pytest

from aidd_agent.chat_agent import ChatAgent
from aidd_agent.chat_web import make_server
from aidd_agent.domain_tools import DomainTools
from aidd_agent.workflow_results import handle, summarize


def kras_report():
    report = dict(status='complete', scope='existing_pool_grid_docking_handoff', targets={}, exports={})
    for name in ('KRAS_G12D', 'KRAS_G12V'):
        report['targets'][name] = dict(pool=100, processed=100, strict=2, balanced=20, lenient=80, unresolved=20)
        report['exports'][name] = dict(status='complete', target=name, requested_cap=15,
            selected=15, exported=13, quarantined=2, folder='/remote/' + name,
            batches=[dict(path='/remote/' + name + '/batch_001.mol2', molecules=13, sha256='a'*64)])
    return report


@pytest.fixture
def setup(tmp_path):
    imports = tmp_path / 'imports'; imports.mkdir()
    app = ChatAgent(tmp_path / 'storage', start=False, result_import_roots=[imports])
    path = imports / 'report.json'; path.write_text(json.dumps(kras_report()))
    yield app, app.new_session(), path
    app.close()


def test_import_snapshot_dedup_freshness_and_no_compute(setup):
    app, sid, path = setup
    command = '/workflow_results ' + json.dumps(dict(operation='attach', kind='kras_export', report=str(path)))
    response = json.loads(app.ask(sid, command, 'deepseek'))
    assert response['computation_launched'] is False
    assert response['payload_verification'] == 'not_performed'
    assert not app.jobs(sid)
    repeated = json.loads(app.ask(sid, command, 'deepseek'))
    assert repeated['attachment_id'] == response['attachment_id']
    rows = app.snapshot(sid)['workflow_results']
    assert len(rows) == 1 and 'report_snapshot' not in rows[0]
    assert rows[0]['summary']['targets']['KRAS_G12D']['exported'] == 13
    args = dict(operation='check', attachment_id=response['attachment_id'])
    assert handle(app, sid, args)['source_status'] == 'unchanged'
    path.write_text(json.dumps(kras_report(), indent=2))
    assert handle(app, sid, args)['source_status'] == 'changed_reimport_to_update'
    path.unlink()
    assert handle(app, sid, args)['source_status'] == 'unavailable_or_invalid'
    content = handle(app, sid, dict(operation='read', attachment_id=response['attachment_id'], pointer='/exports/KRAS_G12D/exported'))
    assert content['content']['value'] == 13


def test_sessions_paths_and_model_invented_paths_rejected(setup, tmp_path):
    app, sid, path = setup
    result = handle(app, sid, dict(operation='attach', kind='kras_export', report=str(path)))
    with pytest.raises(ValueError, match='session'):
        handle(app, app.new_session(), dict(operation='read', attachment_id=result['attachment_id']))
    foreign = tmp_path / 'foreign.json'; foreign.write_text(path.read_text())
    with pytest.raises(ValueError, match='outside'):
        handle(app, sid, dict(operation='attach', kind='kras_export', report=str(foreign)))
    with pytest.raises(ValueError, match='outside'):
        handle(app, sid, dict(operation='attach', kind='kras_export', report=str(path.parent / '..' / 'foreign.json')))
    tools = DomainTools(app, sid, 'Show my results', 'deepseek')
    with pytest.raises(ValueError, match='exact report path'):
        tools.call('workflow_results', dict(operation='attach', kind='kras_export', report=str(path)))
    assert len(tools.call('workflow_results', {})['attachments']) == 1


@pytest.mark.parametrize('mutation', ['sum', 'pool', 'quota', 'hash', 'duplicate', 'boolean', 'scope'])
def test_inconsistent_kras_reports_rejected(mutation):
    report = kras_report(); target = report['targets']['KRAS_G12D']; export = report['exports']['KRAS_G12D']
    if mutation == 'sum': export['exported'] += 1
    if mutation == 'pool': target['processed'] -= 1
    if mutation == 'quota': export['requested_cap'] = 9
    if mutation == 'hash': export['batches'][0]['sha256'] = 'invented'
    if mutation == 'duplicate': export['batches'] *= 2
    if mutation == 'boolean': target['strict'] = True
    if mutation == 'scope': report['scope'] = 'sampled_pool_grid_pilot'
    with pytest.raises(ValueError): summarize('kras_export', report)


def af3_report():
    return dict(mode='prediction', status='complete', jobs=[dict(id='ligand1', status='model_generated',
        input_sha256='b'*64, model='/remote/ligand1_model.cif', model_sha256='c'*64, confidence=dict(iptm=.7))])


def test_af3_host_checks_are_not_models_and_failures_preserved():
    report = af3_report()
    assert summarize('af3_batch', report)['generated_models'] == 1
    checked = dict(mode='host_check', status='host_checks_passed_not_inference',
        jobs=[dict(id='ligand1', status='host_checks_passed', input_sha256='b'*64)])
    assert summarize('af3_batch', checked)['generated_models'] == 0
    checked['status'] = 'complete'
    with pytest.raises(ValueError): summarize('af3_batch', checked)
    report['jobs'][0]['status'] = 'failed'
    with pytest.raises(ValueError): summarize('af3_batch', report)
    report['status'] = 'failed'
    assert summarize('af3_batch', report)['failed_jobs'] == 1
    report = af3_report(); report['jobs'][0].pop('model_sha256')
    with pytest.raises(ValueError): summarize('af3_batch', report)


def test_af3_import_and_paged_confidence(setup):
    app, sid, path = setup
    path.write_text(json.dumps(af3_report()))
    result = handle(app, sid, dict(operation='attach', kind='af3_batch', report=str(path)))
    item = handle(app, sid, dict(operation='read', attachment_id=result['attachment_id'], pointer='/jobs/0/confidence'))
    assert item['content']['value'] == {'iptm': .7}
    assert not app.jobs(sid)


def test_af3_batch_totals_preserve_pending_jobs():
    report = af3_report(); report.update(status='running', total_jobs=3)
    assert summarize('af3_batch', report)['pending_jobs'] == 2
    report['status'] = 'complete'
    with pytest.raises(ValueError, match='total_jobs'): summarize('af3_batch', report)
    report['total_jobs'] = 1
    assert summarize('af3_batch', report)['pending_jobs'] == 0


def test_http_result_attachment_preserves_auth_and_state(setup):
    app, sid, path = setup
    server = make_server(app, port=0, token='fixture-token')
    thread = threading.Thread(target=server.serve_forever, daemon=True); thread.start()
    url = 'http://127.0.0.1:' + str(server.server_port) + '/api/workflow-results'
    body = json.dumps(dict(session=sid, arguments=dict(operation='attach', kind='kras_export', report=str(path)))).encode()
    try:
        with pytest.raises(HTTPError) as error: urlopen(Request(url, data=body))
        assert error.value.code == 401
        with urlopen(Request(url, data=body, headers={'Authorization': 'Bearer fixture-token'})) as response:
            assert json.load(response)['status'] == 'attached'
        assert app.snapshot(sid)['workflow_results'][0]['summary']['status'] == 'complete'
        assert not app.jobs(sid)
    finally:
        server.shutdown(); server.server_close(); thread.join()


def test_af3_prediction_manifest_keeps_provenance():
    manifest = dict(backend='alphafold3', construct_name='fixture', structure_path='/remote/model.cif',
                    structure_sha256='d'*64, chain_ids=['A', 'B'], confidence={'iptm': .6})
    summary = summarize('af3_prediction', manifest)
    assert summary['chain_ids'] == ['A', 'B'] and summary['status'] == 'model_reported'
