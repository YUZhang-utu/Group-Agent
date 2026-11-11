"""Small real geometry/export checks; no live provider or production library."""
import hashlib
import io
import json
from pathlib import Path
import sqlite3

import numpy as np
import pytest

from aidd_agent import ligand_library as search
from aidd_agent.chat_agent import ChatAgent
from aidd_agent.chat_uploads import upload
from aidd_agent.prompt_workflow import run_plan
from aidd_agent.prompt_smoke import OfflineServices
from test_domain_agent import scripted
from test_mol2_rdkit import _block


def tiny_library(root):
    from aidd_agent.conformer_artifacts import build_shard
    from aidd_agent.chemical_companion import build_chemical_companion_catalog
    from aidd_agent.mol2 import iter_mol2_blocks
    batch = root/'batch'; batch.mkdir()
    raw = _block('1 C1 -1.4 0 0 C.3 1 M 0\n2 C2 0 0 0 C.2 1 M 0\n3 O 0 1.2 0 O.2 1 M 0\n4 C3 1.4 -.4 .2 C.3 1 M 0\n',
                 '1 1 2 1\n2 2 3 2\n3 2 4 1\n', 4, 3)
    source = root/'molecules.mol2'
    source.write_text(raw + raw.replace('\ntest\n', '\nsecond\n'), newline='')
    with sqlite3.connect(batch/'registry.sqlite3') as db:
        db.execute('CREATE TABLE molecule(id,library_id,source_name)')
        db.execute('CREATE TABLE conformer(id,molecule_id,source_path,source_record_index,content_sha256,source_record_name)')
        for i, (_, text) in enumerate(iter_mol2_blocks(source)):
            db.execute('INSERT INTO molecule VALUES(?,?,?)', (f'M{i}', 'LIB', 'original '+str(i)))
            db.execute('INSERT INTO conformer VALUES(?,?,?,?,?,?)',
                       (f'C{i}', f'M{i}', str(source.resolve()), i, hashlib.sha256(text.encode()).hexdigest(), 'record '+str(i)))
    artifacts = batch/'artifacts'
    manifest = build_shard(batch/'registry.sqlite3', 'LIB', source, artifacts, 0)
    catalog = dict(format='aidd-conformer-artifact-catalog', version=1, library_id='LIB', conformers=2,
        shards=[dict(name=source.stem, path=str(artifacts/source.stem), global_id_start=0, conformers=2,
                     manifest_sha256=hashlib.sha256((artifacts/source.stem/'manifest.json').read_bytes()).hexdigest())])
    (artifacts/'catalog.json').write_text(json.dumps(catalog))
    build_chemical_companion_catalog(None, 'LIB', artifacts/'catalog.json', batch/'chemical', workers=1)
    return batch, raw


@pytest.mark.parametrize('source_mode', ['upload', 'prepared_cocrystal_query'])
def test_uploaded_query_full_scan_and_real_mol2_export(tmp_path, source_mode):
    batch, raw = tiny_library(tmp_path)
    app = ChatAgent(tmp_path/'chat', start=False)
    try:
        sid = app.new_session()
        attachment = upload(app, sid, 'query.mol2', io.BytesIO(raw.encode()), len(raw.encode()))
        ctx = app.context
        if source_mode == 'upload':
            args = dict(attachment_id=attachment['id'])
        else:
            from aidd_agent.prompt_workflow import project_root
            from aidd_agent.expanded_wee1 import fingerprint
            project = project_root(Path(ctx['db']), ctx['user_id'], ctx['project_id'])
            query = project/'native.npz'; search.prepare_ligand(attachment['path'], query, 'TEST:LIG:A:1')
            child = project/'consensus.json'
            child.write_text(json.dumps(dict(prepared_complexes=[dict(query_id='TEST:LIG:A:1', query_npz=str(query))],
                sources=fingerprint([query, query.with_suffix('.manifest.json')]))))
            report_path = project/'prior.json'
            report_path.write_text(json.dumps(dict(steps=dict(prepare=dict(action='structure_consensus', status='complete', result=dict(report=str(child)))))))
            with app.connect() as db:
                db.execute('INSERT INTO jobs(id,session,status,report,request,created) VALUES(?,?,?,?,?,?)',
                           ('prior-task', sid, 'complete', str(report_path), 'Prepared query fixture, not a live cocrystal discovery', 1))
            found = search.handle(app, sid, dict(operation='inspect', task_id='prior-task'))
            assert found['queries'][0]['query_id'] == 'TEST:LIG:A:1'
            args = dict(task_id='prior-task', query_id='TEST:LIG:A:1')
        receipt = search.handle(app, sid, dict(operation='start', export_molecules=1, **args))
        runtime = tmp_path/'runtime.json'
        runtime.write_text(json.dumps(dict(search=dict(batch=str(batch), workers=1, refine_chunk=1))))
        job = app.task(sid, receipt['task_id'])
        if source_mode == 'upload':
            app.runtime = runtime; app.allow_compute = True
            app.execute(job)
            assert app.task(sid, job['id'])['status'] == 'complete', app.task(sid, job['id']).get('error')
            report = json.loads(Path(job['report']).read_text())
        else:
            report, _ = run_plan(Path(ctx['db']), ctx['user_id'], ctx['project_id'], Path(job['plan']),
                runtime=runtime, services=OfflineServices(), allow_compute=True)
        assert report['status'] == 'complete', report
        result = report['steps']['search_export']['result']
        assert result['full_coverage'] and result['processed_conformers'] == 2
        assert result['ranked_molecules'] == 2 and result['exported_molecules'] == 1
        original = list(Path(result['original_mol2_directory']).glob('*.mol2'))
        posed = list(Path(result['posed_mol2_directory']).glob('*.mol2'))
        assert len(original) == len(posed) == 1 and original[0].read_text() == raw
        assert 'original 0' in Path(result['names']).read_text(encoding='utf-8-sig')
        assert Path(result['ranking']).exists()
        app.update(job['id'], status='complete')
        reply = app.ask(sid, '/results '+job['id'], 'deepseek')
        assert result['original_mol2_directory'] in reply and 'processed_conformers: 2' in reply
        with sqlite3.connect(result['ranking']) as db:
            assert db.execute('SELECT COUNT(*) FROM ranking').fetchone()[0] == 2
        # Direct adapter resume verifies/reuses chunks and the sealed export.
        step = json.loads(Path(job['plan']).read_text())['plan']['steps'][0]
        execution = Path(job['plan']).parent/'execution'
        second = search.execute(step['params'], execution/'search_export', execution,
                                json.loads(runtime.read_text()), True)
        assert second['stage']['chunks_reused'] == 2 and second['exported_molecules'] == 1
        if source_mode == 'prepared_cocrystal_query':
            query.write_bytes(b'changed')
            with pytest.raises(ValueError):
                search.handle(app, sid, dict(operation='start', export_molecules=1, **args))
        else:
            chemistry_path = batch/'chemical/catalog.json'
            chemistry = json.loads(chemistry_path.read_text()); chemistry['library_id'] = 'WRONG'
            chemistry_path.write_text(json.dumps(chemistry))
            with pytest.raises(ValueError, match='different libraries'):
                search.execute(step['params'], execution/'search_export', execution,
                               json.loads(runtime.read_text()), True)
    finally: app.close()


def test_ranking_deduplicates_molecules_with_stable_pose_ties(tmp_path):
    stage = tmp_path/'scores'; (stage/'chunks').mkdir(parents=True)
    arrays = dict(global_ids=np.arange(3), molecule_ids=np.array(['M0','M0','M1']))
    arrays[search.OBJECTIVE+'__objective'] = np.array([.5,.8,.8])
    arrays[search.OBJECTIVE+'__transform'] = np.tile(np.eye(4).reshape(1,16), (3,1))
    np.savez(stage/'chunks/chunk-000000.npz', **arrays)
    ranking = tmp_path/'ranking.sqlite'
    assert search.rank_chunks(stage, 3, 3, ranking, 'query') == 2
    with sqlite3.connect(ranking) as db:
        assert db.execute('SELECT mid FROM ranking ORDER BY rank').fetchall() == [('M0',), ('M1',)]
        assert json.loads(db.execute("SELECT payload FROM poses WHERE mid='M0'").fetchone()[0])['global_id'] == 1


def test_attachment_ownership_and_compute_gate(tmp_path):
    from aidd_agent.workflow_errors import Blocked
    app = ChatAgent(tmp_path/'chat', start=False)
    try:
        from test_research_projects import create
        sid = app.new_session(); other = create(app,sid,'Different target')['session']
        file = upload(app, sid, 'ligand.sdf', io.BytesIO(b'not yet parsed'), 14)
        with pytest.raises(ValueError, match='unavailable'):
            search.handle(app, other, dict(operation='start', attachment_id=file['id'], export_molecules=2))
        receipt = search.handle(app, sid, dict(operation='start', attachment_id=file['id'], export_molecules=2))
        job = app.task(sid, receipt['task_id']); envelope = json.loads(Path(job['plan']).read_text())
        with pytest.raises(Blocked):
            search.execute(envelope['plan']['steps'][0]['params'], tmp_path, Path(job['plan']).parent/'execution', {}, False)
    finally: app.close()


def test_uploaded_two_dimensional_query_is_not_silently_embedded(tmp_path):
    from rdkit import Chem
    from rdkit.Chem import AllChem
    mol = Chem.MolFromSmiles('CCCO'); AllChem.Compute2DCoords(mol)
    source = tmp_path/'flat.sdf'
    with Chem.SDWriter(str(source)) as writer: writer.write(mol)
    with pytest.raises(ValueError, match='3D'):
        search.prepare_ligand(source, tmp_path/'query.npz', 'flat')


def test_chat_submission_returns_receipt_without_second_model_call(tmp_path):
    app = ChatAgent(tmp_path/'chat', start=False)
    try:
        sid = app.new_session()
        file = upload(app, sid, 'query.sdf', io.BytesIO(b'parsed by worker'), 16)
        app.domain_planner = scripted([dict(kind='tool', tool='ligand_library', purpose='Search and export',
            arguments=dict(operation='start', attachment_id=file['id'], export_molecules=10))])
        reply = app.ask(sid, 'Search the library with my uploaded ligand and export 10 MOL2 candidates.', 'deepseek')
        assert len(app.jobs(sid)) == 1 and app.task(sid)['id'] in reply
    finally: app.close()
