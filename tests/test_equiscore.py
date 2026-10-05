"""Engineering fixtures only; no EquiScore GPU accuracy claim."""
import csv
import importlib.util
import json
from pathlib import Path
import sqlite3
from types import SimpleNamespace

import pytest

from aidd_agent import equiscore as es
from aidd_agent.block_ranking import inspect
from aidd_agent.final_work_blocks import sha
from aidd_agent.joint_spatial_profiles import save


def seal(root, **values):
    save(root / 'report.json', dict(status='complete', **values,
         output_hashes={p.relative_to(root).as_posix(): sha(p) for p in root.rglob('*') if p.is_file() and p.name != 'report.json'}))
    return root / 'report.json'


def source_panel(root):
    sample, prep, docking = [root / name for name in ('sample', 'prep', 'docking')]
    for directory in (sample, prep, docking): directory.mkdir()
    with sqlite3.connect(sample / 'samples.sqlite') as db:
        db.executescript('CREATE TABLE selected(cid TEXT, mid TEXT); CREATE TABLE sample(scheme TEXT,block_id TEXT,cid TEXT); CREATE TABLE population(scheme TEXT,block_id TEXT,n INTEGER);')
        for bid in ('A', 'B'):
            db.execute('INSERT INTO population VALUES(?,?,?)', ('E095', bid, 1000))
            for i in range(20):
                db.execute('INSERT INTO selected VALUES(?,?)', (f'{bid}{i}', f'{bid}M{i // 2}'))
                db.execute('INSERT INTO sample VALUES(?,?,?)', ('E095', bid, f'{bid}{i}'))
    sample_report = seal(sample, kind='block_sample')
    (prep / 'receptors').mkdir()
    (prep / 'receptors/R1.mol2').write_text('Fixture receptor')
    (prep / 'protein-input.pdb').write_text('Fixture original PDB')
    save(prep / 'profile.json', dict(receptors=[dict(id='R1', mol2=str(prep / 'protein.mol2'))]))
    save(prep / 'signature.json', dict(sample_sha256=sha(sample_report)))
    prepared_report = seal(prep, kind='block_plants_prepare', sample_report=str(sample_report))
    save(docking / 'signature.json', dict(prepared_sha256=sha(prepared_report)))
    rows = [dict(cid=f'{bid}{i}', alias=f'{bid}{i}', receptor='R1', job='fixture', status='ok',
                 score=-100 if bid == 'A' else -50, pose_name=f'{bid}{i}', pose_file='saved.mol2', pose_index=i)
            for bid in ('A', 'B') for i in range(20)]
    es.write_csv(docking / 'scores.csv', rows)
    return seal(docking, kind='block_plants_run', score='ChemPLP', prepared_report=str(prepared_report), sample_report=str(sample_report))


@pytest.fixture
def engine(tmp_path, monkeypatch):
    source = source_panel(tmp_path)
    profile = tmp_path / 'profile.json'; profile.write_text('{}')
    worker = tmp_path / 'worker.py'; worker.write_text('# Fixture worker')
    monkeypatch.setattr(es, 'load_profile', lambda _: (dict(python='fixture-python'), worker, {str(profile): sha(profile)}))
    calls = []
    def execute(command, **kwargs):
        job_path = Path(command[-1]); job = json.loads(job_path.read_text()); calls.append(job)
        with sqlite3.connect(job_path.parent / 'predictions.sqlite') as db:
            db.execute('CREATE TABLE IF NOT EXISTS scores(cid TEXT,receptor TEXT,status TEXT,score REAL,error TEXT,PRIMARY KEY(cid,receptor))')
            with open(job['scores'], newline='') as stream:
                for row in csv.DictReader(stream):
                    if job['selected_cids'] is not None and row['cid'] not in job['selected_cids']: continue
                    # B wins EquiScore; its odd-index conformers beat the even ones.
                    value = (0.7 if row['cid'][0] == 'B' else 0.1) + int(row['cid'][1:]) / 100
                    db.execute('INSERT OR IGNORE INTO scores VALUES(?,?,?,?,?)', (row['cid'], row['receptor'], 'ok', value, ''))
        save(job_path.parent / 'environment.json', dict(fixture=True))
        return SimpleNamespace(returncode=0)
    monkeypatch.setattr(es.subprocess, 'run', execute)
    return source, profile, calls


def test_pilot_full_rank_reversal_and_molecule_conformers(engine, tmp_path):
    source, profile, calls = engine
    pilot = tmp_path / 'pilot'; full = tmp_path / 'full'; analysis = tmp_path / 'analysis'
    with pytest.raises(ValueError, match='successful technical pilot'):
        es.run(source, profile, full, scope='full')
    assert not calls
    result = es.run(source, profile, pilot)
    assert result['pairs'] == result['scored'] == 32
    assert result['scope'] == 'pilot'
    with pytest.raises(ValueError, match='Only a full'):
        es.analyze(pilot / 'report.json', analysis)
    es.run(source, profile, pilot)
    assert len(calls) == 1  # Completed receipt reuse never submits inference again.
    result = es.run(source, profile, full, scope='full', pilot_report=pilot / 'report.json')
    assert result['scored'] == 40
    assert calls[-1]['selected_cids'] is None
    assert calls[-1]['pilot_environment'] == {'fixture': True}
    result = es.analyze(full / 'report.json', analysis)
    assert result['ranking']['better'] == 'higher'
    ranks = inspect(analysis / 'report.json', scheme='E095', receptor='R1')
    assert [row['block_id'] for row in ranks['blocks']] == ['B', 'A']
    assert ranks['blocks'][0]['top_n_mean'] == pytest.approx(0.8)
    top = inspect(analysis / 'report.json', operation='top', scheme='E095', receptor='R1', block_id='B')
    assert [row['cid'] for row in top['candidates']] == ['B19', 'B17', 'B15', 'B13', 'B11']
    assert top['score'] == 'EquiScore'
    with (analysis / 'rank_comparison.csv').open() as stream: comparison = list(csv.DictReader(stream))
    assert comparison[0]['chemplp_rank'] == '1' and comparison[0]['equiscore_rank'] == '2'
    assert comparison[1]['rank_improvement'] == '1'


def test_changed_pilot_inputs_block_full(engine, tmp_path):
    source, profile, calls = engine
    es.run(source, profile, tmp_path / 'pilot')
    profile.write_text('{"changed": true}')
    with pytest.raises(ValueError, match='different inputs'):
        es.run(source, profile, tmp_path / 'full', scope='full', pilot_report=tmp_path / 'pilot/report.json')
    assert len(calls) == 1


def test_startup_failure_displays_underlying_exception(engine, tmp_path, monkeypatch):
    source, profile, _ = engine
    def fail(command, **kwargs):
        kwargs['stdout'].write('Traceback (most recent call last):\nValueError: fixture receptor mismatch\n')
        kwargs['stdout'].flush()
        return SimpleNamespace(returncode=1)
    monkeypatch.setattr(es.subprocess, 'run', fail)
    with pytest.raises(RuntimeError, match='worker exit code 1') as caught:
        es.run(source, profile, tmp_path / 'pilot')
    assert 'ValueError: fixture receptor mismatch' in str(caught.value)
    assert not (tmp_path / 'pilot/report.json').exists()


def test_worker_log_tail_excludes_previous_attempt_and_is_bounded(tmp_path):
    log = tmp_path / 'worker.log'
    log.write_bytes(b'Previous unrelated error\n')
    start = log.stat().st_size
    with log.open('ab') as stream:
        stream.write(b'Current failure\n')
    error = es.worker_failure(tmp_path, 'Failed', 2, start)
    assert 'Previous unrelated' not in str(error) and 'Current failure' in str(error)
    with log.open('ab') as stream:
        stream.write(b'x' * 15000 + b'\nRuntimeError: final cause\n')
    assert len(str(es.worker_failure(tmp_path, 'Failed', 2, start))) < 12500
    assert 'RuntimeError: final cause' in str(es.worker_failure(tmp_path, 'Failed', 2, start))


def test_partial_receipt_resumes_only_missing_pairs(engine, tmp_path):
    source, profile, calls = engine
    out = tmp_path / 'pilot'
    es.run(source, profile, out)
    with sqlite3.connect(out / 'predictions.sqlite') as db:
        db.execute('DELETE FROM scores WHERE rowid=(SELECT min(rowid) FROM scores)')
    report = json.loads((out / 'report.json').read_text())
    report.update(status='partial', failed_pairs=1)
    report['output_hashes']['predictions.sqlite'] = sha(out / 'predictions.sqlite')
    save(out / 'report.json', report)
    assert es.run(source, profile, out)['status'] == 'complete'
    assert len(calls) == 2


@pytest.mark.parametrize('status,score,error', [('ok', float('inf'), 'Invalid'), ('ok', -1, 'Invalid'), ('model_failed', 0.0, 'must be null'), ('invented', None, 'Unknown')])
def test_invalid_predictions_rejected(tmp_path, status, score, error):
    es.write_csv(tmp_path / 'scores.csv', [dict(cid='C', receptor='R', status='ok', score=-10)])
    with sqlite3.connect(tmp_path / 'predictions.sqlite') as db:
        db.execute('CREATE TABLE scores(cid,receptor,status,score,error)')
        db.execute('INSERT INTO scores VALUES(?,?,?,?,?)', ('C', 'R', status, score, ''))
    with pytest.raises(ValueError, match=error): es.collect(tmp_path / 'scores.csv', tmp_path / 'predictions.sqlite')


def test_missing_unknown_and_duplicate_predictions(tmp_path):
    es.write_csv(tmp_path / 'scores.csv', [dict(cid='C', receptor='R', status='ok', score=-10)])
    path = tmp_path / 'predictions.sqlite'
    with sqlite3.connect(path) as db: db.execute('CREATE TABLE scores(cid,receptor,status,score,error)')
    missing = es.collect(tmp_path / 'scores.csv', path)[0]
    assert missing['status'] == 'not_run' and missing['score'] is None
    with sqlite3.connect(path) as db: db.execute("INSERT INTO scores VALUES('unknown','R','ok',0.8,'')")
    with pytest.raises(ValueError, match='Unknown/unrequested'): es.collect(tmp_path / 'scores.csv', path)
    with sqlite3.connect(path) as db: db.execute("INSERT INTO scores VALUES('unknown','R','ok',0.8,'')")
    with pytest.raises(ValueError, match='Duplicate prediction'): es.collect(tmp_path / 'scores.csv', path)


def test_worker_rejects_unsealed_pose_and_changed_input(tmp_path):
    filename = Path(__file__).resolve().parents[1] / 'scripts/equiscore_worker.py'
    spec = importlib.util.spec_from_file_location('equiscore_worker', filename)
    worker = importlib.util.module_from_spec(spec); spec.loader.exec_module(worker)
    pose = tmp_path / 'attempt/poses/docked_ligands.mol2'; pose.parent.mkdir(parents=True)
    pose.write_text('fixture pose')
    ligand = pose.parent.parent / 'ligands.mol2'; ligand.write_text('fixture ligand')
    row = dict(pose_file=str(pose))
    with pytest.raises(ValueError, match='not sealed'): worker.verify_pose_inputs(row, tmp_path, {}, set())
    hashes = {p.relative_to(tmp_path).as_posix(): sha(p) for p in (pose, ligand)}
    worker.verify_pose_inputs(row, tmp_path, hashes, set())
    ligand.write_text('changed ligand')
    with pytest.raises(ValueError, match='changed'): worker.verify_pose_inputs(row, tmp_path, hashes, set())


@pytest.fixture
def worker_module():
    filename = Path(__file__).resolve().parents[1] / 'scripts/equiscore_worker.py'
    spec = importlib.util.spec_from_file_location('equiscore_worker', filename)
    worker = importlib.util.module_from_spec(spec); spec.loader.exec_module(worker)
    return worker


def test_screening_checkpoint_ignores_only_unused_known_keys(worker_module):
    from unittest.mock import Mock
    model = Mock()
    model.state_dict.return_value = {'layer.weight': object()}
    weight = object()
    state = {'layer.weight': weight, 'mu': 0, 'dev': 1}
    assert worker_module.load_screening_state(model, {'model': state}) == ['dev', 'mu']
    model.load_state_dict.assert_called_once_with({'layer.weight': weight}, strict=True)
    assert set(state) == {'layer.weight', 'mu', 'dev'}  # Never mutate the checkpoint.


def test_screening_checkpoint_unknown_extra_rejected(worker_module):
    from unittest.mock import Mock
    model = Mock(); model.state_dict.return_value = {'weight': object()}
    with pytest.raises(ValueError, match='unknown'):
        worker_module.load_screening_state(model, {'model': {'weight': 1, 'mu': 0, 'unknown': 2}})
    model.load_state_dict.assert_not_called()


@pytest.mark.parametrize('failure', ['Missing key: layer.weight', 'size mismatch for layer.weight'])
def test_screening_checkpoint_strict_errors_propagate(worker_module, failure):
    from unittest.mock import Mock
    model = Mock(); model.state_dict.return_value = {'layer.weight': object()}
    model.load_state_dict.side_effect = RuntimeError(failure)
    with pytest.raises(RuntimeError, match=failure):
        worker_module.load_screening_state(model, {'model': {'mu': 0, 'dev': 1}})
    assert model.load_state_dict.call_args.kwargs == {'strict': True}


def test_screening_checkpoint_keeps_keys_expected_by_model(worker_module):
    from unittest.mock import Mock
    model = Mock(); model.state_dict.return_value = {'mu': 0, 'dev': 1}
    assert worker_module.load_screening_state(model, {'model': {'mu': 7, 'dev': 8}}) == []
    model.load_state_dict.assert_called_once_with({'mu': 7, 'dev': 8}, strict=True)


@pytest.mark.parametrize('modern', [True, False])
def test_pinned_checkpoint_loading_across_torch_defaults(worker_module, tmp_path, monkeypatch, modern):
    path = tmp_path / 'weights.pt'; path.write_bytes(b'fixture checkpoint')
    monkeypatch.setattr(worker_module, 'WEIGHT_SHA', sha(path))
    calls = []
    def new_load(filename, map_location=None, weights_only=True):
        calls.append((filename, map_location, weights_only)); return {'model': {}}
    def old_load(filename, map_location=None):
        calls.append((filename, map_location)); return {'model': {}}
    torch = SimpleNamespace(load=new_load if modern else old_load)
    assert worker_module.read_checkpoint(torch, path) == {'model': {}}
    assert calls == [(str(path), 'cpu', False)] if modern else calls == [(str(path), 'cpu')]
    path.write_bytes(b'changed checkpoint')
    with pytest.raises(ValueError, match='hash mismatch'): worker_module.read_checkpoint(torch, path)
    assert len(calls) == 1


def test_profile_preserves_python_symlink_path(tmp_path, monkeypatch):
    repo = tmp_path / 'repo'; repo.mkdir()
    weights = repo / 'workdir/official_weight/save_model_screen.pt'
    weights.parent.mkdir(parents=True); weights.write_bytes(b'fixture checkpoint')
    monkeypatch.setattr(es, 'WEIGHT_SHA', sha(weights))
    monkeypatch.setattr(es.subprocess, 'run', lambda *a, **k: SimpleNamespace(stdout=es.COMMIT, returncode=0))
    interpreter = tmp_path / 'overlay/bin/python'
    original_resolve = Path.resolve
    def resolve(path, *args, **kwargs):
        if path == interpreter: return tmp_path / 'base/bin/python'
        return original_resolve(path, *args, **kwargs)
    monkeypatch.setattr(Path, 'resolve', resolve)
    profile = tmp_path / 'profile.json'
    save(profile, dict(python=str(interpreter), repository=str(repo)))
    config, _, _ = es.load_profile(profile)
    assert config['python'] == str(interpreter)


def test_existing_environment_setup_stops_before_install_on_gpu_failure(tmp_path, monkeypatch):
    import os
    filename = Path(__file__).resolve().parents[1] / 'scripts/setup_equiscore_existing.py'
    spec = importlib.util.spec_from_file_location('setup_equiscore_existing', filename)
    setup = importlib.util.module_from_spec(spec); spec.loader.exec_module(setup)
    monkeypatch.setattr(setup, 'os', SimpleNamespace(name='posix', environ=os.environ, path=os.path))
    monkeypatch.setattr(setup.sys, 'argv', ['setup', '--prefix', str(tmp_path / 'overlay')])
    monkeypatch.setattr(setup.sys, 'version_info', (3, 9, 23))
    calls = []
    def fail(command, **kwargs):
        calls.append(command)
        raise RuntimeError('DGL CUDA fixture failure')
    monkeypatch.setattr(setup.subprocess, 'run', fail)
    with pytest.raises(RuntimeError, match='DGL CUDA fixture'):
        setup.main()
    assert len(calls) == 1 and calls[0][-1] == '--gpu-check'
    assert not (tmp_path / 'overlay').exists()


def test_overlay_constraints_fix_numpy_abi_without_changing_gpu_stack(monkeypatch):
    filename = Path(__file__).resolve().parents[1] / 'scripts/setup_equiscore_existing.py'
    spec = importlib.util.spec_from_file_location('setup_equiscore_existing', filename)
    setup = importlib.util.module_from_spec(spec); spec.loader.exec_module(setup)
    versions = {'numpy': '2.0.2', 'MDAnalysis': '2.7.0', 'torch': '2.7.0+cu128',
                'dgl': '2.5.0+cu121', 'rdkit': '2025.9.2', 'pandas': '2.2.3'}
    monkeypatch.setattr(setup.importlib.metadata, 'distributions', lambda: [
        SimpleNamespace(metadata={'Name': key}, version=value) for key, value in versions.items()])
    constraints = set(setup.package_constraints())
    assert {'numpy==1.26.4', 'mdanalysis==2.7.0', 'prolif==1.1.0'} <= constraints
    assert {'torch==2.7.0+cu128', 'dgl==2.5.0+cu121', 'rdkit==2025.9.2', 'pandas==2.2.3'} <= constraints
    assert 'numpy==2.0.2' not in constraints
