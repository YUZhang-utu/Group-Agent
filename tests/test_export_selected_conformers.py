import importlib.util
import json
from pathlib import Path
import sqlite3

import numpy as np
import pytest


spec = importlib.util.spec_from_file_location('selected_export', Path(__file__).parents[1] / 'scripts/export_selected_conformers.py')
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


def fixture(tmp):
    batch, search, out = tmp / 'batch', tmp / 'search', tmp / 'out'
    shard = batch / 'artifacts/shard'
    shard.mkdir(parents=True)
    search.mkdir()
    source = tmp / 'source.mol2'
    records = [f'@<TRIPOS>MOLECULE\nname_conf{i}\n0 0 0 0 0\nSMALL\nNO_CHARGES\n' for i in range(4)]
    source.write_text(''.join(records), encoding='utf-8')
    cids, mids = ['c0', 'c1', 'c2', 'c3'], ['m0', 'm0', 'm1', 'm2']
    for name, data in [('conformer_ids.bin', cids), ('molecule_ids.bin', mids)]:
        np.asarray(data, dtype='S16').tofile(shard / name)
    module.save(shard / 'manifest.json', {'files': {p.name: {'sha256': module.sha(p)} for p in shard.glob('*.bin')}})
    module.save(batch / 'artifacts/catalog.json', {'shards': [{'path': str(shard), 'name': 'shard', 'global_id_start': 0, 'conformers': 4}]})
    with sqlite3.connect(batch / 'registry.sqlite3') as db:
        db.executescript('CREATE TABLE molecule(id TEXT PRIMARY KEY,source_name TEXT);'
                         'CREATE TABLE conformer(id TEXT PRIMARY KEY,molecule_id TEXT,source_path TEXT,source_record_index INT,source_record_name TEXT,content_sha256 TEXT);')
        db.executemany('INSERT INTO molecule VALUES (?,?)', [(m, m) for m in sorted(set(mids))])
        db.executemany('INSERT INTO conformer VALUES (?,?,?,?,?,?)',
                       [(cids[i], mids[i], str(source), i, f'name_conf{i}', module.hashlib.sha256(r.encode()).hexdigest()) for i, r in enumerate(records)])
    with sqlite3.connect(search / 'ranking.sqlite') as db:
        db.execute('CREATE TABLE ranking(mid TEXT PRIMARY KEY,rank INT)')
        db.executemany('INSERT INTO ranking VALUES (?,?)', [('m0', 2), ('m1', 1)])
    np.save(search / 'selected-molecules.npy', np.asarray(['m0', 'm1'], dtype='S16'))
    np.save(search / 'selected-conformers.npy', np.asarray([0, 1, 2]))
    seal(search)
    return batch, search, out, source, records


def seal(search):
    module.save(search / 'report.json', {'status': 'complete', 'retrieval': {'outputs': {
        str(p): module.sha(p) for p in search.glob('*.npy')}},
        'outputs': {str(search / 'ranking.sqlite'): module.sha(search / 'ranking.sqlite')}})


def test_all_conformers_and_resume(tmp_path):
    batch, search, out, source, records = fixture(tmp_path)
    result = module.run(batch, search, out)
    assert result['molecules'] == 2 and result['conformers'] == 3
    assert (out / 'part-00001.mol2').read_text() == ''.join(records[:3])
    digest = module.sha(out / 'part-00001.mol2')
    module.run(batch, search, out, resume=True)
    assert module.sha(out / 'part-00001.mol2') == digest
    (out / 'part-00001.mol2').write_text('changed')
    with pytest.raises(ValueError, match='part hash'):
        module.run(batch, search, out, resume=True)


def test_missing_conformer_rejected(tmp_path):
    batch, search, out, _, _ = fixture(tmp_path)
    np.save(search / 'selected-conformers.npy', np.asarray([0, 2]))
    seal(search)
    with pytest.raises(ValueError, match='all registered conformers'):
        module.run(batch, search, out)
    assert not (out / 'report.json').exists()


def test_changed_raw_record_rejected_and_resumed(tmp_path):
    batch, search, out, source, records = fixture(tmp_path)
    source.write_text(''.join(records).replace('name_conf1', 'bad_name'))
    with pytest.raises(ValueError, match='Source record changed'):
        module.run(batch, search, out)
    assert not (out / 'report.json').exists()
    source.write_text(''.join(records))
    assert module.run(batch, search, out, resume=True)['conformers'] == 3
