import importlib.util
import json
from pathlib import Path
import sqlite3

import numpy as np
import pytest

spec = importlib.util.spec_from_file_location(
    'bulk_precompute', Path(__file__).resolve().parents[1] / 'scripts' / 'precompute_library_batch.py')
bulk = importlib.util.module_from_spec(spec)
spec.loader.exec_module(bulk)


def block(name, offset=0, atom_type='C.3'):
    return f'''@<TRIPOS>MOLECULE
{name}
2 1 0 0 0
SMALL
USER_CHARGES
@<TRIPOS>ATOM
1 C1 {offset} 0 0 {atom_type} 1 MOL 0.0
2 N1 {offset+1} 0 0 N.3 1 MOL 0.0
@<TRIPOS>BOND
1 1 2 1
'''


def test_snapshot_detects_source_change(tmp_path):
    p = tmp_path/'one.mol2'; p.write_text(block('A'))
    item = bulk.snapshot(p); bulk.check_source(item)
    p.write_text(block('A')+'\n')
    with pytest.raises(ValueError, match='Source changed'):
        bulk.check_source(item)


def test_old_ids_restored_and_new_ids_do_not_replace_them(tmp_path):
    source = tmp_path/'old.mol2'; source.write_text(block('A_conf0')+block('A_conf1', 2))
    shard = tmp_path/'old'; shard.mkdir()
    np.asarray([b'MOL-OLD', b'MOL-OLD'], dtype='S16').tofile(shard/'molecule_ids.bin')
    np.asarray([b'CNF-OLD0', b'CNF-OLD1'], dtype='S16').tofile(shard/'conformer_ids.bin')
    db = tmp_path/'build.sqlite3'; bulk.private_registry(db, 'LIB-OLD')
    report = bulk.register_source(db, 'LIB-OLD', source.resolve(), bulk.digest(source), shard)
    assert report['inserted'] == 2
    assert bulk.register_source(db, 'LIB-OLD', source.resolve(), bulk.digest(source), shard)['inserted'] == 2
    new = tmp_path/'new.mol2'; new.write_text(block('B_conf0'))
    bulk.register_source(db, 'LIB-OLD', new.resolve(), bulk.digest(new))
    with sqlite3.connect(db) as connection:
        assert connection.execute('SELECT id FROM molecule WHERE source_name="A"').fetchone()[0] == 'MOL-OLD'
        assert connection.execute('SELECT COUNT(*) FROM conformer').fetchone()[0] == 3
        assert connection.execute('SELECT COUNT(*) FROM batch_source').fetchone()[0] == 2


def test_conflicting_names_rollback_entire_new_source(tmp_path):
    db = tmp_path/'db'; bulk.private_registry(db, 'L')
    p = tmp_path/'bad.mol2'; p.write_text(block('A_conf0')+block('A_conf0', 3))
    with pytest.raises(ValueError, match='Conformer index conflict'):
        bulk.register_source(db, 'L', p.resolve(), bulk.digest(p))
    with sqlite3.connect(db) as c:
        assert c.execute('SELECT COUNT(*) FROM conformer').fetchone()[0] == 0
        assert c.execute('SELECT COUNT(*) FROM molecule').fetchone()[0] == 0
        assert c.execute('SELECT COUNT(*) FROM batch_source').fetchone()[0] == 0


def test_ordered_topology_change_is_not_silently_template_cached(tmp_path):
    db = tmp_path/'db'; bulk.private_registry(db, 'L')
    p = tmp_path/'bad.mol2'; p.write_text(block('A_conf0')+block('A_conf1', atom_type='O.3'))
    with pytest.raises(ValueError, match='topology conflict'):
        bulk.register_source(db, 'L', p.resolve(), bulk.digest(p))


def test_exact_duplicate_record_retains_original_source_mapping(tmp_path):
    db = tmp_path/'db'; bulk.private_registry(db, 'L')
    a = tmp_path/'a.mol2'; b = tmp_path/'b.mol2'
    a.write_text(block('A')); b.write_text(block('A')+block('B'))
    bulk.register_source(db, 'L', a.resolve(), bulk.digest(a))
    result = bulk.register_source(db, 'L', b.resolve(), bulk.digest(b))
    assert (result['inserted'], result['duplicates']) == (1, 1)
    with sqlite3.connect(db) as c:
        assert c.execute('SELECT source_path FROM conformer WHERE source_record_name="A"').fetchone()[0] == str(a.resolve())


def test_old_source_count_mismatch_rolls_back(tmp_path):
    db = tmp_path/'db'; bulk.private_registry(db, 'L')
    p = tmp_path/'old.mol2'; p.write_text(block('A'))
    shard = tmp_path/'old'; shard.mkdir()
    np.asarray([b'M1', b'M2'], dtype='S16').tofile(shard/'molecule_ids.bin')
    np.asarray([b'C1', b'C2'], dtype='S16').tofile(shard/'conformer_ids.bin')
    with pytest.raises(ValueError, match='count mismatch'):
        bulk.register_source(db, 'L', p.resolve(), bulk.digest(p), shard)
    with sqlite3.connect(db) as c:
        assert c.execute('SELECT COUNT(*) FROM conformer').fetchone()[0] == 0


def test_integrity_rejects_corrupted_completed_shard(tmp_path):
    p = tmp_path/'array.bin'; p.write_bytes(b'abc')
    bulk.write(tmp_path/'manifest.json', {'files': {'array.bin': {'bytes': 3, 'sha256': bulk.digest(p)}}})
    bulk.verify(tmp_path)
    p.write_bytes(b'xyz')
    with pytest.raises(ValueError, match='integrity'):
        bulk.verify(tmp_path)


def test_partial_recovery_preserves_bytes_and_rejects_external_move(tmp_path):
    output = tmp_path/'run'; output.mkdir()
    p = output/'.x.partial'; p.mkdir(); (p/'raw.bin').write_bytes(b'partial data')
    bulk.recover_partial(output, p)
    assert not p.exists()
    assert next((output/'recovery').rglob('raw.bin')).read_bytes() == b'partial data'
    external = tmp_path/'external'; external.mkdir()
    with pytest.raises(ValueError, match='external partial'):
        bulk.recover_partial(output, external)
    assert external.exists()


def test_manifest_rejects_directory_traversal(tmp_path):
    directory = tmp_path/'shard'; directory.mkdir()
    p = tmp_path/'outside'; p.write_bytes(b'x')
    bulk.write(directory/'manifest.json', {'files': {'../outside': {'bytes': 1, 'sha256': bulk.digest(p)}}})
    with pytest.raises(ValueError, match='escapes'):
        bulk.verify(directory)
