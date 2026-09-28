import json
import sqlite3

import numpy as np
import pytest

from aidd_agent.conformer_block_store import fit
from aidd_agent.conformer_block_validate import validate_model, independent_omega, validate
from aidd_agent.macrocycle_identity_audit import sha
from aidd_agent.macrocycle_blocks import angle


def fixture(tmp_path):
    rows = [dict(conformer_id=f'C{i}', molecule_id=f'M{i//3}', hard_group='g',
                 descriptor=[float(i), float(i % 2)]) for i in range(9)]
    descriptor = tmp_path / 'descriptors.sqlite'
    with sqlite3.connect(descriptor) as db:
        db.execute('CREATE TABLE descriptor(cid TEXT PRIMARY KEY,payload TEXT)')
        db.executemany('INSERT INTO descriptor VALUES(?,?)', ((r['conformer_id'], json.dumps(r)) for r in rows))
    model = tmp_path / 'model'
    fit(rows, model, 4, {'v': 1})
    return model, descriptor


def refresh_hashes(model):
    receipt = json.loads((model / 'report.json').read_text())
    receipt['output_hashes'] = {name: sha(model / name) for name in ('blocks.sqlite', 'memberships.jsonl')}
    (model / 'report.json').write_text(json.dumps(receipt))


def test_exhaustive_gate_accepts_valid_model(tmp_path):
    model, descriptor = fixture(tmp_path)
    result = validate_model(model, descriptor, {'v': 1})
    assert result['conformers_checked'] == 9
    assert result['structural_gate'] == 'passed'


def test_wrong_leaf_rejected_even_with_updated_checksums(tmp_path):
    model, descriptor = fixture(tmp_path)
    with sqlite3.connect(model / 'blocks.sqlite') as db:
        a = db.execute('SELECT node FROM point WHERE cid="C0"').fetchone()[0]
        b = db.execute('SELECT node FROM point WHERE cid="C8"').fetchone()[0]
        assert a != b
        db.execute('UPDATE point SET node=? WHERE cid="C0"', (b,))
        db.execute('UPDATE point SET node=? WHERE cid="C8"', (a,))
        with (model / 'memberships.jsonl').open('w') as stream:
            for cid, mid, node, provenance in db.execute('SELECT cid,mid,node,provenance FROM point ORDER BY cid'):
                stream.write(json.dumps(dict(conformer_id=cid, molecule_id=mid, block_id=node,
                                            provenance=json.loads(provenance)))+'\n')
    refresh_hashes(model)
    with pytest.raises(ValueError, match='wrong frozen-tree leaf'):
        validate_model(model, descriptor, {'v': 1})


@pytest.mark.parametrize('mutation,reason', [
    ("UPDATE point SET mid='WRONG' WHERE cid='C0'", 'identity'),
    ("UPDATE point SET group_id='WRONG' WHERE cid='C0'", 'stratum'),
    ("UPDATE tree SET n=100 WHERE axis IS NULL", 'capacity|counts'),
    ("UPDATE tree SET radius=999 WHERE axis IS NULL", 'radius'),
])
def test_semantic_corruption_not_just_checksums(tmp_path, mutation, reason):
    model, descriptor = fixture(tmp_path)
    with sqlite3.connect(model / 'blocks.sqlite') as db:
        db.execute(mutation)
    refresh_hashes(model)
    with pytest.raises(ValueError, match=reason):
        validate_model(model, descriptor, {'v': 1})


def test_export_missing_row_detected_after_rehash(tmp_path):
    model, descriptor = fixture(tmp_path)
    path = model / 'memberships.jsonl'
    path.write_text(''.join(path.read_text().splitlines(keepends=True)[:-1]))
    refresh_hashes(model)
    with pytest.raises(ValueError, match='coverage'):
        validate_model(model, descriptor, {'v': 1})


def test_capacity_gate_even_when_metadata_agrees(tmp_path):
    model, descriptor = fixture(tmp_path)
    with sqlite3.connect(model / 'blocks.sqlite') as db:
        db.execute("UPDATE metadata SET value='1' WHERE key='capacity'")
    receipt = json.loads((model / 'report.json').read_text())
    receipt['capacity'] = 1
    (model / 'report.json').write_text(json.dumps(receipt))
    refresh_hashes(model)
    with pytest.raises(ValueError, match='capacity'):
        validate_model(model, descriptor, {'v': 1})


def test_independent_angle_formula():
    rng = np.random.default_rng(12)
    for _ in range(100):
        points = rng.normal(size=(4, 3))
        difference = (independent_omega(points)-np.degrees(angle(points))+180) % 360-180
        assert abs(difference) < 1e-10
    with pytest.raises(ValueError, match='Degenerate'):
        independent_omega(np.zeros((4, 3)))


def test_failed_build_emits_explicit_hold_receipt(tmp_path):
    build = tmp_path / 'build'; build.mkdir()
    (build / 'report.json').write_text('{"status":"running"}')
    result = validate(build, tmp_path / 'validation')
    assert result['status'] == 'failed' and result['release_gate'] == 'hold'
    assert result['errors']
    assert (tmp_path / 'validation/report.json').is_file()
