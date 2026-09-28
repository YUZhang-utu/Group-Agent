import json
import sqlite3

import pytest

from aidd_agent.macrocycle_identity_audit import sha
from aidd_agent.macrocycle_preflight import inventory, block_readiness


def test_inventory_coverage_and_output_guard(tmp_path):
    source = tmp_path / 'sources'
    source.mkdir()
    (source / 'a.mol2').write_text('source a')
    (source / 'b.mol2').write_text('source b')
    (source / 'a.csv').write_text('Name,SMILES\n')
    (source / 'orphan.csv').write_text('Name,SMILES\n')
    with pytest.raises(ValueError, match='outside'):
        inventory(source, source / 'output')
    result = inventory(source, tmp_path / 'inventory')
    assert result['mol2_files'] == 2 and result['csv_files'] == 1
    assert result['mol2_without_csv'] == [str(source / 'b.mol2')]
    assert result['undiscovered_csv'] == [str(source / 'orphan.csv')]
    assert all(len(row['sha256']) == 64 for row in result['files'])
    with pytest.raises(FileExistsError):
        inventory(source, tmp_path / 'inventory')


def fixture_audit(tmp_path):
    audit = tmp_path / 'audit'
    audit.mkdir()
    rows = [
        dict(name='A', conformer_name='A_conf1', mol2_canonical_smiles='one'),
        dict(name='A', conformer_name='A_conf2', mol2_canonical_smiles='two'),
        dict(name='B', conformer_name='B_conf1', mol2_canonical_smiles='three'),
        dict(name='B', conformer_name='B_conf1', mol2_canonical_smiles='three'),
        dict(name='C', conformer_name='C_conf1', name_mapping='unverified',
             verification_source='mol2_only', omega_states=['boundary'], name_rotations=None),
        dict(name='D', conformer_name='D_conf1', status='review', reason='chemical_graph_mismatch'),
    ]
    with sqlite3.connect(audit / 'audit.sqlite') as db:
        db.execute('CREATE TABLE conformer(name TEXT, conf_name TEXT, payload TEXT)')
        for item in rows:
            row = dict(status='ok', ring_source_atom_ids=[1, 2, 3], name_mapping='verified_against_mol2')
            row.update(item)
            db.execute('INSERT INTO conformer VALUES(?,?,?)', (row['name'], row['conformer_name'], json.dumps(row)))
    (audit / 'issues.jsonl').write_text('')
    report = dict(status='complete', output_hashes={p: sha(audit / p) for p in ('audit.sqlite', 'issues.jsonl')})
    (audit / 'report.json').write_text(json.dumps(report))
    return audit


def test_readiness_cross_record_conflicts_and_flag_retention(tmp_path):
    audit = fixture_audit(tmp_path)
    before = sha(audit / 'audit.sqlite')
    result = block_readiness(audit, tmp_path / 'ready')
    assert result['counts'] == {'review': 5, 'geometry_candidate': 1}
    rows = [json.loads(line) for line in (tmp_path / 'ready/conformer_readiness.jsonl').read_text().splitlines()]
    assert len(rows) == 6
    assert 'same_name_chemistry_conflict' in rows[0]['review_reasons']
    assert 'duplicate_conformer_name' in rows[2]['review_reasons']
    assert rows[4]['state'] == 'geometry_candidate'
    assert len(rows[4]['flags']) == 2
    assert all(not row['production_block_assignment'] for row in rows)
    assert sha(audit / 'audit.sqlite') == before
    assert all(row['minimum_blocks_ignoring_strata'] == 1 for row in result['capacity_planning'])


@pytest.mark.parametrize('change', ['hash', 'status'])
def test_reject_stale_or_incomplete_audit(tmp_path, change):
    audit = fixture_audit(tmp_path)
    if change == 'hash':
        (audit / 'issues.jsonl').write_text('modified')
    else:
        report = json.loads((audit / 'report.json').read_text())
        report['status'] = 'running'
        (audit / 'report.json').write_text(json.dumps(report))
    with pytest.raises(ValueError):
        block_readiness(audit, tmp_path / 'ready')
    assert not (tmp_path / 'ready').exists()
