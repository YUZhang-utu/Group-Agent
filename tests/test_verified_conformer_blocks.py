import json
import sqlite3
from types import SimpleNamespace

import numpy as np
import pytest
from rdkit import Chem

from aidd_agent.macrocycle_descriptors import describe_peptide
from aidd_agent.conformer_block_store import fit, open_model, route, propose, sample
from aidd_agent.verified_macrocycle_blocks import registry_identity


def molecule():
    mol = Chem.MolFromSmiles('N1[C@@H](C)C(=O)N(C)[C@@H](C)C(=O)N2CCC[C@H]2C1=O')
    coordinates = np.random.default_rng(81).normal(size=(mol.GetNumAtoms(), 3))
    conf = Chem.Conformer(mol.GetNumAtoms())
    for i, point in enumerate(coordinates):
        conf.SetAtomPosition(i, point)
    mol.AddConformer(conf)
    return mol


@pytest.mark.parametrize('variant', ['backbone', 'chemistry', 'typed'])
def test_peptide_proline_transform_and_reindexing(variant):
    mol = molecule()
    expected = describe_peptide(mol, 'c--A-Anme-P-c', variant)
    perm = list(map(int, np.random.default_rng(2).permutation(mol.GetNumAtoms())))
    changed = Chem.RenumberAtoms(mol, perm)
    rotation = np.array([[0, -1, 0], [1, 0, 0], [0, 0, 1]])
    coordinates = np.array(changed.GetConformer().GetPositions()) @ rotation + [7, -9, 2]
    for i, point in enumerate(coordinates):
        changed.GetConformer().SetAtomPosition(i, point)
    actual = describe_peptide(changed, 'c--A-Anme-P-c', variant)
    assert actual['hard_group'] == expected['hard_group']
    assert len(actual['ring_atoms']) == 9
    np.testing.assert_allclose(actual['descriptor'], expected['descriptor'], atol=1e-10)
    assert len(actual['descriptor']) == {'backbone': 18, 'chemistry': 51, 'typed': 177}[variant]


def test_graph_only_and_equivalent_rotations():
    mol = Chem.MolFromSmiles('N1CC(=O)NCC(=O)NCC1=O')
    conf = Chem.Conformer(mol.GetNumAtoms())
    for i, point in enumerate(np.random.default_rng(6).normal(size=(mol.GetNumAtoms(), 3))):
        conf.SetAtomPosition(i, point)
    mol.AddConformer(conf)
    a = describe_peptide(mol, 'c--G-G-G-c')
    b = describe_peptide(Chem.RenumberAtoms(mol, list(reversed(range(mol.GetNumAtoms())))))
    assert a['equivalent_rotations'] == 3
    assert b['name_mapping_verified'] is False
    np.testing.assert_allclose(a['descriptor'], b['descriptor'], atol=1e-10)


def test_sidechain_donor_acceptor_channels_are_separate():
    mol = Chem.MolFromSmiles('N1[C@@H](CCO)C(=O)N(C)[C@@H](C)C(=O)N2CCC[C@H]2C1=O')
    conf = Chem.Conformer(mol.GetNumAtoms())
    for i, point in enumerate(np.random.default_rng(1).normal(size=(mol.GetNumAtoms(), 3))):
        conf.SetAtomPosition(i, point)
    mol.AddConformer(conf)
    result = describe_peptide(mol)
    from aidd_agent.macrocycle_descriptors import FAMILIES
    vectors = np.array(result['descriptor']).reshape(3, -1)
    donor = FAMILIES.index('Donor'); acceptor = FAMILIES.index('Acceptor')
    for i, features in enumerate(result['typed_feature_atoms']):
        if any(family == 'Donor' for family, _ in features):
            assert vectors[i, 10+donor] > 0 and vectors[i, 10+acceptor] > 0
            assert np.linalg.norm(vectors[i, 17+6*donor:20+6*donor]) > 0
            assert np.linalg.norm(vectors[i, 17+6*acceptor:20+6*acceptor]) > 0
            break
    else:
        pytest.fail('Expected separately typed alcohol donor and acceptor features')


def test_disk_statistics_cross_chunk_boundary(tmp_path):
    result = fit(rows(2051), tmp_path / 'large', 80, {'v': 1})
    assert result['counts']['conformers'] == 2051
    assert max(result['block_size_histogram']) <= 80


def rows(n=19):
    return [dict(conformer_id=f'C{i:03}', molecule_id=f'M{i//3}', hard_group='a' if i % 2 else 'b',
                 descriptor=[float(i), float(i % 3)]) for i in range(n)]


def test_disk_fit_capacity_determinism_and_frozen_routing(tmp_path):
    data = rows()
    a = fit(iter(data), tmp_path / 'a', 4, {'version': 'test'})
    b = fit(iter(reversed(data)), tmp_path / 'b', 4, {'version': 'test'})
    assert a['model_id'] == b['model_id']
    assert a['counts']['conformers'] == 19
    assert max(a['block_size_histogram']) <= 4
    assert (tmp_path / 'a/memberships.jsonl').read_bytes() == (tmp_path / 'b/memberships.jsonl').read_bytes()
    db, _ = open_model(tmp_path / 'a')
    try:
        for row in data:
            node, reason = route(db, row)
            assert reason is None
            assert node == db.execute('SELECT node FROM point WHERE cid=?', (row['conformer_id'],)).fetchone()[0]
        assert db.execute('SELECT count(DISTINCT node) FROM point WHERE mid="M0"').fetchone()[0] > 1
    finally:
        db.close()


def test_incremental_capacity_unknown_and_schema(tmp_path):
    data = [dict(conformer_id=f'A{i}', molecule_id=f'M{i}', hard_group='g', descriptor=[0.]) for i in range(5)]
    fit(data, tmp_path / 'model', 4, {'v': 1})
    new = [dict(conformer_id=f'Z{i}', molecule_id=f'NEW{i}', hard_group='g', descriptor=[0.]) for i in range(3)]
    new.append(dict(conformer_id='U', molecule_id='U', hard_group='unknown', descriptor=[0.]))
    new.append(dict(conformer_id='Y', molecule_id='Y', hard_group='g', descriptor=[100.]))
    new.append(data[0])
    result = propose(new, tmp_path / 'model', tmp_path / 'proposal', {'v': 1})
    assert result['counts'] == dict(proposed=1, capacity_overflow=2, unknown_hard_group=1,
                                    outside_frozen_radius=1, existing_conformer_id=1)
    with pytest.raises(ValueError, match='schema'):
        propose([], tmp_path / 'model', tmp_path / 'wrong', {'v': 2})
    db, _ = open_model(tmp_path / 'model')
    assert db.execute('SELECT count(*) FROM point').fetchone()[0] == 5
    db.close()


def test_samples_unique_molecules_and_all_in_block_conformers(tmp_path):
    fit(rows(), tmp_path / 'model', 10, {'v': 1})
    result = sample(tmp_path / 'model', tmp_path / 'samples', (1, 100), seed=7)
    one = [json.loads(line) for line in (tmp_path / 'samples/sample-1.jsonl').read_text().splitlines()]
    assert all(len({r['molecule_id'] for r in one if r['block_id'] == block}) == 1
               for block in {r['block_id'] for r in one})
    assert result['samples'][1]['conformer_evaluations'] == 19
    assert result['samples'][1]['unique_molecules'] == 7
    assert result['samples'][1]['block_molecule_pairs'] > 7
    db, _ = open_model(tmp_path / 'model')
    for block, mid in {(r['block_id'], r['molecule_id']) for r in one}:
        expected = {r[0] for r in db.execute('SELECT cid FROM point WHERE node=? AND mid=?', (block, mid))}
        assert expected == {r['conformer_id'] for r in one if (r['block_id'], r['molecule_id']) == (block, mid)}
    db.close()


def test_registry_requires_exact_provenance(tmp_path):
    db = sqlite3.connect(':memory:')
    db.executescript('''CREATE TABLE molecule(id,library_id,source_name);
        CREATE TABLE conformer(id,molecule_id,source_record_name,conformer_index,
        content_sha256,topology_sha256,source_path,source_record_index,atom_count,bond_count);''')
    source = tmp_path / 'a.mol2'
    db.execute('INSERT INTO molecule VALUES(?,?,?)', ('M', 'L', 'name'))
    db.execute('INSERT INTO conformer VALUES(?,?,?,?,?,?,?,?,?,?)',
               ('C', 'M', 'name_conf1', 1, 'hash', 'topology', str(source), 0, 20, 21))
    record = SimpleNamespace(molecule_name='name', name='name_conf1', conformer_index=1,
        content_sha256='hash', topology_sha256='topology', source_path=source, record_index=0,
        atom_count=20, bond_count=21)
    assert registry_identity(db, 'L', record) == ('C', 'M')
    record.content_sha256 = 'different'
    with pytest.raises(ValueError, match='provenance'):
        registry_identity(db, 'L', record)
    db.close()


def test_invalid_and_duplicate_descriptor_identity(tmp_path):
    data = rows(1)
    with pytest.raises(sqlite3.IntegrityError):
        fit(data+data, tmp_path / 'duplicate', 2, {})
    data[0]['descriptor'] = [float('nan')]
    with pytest.raises(ValueError, match='descriptor'):
        fit(data, tmp_path / 'invalid', 2, {})


def test_overlapping_molecule_replay_equal_conformer_cost(tmp_path):
    import csv
    from aidd_agent.conformer_block_replay import replay
    data = rows(39)
    fit(data, tmp_path / 'model', 20, {})
    path = tmp_path / 'scores.csv'
    with path.open('w', newline='') as stream:
        writer = csv.writer(stream)
        writer.writerow(['conformer_id', 'score', 'status'])
        for i, row in enumerate(data):
            writer.writerow([row['conformer_id'], i if i else '', '' if i else 'no_surviving_pose'])
    result = replay(tmp_path / 'model', path, initial=1, fraction=.5, top=4)
    assert result['adaptive']['scored_conformers'] == result['uniform_conformer']['scored_conformers']
    assert result['adaptive']['scored_conformers'] <= result['conformer_budget']
    assert all(row['effective_k'] == max(1, int(np.ceil(row['sampled_molecules']*.05))) for row in result['block_tails'])
    full = replay(tmp_path / 'model', path, initial=1, fraction=1, top=4)
    assert full['adaptive']['molecule_ranking_recall'] == 1
    assert full['adaptive']['evaluated_unique_molecules'] == 13
    with pytest.raises(ValueError, match='exceed'):
        replay(tmp_path / 'model', path, initial=100, fraction=.1)
    with pytest.raises(ValueError, match='memory limit'):
        replay(tmp_path / 'model', path, max_conformers=1)
    path.write_text('conformer_id,score\nC000,1\n')
    with pytest.raises(ValueError, match='every'):
        replay(tmp_path / 'model', path)
