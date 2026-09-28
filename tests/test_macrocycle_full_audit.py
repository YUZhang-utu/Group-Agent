import json
from pathlib import Path

import pytest
from rdkit import Chem

from aidd_agent.macrocycle_peptide import peptide, tokens
from aidd_agent.macrocycle_full_audit import csv_one, mol_one
from aidd_agent.mol2 import parse_mol2_block


SMILES = 'N1[C@@H](C)C(=O)N(C)[C@@H](C)C(=O)N2CCC[C@H]2C1=O'
NAME = 'c--A-Anme-P-c'


def test_tokens_and_proline_main_ring():
    assert tokens('c--A-W-Lnme-AdFnme-dWnme-c') == ['A', 'W', 'Lnme', 'A', 'dFnme', 'dWnme']
    mol = Chem.MolFromSmiles(SMILES)
    result = peptide(mol, NAME)
    assert len(result['units']) == 3
    assert len(set(result['ring_atoms'])) == 9
    assert result['tokens'] == ['A', 'Anme', 'P']
    assert result['name_rotations'] == 1
    permuted = Chem.RenumberAtoms(mol, list(reversed(range(mol.GetNumAtoms()))))
    other = peptide(permuted, NAME)
    assert other['tokens'] == result['tokens']
    explicit = Chem.RemoveHs(Chem.AddHs(mol))
    for atom in explicit.GetAtoms():
        atom.SetNumExplicitHs(atom.GetTotalNumHs())
        atom.SetNoImplicit(True)
    assert peptide(explicit, NAME)['tokens'] == result['tokens']
    with pytest.raises(ValueError, match='name_structure_residue_mismatch'):
        peptide(mol, 'c--dA-Anme-P-c')
    with pytest.raises(ValueError, match='name_structure_residue_mismatch'):
        peptide(mol, 'c--A-A-P-c')


def test_repeated_residues_and_unknown_names():
    result = peptide(Chem.MolFromSmiles('N1CC(=O)NCC(=O)NCC1=O'), 'c--G-G-G-c')
    assert result['name_rotations'] == 3
    with pytest.raises(ValueError, match='unsupported_residue_name'):
        tokens('c--A-Z-P-c')
    bad = csv_one(('source.csv', 2, 'c--dA-Anme-P-c', SMILES))
    assert bad['status'] == 'review'
    assert bad['reason'] == 'name_structure_residue_mismatch'


def test_mol_worker_missing_and_duplicate_name():
    raw = '@<TRIPOS>MOLECULE\nc--G-G-G-c_conf1\n1 0 0 0 0\nSMALL\nNO_CHARGES\n@<TRIPOS>ATOM\n1 C1 0 0 0 C.3 1 UNK 0\n@<TRIPOS>BOND\n'
    record = parse_mol2_block(Path('source.mol2'), 0, raw)
    assert mol_one((record, None))['reason'] == 'name_missing_from_csv'
    assert mol_one((record, {'duplicate': True}))['reason'] == 'duplicate_csv_name'


def test_full_scan_reports_missing_molecules_and_suffix_shortfall(tmp_path):
    import csv
    from aidd_agent.macrocycle_full_audit import run
    source = tmp_path / 'source.mol2'
    source.write_text('@<TRIPOS>MOLECULE\nc--G-G-G-c_conf1\n1 0 0 0 0\nSMALL\nNO_CHARGES\n@<TRIPOS>ATOM\n1 C1 0 0 0 C.3 1 UNK 0\n@<TRIPOS>BOND\n')
    mapping = tmp_path / 'source.csv'
    with mapping.open('w', newline='') as stream:
        writer = csv.writer(stream)
        writer.writerow(['SMILES', 'Name'])
        writer.writerow([SMILES, NAME])
    output = tmp_path / 'audit'
    report = run([mapping], [source], output, workers=1)
    assert report['status'] == 'complete'
    assert report['conformer_counts'] == {'name_missing_from_csv': 1}
    assert report['csv_names_without_mol2'] == 1
    assert report['non_123_suffix_groups'] == 1
    assert report['conformer_count_histogram'] == {1: 1}
    issues = [json.loads(line) for line in (output / 'issues.jsonl').read_text().splitlines()]
    assert {row['kind'] for row in issues} == {'mol2', 'non_123_suffix_group', 'csv_name_without_mol2'}
    with pytest.raises(FileExistsError):
        run([mapping], [source], output, workers=1)


def test_discovery_optional_same_stem_csv(tmp_path):
    from aidd_agent.macrocycle_full_audit import discover
    (tmp_path / 'nested').mkdir()
    for name in ('a.mol2', 'a.csv', 'nested/b.mol2', 'unrelated.csv'):
        (tmp_path / name).write_text('')
    csv_paths, mol_paths = discover(tmp_path)
    assert [p.name for p in csv_paths] == ['a.csv']
    assert [p.name for p in mol_paths] == ['a.mol2', 'b.mol2']


def test_mol2_only_geometry_and_no_csv_conflict_fallback(monkeypatch):
    import numpy as np
    import aidd_agent.macrocycle_full_audit as module
    mol = Chem.MolFromSmiles(SMILES)
    conf = Chem.Conformer(mol.GetNumAtoms())
    for i, point in enumerate(np.random.default_rng(8).normal(size=(mol.GetNumAtoms(), 3))):
        conf.SetAtomPosition(i, point)
    mol.AddConformer(conf)
    atoms = ''.join(f'{i+1} X 0 0 0 C.3 1 UNK 0\n' for i in range(mol.GetNumAtoms()))
    from types import SimpleNamespace
    record = SimpleNamespace(molecule_name=NAME, name=NAME+'_conf1', conformer_index=1,
        source_path=Path('source.mol2'), record_index=0, content_sha256='sourcehash',
        raw_text='@<TRIPOS>ATOM\n'+atoms+'@<TRIPOS>BOND\n')
    monkeypatch.setattr(module, 'load_rdkit_mol2', lambda *args: (Chem.Mol(mol), 'strict'))
    row = mol_one((record, None, True))
    assert row['status'] == 'ok'
    assert row['verification_source'] == 'mol2_only'
    assert row['csv_identity_verified'] is False
    assert len(row['omega_degrees']) == 3
    record.molecule_name = 'unknown_name'
    row = mol_one((record, None, True))
    assert row['status'] == 'ok' and row['name_mapping'] == 'unverified'
    assert row['name_rotations'] is None
    assert mol_one((record, {'duplicate': True}, True))['reason'] == 'duplicate_csv_name'
    bad = dict(status='ok', plain='wrong', iso='wrong')
    assert mol_one((record, bad, True))['reason'] == 'chemical_graph_mismatch'


def test_summary_histograms_match_sql_without_global_sort():
    import io
    import sqlite3
    from aidd_agent.macrocycle_full_audit import summarize
    db = sqlite3.connect(':memory:')
    db.executescript('''CREATE TABLE molecule(name TEXT, occurrences INTEGER, payload TEXT);
        CREATE TABLE conformer(name TEXT, conf_index INTEGER, conf_name TEXT,
            status TEXT, reason TEXT, payload TEXT);
        CREATE INDEX conf_name ON conformer(name);''')
    for i in range(15):
        payload = dict(ring_size=9 if i % 2 else 12, omega_states=['trans', 'cis'],
                       verification_source='mol2_only', name_rotations=1,
                       name_mapping='verified_against_mol2', mol2_canonical_smiles='C')
        db.execute('INSERT INTO conformer VALUES(?,?,?,?,?,?)',
                   (str(i // 3), i % 3 + 1, str(i), 'ok', None, json.dumps(payload)))
    expected = {}
    for field, key in [('ring_size_histogram', 'ring_size'),
                       ('omega_pattern_histogram', 'omega_states'),
                       ('verification_source_counts', 'verification_source')]:
        expected[field] = dict(db.execute(
            f"SELECT json_extract(payload,'$.{key}'),count(*) FROM conformer GROUP BY 1"))
    statements = []
    db.set_trace_callback(statements.append)
    report = {}
    summarize(db, report, io.StringIO())
    assert all(report[key] == value for key, value in expected.items())
    assert not any('GROUP BY 1' in query for query in statements)
    db.close()


def test_failed_summary_recovery_and_guards(tmp_path, monkeypatch):
    import sqlite3
    import aidd_agent.macrocycle_full_audit as module
    source = tmp_path / 'source.mol2'
    source.write_text('@<TRIPOS>MOLECULE\nc--G-G-G-c_conf1\n1 0 0 0 0\nSMALL\nNO_CHARGES\n@<TRIPOS>ATOM\n1 C1 0 0 0 C.3 1 UNK 0\n@<TRIPOS>BOND\n')
    output = tmp_path / 'audit'
    real = module.summarize

    def fail(db, report, issues):
        real(db, report, issues)
        raise sqlite3.OperationalError('database or disk is full')

    monkeypatch.setattr(module, 'summarize', fail)
    with pytest.raises(sqlite3.OperationalError):
        module.run([], [source], output, workers=1)
    original = (output / 'report.json').read_text()
    db_hash = module.sha(output / 'audit.sqlite')
    monkeypatch.setattr(module, 'summarize', real)
    saved = source.read_bytes()
    source.write_bytes(saved + b'\n')
    with pytest.raises(ValueError, match='hash changed'):
        module.recover_summary(output)
    source.write_bytes(saved)
    modified = json.loads(original)
    modified.pop('failed_stage', None)
    modified.pop('unique_csv_names', None)
    (output / 'report.json').write_text(json.dumps(modified))
    with pytest.raises(ValueError, match='No evidence ingestion reached summary'):
        module.recover_summary(output)
    modified = json.loads(original)
    modified['mol2_sources'][0]['records'] += 1
    (output / 'report.json').write_text(json.dumps(modified))
    with pytest.raises(ValueError, match='coverage mismatch'):
        module.recover_summary(output)
    modified['mol2_sources'] = []
    (output / 'report.json').write_text(json.dumps(modified))
    with pytest.raises(ValueError, match='Incomplete source manifest'):
        module.recover_summary(output)
    # Legacy reports have no failed_stage, and progress can add a zero ok counter.
    legacy = json.loads(original)
    legacy.pop('failed_stage', None)
    legacy['conformer_counts']['ok'] = 0
    original = json.dumps(legacy)
    (output / 'report.json').write_text(original)
    # A failed recovery can be retried without duplicating derived issues.
    monkeypatch.setattr(module, 'summarize', fail)
    with pytest.raises(sqlite3.OperationalError):
        module.recover_summary(output)
    assert (output / 'report.json').read_text() == original
    monkeypatch.setattr(module, 'summarize', real)
    report = module.recover_summary(output)
    assert report['status'] == 'complete'
    assert module.sha(output / 'audit.sqlite') == db_hash
    assert report['code_hashes'] == json.loads(original)['code_hashes']
    issues = [json.loads(line) for line in (output / 'issues.jsonl').read_text().splitlines()]
    assert sum(row['kind'] == 'non_123_suffix_group' for row in issues) == 1
    assert Path(report['summary_recovery']['backup'], 'report.json').read_text() == original
    from aidd_agent.macrocycle_preflight import block_readiness
    ready = block_readiness(output, tmp_path / 'readiness')
    assert ready['status'] == 'preparation_complete'
    with pytest.raises(ValueError, match='requires a failed audit'):
        module.recover_summary(output)
