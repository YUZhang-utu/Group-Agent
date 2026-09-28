import csv

from rdkit import Chem

from aidd_agent.macrocycle_identity_audit import compare, read_names


def test_chemical_identity_and_stereo_are_separate():
    reference = 'N[C@@H](C)C(=O)O'
    observed = Chem.AddHs(Chem.MolFromSmiles(reference))
    result = compare(reference, observed)
    assert result['status'] == 'identity_match'
    assert len(result['csv_atom_to_mol2_heavy_index']) == 6
    assert compare('N[C@H](C)C(=O)O', observed)['status'] == 'stereo_mismatch'
    assert compare('NCC(=O)O', observed)['status'] == 'graph_mismatch'
    assert compare('CC(O)C(=O)O', Chem.MolFromSmiles('CC(O)C(=O)O'))['status'] == 'stereo_incomplete'
    assert compare('[NH3+]CC(=O)O', Chem.MolFromSmiles('NCC(=O)O'))['status'] == 'graph_mismatch'
    assert compare('[13CH3]O', Chem.MolFromSmiles('CO'))['status'] == 'graph_mismatch'


def test_aromatic_hydrogen_and_atom_order():
    smiles = 'c1cc[nH]c1'
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    mol = Chem.RenumberAtoms(mol, list(reversed(range(mol.GetNumAtoms()))))
    result = compare(smiles, mol)
    assert result['status'] == 'identity_match'
    assert len(set(result['csv_atom_to_mol2_heavy_index'])) == 5


def test_csv_join_uses_name_not_position_and_keeps_duplicates(tmp_path):
    paths = [tmp_path / 'one.csv', tmp_path / 'two.csv']
    for path, rows in zip(paths, [[('CO', 'B'), ('CC', 'A')], [('CN', 'A')]]):
        with path.open('w', newline='', encoding='utf-8') as stream:
            writer = csv.writer(stream)
            writer.writerow(['SMILES', 'Name'])
            writer.writerows(rows)
    entries, metadata = read_names(paths, {'A'})
    assert set(entries) == {'A'}
    assert [r['smiles'] for r in entries['A']] == ['CC', 'CN']
    assert [r['csv_row'] for r in entries['A']] == [3, 2]
    assert [r['rows'] for r in metadata] == [2, 1]
