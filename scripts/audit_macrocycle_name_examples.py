"""Check supplied naming examples against backbone N methylation and alpha stereo."""
import argparse
import csv
import json
from pathlib import Path

from rdkit import Chem
from aidd_agent.mol2 import iter_mol2_records, load_rdkit_mol2


def units(mol):
    Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
    query = Chem.MolFromSmarts('[N]-[C;X4]-[C](=[O])')
    result = []
    for n, ca, c, _ in mol.GetSubstructMatches(query):
        atom = mol.GetAtomWithIdx(n)
        methyl = [a.GetIdx() for a in atom.GetNeighbors()
                  if a.GetAtomicNum() == 6 and a.GetTotalNumHs(includeNeighbors=True) == 3
                  and sum(b.GetAtomicNum() > 1 for b in a.GetNeighbors()) == 1]
        alpha = mol.GetAtomWithIdx(ca)
        result.append(dict(N=n, CA=ca, C=c,
            CA_CIP=alpha.GetProp('_CIPCode') if alpha.HasProp('_CIPCode') else None,
            N_H=atom.GetTotalNumHs(includeNeighbors=True), N_methyl_atoms=methyl))
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise ValueError('Use a new output file')
    names = {'c--A-W-Lnme-A-dLnme-dVNMe-c', 'c--A-W-Lnme-A-dLnme-VNMe-c',
             'c--A-W-Lnme-AdFnme-Wnme-c', 'c--A-W-Lnme-AdFnme-dWnme-c'}
    output = []
    for path in sorted(args.source.glob('split_000*.csv')):
        with path.open(encoding='utf-8-sig', newline='') as stream:
            for row_number, row in enumerate(csv.DictReader(stream), 2):
                if row['Name'] in names:
                    output.append(dict(name=row['Name'], source=str(path), csv_row=row_number,
                        smiles=row['SMILES'], units=units(Chem.MolFromSmiles(row['SMILES']))))
    record = next(iter_mol2_records(args.source / 'split_0001.mol2'))
    mol, mode = load_rdkit_mol2(record.raw_text, record.name)
    atoms = record.raw_text.split('@<TRIPOS>ATOM')[1].split('@<TRIPOS>')[0]
    ids = [int(line.split()[0]) for line in atoms.splitlines() if line.strip()]
    observed = units(mol)
    for row in observed:
        for key in ('N', 'CA', 'C'):
            row[key + '_mol2_id'] = ids[row[key]]
        row['methyl_mol2_ids'] = [ids[i] for i in row['N_methyl_atoms']]
    output.append(dict(name=record.name, source=str(record.source_path),
        content_sha256=record.content_sha256, sanitization=mode, units=observed))
    args.output.write_text(json.dumps(output, indent=2), encoding='utf-8')
    print(json.dumps(output, indent=2))


if __name__ == '__main__':
    main()
