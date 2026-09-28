"""Join source SMILES to bounded MOL2 prefixes by exact name, with identity checks."""
import argparse
from collections import Counter
import csv
import hashlib
import itertools
import json
from pathlib import Path

from .mol2 import iter_mol2_records, load_rdkit_mol2


def sha(path):
    value = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            value.update(block)
    return value.hexdigest()


def read_names(paths, wanted):
    """Keep all matching rows, including duplicates; do not assume shard ordering."""
    entries = {}; metadata = []
    for path in map(Path, paths):
        count = 0
        with path.open(encoding='utf-8-sig', newline='') as stream:
            reader = csv.DictReader(stream)
            if not {'SMILES', 'Name'} <= set(reader.fieldnames or []):
                raise ValueError(f'Missing SMILES/Name columns: {path}')
            for row_number, row in enumerate(reader, 2):
                count += 1
                if row['Name'] in wanted:
                    entries.setdefault(row['Name'], []).append(dict(
                        csv_path=str(path.resolve()), csv_row=row_number, smiles=row['SMILES']))
        metadata.append(dict(path=str(path.resolve()), rows=count, sha256=sha(path)))
    return entries, metadata


def compare(smiles, mol):
    from rdkit import Chem
    reference = Chem.MolFromSmiles(smiles)
    if reference is None:
        return dict(status='invalid_smiles')
    reference = Chem.RemoveHs(reference)
    observed = Chem.RemoveHs(Chem.Mol(mol))
    # Atom-map numbers are provenance annotations, not chemical identity.
    for graph in (reference, observed):
        for atom in graph.GetAtoms():
            atom.SetAtomMapNum(0)
    def canonical(graph, stereo=True):
        return Chem.MolToSmiles(graph, canonical=True, isomericSmiles=stereo)
    # Removing stereo separately preserves isotope identity in the graph check.
    plain_reference = Chem.Mol(reference); plain_observed = Chem.Mol(observed)
    Chem.RemoveStereochemistry(plain_reference); Chem.RemoveStereochemistry(plain_observed)
    graph_equal = canonical(plain_reference) == canonical(plain_observed)
    stereo_equal = canonical(reference) == canonical(observed)
    def unspecified(graph):
        return sum(str(info.specified) == 'Unspecified' for info in Chem.FindPotentialStereo(graph))
    unknown_reference = unspecified(reference); unknown_observed = unspecified(observed)
    result = dict(status='graph_mismatch' if not graph_equal else
        ('stereo_mismatch' if not stereo_equal else
         ('stereo_incomplete' if unknown_reference or unknown_observed else 'identity_match')),
        heavy_graph_equal=graph_equal, isomeric_smiles_equal=stereo_equal,
        csv_canonical_smiles=canonical(reference), mol2_canonical_smiles=canonical(observed),
        csv_unspecified_stereo=unknown_reference, mol2_unspecified_stereo=unknown_observed)
    # Equal whole graphs allow a full atom correspondence. Symmetry may yield
    # alternatives; this is one valid correspondence, not a unique BB assignment.
    if graph_equal and stereo_equal:
        match = observed.GetSubstructMatch(reference, useChirality=True)
        if len(match) == reference.GetNumAtoms() == observed.GetNumAtoms():
            result['csv_atom_to_mol2_heavy_index'] = list(match)
            result['mapping_scope'] = 'One graph correspondence; symmetry uniqueness not established'
    return result


def audit(mol2_paths, csv_paths, output, per_file=1000):
    if not 1 <= per_file <= 10000:
        raise ValueError('per_file must be 1..10000')
    output = Path(output)
    if output.exists():
        raise ValueError('Use a new output directory')
    records = [record for path in map(Path, mol2_paths)
               for record in itertools.islice(iter_mol2_records(path), per_file)]
    entries, csv_metadata = read_names(csv_paths, {r.molecule_name for r in records})
    output.mkdir(parents=True)
    counts = Counter(); examples = {}; matched_names = set()
    with (output / 'identities.jsonl').open('w', encoding='utf-8') as stream:
        for record in records:
            row = dict(source_name=record.molecule_name, conformer_name=record.name,
                source_path=str(record.source_path), source_record_index=record.record_index,
                content_sha256=record.content_sha256)
            matches = entries.get(record.molecule_name, [])
            row['csv_matches'] = matches
            if not matches:
                row['status'] = 'name_not_found'
            elif len(matches) != 1:
                row['status'] = 'duplicate_name_requires_review'
            else:
                matched_names.add(record.molecule_name)
                try:
                    mol, mode = load_rdkit_mol2(record.raw_text, record.name)
                    row['sanitization'] = mode
                    if mode != 'strict':
                        row['status'] = 'non_strict_mol2_requires_review'
                    else:
                        row.update(compare(matches[0]['smiles'], mol))
                        atoms = record.raw_text.split('@<TRIPOS>ATOM', 1)[1].split('@<TRIPOS>', 1)[0]
                        ids = [int(line.split()[0]) for line in atoms.splitlines() if line.strip()]
                        row['mol2_heavy_source_atom_ids'] = [ids[a.GetIdx()] for a in mol.GetAtoms() if a.GetAtomicNum() > 1]
                except ValueError as exc:
                    row.update(status='parse_or_comparison_error', error=str(exc))
            counts[row['status']] += 1
            examples.setdefault(row['status'], [])
            if len(examples[row['status']]) < 5:
                examples[row['status']].append(record.name)
            stream.write(json.dumps(row) + '\n')
    report = dict(status='complete', scope='MOL2 prefixes joined against all supplied CSV rows by exact name',
        conformers=len(records), molecule_names=len({r.molecule_name for r in records}),
        uniquely_joined_names=len(matched_names), counts=dict(counts), examples=examples,
        csv_sources=csv_metadata, per_file=per_file,
        mol2_sources=[dict(path=str(Path(p).resolve()), bytes=Path(p).stat().st_size,
                           mtime_ns=Path(p).stat().st_mtime_ns) for p in mol2_paths],
        code_sha256=sha(__file__), output_sha256=sha(output / 'identities.jsonl'),
        limitations=['No source modification, building-block label assignment or main-ring recovery',
                    'Not a uniform sample or full-library validation; MOL2 hashes are per sampled record',
                    'SMILES alone does not define the geometric cis/trans state of every amide'])
    (output / 'report.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mol2', nargs='+', type=Path, required=True)
    parser.add_argument('--csv', nargs='+', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--per-file', type=int, default=1000)
    args = parser.parse_args()
    print(json.dumps(audit(args.mol2, args.csv, args.output, args.per_file), indent=2))


if __name__ == '__main__':
    main()
