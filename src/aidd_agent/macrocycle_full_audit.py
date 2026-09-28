"""Disk-backed full-source audit; all supplied CSV rows and MOL2 conformers."""
import argparse
from collections import Counter, deque
from concurrent.futures import ProcessPoolExecutor
import csv
import json
from pathlib import Path
import sqlite3
import time

import numpy as np
from rdkit import Chem, rdBase
from .macrocycle_identity_audit import sha
from .macrocycle_peptide import peptide
from .macrocycle_blocks import angle, amide_state
from .mol2 import iter_mol2_records, load_rdkit_mol2


def canonical(mol):
    graph = Chem.RemoveHs(Chem.Mol(mol))
    for atom in graph.GetAtoms():
        atom.SetAtomMapNum(0)
    iso = Chem.MolToSmiles(graph, isomericSmiles=True)
    unknown = sum(str(i.specified) == 'Unspecified' for i in Chem.FindPotentialStereo(graph))
    Chem.RemoveStereochemistry(graph)
    return iso, Chem.MolToSmiles(graph, isomericSmiles=True), unknown


def csv_one(item):
    path, row, name, smiles = item
    out = dict(path=path, row=row, name=name, smiles=smiles)
    try:
        with rdBase.BlockLogs():
            mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError('invalid_smiles')
        iso, plain, unknown = canonical(mol)
        out.update(iso=iso, plain=plain, unspecified_stereo=unknown)
        out.update(peptide(mol, name))
        out['status'] = 'ok' if not unknown else 'unspecified_stereo'
    except (ValueError, RuntimeError) as exc:
        out.update(status='review', reason=str(exc))
    return out


def mol_one(item):
    record, reference = item[:2]
    allow_missing = item[2] if len(item) > 2 else False
    out = dict(name=record.molecule_name, conformer_name=record.name,
        conformer_index=record.conformer_index, path=str(record.source_path),
        record_index=record.record_index, content_sha256=record.content_sha256)
    out['verification_source'] = 'csv_and_mol2' if reference is not None else 'mol2_only'
    try:
        if reference is None and not allow_missing:
            raise ValueError('name_missing_from_csv')
        if reference is not None and reference.get('duplicate'):
            raise ValueError('duplicate_csv_name')
        if reference is not None and reference['status'] != 'ok':
            raise ValueError('csv_requires_review:' + reference.get('reason', reference['status']))
        mol, mode = load_rdkit_mol2(record.raw_text, record.name)
        if mode != 'strict':
            raise ValueError('non_strict_mol2')
        iso, plain, unknown = canonical(mol)
        out.update(mol2_canonical_smiles=iso, csv_identity_verified=False)
        if reference is not None and plain != reference['plain']:
            raise ValueError('chemical_graph_mismatch')
        if reference is not None and iso != reference['iso']:
            raise ValueError('stereochemistry_mismatch')
        if unknown:
            raise ValueError('unspecified_stereo')
        out['csv_identity_verified'] = reference is not None
        try:
            mapped = peptide(mol, record.molecule_name)
            out['name_mapping'] = 'verified_against_mol2'
        except ValueError as exc:
            if reference is not None:
                raise
            mapped = peptide(mol)
            out.update(name_mapping='unverified', name_mapping_reason=str(exc),
                       alignment_scope='Graph-only cycle; cross-conformer residue alignment not established')
        heavy = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() > 1]
        points = np.asarray(mol.GetConformer().GetPositions())[heavy]
        if not np.isfinite(points).all():
            raise ValueError('nonfinite_coordinates')
        atom_lines = record.raw_text.split('@<TRIPOS>ATOM', 1)[1].split('@<TRIPOS>', 1)[0]
        source_ids = [int(line.split()[0]) for line in atom_lines.splitlines() if line.strip()]
        if len(source_ids) != mol.GetNumAtoms() or len(set(source_ids)) != len(source_ids):
            raise ValueError('invalid_source_atom_ids')
        units = mapped['units']; omega = []
        for i, (_, ca, carbon, _) in enumerate(units):
            nxt = units[(i + 1) % len(units)]
            omega.append(float(np.degrees(angle(points[[ca, carbon, nxt[0], nxt[1]]]))))
        out.update(status='ok', ring_size=len(mapped['ring_atoms']),
            name_rotations=mapped['name_rotations'], omega_degrees=omega,
            omega_states=[amide_state(a) for a in omega],
            ring_source_atom_ids=[source_ids[heavy[i]] for i in mapped['ring_atoms']])
    except (ValueError, RuntimeError) as exc:
        out.update(status='review', reason=str(exc))
    return out


def batch_worker(payload):
    kind, rows = payload
    function = csv_one if kind == 'csv' else mol_one
    return [function(row) for row in rows]


def bounded(pool, kind, iterable, workers):
    pending = deque(); batch = []
    for row in iterable:
        batch.append(row)
        if len(batch) == 100:
            pending.append(pool.submit(batch_worker, (kind, batch))); batch = []
            if len(pending) >= workers * 2:
                yield from pending.popleft().result()
    if batch:
        pending.append(pool.submit(batch_worker, (kind, batch)))
    while pending:
        yield from pending.popleft().result()


def discover(source):
    source = Path(source).resolve()
    if not source.is_dir():
        raise ValueError('Source directory does not exist')
    mol2 = sorted(p.resolve() for p in source.rglob('*') if p.is_file() and p.suffix.lower() == '.mol2')
    if not mol2:
        raise ValueError('No MOL2 files found')
    csv_paths = sorted({p.with_suffix('.csv') for p in mol2 if p.with_suffix('.csv').is_file()})
    return csv_paths, mol2


def run(csv_paths, mol2_paths, output, workers=4, allow_missing_csv=False):
    if not 1 <= workers <= 32:
        raise ValueError('workers must be 1..32')
    csv_paths = sorted({Path(p).resolve() for p in csv_paths})
    mol2_paths = sorted({Path(p).resolve() for p in mol2_paths})
    if not mol2_paths or any(not p.is_file() for p in csv_paths + mol2_paths):
        raise ValueError('No MOL2 inputs or missing input file')
    output = Path(output)
    output.mkdir(parents=True, exist_ok=False)
    started = time.monotonic(); counts = Counter(); csv_counts = Counter(); omega_counts = Counter()
    report = dict(status='running', scope='All records in supplied sources, not unseen library shards',
                  workers=workers, csv_sources=[], mol2_sources=[], allow_missing_csv=allow_missing_csv,
                  input_manifest=dict(csv=[str(p) for p in csv_paths], mol2=[str(p) for p in mol2_paths],
                    mol2_without_same_stem_csv=[str(p) for p in mol2_paths if p.with_suffix('.csv') not in csv_paths]),
                  code_hashes={p.name:sha(p) for p in
                    [Path(__file__), Path(__file__).with_name('macrocycle_peptide.py'),
                     Path(__file__).with_name('macrocycle_blocks.py'), Path(__file__).with_name('mol2.py')]})
    def checkpoint(stage):
        report.update(stage=stage, csv_counts=dict(csv_counts), conformer_counts=dict(counts),
                      wall_seconds=time.monotonic() - started)
        temp = output / 'report.tmp'
        temp.write_text(json.dumps(report, indent=2), encoding='utf-8')
        temp.replace(output / 'report.json')
    checkpoint('csv')
    db = sqlite3.connect(output / 'audit.sqlite')
    db.executescript('''CREATE TABLE molecule(name TEXT PRIMARY KEY, occurrences INTEGER, payload TEXT);
        CREATE TABLE conformer(name TEXT, conf_index INTEGER, conf_name TEXT, status TEXT, reason TEXT, payload TEXT);
        CREATE INDEX conf_name ON conformer(name);''')
    try:
        with ProcessPoolExecutor(max_workers=workers) as pool, (output / 'issues.jsonl').open('w', encoding='utf-8') as issues:
            def csv_rows():
                for path in map(Path, csv_paths):
                    count = 0
                    with path.open(encoding='utf-8-sig', newline='') as stream:
                        reader = csv.DictReader(stream)
                        if not {'Name', 'SMILES'} <= set(reader.fieldnames or []):
                            raise ValueError('missing_csv_columns')
                        for rownum, row in enumerate(reader, 2):
                            count += 1
                            yield str(path.resolve()), rownum, row['Name'], row['SMILES']
                    report['csv_sources'].append(dict(path=str(path.resolve()), rows=count, sha256=sha(path)))
            for i, row in enumerate(bounded(pool, 'csv', csv_rows(), workers), 1):
                csv_counts[row['status'] if row['status'] == 'ok' else row.get('reason', row['status'])] += 1
                previous = db.execute('SELECT occurrences FROM molecule WHERE name=?', (row['name'],)).fetchone()
                if previous:
                    db.execute('UPDATE molecule SET occurrences=occurrences+1 WHERE name=?', (row['name'],))
                    issues.write(json.dumps(dict(kind='duplicate_csv_name', **row)) + '\n')
                else:
                    db.execute('INSERT INTO molecule VALUES(?,1,?)', (row['name'], json.dumps(row)))
                if row['status'] != 'ok':
                    issues.write(json.dumps(dict(kind='csv', **row)) + '\n')
                if i % 10000 == 0:
                    db.commit(); checkpoint('csv'); print(f'CSV rows {i}', flush=True)
            db.commit(); checkpoint('mol2')
            def mol_rows():
                last_name = None; cached = None
                for path in map(Path, mol2_paths):
                    count = 0
                    for record in iter_mol2_records(path):
                        count += 1
                        if record.molecule_name != last_name:
                            entry = db.execute('SELECT occurrences,payload FROM molecule WHERE name=?', (record.molecule_name,)).fetchone()
                            cached = json.loads(entry[1]) if entry else None
                            if entry and entry[0] != 1:
                                cached['duplicate'] = True
                            last_name = record.molecule_name
                        yield record, cached, allow_missing_csv
                    report['mol2_sources'].append(dict(path=str(path.resolve()), records=count, sha256=sha(path)))
            for i, row in enumerate(bounded(pool, 'mol2', mol_rows(), workers), 1):
                counts[row.get('reason', row['status'])] += 1
                omega_counts.update(row.get('omega_states', []))
                db.execute('INSERT INTO conformer VALUES(?,?,?,?,?,?)',
                    (row['name'], row['conformer_index'], row['conformer_name'], row['status'], row.get('reason'), json.dumps(row)))
                if row['status'] != 'ok':
                    issues.write(json.dumps(dict(kind='mol2', **row)) + '\n')
                if i % 10000 == 0:
                    db.commit(); checkpoint('mol2'); print(f'MOL2 records {i}; issues {i-counts["ok"]}', flush=True)
            db.commit()
            report['unique_csv_names'] = db.execute('SELECT count(*) FROM molecule').fetchone()[0]
            report['duplicate_csv_names'] = db.execute('SELECT count(*) FROM molecule WHERE occurrences>1').fetchone()[0]
            report['csv_names_without_mol2'] = db.execute('SELECT count(*) FROM molecule m WHERE NOT EXISTS(SELECT 1 FROM conformer c WHERE c.name=m.name)').fetchone()[0]
            report['conformer_count_histogram'] = dict(db.execute('SELECT n,count(*) FROM (SELECT count(*) n FROM conformer GROUP BY name) GROUP BY n').fetchall())
            report['duplicate_conformer_name_groups'] = db.execute('SELECT count(*) FROM (SELECT conf_name FROM conformer GROUP BY conf_name HAVING count(*)>1)').fetchone()[0]
            report['non_123_suffix_groups'] = db.execute('SELECT count(*) FROM (SELECT name FROM conformer GROUP BY name HAVING count(*)!=3 OR count(DISTINCT conf_index)!=3 OR min(conf_index)!=1 OR max(conf_index)!=3)').fetchone()[0]
            for name, count in db.execute('SELECT name,count(*) FROM conformer GROUP BY name HAVING count(*)!=3 OR count(DISTINCT conf_index)!=3 OR min(conf_index)!=1 OR max(conf_index)!=3'):
                issues.write(json.dumps(dict(kind='non_123_suffix_group', name=name, records=count)) + '\n')
            report['ring_size_histogram'] = dict(db.execute("SELECT json_extract(payload,'$.ring_size'),count(*) FROM conformer WHERE status='ok' GROUP BY 1").fetchall())
            report['ambiguous_name_rotation_conformers'] = db.execute("SELECT count(*) FROM conformer WHERE status='ok' AND json_extract(payload,'$.name_rotations')>1").fetchone()[0]
            report['names_with_multiple_omega_patterns'] = db.execute("SELECT count(*) FROM (SELECT name FROM conformer WHERE status='ok' AND json_extract(payload,'$.name_mapping')='verified_against_mol2' GROUP BY name HAVING count(DISTINCT json_extract(payload,'$.omega_states'))>1)").fetchone()[0]
            report['omega_pattern_histogram'] = dict(db.execute("SELECT json_extract(payload,'$.omega_states'),count(*) FROM conformer WHERE status='ok' GROUP BY 1").fetchall())
            for (name,) in db.execute('SELECT name FROM molecule m WHERE NOT EXISTS(SELECT 1 FROM conformer c WHERE c.name=m.name)'):
                issues.write(json.dumps(dict(kind='csv_name_without_mol2', name=name)) + '\n')
            report['omega_counts'] = dict(omega_counts)
            report['verification_source_counts'] = dict(db.execute("SELECT json_extract(payload,'$.verification_source'),count(*) FROM conformer GROUP BY 1").fetchall())
            report['mol2_only_passed'] = db.execute("SELECT count(*) FROM conformer WHERE status='ok' AND json_extract(payload,'$.verification_source')='mol2_only'").fetchone()[0]
            report['unverified_name_mapping_conformers'] = db.execute("SELECT count(*) FROM conformer WHERE json_extract(payload,'$.name_mapping')='unverified'").fetchone()[0]
            report['same_name_multiple_mol2_chemistries'] = db.execute("SELECT count(*) FROM (SELECT name FROM conformer GROUP BY name HAVING count(DISTINCT json_extract(payload,'$.mol2_canonical_smiles'))>1)").fetchone()[0]
            for (name,) in db.execute("SELECT name FROM conformer GROUP BY name HAVING count(DISTINCT json_extract(payload,'$.mol2_canonical_smiles'))>1"):
                issues.write(json.dumps(dict(kind='same_name_multiple_mol2_chemistries', name=name)) + '\n')
        report.update(status='complete',
            limitations=['No production registry changes or full-library coverage beyond supplied files',
                'Name rotation may be ambiguous for repeated identical residues',
                'Geometric omega bins are exploratory, not energy barriers; block clustering not run'])
        checkpoint('complete')
    except BaseException as exc:
        db.commit(); report.update(status='failed', error=str(exc)); checkpoint('failed'); raise
    finally:
        db.close()
    report['output_hashes'] = {p.name:sha(p) for p in (output / 'audit.sqlite', output / 'issues.jsonl')}
    checkpoint('complete')
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source-dir', type=Path, help='Recursively discover MOL2 and optional same-stem CSV files')
    parser.add_argument('--csv', nargs='+', type=Path, default=[])
    parser.add_argument('--mol2', nargs='+', type=Path, default=[])
    parser.add_argument('--allow-missing-csv', action='store_true')
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--workers', type=int, default=4)
    args = parser.parse_args()
    if args.source_dir:
        if args.csv or args.mol2:
            parser.error('Use --source-dir or explicit file lists, not both')
        args.csv, args.mol2 = discover(args.source_dir)
    print(json.dumps(run(args.csv, args.mol2, args.output, args.workers,
                        args.allow_missing_csv or bool(args.source_dir)), indent=2))


if __name__ == '__main__':
    main()
