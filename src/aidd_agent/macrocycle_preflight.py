"""Freeze source inventory and prepare audit-backed conformer block admission."""
import argparse
from collections import Counter
import csv
import json
from pathlib import Path
import sqlite3

from .macrocycle_full_audit import discover
from .macrocycle_identity_audit import sha


def fresh_output(output, source):
    output, source = Path(output).resolve(), Path(source).resolve()
    if output == source or source in output.parents:
        raise ValueError('Output must be outside the input tree')
    if output.exists():
        raise FileExistsError('Use a new output directory')
    return output


def write_report(output, report):
    report['implementation_sha256'] = sha(Path(__file__))
    (output / 'report.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
    return report


def inventory(source, output):
    source = Path(source).resolve()
    output = fresh_output(output, source)
    csv_paths, mol_paths = discover(source)
    records = []
    for path in sorted(csv_paths + mol_paths):
        before = path.stat()
        checksum = sha(path)
        after = path.stat()
        if (before.st_size, before.st_mtime_ns) != (after.st_size, after.st_mtime_ns):
            raise ValueError('Source changed during inventory: ' + str(path))
        row = dict(path=str(path), relative_path=path.relative_to(source).as_posix(),
                   bytes=after.st_size, mtime_ns=after.st_mtime_ns, sha256=checksum)
        if path in csv_paths:
            with path.open(encoding='utf-8-sig', newline='') as stream:
                fields = csv.DictReader(stream).fieldnames or []
            row['required_columns_present'] = {'Name', 'SMILES'} <= set(fields)
        records.append(row)
    output.mkdir(parents=True)
    return write_report(output, dict(status='inventory_complete', source_dir=str(source),
        scope='Discovered files only; no chemical audit or claim of complete workstation coverage',
        files=records, mol2_files=len(mol_paths), csv_files=len(csv_paths),
        mol2_without_csv=[str(p) for p in mol_paths if p.with_suffix('.csv') not in csv_paths],
        undiscovered_csv=[str(p.resolve()) for p in sorted(source.rglob('*'))
                          if p.is_file() and p.suffix.lower() == '.csv' and p.resolve() not in csv_paths],
        source_bytes=sum(r['bytes'] for r in records),
        invalid_csv_headers=[r['path'] for r in records if r.get('required_columns_present') is False],
        limits=['File hashing performs one full byte read; conformer counts are not estimated',
                'Inventory does not lock sources; compare source hashes with the completed audit']))


def block_readiness(audit, output):
    audit = Path(audit).resolve()
    output = fresh_output(output, audit)
    receipt = json.loads((audit / 'report.json').read_text(encoding='utf-8'))
    if receipt.get('status') != 'complete':
        raise ValueError('Audit is not complete')
    for name in ('audit.sqlite', 'issues.jsonl'):
        expected = receipt.get('output_hashes', {}).get(name)
        if not expected or sha(audit / name) != expected:
            raise ValueError('Audit artifact hash mismatch: ' + name)
    counts = Counter(); flags = Counter()
    with sqlite3.connect((audit / 'audit.sqlite').as_uri() + '?mode=ro', uri=True) as db:
        # Disk-backed temporary sets avoid collecting every conflicting name in Python.
        db.execute('PRAGMA temp_store=FILE')
        db.execute('CREATE TEMP TABLE duplicates AS SELECT conf_name FROM conformer GROUP BY conf_name HAVING count(*)>1')
        db.execute('CREATE UNIQUE INDEX duplicate_key ON duplicates(conf_name)')
        db.execute("CREATE TEMP TABLE conflicts AS SELECT name FROM conformer GROUP BY name HAVING count(DISTINCT json_extract(payload,'$.mol2_canonical_smiles'))>1")
        db.execute('CREATE UNIQUE INDEX conflict_key ON conflicts(name)')
        output.mkdir(parents=True)
        with (output / 'conformer_readiness.jsonl').open('w', encoding='utf-8') as stream:
            query = '''SELECT c.payload, d.conf_name IS NOT NULL, x.name IS NOT NULL
                FROM conformer c LEFT JOIN duplicates d ON d.conf_name=c.conf_name
                LEFT JOIN conflicts x ON x.name=c.name ORDER BY c.rowid'''
            for payload, duplicate, conflict in db.execute(query):
                row = json.loads(payload)
                reasons = []
                if row['status'] != 'ok':
                    reasons.append(row.get('reason', 'audit_review'))
                if duplicate:
                    reasons.append('duplicate_conformer_name')
                if conflict:
                    reasons.append('same_name_chemistry_conflict')
                if not row.get('ring_source_atom_ids'):
                    reasons.append('missing_peptide_atom_map')
                state = 'review' if reasons else 'geometry_candidate'
                row_flags = []
                if row.get('name_mapping') != 'verified_against_mol2':
                    row_flags.append('cross_conformer_alignment_unverified')
                if (row.get('name_rotations') or 0) > 1:
                    row_flags.append('equivalent_cyclic_rotations')
                if 'boundary' in row.get('omega_states', []):
                    row_flags.append('boundary_omega')
                counts[state] += 1
                flags.update(row_flags)
                entry = {key: row.get(key) for key in (
                    'name', 'conformer_name', 'conformer_index', 'path', 'record_index',
                    'content_sha256', 'verification_source', 'csv_identity_verified',
                    'ring_source_atom_ids', 'name_rotations', 'omega_states')}
                entry.update(state=state, review_reasons=reasons, flags=row_flags,
                             registry_join='pending', production_block_assignment=False)
                stream.write(json.dumps(entry) + '\n')
    eligible = counts['geometry_candidate']
    return write_report(output, dict(status='preparation_complete',
        audit_report_sha256=sha(audit / 'report.json'), audit_dir=str(audit),
        audit_artifact_hashes=receipt['output_hashes'], counts=dict(counts), flags=dict(flags),
        capacity_planning=[dict(capacity=n, minimum_blocks_ignoring_strata=(eligible+n-1)//n)
                           for n in (10000, 20000, 30000)],
        manifest_sha256=sha(output / 'conformer_readiness.jsonl'),
        production_ready=False, limitations=[
            'No registry identities invented; source-to-registry hash and atom-order join pending',
            'No descriptor extraction or clustering performed; capacity counts are lower bounds',
            'Flags preserve geometry but do not establish unique residue alignment',
            'Missing conformer suffixes do not discard existing conformers']))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--source-dir', type=Path)
    group.add_argument('--audit', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    result = inventory(args.source_dir, args.output) if args.source_dir else block_readiness(args.audit, args.output)
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
