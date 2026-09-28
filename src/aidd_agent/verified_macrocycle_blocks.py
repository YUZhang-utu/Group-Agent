"""Audit-backed source extraction, frozen capacity blocks and sampling preparation."""
import argparse
from collections import Counter
import json
from pathlib import Path
import sqlite3
import time

from rdkit import RDConfig, rdBase

from .conformer_block_store import fit, sample
from .macrocycle_blocks import digest
from .macrocycle_descriptors import describe_peptide, VERSION, VARIANTS, FAMILIES
from .macrocycle_identity_audit import sha
from .macrocycle_preflight import block_readiness, fresh_output
from .mol2 import iter_mol2_records, load_rdkit_mol2


def registry_identity(db, library, record):
    matches = db.execute('''SELECT c.id,c.molecule_id FROM conformer c
        JOIN molecule m ON m.id=c.molecule_id WHERE m.library_id=? AND
        m.source_name=? AND c.source_record_name=? AND c.conformer_index=? AND
        c.content_sha256=? AND c.topology_sha256=? AND c.source_path=? AND
        c.source_record_index=? AND c.atom_count=? AND c.bond_count=?''',
        (library, record.molecule_name, record.name, record.conformer_index,
         record.content_sha256, record.topology_sha256, str(record.source_path.resolve()),
         record.record_index, record.atom_count, record.bond_count)).fetchall()
    if len(matches) != 1:
        raise ValueError('Registry provenance match is absent or ambiguous')
    return matches[0]


def run(audit, output, variant='typed', capacities=(10000, 20000, 30000),
        registry=None, library=None, sample_sizes=(100, 500)):
    if variant not in VARIANTS or not capacities or len(set(capacities)) != len(capacities):
        raise ValueError('Invalid variant or repeated/empty capacities')
    if any(type(n) is not int or n < 1 for n in capacities):
        raise ValueError('Positive integer capacities required')
    if bool(registry) != bool(library):
        raise ValueError('Registry and library ID must be supplied together')
    audit = Path(audit).resolve(); output = fresh_output(output, audit)
    receipt = json.loads((audit / 'report.json').read_text(encoding='utf-8'))
    # Raw source byte validation is independent of stored audit database validation.
    sources = receipt.get('mol2_sources', [])
    if not sources:
        raise ValueError('Audit contains no source files')
    for source in sources:
        path = Path(source['path']).resolve()
        if output == path.parent or path.parent in output.parents:
            raise ValueError('Output must be outside source trees')
        if sha(path) != source['sha256']:
            raise ValueError('MOL2 source hash differs from audit: ' + str(path))
    output.mkdir(parents=True)
    block_readiness(audit, output / 'admission')
    schema = dict(version=VERSION, variant=variant, rdkit=rdBase.rdkitVersion,
                  feature_families=FAMILIES,
                  feature_definition_sha256=sha(Path(RDConfig.RDDataDir) / 'BaseFeatures.fdef'),
                  scales=dict(side_atom_count=10, feature_count=4, local_coordinates_angstrom=5))
    started = time.monotonic(); counts = Counter()
    identity_scope = 'production_registry' if registry else 'source_only_not_production_ids'
    registry_db = sqlite3.connect(Path(registry).resolve().as_uri()+'?mode=ro', uri=True) if registry else None
    try:
        with sqlite3.connect(output / 'descriptors.sqlite') as db:
            db.executescript('''CREATE TABLE admission(path TEXT, record_index INTEGER,
                payload TEXT, visited INTEGER DEFAULT 0, PRIMARY KEY(path,record_index));
                CREATE TABLE descriptor(cid TEXT PRIMARY KEY, payload TEXT);''')
            with (output / 'admission/conformer_readiness.jsonl').open(encoding='utf-8') as stream:
                for line in stream:
                    row = json.loads(line)
                    db.execute('INSERT INTO admission(path,record_index,payload) VALUES(?,?,?)',
                               (row['path'], row['record_index'], line))
            db.commit()
            with (output / 'review.jsonl').open('w', encoding='utf-8') as review:
                for source in sources:
                    source_count = 0
                    for record in iter_mol2_records(Path(source['path'])):
                        source_count += 1; counts['source_conformers'] += 1
                        key = (str(record.source_path.resolve()), record.record_index)
                        found = db.execute('SELECT payload,visited FROM admission WHERE path=? AND record_index=?', key).fetchone()
                        if found is None or found[1]:
                            raise ValueError('Source record missing from audit or visited twice')
                        admitted = json.loads(found[0])
                        if (admitted['content_sha256'], admitted['conformer_name'], admitted['name']) != (
                                record.content_sha256, record.name, record.molecule_name):
                            raise ValueError('Source record identity/hash mismatch')
                        db.execute('UPDATE admission SET visited=1 WHERE path=? AND record_index=?', key)
                        try:
                            if admitted['state'] != 'geometry_candidate':
                                raise ValueError('audit_review:' + ','.join(admitted['review_reasons']))
                            mol, mode = load_rdkit_mol2(record.raw_text, record.name)
                            if mode != 'strict':
                                raise ValueError('Strict chemistry required')
                            name = None if 'cross_conformer_alignment_unverified' in admitted['flags'] else record.molecule_name
                            description = describe_peptide(mol, name, variant)
                            atoms = record.raw_text.split('@<TRIPOS>ATOM', 1)[1].split('@<TRIPOS>', 1)[0]
                            ids = [int(line.split()[0]) for line in atoms.splitlines() if line.strip()]
                            heavy_ids = [ids[a.GetIdx()] for a in mol.GetAtoms() if a.GetAtomicNum() > 1]
                            ring_ids = [heavy_ids[i] for i in description['ring_atoms']]
                            if set(ring_ids) != set(admitted['ring_source_atom_ids']):
                                raise ValueError('Peptide atom map differs from audited map')
                            if registry_db:
                                cid, mid = registry_identity(registry_db, library, record)
                            else:
                                cid = 'SOURCE-C-' + digest([*key, record.content_sha256])
                                mid = 'SOURCE-M-' + digest(record.molecule_name)
                            provenance = dict(source_name=record.molecule_name, source_record_name=record.name,
                                source_path=key[0], source_record_index=key[1], content_sha256=record.content_sha256,
                                identity_scope=identity_scope, library_id=library, flags=admitted['flags'],
                                csv_identity_verified=admitted['csv_identity_verified'],
                                verification_source=admitted['verification_source'],
                                heavy_atom_source_ids=heavy_ids, ring_source_atom_ids=ring_ids,
                                units=description['units'], typed_feature_atoms=description['typed_feature_atoms'],
                                omega_states=description['omega_states'], equivalent_rotations=description['equivalent_rotations'])
                            row = dict(conformer_id=cid, molecule_id=mid, descriptor=description['descriptor'],
                                       hard_group=description['hard_group'], provenance=provenance)
                            db.execute('INSERT INTO descriptor VALUES(?,?)', (cid, json.dumps(row)))
                            counts['described_conformers'] += 1
                        except ValueError as exc:
                            counts['review_conformers'] += 1
                            review.write(json.dumps(dict(path=key[0], record_index=key[1],
                                conformer_name=record.name, reason=str(exc)))+'\n')
                        if counts['source_conformers'] % 1000 == 0:
                            db.commit()
                            print(f"Extracted {counts['source_conformers']} source conformers", flush=True)
                    if source_count != source['records'] or sha(Path(source['path'])) != source['sha256']:
                        raise ValueError('Source changed while extracting descriptors')
                if db.execute('SELECT count(*) FROM admission WHERE visited=0').fetchone()[0]:
                    raise ValueError('Audit contains records absent from source traversal')
            db.commit()
            reports = []
            for capacity in capacities:
                rows = (json.loads(row[0]) for row in db.execute('SELECT payload FROM descriptor ORDER BY cid'))
                model = output / f'capacity-{capacity}'
                reports.append(fit(rows, model, capacity, schema))
                sample(model, output / f'samples-{capacity}', sizes=sample_sizes)
    finally:
        if registry_db:
            registry_db.close()
    result = dict(status='complete', schema=schema, counts=dict(counts), identity_scope=identity_scope,
        audit_report_sha256=sha(audit / 'report.json'), capacities=reports, wall_seconds=time.monotonic()-started,
        output_hashes={p: sha(output / p) for p in ('descriptors.sqlite', 'review.jsonl')},
        code_hashes={p: sha(Path(__file__).with_name(p)) for p in (
            'verified_macrocycle_blocks.py', 'macrocycle_descriptors.py', 'macrocycle_peptide.py',
            'conformer_block_store.py', 'macrocycle_preflight.py')},
        limitations=['Offline source-backed blocks, no production search dispatch or block rejection',
                    'Source descriptors use exact MOL2 coordinates; no quantized artifact-order assumption',
                    'Full-library speed, recall and feature weights require workstation benchmarking',
                    'Interrupted builds require a fresh output; incremental updates are separate proposals'])
    (output / 'report.json').write_text(json.dumps(result, indent=2), encoding='utf-8')
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--audit', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--variant', choices=VARIANTS, default='typed')
    parser.add_argument('--capacities', type=int, nargs='+', default=[10000, 20000, 30000])
    parser.add_argument('--registry', type=Path)
    parser.add_argument('--library-id')
    args = parser.parse_args()
    print(json.dumps(run(args.audit, args.output, args.variant, tuple(args.capacities),
                         args.registry, args.library_id), indent=2))


if __name__ == '__main__':
    main()
