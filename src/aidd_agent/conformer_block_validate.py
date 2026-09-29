"""Exhaustive membership checks and sampled raw-source geometry verification."""
import argparse
from collections import Counter
import hashlib
from itertools import zip_longest
import json
from pathlib import Path
import sqlite3

import numpy as np
from rdkit import Chem

from .conformer_block_store import open_model
from .macrocycle_descriptors import describe_peptide
from .macrocycle_identity_audit import sha
from .mol2 import iter_mol2_records, load_rdkit_mol2


def require(condition, message):
    if not condition:
        raise ValueError(message)


def readonly(path):
    return sqlite3.connect(Path(path).resolve().as_uri()+'?mode=ro', uri=True)


def validate_model(model, descriptors, schema):
    model = Path(model)
    db, report = open_model(model)
    try:
        require(report['schema'] == schema, 'Model/build descriptor schema mismatch')
        require(sha(model / 'memberships.jsonl') == report['output_hashes']['memberships.jsonl'],
                'Membership export hash mismatch')
        db.execute('ATTACH DATABASE ? AS input', (Path(descriptors).resolve().as_uri()+'?mode=ro',))
        metadata = {key: json.loads(value) for key, value in db.execute('SELECT key,value FROM metadata')}
        require(all(metadata[key] == report[key] for key in ('schema', 'model_id', 'capacity', 'version')),
                'Model metadata differs from receipt')
        trees = {}
        for row in db.execute('SELECT node,group_id,path,n,axis,cut,cut_id,left_node,right_node,centroid,radius FROM tree'):
            node, group, path, n, axis, cut, cut_id, left, right, center, radius = row
            trees[node] = dict(group=group, path=path, n=n, axis=axis, cut=cut,
                               cut_id=cut_id, left=left, right=right,
                               center=json.loads(center) if center else None, radius=radius)
        roots = {row['group']: node for node, row in trees.items() if row['path'] == ''}
        leaves = {node for node, row in trees.items() if row['axis'] is None}
        for node, row in trees.items():
            require(row['n'] > 0, 'Empty tree node')
            if node in leaves:
                require(row['n'] <= report['capacity'], 'Leaf exceeds capacity')
                require(row['radius'] is not None and np.isfinite(row['radius']) and row['radius'] >= 0,
                        'Invalid leaf radius')
            else:
                for key, bit in (('left', '0'), ('right', '1')):
                    child = trees.get(row[key])
                    require(child is not None and child['group'] == row['group'] and
                            child['path'] == row['path']+bit, 'Broken tree edge or cross-stratum child')
                require(trees[row['left']]['n']+trees[row['right']]['n'] == row['n'], 'Tree counts do not add up')
        node_counts = Counter(); fingerprint = hashlib.sha256(); count = 0
        query = '''SELECT p.cid,p.mid,p.group_id,p.vector,p.provenance,p.node,d.payload
            FROM point p LEFT JOIN input.descriptor d ON d.cid=p.cid ORDER BY p.cid'''
        with (model / 'memberships.jsonl').open(encoding='utf-8') as export:
            for pair, line in zip_longest(db.execute(query), export):
                require(pair is not None and line is not None, 'Membership export coverage mismatch')
                cid, mid, group, encoded, provenance, assigned, payload = pair
                require(payload is not None, 'Block contains a conformer absent from extraction')
                original = json.loads(payload); vector = np.asarray(json.loads(encoded))
                require((cid, mid, group) == (original['conformer_id'], original['molecule_id'], original['hard_group']),
                        'Descriptor identity or hard stratum changed')
                require(vector.ndim == 1 and len(vector) > 0 and np.isfinite(vector).all() and
                        np.array_equal(vector, original['descriptor']), 'Descriptor values differ from extraction')
                require(json.loads(provenance) == original.get('provenance', {}), 'Source provenance changed')
                expected_export = dict(conformer_id=cid, molecule_id=mid, block_id=assigned,
                                       provenance=json.loads(provenance))
                require(json.loads(line) == expected_export, 'Export differs from stored membership')
                require(assigned in leaves and trees[assigned]['group'] == group, 'Mixed stratum or non-leaf assignment')
                require(group in roots, 'Missing tree root')
                node = roots[group]
                for _ in range(len(trees)+1):
                    node_counts[node] += 1
                    branch = trees[node]
                    if branch['axis'] is None:
                        break
                    axis = branch['axis']
                    require(type(axis) is int and 0 <= axis < len(vector), 'Invalid split axis')
                    node = branch['left'] if (vector[axis], cid) < (branch['cut'], branch['cut_id']) else branch['right']
                else:
                    raise ValueError('Tree traversal cycle')
                require(node == assigned, 'Conformer assigned to the wrong frozen-tree leaf')
                fingerprint.update((json.dumps((cid, mid, group, encoded))+'\n').encode())
                count += 1
        require(db.execute('SELECT count(*) FROM input.descriptor').fetchone()[0] == count,
                'Extracted conformers missing from blocks')
        require(count > 0, 'No conformers assigned; an empty model cannot pass')
        require(all(node_counts[node] == row['n'] for node, row in trees.items()), 'Stored node population is incorrect')
        require(fingerprint.hexdigest() == report['descriptor_content_sha256'], 'Descriptor fingerprint mismatch')
        require(report['counts']['conformers'] == count and report['counts']['blocks'] == len(leaves),
                'Receipt counts differ from memberships')
        require(report['counts']['molecules'] == db.execute('SELECT count(DISTINCT mid) FROM point').fetchone()[0],
                'Receipt molecule count mismatch')
        maximum_radius = 0.0
        for node in sorted(leaves):
            n = 0; mean = None; farthest = 0.0
            center = np.asarray(trees[node]['center'])
            cursor = db.execute('SELECT vector FROM point WHERE node=? ORDER BY cid', (node,))
            while batch := cursor.fetchmany(1024):
                x = np.asarray([json.loads(row[0]) for row in batch])
                require(x.shape[1:] == center.shape and np.isfinite(center).all(), 'Invalid centroid dimension')
                avg = x.mean(axis=0)
                mean = avg if mean is None else mean+(avg-mean)*len(x)/(n+len(x))
                n += len(x)
                farthest = max(farthest, float(np.linalg.norm(x-center, axis=1).max()))
            require(np.allclose(mean, center, atol=1e-9, rtol=1e-9), 'Leaf centroid differs from members')
            require(np.isclose(farthest, trees[node]['radius'], atol=1e-8, rtol=1e-9), 'Leaf radius differs from members')
            maximum_radius = max(maximum_radius, farthest)
        return dict(capacity=report['capacity'], conformers_checked=count, blocks_checked=len(leaves),
                    max_block_radius=maximum_radius, structural_gate='passed')
    finally:
        db.close()


def independent_omega(points):
    """Plane-normal formula, independent of the descriptor's projected-vector formula."""
    p = np.asarray(points, dtype=float)
    first, bond, last = p[1]-p[0], p[2]-p[1], p[3]-p[2]
    n1, n2 = np.cross(first, bond), np.cross(bond, last)
    require(min(np.linalg.norm(n1), np.linalg.norm(n2), np.linalg.norm(bond)) > 1e-10,
            'Degenerate sampled omega geometry')
    return float(np.degrees(np.arctan2(np.dot(np.cross(n1, n2), bond/np.linalg.norm(bond)), np.dot(n1, n2))))


def validate_source_sample(build, schema, output, limit):
    require(type(limit) is int and limit > 0, 'Positive source sample limit required')
    with readonly(build / 'descriptors.sqlite') as db:
        # SOURCE IDs are hash-based; registry IDs are not necessarily random. This
        # is a deterministic inspection panel, not an unbiased population estimate.
        selected = [json.loads(row[0]) for row in db.execute('SELECT payload FROM descriptor ORDER BY cid LIMIT ?', (limit,))]
    wanted = {(row['provenance']['source_path'], row['provenance']['source_record_index']): row for row in selected}
    require(len(wanted) == len(selected) and bool(wanted), 'Invalid sampled source provenance')
    remaining = set(wanted); inspected = 0
    with (output / 'source_inspection.jsonl').open('w', encoding='utf-8') as stream:
        for source in sorted({path for path, _ in wanted}):
            for record in iter_mol2_records(Path(source)):
                key = (str(record.source_path.resolve()), record.record_index)
                if key not in wanted:
                    continue
                row = wanted[key]; provenance = row['provenance']
                require((record.content_sha256, record.name, record.molecule_name) == (
                    provenance['content_sha256'], provenance['source_record_name'], provenance['source_name']),
                    'Sampled raw-source identity mismatch')
                mol, mode = load_rdkit_mol2(record.raw_text, record.name)
                require(mode == 'strict', 'Sampled MOL2 is not strict chemistry')
                name = None if 'cross_conformer_alignment_unverified' in provenance['flags'] else record.molecule_name
                described = describe_peptide(mol, name, schema['variant'])
                atoms = record.raw_text.split('@<TRIPOS>ATOM', 1)[1].split('@<TRIPOS>', 1)[0]
                ids = [int(line.split()[0]) for line in atoms.splitlines() if line.strip()]
                heavy_ids = [ids[a.GetIdx()] for a in mol.GetAtoms() if a.GetAtomicNum() > 1]
                require(heavy_ids == provenance['heavy_atom_source_ids'] and
                        [heavy_ids[i] for i in described['ring_atoms']] == provenance['ring_source_atom_ids'],
                        'Sampled source atom IDs differ from descriptor provenance')
                require(json.loads(json.dumps(described['typed_feature_atoms'])) == provenance['typed_feature_atoms'],
                        'Sampled typed feature atom mapping differs')
                require(described['hard_group'] == row['hard_group'] and np.allclose(
                    described['descriptor'], row['descriptor'], atol=1e-9, rtol=1e-9),
                    'Sampled source descriptor does not reproduce')
                heavy = Chem.RemoveHs(Chem.Mol(mol)); xyz = np.asarray(heavy.GetConformer().GetPositions())
                units = provenance['units']; degrees = []
                require(units == described['units'], 'Sampled residue mapping differs')
                for i, (_, ca, carbon, _) in enumerate(units):
                    nxt = units[(i+1) % len(units)]
                    degrees.append(independent_omega(xyz[[ca, carbon, nxt[0], nxt[1]]]))
                states = ['cis' if abs(a) <= 30 else 'trans' if abs(a) >= 150 else 'boundary' for a in degrees]
                require(states == provenance['omega_states'], 'Independent omega states disagree')
                expected = np.degrees(np.asarray(described['torsions_radians'])[:, 2])
                require(np.max(np.abs((np.asarray(degrees)-expected+180) % 360-180)) < 1e-7,
                        'Independent omega angles disagree')
                stream.write(json.dumps(dict(conformer_id=row['conformer_id'], molecule_id=row['molecule_id'],
                    source_name=record.molecule_name, source_path=source, source_record_index=record.record_index,
                    content_sha256=record.content_sha256, ring_source_atom_ids=provenance['ring_source_atom_ids'],
                    omega_degrees=degrees, omega_states=states, flags=provenance['flags']))+'\n')
                remaining.remove(key); inspected += 1
                if not any(path == source for path, _ in remaining):
                    break
    require(not remaining, 'Sampled source records were not found')
    return inspected


def validate(build, output, source_samples=200):
    build, output = Path(build).resolve(), Path(output).resolve()
    require(output != build and build not in output.parents, 'Validation output must be outside the build tree')
    output.mkdir(parents=True, exist_ok=False)
    result = dict(status='failed', structural_gate='failed', release_gate='hold',
                  scientific_quality='not_established_by_integrity_checks', errors=[])
    try:
        receipt = json.loads((build / 'report.json').read_text(encoding='utf-8'))
        require(receipt.get('status') == 'complete', 'Block build is incomplete')
        for name in ('descriptors.sqlite', 'review.jsonl'):
            require(sha(build / name) == receipt['output_hashes'][name], 'Build artifact hash mismatch: '+name)
        with readonly(build / 'descriptors.sqlite') as db:
            described = db.execute('SELECT count(*) FROM descriptor').fetchone()[0]
            admitted = db.execute('SELECT count(*) FROM admission').fetchone()[0]
            require(db.execute('SELECT count(*) FROM admission WHERE visited!=1').fetchone()[0] == 0,
                    'Not all audit records visited exactly once')
            db.execute('PRAGMA temp_store=FILE')
            db.execute('CREATE TEMP TABLE delivered(path TEXT,idx INTEGER,PRIMARY KEY(path,idx))')
            flags = Counter(); verification_sources = Counter()
            for (payload,) in db.execute('SELECT payload FROM descriptor'):
                descriptor = json.loads(payload); provenance = descriptor['provenance']
                flags.update(provenance.get('flags', []))
                verification_sources[provenance.get('verification_source', 'unspecified')] += 1
                key = (provenance['source_path'], provenance['source_record_index'])
                admitted_row = db.execute('SELECT payload FROM admission WHERE path=? AND record_index=?', key).fetchone()
                require(admitted_row is not None, 'Descriptor source missing from admission')
                admission = json.loads(admitted_row[0])
                require(admission['state'] == 'geometry_candidate' and
                        (admission['content_sha256'], admission['conformer_name'], admission['name']) ==
                        (provenance['content_sha256'], provenance['source_record_name'], provenance['source_name']),
                        'Descriptor does not match admitted source identity')
                db.execute('INSERT INTO delivered VALUES(?,?)', key)
            review_count = 0
            with (build / 'review.jsonl').open(encoding='utf-8') as stream:
                for line in stream:
                    row = json.loads(line); key = (row['path'], row['record_index'])
                    admitted_row = db.execute('SELECT payload FROM admission WHERE path=? AND record_index=?', key).fetchone()
                    require(admitted_row is not None and json.loads(admitted_row[0])['conformer_name'] == row['conformer_name'],
                            'Review record does not match admitted source identity')
                    db.execute('INSERT INTO delivered VALUES(?,?)', key)
                    review_count += 1
            require(db.execute('SELECT count(*) FROM delivered').fetchone()[0] == admitted,
                    'Audit records missing or duplicated between descriptors and review')
        counts = receipt['counts']
        require(described+review_count == admitted == counts['source_conformers'], 'Source accounting does not balance')
        require(described == counts.get('described_conformers', 0) and review_count == counts.get('review_conformers', 0),
                'Build receipt count mismatch')
        capacities = [row['capacity'] for row in receipt['capacities']]
        require(bool(capacities) and len(set(capacities)) == len(capacities), 'Invalid capacity comparison manifest')
        models = []
        for embedded in receipt['capacities']:
            model = build / f"capacity-{embedded['capacity']}"
            require(json.loads((model / 'report.json').read_text(encoding='utf-8')) == embedded,
                    'Model receipt differs from build receipt')
            models.append(validate_model(model, build / 'descriptors.sqlite', receipt['schema']))
        inspected = validate_source_sample(build, receipt['schema'], output, source_samples)
        with (output / 'source_inspection.jsonl').open(encoding='utf-8') as stream:
            sampled_ids = [(json.loads(line)['conformer_id'],) for line in stream]
        for model in models:
            with readonly(build / f"capacity-{model['capacity']}" / 'blocks.sqlite') as db:
                db.execute('CREATE TEMP TABLE inspected(cid TEXT PRIMARY KEY)')
                db.executemany('INSERT INTO inspected VALUES(?)', sampled_ids)
                model['raw_source_sampled_blocks'] = db.execute(
                    'SELECT count(DISTINCT node) FROM point JOIN inspected USING(cid)').fetchone()[0]
        result.update(status='complete', structural_gate='passed',
            release_gate='review_required' if review_count else 'passed_for_offline_use',
            input_report_sha256=sha(build / 'report.json'), source_conformers=admitted,
            described_conformers=described, review_conformers=review_count,
            descriptor_flags=dict(flags), verification_sources=dict(verification_sources),
            source_samples_checked=inspected, source_sample_scope='Deterministic bounded inspection, not full raw-source recomputation',
            models=models, source_inspection_sha256=sha(output / 'source_inspection.jsonl'),
            limitations=['No claim of optimal clusters, energy basins, activity or screening recall',
                        'Do not reject unsampled blocks before held-out reference-score evaluation',
                        'Source-only IDs still require explicit registry joining for production'])
    except (ValueError, KeyError, TypeError, OSError, sqlite3.Error) as exc:
        result['errors'].append(str(exc))
    result['validator_sha256'] = sha(Path(__file__))
    (output / 'report.json').write_text(json.dumps(result, indent=2), encoding='utf-8')
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--build', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--source-samples', type=int, default=200)
    parser.add_argument('--allow-review-for-offline', action='store_true',
                        help='Exit successfully after complete structural validation while preserving review_required')
    args = parser.parse_args()
    result = validate(args.build, args.output, args.source_samples)
    print(json.dumps(result, indent=2))
    reviewed_offline = (args.allow_review_for_offline and result.get('status') == 'complete'
                        and result.get('structural_gate') == 'passed' and not result.get('errors')
                        and result.get('release_gate') == 'review_required')
    raise SystemExit(0 if result['release_gate'] == 'passed_for_offline_use' or reviewed_offline else 2)


if __name__ == '__main__':
    main()
