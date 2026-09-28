"""Disk-backed deterministic capacity trees and molecule-aware sample manifests."""
from collections import Counter
import hashlib
import json
from pathlib import Path
import sqlite3

import numpy as np

from .macrocycle_blocks import digest
from .macrocycle_identity_audit import sha

VERSION = 'disk-conformer-median-tree-v1'


def moments(db, node):
    count = 0; mean = None; m2 = None
    cursor = db.execute('SELECT vector FROM point WHERE node=? ORDER BY cid', (node,))
    while rows := cursor.fetchmany(1024):
        x = np.asarray([json.loads(row[0]) for row in rows], dtype=float)
        if x.ndim != 2 or not x.shape[1] or not np.isfinite(x).all():
            raise ValueError('Invalid descriptor matrix')
        avg = x.mean(axis=0); batch_m2 = ((x-avg)**2).sum(axis=0)
        if mean is None:
            mean, m2, count = avg, batch_m2, len(x)
        else:
            delta = avg-mean; new = count+len(x)
            m2 += batch_m2 + delta**2 * count*len(x)/new
            mean += delta * len(x)/new; count = new
    return count, mean, m2


def fit(rows, output, capacity, schema):
    if type(capacity) is not int or capacity < 1:
        raise ValueError('Positive integer capacity required')
    schema = json.loads(json.dumps(schema))
    output = Path(output)
    output.mkdir(parents=True, exist_ok=False)
    counts = Counter()
    with sqlite3.connect(output / 'blocks.sqlite') as db:
        db.execute('PRAGMA temp_store=FILE')
        db.executescript('''CREATE TABLE point(cid TEXT PRIMARY KEY, mid TEXT NOT NULL,
            group_id TEXT NOT NULL, vector TEXT NOT NULL, provenance TEXT, node TEXT);
            CREATE INDEX node_index ON point(node);
            CREATE INDEX group_index ON point(group_id);
            CREATE TABLE tree(node TEXT PRIMARY KEY, group_id TEXT, path TEXT, n INTEGER,
            axis INTEGER, cut REAL, cut_id TEXT, left_node TEXT, right_node TEXT,
            centroid TEXT, radius REAL);
            CREATE UNIQUE INDEX tree_group_path ON tree(group_id,path);
            CREATE TABLE metadata(key TEXT PRIMARY KEY,value TEXT);''')
        for row in rows:
            vector = np.asarray(row['descriptor'], dtype=float)
            if vector.ndim != 1 or not len(vector) or not np.isfinite(vector).all():
                raise ValueError('Invalid descriptor')
            if not row['conformer_id'] or not row['molecule_id'] or not row['hard_group']:
                raise ValueError('Missing conformer/molecule/group identity')
            db.execute('INSERT INTO point VALUES(?,?,?,?,?,NULL)',
                       (row['conformer_id'], row['molecule_id'], row['hard_group'],
                        json.dumps(vector.tolist()), json.dumps(row.get('provenance', {}))))
        db.commit()
        fingerprint = hashlib.sha256()
        for row in db.execute('SELECT cid,mid,group_id,vector FROM point ORDER BY cid'):
            fingerprint.update((json.dumps(row)+'\n').encode())
        model_id = digest([VERSION, schema, capacity, fingerprint.hexdigest()])
        for key, value in dict(version=VERSION, schema=schema, capacity=capacity,
                               model_id=model_id).items():
            db.execute('INSERT INTO metadata VALUES(?,?)', (key, json.dumps(value)))
        visited_nodes = 0
        groups = db.execute('SELECT DISTINCT group_id FROM point ORDER BY group_id')
        for (group,) in groups:
            root = digest([model_id, group, ''])
            db.execute('UPDATE point SET node=? WHERE group_id=?', (root, group))
            pending = [(root, '')]
            while pending:
                node, path = pending.pop()
                n, mean, m2 = moments(db, node)
                if n > capacity:
                    axis = int(np.argmax(m2/n))
                    expression = f"json_extract(vector,'$[{axis}]')"
                    cut, cut_id = db.execute(f'SELECT {expression},cid FROM point WHERE node=? ORDER BY {expression},cid LIMIT 1 OFFSET ?', (node, n//2)).fetchone()
                    left, right = (digest([model_id, group, path+bit]) for bit in '01')
                    db.execute(f'''UPDATE point SET node=CASE WHEN {expression} < ? OR
                        ({expression} = ? AND cid < ?) THEN ? ELSE ? END WHERE node=?''',
                        (cut, cut, cut_id, left, right, node))
                    db.execute('INSERT INTO tree VALUES(?,?,?,?,?,?,?,?,?,?,?)',
                               (node, group, path, n, axis, cut, cut_id, left, right, None, None))
                    pending.extend(((right, path+'1'), (left, path+'0')))
                else:
                    radius = 0.0
                    cursor = db.execute('SELECT vector FROM point WHERE node=?', (node,))
                    while batch := cursor.fetchmany(1024):
                        x = np.asarray([json.loads(row[0]) for row in batch])
                        radius = max(radius, float(np.linalg.norm(x-mean, axis=1).max()))
                    db.execute('INSERT INTO tree VALUES(?,?,?,?,?,?,?,?,?,?,?)',
                               (node, group, path, n, None, None, None, None, None,
                                json.dumps(mean.tolist()), radius))
                    counts['blocks'] += 1; counts['conformers'] += n
                visited_nodes += 1
                if visited_nodes % 128 == 0:
                    db.commit()
        db.execute('CREATE INDEX block_molecule ON point(node,mid)')
        counts['molecules'] = db.execute('SELECT count(DISTINCT mid) FROM point').fetchone()[0]
        sizes = dict(db.execute('SELECT n,count(*) FROM tree WHERE axis IS NULL GROUP BY n'))
        with (output / 'memberships.jsonl').open('w', encoding='utf-8') as stream:
            for cid, mid, node, provenance in db.execute('SELECT cid,mid,node,provenance FROM point ORDER BY cid'):
                stream.write(json.dumps(dict(conformer_id=cid, molecule_id=mid,
                    block_id=node, provenance=json.loads(provenance)))+'\n')
        db.commit()
    report = dict(status='complete', version=VERSION, model_id=model_id, schema=schema,
        capacity=capacity, counts=dict(counts), block_size_histogram=sizes,
        descriptor_content_sha256=fingerprint.hexdigest(),
        output_hashes={p: sha(output / p) for p in ('blocks.sqlite', 'memberships.jsonl')},
        implementation_sha256=sha(Path(__file__)),
        scope='Offline fitted conformer blocks; no screening rejection or recall claim')
    (output / 'report.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
    return report


def open_model(directory):
    directory = Path(directory).resolve()
    report = json.loads((directory / 'report.json').read_text(encoding='utf-8'))
    if report.get('status') != 'complete' or sha(directory / 'blocks.sqlite') != report['output_hashes']['blocks.sqlite']:
        raise ValueError('Incomplete or modified block model')
    db = sqlite3.connect((directory / 'blocks.sqlite').as_uri()+'?mode=ro', uri=True)
    return db, report


def route(db, row):
    node = db.execute("SELECT node FROM tree WHERE group_id=? AND path=''", (row['hard_group'],)).fetchone()
    if node is None:
        return None, 'unknown_hard_group'
    while True:
        state = db.execute('SELECT axis,cut,cut_id,left_node,right_node,centroid,radius FROM tree WHERE node=?', node).fetchone()
        axis, cut, cut_id, left, right, centroid, radius = state
        vector = np.asarray(row['descriptor'], dtype=float)
        if vector.ndim != 1 or not np.isfinite(vector).all():
            return None, 'invalid_descriptor'
        if axis is None:
            center = np.asarray(json.loads(centroid))
            if vector.shape != center.shape:
                return None, 'descriptor_dimension_mismatch'
            if np.linalg.norm(vector-center) > radius+1e-8:
                return None, 'outside_frozen_radius'
            return node[0], None
        if axis >= len(vector):
            return None, 'descriptor_dimension_mismatch'
        node = (left if (vector[axis], row['conformer_id']) < (cut, cut_id) else right,)


def propose(rows, model, output, schema):
    """Write deterministic incremental proposals; never change the frozen model."""
    db, report = open_model(model)
    try:
        if json.loads(json.dumps(schema)) != report['schema']:
            raise ValueError('Descriptor schema mismatch')
        output = Path(output); output.mkdir(parents=True, exist_ok=False)
        # Staging new records on disk makes ordering and duplicate rejection explicit.
        with sqlite3.connect(output / 'proposals.sqlite') as staging:
            staging.execute('CREATE TABLE incoming(cid TEXT PRIMARY KEY,payload TEXT)')
            staging.executemany('INSERT INTO incoming VALUES(?,?)',
                                ((r['conformer_id'], json.dumps(r)) for r in rows))
            counts = Counter(); occupancy = Counter()
            with (output / 'proposals.jsonl').open('w', encoding='utf-8') as stream:
                for cid, payload in staging.execute('SELECT cid,payload FROM incoming ORDER BY cid'):
                    row = json.loads(payload)
                    if db.execute('SELECT 1 FROM point WHERE cid=?', (cid,)).fetchone():
                        node, reason = None, 'existing_conformer_id'
                    else:
                        node, reason = route(db, row)
                    if node:
                        n = db.execute('SELECT n FROM tree WHERE node=?', (node,)).fetchone()[0]
                        if n+occupancy[node] >= report['capacity']:
                            reason = 'capacity_overflow'; node = None
                        else:
                            occupancy[node] += 1
                    state = reason or 'proposed'
                    counts[state] += 1
                    stream.write(json.dumps(dict(conformer_id=cid, molecule_id=row['molecule_id'],
                                                block_id=node, state=state))+'\n')
        result = dict(status='proposal_only', model_id=report['model_id'], counts=dict(counts),
                      limitation='Do not combine separate proposal batches without rechecking occupancy')
        (output / 'report.json').write_text(json.dumps(result, indent=2), encoding='utf-8')
        return result
    finally:
        db.close()


def main():
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--model', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--descriptor-db', type=Path, help='Propose rows from an E080 descriptors.sqlite')
    parser.add_argument('--schema-report', type=Path, help='E080 report defining incoming descriptor schema')
    parser.add_argument('--sample-sizes', type=int, nargs='+', default=[100, 500])
    parser.add_argument('--seed', type=int, default=20260928)
    args = parser.parse_args()
    if args.descriptor_db:
        if not args.schema_report:
            parser.error('Incremental proposals require --schema-report')
        receipt = json.loads(args.schema_report.read_text(encoding='utf-8'))
        if receipt.get('status') != 'complete' or sha(args.descriptor_db) != receipt.get('output_hashes', {}).get('descriptors.sqlite'):
            raise ValueError('Incoming descriptor receipt/hash mismatch')
        with sqlite3.connect(args.descriptor_db.resolve().as_uri()+'?mode=ro', uri=True) as db:
            result = propose((json.loads(row[0]) for row in db.execute('SELECT payload FROM descriptor')),
                             args.model, args.output, receipt['schema'])
    else:
        if args.schema_report:
            parser.error('--schema-report requires --descriptor-db')
        result = sample(args.model, args.output, tuple(args.sample_sizes), args.seed)
    print(json.dumps(result, indent=2))


def sample(model, output, sizes=(100, 500), seed=20260928):
    if any(type(n) is not int or n < 1 for n in sizes) or len(set(sizes)) != len(sizes):
        raise ValueError('Unique positive sample sizes required')
    db, report = open_model(model)
    try:
        output = Path(output); output.mkdir(parents=True, exist_ok=False)
        db.execute('PRAGMA temp_store=FILE')
        db.create_function('sample_key', 2, lambda block, mid: digest([seed, block, mid]))
        db.execute('''CREATE TEMP TABLE selection AS SELECT node,mid,
            row_number() OVER(PARTITION BY node ORDER BY sample_key(node,mid),mid) AS rank
            FROM (SELECT DISTINCT node,mid FROM point)''')
        db.execute('CREATE INDEX selection_lookup ON selection(node,mid)')
        results = []
        for size in sizes:
            path = output / f'sample-{size}.jsonl'
            count = 0
            with path.open('w', encoding='utf-8') as stream:
                for block, mid, cid in db.execute('''SELECT p.node,p.mid,p.cid FROM point p
                    JOIN selection s ON s.node=p.node AND s.mid=p.mid WHERE s.rank<=?
                    ORDER BY p.node,p.mid,p.cid''', (size,)):
                    stream.write(json.dumps(dict(block_id=block, molecule_id=mid, conformer_id=cid))+'\n')
                    count += 1
            pairs, unique = db.execute('SELECT count(*),count(DISTINCT mid) FROM selection WHERE rank<=?', (size,)).fetchone()
            results.append(dict(requested_per_block=size, block_molecule_pairs=pairs,
                                unique_molecules=unique, conformer_evaluations=count, sha256=sha(path)))
        result = dict(model_id=report['model_id'], seed=seed, samples=results,
            scope='All in-block conformers of sampled molecules; reuse scores by conformer ID across budgets')
        (output / 'report.json').write_text(json.dumps(result, indent=2), encoding='utf-8')
        return result
    finally:
        db.close()


if __name__ == '__main__':
    main()
