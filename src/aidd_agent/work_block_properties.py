"""Resumable source profiles and population-floored side-chain work-block refinement."""
import argparse
from concurrent.futures import ProcessPoolExecutor
import json
from pathlib import Path
import sqlite3
import time

import numpy as np
from rdkit import rdBase, RDConfig

from .boundary_pair_review import prepare_record, readonly
from .calibrated_boundary_routing import feature_matrix
from .final_work_blocks import read, sealed, sha, identifier, write_csv
from .mol2 import iter_mol2_blocks, parse_mol2_block

VERSION = 'work-block-property-summary-v1'
WIDTH = 69


def summarize(matrix):
    """Order-invariant broad chemistry and actual local sterics, not sequence identity."""
    x = np.asarray(matrix, dtype=float)
    if x.ndim != 2 or x.shape[1] != 29 or not len(x) or not np.isfinite(x).all():
        raise ValueError('Invalid source profile')
    x = x[:, :23]
    return np.concatenate((x.mean(0), x.std(0), x.max(0))).astype('<f4')


def compute_profile(task):
    path, index, text, payload, variant = task
    row = json.loads(payload)
    record = parse_mol2_block(Path(path), index, text)
    profile = prepare_record(record, row, variant)
    return row['conformer_id'], summarize(feature_matrix(profile, row, variant)).tobytes()


def profile_sources():
    names = ['work_block_properties.py', 'boundary_pair_review.py', 'calibrated_boundary_routing.py',
             'macrocycle_property_profiles.py', 'macrocycle_descriptors.py', 'macrocycle_peptide.py', 'mol2.py']
    return {n: sha(Path(__file__).with_name(n)) for n in names}


def prepare(blocks, build, output, workers=4, resume=False):
    """Commit bounded batches; incomplete runs never expose a completed receipt."""
    started = time.monotonic()
    if workers < 1 or workers > 64:
        raise ValueError('Workers must be between 1 and 64')
    blocks, build, output = [Path(p).resolve() for p in (blocks, build, output)]
    br = sealed(blocks, ['work_blocks.sqlite'])
    routing = Path(br['routing']); model = Path(br['model'])
    if any(output == p or p in output.parents or output in p.parents for p in (blocks, build, routing, model)):
        raise ValueError('Use a separate profile output')
    rr = sealed(routing, ['features.sqlite'])
    sealed(model, ['blocks.sqlite'])
    for path, digest in br['sources'].items():
        if sha(path) != digest:
            raise ValueError('Parent work-block provenance changed')
    for name in ['calibrated_boundary_routing.py', 'boundary_pair_review.py', 'macrocycle_property_profiles.py', 'macrocycle_descriptors.py']:
        if sha(Path(__file__).with_name(name)) != rr['code_hashes'][name]:
            raise ValueError('Do not mix E091 features with a changed profile implementation: ' + name)
    if rr['sources'].get(str(build / 'report.json')) != sha(build / 'report.json'):
        raise ValueError('Descriptor build does not match E091')
    source_report = read(build / 'report.json')
    if source_report['schema']['rdkit'] != rdBase.rdkitVersion:
        raise ValueError('Use the RDKit version recorded in the source build')
    # This full sequential hash prevents resuming against modified descriptors.
    if sha(build / 'descriptors.sqlite') != source_report['output_hashes']['descriptors.sqlite']:
        raise ValueError('Descriptor database changed')
    signature = dict(version=VERSION, blocks_report_sha256=sha(blocks / 'report.json'),
                     build_report_sha256=sha(build / 'report.json'), code_hashes=profile_sources(),
                     rdkit=rdBase.rdkitVersion, numpy=np.__version__,
                     feature_definition_sha256=sha(Path(RDConfig.RDDataDir) / 'BaseFeatures.fdef'))
    if signature['feature_definition_sha256'] != source_report['schema']['feature_definition_sha256']:
        raise ValueError('Chemical feature definitions differ from the source build')
    if output.exists():
        if not resume or read(output / 'signature.json') != signature:
            raise ValueError('Resume requires unchanged inputs/code and --resume')
        if (output / 'report.json').exists():
            sealed(output, ['profiles.sqlite'])
            return read(output / 'report.json')
    else:
        output.mkdir()
        (output / 'signature.json').write_text(json.dumps(signature, indent=2), encoding='utf-8')
    with sqlite3.connect(output / 'profiles.sqlite') as db, readonly(model / 'blocks.sqlite') as original, readonly(build / 'descriptors.sqlite') as descriptors:
        db.executescript('''CREATE TABLE IF NOT EXISTS files(id INTEGER PRIMARY KEY,path TEXT UNIQUE);
            CREATE TABLE IF NOT EXISTS item(cid TEXT PRIMARY KEY,parent TEXT NOT NULL,file_id INTEGER,idx INTEGER,vector BLOB,origin TEXT);
            CREATE TABLE IF NOT EXISTS meta(key TEXT PRIMARY KEY,value TEXT);''')
        planned = db.execute("SELECT value FROM meta WHERE key='planned'").fetchone()
        if not planned:
            original.execute('ATTACH DATABASE ? AS work', ((blocks / 'work_blocks.sqlite').as_uri() + '?mode=ro',))
            original.execute('ATTACH DATABASE ? AS cache', ((routing / 'features.sqlite').as_uri() + '?mode=ro',))
            sql = '''SELECT p.cid,COALESCE(o.block_id,c.block_id),
                json_extract(p.provenance,'$.source_path'),json_extract(p.provenance,'$.source_record_index'),f.n,f.vector
                FROM point p LEFT JOIN work.class_map c ON c.group_id=p.group_id
                LEFT JOIN work.override o ON o.cid=p.cid LEFT JOIN cache.feature f ON f.cid=p.cid
                WHERE COALESCE(o.block_id,c.block_id) IS NOT NULL AND COALESCE(o.block_id,c.block_id)!='special-exhaustive' '''
            files = {path: fid for fid, path in db.execute('SELECT id,path FROM files')}
            for number, (cid, parent, path, idx, n, vector) in enumerate(original.execute(sql), 1):
                if path is None or idx is None:
                    raise ValueError('Missing source locator')
                if path not in files:
                    files[path] = db.execute('INSERT INTO files(path) VALUES(?)', (path,)).lastrowid
                cached = summarize(np.frombuffer(vector, dtype='<f8').reshape(n, 29)).tobytes() if vector is not None else None
                db.execute('INSERT OR IGNORE INTO item VALUES(?,?,?,?,?,?)',
                    (cid, parent, files[path], idx, cached, 'E091' if cached is not None else None))
                if number % 100000 == 0:
                    db.commit(); print('Planned property profiles: ' + str(number), flush=True)
            if db.execute('SELECT count(*) FROM item').fetchone()[0] != br['regular_conformers']:
                raise ValueError('Regular property population mismatch')
            db.execute('CREATE INDEX IF NOT EXISTS profile_file ON item(file_id,idx)')
            db.execute('CREATE INDEX IF NOT EXISTS profile_parent ON item(parent)')
            db.execute("INSERT INTO meta VALUES('planned','true')"); db.commit()
        ready = db.execute('SELECT count(*) FROM item WHERE vector IS NOT NULL').fetchone()[0]
        print(f'Profiles ready: {ready}/{br["regular_conformers"]}; source extraction workers={workers}', flush=True)
        executor = ProcessPoolExecutor(max_workers=workers) if workers > 1 else None
        try:
            for file_id, path in db.execute('SELECT id,path FROM files ORDER BY id').fetchall():
                pending = iter(db.execute('SELECT idx,cid FROM item WHERE file_id=? AND vector IS NULL ORDER BY idx', (file_id,)))
                wanted = next(pending, None)
                if wanted is None:
                    continue
                print('Source profiles: ' + path, flush=True)
                batch = []
                def flush():
                    nonlocal ready
                    results = executor.map(compute_profile, batch, chunksize=4) if executor else map(compute_profile, batch)
                    for cid, vector in results:
                        db.execute("UPDATE item SET vector=?,origin='source_verified' WHERE cid=? AND vector IS NULL", (vector, cid))
                        ready += 1
                    db.commit()
                    (output / 'progress.json').write_text(json.dumps(dict(status='running', profiles_ready=ready,
                        total_regular=br['regular_conformers'], last_source=path)), encoding='utf-8')
                    print(f'Profiles ready: {ready}/{br["regular_conformers"]}', flush=True)
                    batch.clear()
                for index, text in iter_mol2_blocks(Path(path)):
                    if wanted is None:
                        break
                    if index != wanted[0]:
                        continue
                    row = descriptors.execute('SELECT payload FROM descriptor WHERE cid=?', (wanted[1],)).fetchone()
                    if row is None:
                        raise ValueError('Missing source descriptor')
                    batch.append((path, index, text, row[0], source_report['schema']['variant']))
                    wanted = next(pending, None)
                    if len(batch) >= 128:
                        flush()
                if wanted is not None:
                    raise ValueError('Source file ended before all expected profiles')
                if batch:
                    flush()
        finally:
            if executor:
                executor.shutdown(wait=True, cancel_futures=True)
        if db.execute('SELECT count(*) FROM item WHERE vector IS NULL').fetchone()[0]:
            raise ValueError('Incomplete profile coverage')
        origins = {origin: db.execute('SELECT count(*) FROM item WHERE origin=?', (origin,)).fetchone()[0]
                   for origin in ('E091', 'source_verified')}
    result = dict(status='complete', version=VERSION, blocks=str(blocks), total_regular=br['regular_conformers'],
        width=WIDTH, coverage=1.0, origins=origins, wall_seconds_this_invocation=time.monotonic()-started,
        signature=signature, output_hashes={'profiles.sqlite': sha(output / 'profiles.sqlite')},
        limitations=['Order-invariant per-unit mean/std/max are broad work-routing descriptors, not sequence identity',
            'Steric channels use actual local coordinate extents and branching; they are not pocket clash energies',
            'Special conformers are not subdivided or excluded', 'No screening recall or docking validation'])
    (output / 'report.json').write_text(json.dumps(result, indent=2), encoding='utf-8')
    return result


def partition(matrix, minimum=5000, maximum_children=8, minimum_contrast=0.35):
    """Adaptive PCA median work partitions with a hard population floor."""
    x = np.asarray(matrix, dtype=np.float32)
    if x.ndim != 2 or not np.isfinite(x).all() or minimum < 1 or maximum_children < 1 or minimum_contrast <= 0:
        raise ValueError('Invalid property partition inputs')
    leaves = [np.arange(len(x))]; evidence = []
    def proposal(indices):
        if len(indices) < 2 * minimum:
            return None
        values = x[indices]; centered = values - values.mean(0)
        covariance = centered.T @ centered / len(values)
        _, axes = np.linalg.eigh(covariance)
        axis = axes[:, -1]
        pivot = int(np.argmax(abs(axis)))
        if axis[pivot] < 0:
            axis = -axis
        projection = centered @ axis
        ordered = np.argsort(projection, kind='stable')
        sorted_values = projection[ordered]
        options = np.flatnonzero(np.diff(sorted_values) > 1e-7) + 1
        options = options[(options >= minimum) & (options <= len(indices) - minimum)]
        if not len(options):
            return None
        cut = int(options[np.argmin(abs(options - len(indices) / 2))])
        left, right = indices[ordered[:cut]], indices[ordered[cut:]]
        contrast = float(np.sqrt(np.mean((x[left].mean(0) - x[right].mean(0)) ** 2)))
        if contrast < minimum_contrast:
            return None
        return left, right, contrast
    proposals = [proposal(leaves[0])] if maximum_children > 1 else [None]
    while len(leaves) < maximum_children:
        choices = list(enumerate(proposals))
        choices = [(i, p) for i, p in choices if p is not None]
        if not choices:
            break
        index, (left, right, contrast) = max(choices, key=lambda pair: (pair[1][2], -pair[0]))
        leaves[index:index+1] = [left, right]
        proposals[index:index+1] = [proposal(left), proposal(right)] if len(leaves) < maximum_children else [None, None]
        evidence.append(dict(left_count=len(left), right_count=len(right), centroid_contrast=contrast))
    if len(leaves) > 1 and min(map(len, leaves)) < minimum:
        raise AssertionError('Small property leaf')
    return leaves, evidence


def refine(profiles, output, maximum_children=8, minimum_contrast=0.35):
    started = time.monotonic(); profiles, output = Path(profiles).resolve(), Path(output).resolve()
    if not 1 <= maximum_children <= 64 or not math_is_positive(minimum_contrast):
        raise ValueError('Positive finite contrast and 1..64 maximum children required')
    pr = sealed(profiles, ['profiles.sqlite']); blocks = Path(pr['blocks'])
    br = sealed(blocks, ['work_blocks.sqlite'])
    if pr['signature']['blocks_report_sha256'] != sha(blocks / 'report.json') or pr['coverage'] != 1:
        raise ValueError('Complete matching property coverage required')
    if any(output == p or p in output.parents or output in p.parents for p in (profiles, blocks)):
        raise ValueError('Fresh independent property output required')
    output.mkdir(exist_ok=False); rows = []
    with readonly(profiles / 'profiles.sqlite') as source, readonly(blocks / 'work_blocks.sqlite') as work, sqlite3.connect(output / 'property_blocks.sqlite') as dest:
        dest.executescript('''CREATE TABLE membership(cid TEXT PRIMARY KEY,parent_block TEXT,block_id TEXT);
            CREATE TABLE block(block_id TEXT PRIMARY KEY,parent_block TEXT,n INTEGER,metadata TEXT);''')
        for parent, n in work.execute("SELECT block_id,n FROM block WHERE kind!='special_exhaustive' ORDER BY block_id"):
            print('Property refinement: ' + parent + ' / ' + str(n), flush=True)
            scratch = output / 'profile-matrix.tmp'
            matrix = np.memmap(scratch, dtype='<f4', mode='w+', shape=(n, WIDTH))
            cids = []
            for i, (cid, blob) in enumerate(source.execute('SELECT cid,vector FROM item WHERE parent=? ORDER BY cid', (parent,))):
                if i >= n or blob is None:
                    raise ValueError('Property population mismatch')
                matrix[i] = np.frombuffer(blob, dtype='<f4'); cids.append(cid)
            if len(cids) != n:
                raise ValueError('Missing property members')
            leaves, evidence = partition(matrix, br['minimum_regular_block_size'], maximum_children, minimum_contrast)
            for number, indices in enumerate(leaves):
                bid = identifier([VERSION, parent, number, maximum_children, minimum_contrast, sha(profiles / 'report.json')])
                metadata = dict(block_id=bid, parent_block=parent, conformers=len(indices),
                                split=len(leaves) > 1, split_evidence=evidence)
                dest.execute('INSERT INTO block VALUES(?,?,?,?)', (bid, parent, len(indices), json.dumps(metadata)))
                dest.executemany('INSERT INTO membership VALUES(?,?,?)', ((cids[int(i)], parent, bid) for i in indices))
                rows.append(metadata)
            dest.commit(); del matrix; scratch.unlink()
        dest.execute('CREATE INDEX member_block ON membership(block_id)'); dest.commit()
        if dest.execute('SELECT count(*) FROM membership').fetchone()[0] != br['regular_conformers']:
            raise ValueError('Property conservation failed')
    write_csv(output / 'blocks.csv', rows, ['block_id', 'parent_block', 'conformers', 'split', 'split_evidence'])
    result = dict(status='complete', readiness='offline_property_work_blocks_not_validated_search_dispatch',
        parent_blocks=str(blocks), profiles=str(profiles), regular_work_blocks=len(rows),
        regular_conformers=br['regular_conformers'], special_conformers=br['special_conformers'],
        special_pool='special-exhaustive in parent work-block artifact', minimum_regular_block_size=br['minimum_regular_block_size'],
        smallest_regular_block=min((r['conformers'] for r in rows), default=0),
        blocks_below_minimum=sum(r['conformers'] < br['minimum_regular_block_size'] for r in rows),
        maximum_children_per_parent=maximum_children, minimum_centroid_contrast=minimum_contrast,
        logical_capacity_limit=None, wall_seconds=time.monotonic()-started,
        sources={str(p / 'report.json'): sha(p / 'report.json') for p in (blocks, profiles)},
        output_hashes={n: sha(output / n) for n in ['property_blocks.sqlite', 'blocks.csv']},
        implementation_sha256=sha(Path(__file__)),
        limitations=['Adaptive broad property work partitions; no biological-state or optimal-cluster claim',
            'Original identities remain in the parent overlay; property membership never changes E091 acceptance',
            'Centroid contrast and maximum children are engineering defaults requiring held-out evaluation',
            'No old radius reuse or automatic block rejection; small source libraries remain a single block'])
    (output / 'report.json').write_text(json.dumps(result, indent=2), encoding='utf-8')
    return result


def math_is_positive(value):
    return bool(np.isfinite(value) and value > 0)


def main():
    parser = argparse.ArgumentParser(description=__doc__); sub = parser.add_subparsers(dest='command', required=True)
    p = sub.add_parser('prepare')
    for name in ['blocks', 'build', 'output']:
        p.add_argument('--' + name, type=Path, required=True)
    p.add_argument('--workers', type=int, default=4); p.add_argument('--resume', action='store_true')
    p = sub.add_parser('refine')
    p.add_argument('--profiles', type=Path, required=True); p.add_argument('--output', type=Path, required=True)
    p.add_argument('--maximum-children', type=int, default=8); p.add_argument('--minimum-contrast', type=float, default=0.35)
    a = parser.parse_args()
    result = prepare(a.blocks, a.build, a.output, a.workers, a.resume) if a.command == 'prepare' else refine(a.profiles, a.output, a.maximum_children, a.minimum_contrast)
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
