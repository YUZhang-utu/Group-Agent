#!/usr/bin/env python3
"""E032 Linux batch driver; copy this file into an existing Group-Agent checkout.

No original file/catalog is modified. One output directory owns one frozen batch.
Completed sources are reused. A pilot (--max-new-files 1) stops before indexing.
"""
from __future__ import annotations

import argparse
from contextlib import contextmanager
import hashlib
import json
import os
from pathlib import Path
import shutil
import sqlite3
import sys
import time
import uuid

VERSION = 1


def log(message):
    print(time.strftime('%Y-%m-%d %H:%M:%S'), message, flush=True)


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(8 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def read(path):
    return json.loads(Path(path).read_text(encoding='utf-8'))


def write(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + '.tmp')
    with tmp.open('w', encoding='utf-8') as f:
        json.dump(value, f, indent=2)
        f.flush()
        os.fsync(f.fileno())
    os.replace(tmp, path)


def snapshot(path):
    s = path.stat()
    return {'path': str(path.resolve()), 'bytes': s.st_size, 'mtime_ns': s.st_mtime_ns}


def check_source(item):
    if snapshot(Path(item['path'])) != {k: item[k] for k in ('path', 'bytes', 'mtime_ns')}:
        raise ValueError('Source changed since inventory: ' + item['path'])


def verify(directory):
    directory = Path(directory).resolve()
    m = read(directory / 'manifest.json')
    for name, record in m.get('files', {}).items():
        p = (directory / name).resolve()
        if not p.is_relative_to(directory):
            raise ValueError('Manifest path escapes shard')
        if p.stat().st_size != record['bytes'] or digest(p) != record['sha256']:
            raise ValueError('Output integrity failure: ' + str(p))
    for record in m.get('keys', []):
        p = (directory / record['file']).resolve()
        if not p.is_relative_to(directory) or digest(p) != record['sha256']:
            raise ValueError('Posting integrity failure: ' + str(p))
    return m


def resolve_shards(catalog_path):
    catalog_path = Path(catalog_path).resolve()
    c = read(catalog_path)
    expected = 0
    for r in c['shards']:
        p = Path(r['path'])
        if not p.is_dir():
            p = catalog_path.parent / r['name']
        p = p.resolve()
        if digest(p / 'manifest.json') != r['manifest_sha256']:
            raise ValueError('Old manifest mismatch: ' + str(p))
        m = verify(p)
        if m['library_id'] != c['library_id'] or int(m['global_id_start']) != expected:
            raise ValueError('Old library/range mismatch')
        expected += int(m['conformers'])
        r['path'] = str(p)
    if expected != c['conformers']:
        raise ValueError('Old total does not match ranges')
    return c


def recover_partial(output, partial):
    partial = Path(partial)
    if partial.exists():
        if not partial.resolve().is_relative_to(output.resolve()) or partial.is_symlink():
            raise ValueError('Refusing to move an external partial')
        target = output / 'recovery' / (partial.name + '-' + uuid.uuid4().hex[:10])
        target.parent.mkdir(parents=True, exist_ok=True)
        partial.rename(target)
        log('Preserved incomplete output at ' + str(target))


def private_registry(path, library_id):
    from aidd_agent.registry import initialize
    initialize(path)
    with sqlite3.connect(path) as db:
        db.execute('PRAGMA journal_mode=WAL')
        db.execute('INSERT OR IGNORE INTO library VALUES (?,?,?)',
                   (library_id, 'expanded-library-e032', '2026-09-11'))
        db.execute('CREATE TABLE IF NOT EXISTS batch_source '
                   '(path TEXT PRIMARY KEY, sha TEXT, records INTEGER, inserted INTEGER, duplicates INTEGER)')
        db.execute('CREATE INDEX IF NOT EXISTS idx_batch_source_path ON conformer(source_path,source_record_index)')
        db.execute('CREATE INDEX IF NOT EXISTS idx_batch_molecule ON conformer(molecule_id)')


def register_source(dbpath, library_id, source, sha, old_shard=None):
    import numpy as np
    from aidd_agent.mol2 import iter_mol2_records
    from aidd_agent.registry import register_record
    with sqlite3.connect(dbpath) as db:
        db.row_factory = sqlite3.Row
        db.execute('PRAGMA foreign_keys=ON')
        done = db.execute('SELECT * FROM batch_source WHERE path=?', (str(source),)).fetchone()
        if done:
            if done['sha'] != sha:
                raise ValueError('Registered source hash changed')
            return dict(done)
        old_mol = old_conf = None
        if old_shard:
            old_mol = np.memmap(old_shard / 'molecule_ids.bin', dtype='S16', mode='r')
            old_conf = np.memmap(old_shard / 'conformer_ids.bin', dtype='S16', mode='r')
        n = inserted = duplicates = 0
        for record in iter_mol2_records(source):
            # Existing feature templates are reused within a molecule: enforce
            # exactly the same ordered topology before allowing that fast path.
            row = db.execute('SELECT c.topology_sha256 FROM conformer c JOIN molecule m '
                             'ON m.id=c.molecule_id WHERE m.library_id=? AND m.source_name=? LIMIT 1',
                             (library_id, record.molecule_name)).fetchone()
            if row and row[0] != record.topology_sha256:
                raise ValueError(f'Ordered topology conflict: {record.name} in {source}')
            if old_shard:
                if n >= len(old_conf):
                    raise ValueError('Old source contains more records than old artifact')
                mid, cid = old_mol[n].decode(), old_conf[n].decode()
                db.execute('INSERT OR IGNORE INTO molecule VALUES (?,?,?,?)',
                           (mid, library_id, record.molecule_name, '2026-09-11'))
                identity = db.execute('SELECT id FROM molecule WHERE library_id=? AND source_name=?',
                                      (library_id, record.molecule_name)).fetchone()
                if not identity or identity[0] != mid:
                    raise ValueError('Old molecule identity mismatch')
                db.execute('INSERT INTO conformer VALUES (?,?,?,?,?,?,?,?,?,?,?,?)',
                           (cid, mid, record.conformer_index, record.name, record.atom_count,
                            record.bond_count, record.content_sha256, record.topology_sha256,
                            str(source), record.record_index, json.dumps(record.warnings), '2026-09-11'))
                inserted += 1
            else:
                outcome = register_record(db, library_id, record)
                inserted += outcome == 'inserted'
                duplicates += outcome == 'duplicate'
            n += 1
            if n % 50000 == 0:
                log(f'Register {source.name}: {n:,} records')
        if not n or (old_shard and n != len(old_conf)):
            raise ValueError('Empty source or old source/artifact count mismatch')
        db.execute('INSERT INTO batch_source VALUES (?,?,?,?,?)', (str(source), sha, n, inserted, duplicates))
        return dict(path=str(source), sha=sha, records=n, inserted=inserted, duplicates=duplicates)


def artifact_row(directory):
    m = verify(directory)
    return dict(name=directory.name, path=str(directory.resolve()),
                manifest_sha256=digest(directory / 'manifest.json'),
                global_id_start=m['global_id_start'], conformers=m['conformers'],
                heavy_atoms=m['heavy_atoms'], features=m['features'])


def auxiliary(output, row, library_id, dbpath=None, old_chem=None):
    from aidd_agent.chemical_companion import build_chemical_companion_shard
    from aidd_agent.pharmacophore_index import build_shard_pharmacophore_index
    shard = Path(row['path'])
    target = output / 'chemical' / row['name']
    if old_chem:
        target = Path(old_chem['path'])
        if digest(target / 'manifest.json') != old_chem['manifest_sha256']:
            raise ValueError('Old chemistry manifest mismatch')
    if target.exists():
        m = verify(target)
    else:
        recover_partial(output, target.parent / ('.' + target.name + '.partial'))
        m = build_chemical_companion_shard(dbpath, library_id, shard, output / 'chemical',
                                         workers=WORKERS, max_moving_atoms=12)
        verify(target)
    if int(m['global_id_start']) != row['global_id_start'] or m['conformers'] != row['conformers']:
        raise ValueError('Chemistry/v1 range mismatch')
    if m['library_id'] != library_id or m['artifact_v1_manifest_sha256'] != row['manifest_sha256']:
        raise ValueError('Chemistry/v1 provenance mismatch')
    # IDs are checked explicitly even when reusing old chemistry.
    for name in ('conformer_ids.bin', 'molecule_ids.bin'):
        if digest(target / name) != digest(shard / name):
            raise ValueError('Chemistry/v1 identity mismatch')
    chemrow = dict(name=row['name'], path=str(target.resolve()),
                   manifest_sha256=digest(target / 'manifest.json'),
                   global_id_start=row['global_id_start'], conformers=row['conformers'])
    pharma = output / 'pharmacophore' / 'shards' / row['name']
    recover_partial(output, pharma.parent / ('.' + pharma.name + '.partial'))
    build_shard_pharmacophore_index(shard, pharma)
    verify(pharma)
    pharmarow = dict(chemrow, path=str(pharma.resolve()),
                    manifest_sha256=digest(pharma / 'manifest.json'))
    return chemrow, pharmarow


def make_indices(output, catalog, workers):
    """Uniform training across full catalog; per-shard checkpoints, one merge."""
    import numpy as np
    import faiss
    faiss.omp_set_num_threads(workers)
    directory = output / 'faiss'
    directory.mkdir(exist_ok=True)
    stamp = digest(output / 'artifacts' / 'catalog.json')
    total = catalog['conformers']
    trained = directory / 'trained'
    if trained.exists():
        meta = verify(trained)
        if meta['catalog_sha256'] != stamp:
            raise ValueError('FAISS training catalog changed; choose a new batch output')
    else:
        partial = directory / '.trained.partial'
        recover_partial(output, partial)
        partial.mkdir()
        rng = np.random.default_rng(20260911)
        size = min(total, 250000)
        if size < 10000:
            raise ValueError('Full library too small for this IVF-PQ batch preset')
        gids = np.sort(rng.integers(0, total, size=size))
        training = np.empty((size, 60), dtype=np.float32)
        for row in catalog['shards']:
            start = row['global_id_start']; end = start + row['conformers']
            a, b = np.searchsorted(gids, [start, end])
            if b > a:
                raw = np.memmap(Path(row['path']) / 'usrcat.f32.bin', dtype='<f4', mode='r').reshape(-1, 60)
                training[a:b] = raw[gids[a:b] - start]
        if not np.isfinite(training).all():
            raise ValueError('Non-finite USRCAT training data')
        mean, std = training.mean(0), training.std(0)
        std[std < 1e-6] = 1
        training -= mean; training /= std
        nlist = min(4096, 2 ** int(np.log2(max(64, size // 40))))
        index = faiss.IndexIVFPQ(faiss.IndexFlatL2(60), 60, nlist, 20, 8)
        index.cp.seed = 20260911
        index.pq.cp.seed = 20260911
        log(f'Train shared IVF-PQ: {size:,} samples, nlist={nlist}, all-shard sampling')
        t = time.monotonic(); index.train(training)
        faiss.write_index(index, str(partial / 'template.faiss'))
        np.savez(partial / 'transform.npz', mean=mean, std=std)
        write(partial / 'manifest.json', dict(catalog_sha256=stamp, train_seconds=time.monotonic()-t,
              nlist=nlist, train_size=size, seed=20260911, sample_with_replacement=True,
              files={p.name: dict(bytes=p.stat().st_size, sha256=digest(p)) for p in partial.iterdir()}))
        partial.rename(trained)
        del index, training
    meta = verify(trained)
    template_hash = digest(trained / 'template.faiss')
    with np.load(trained / 'transform.npz') as z:
        mean, std = z['mean'], z['std']
    records = []
    for row in catalog['shards']:
        target = directory / 'shards' / row['name']
        if target.exists():
            m = verify(target)
            if m['source_manifest_sha256'] != row['manifest_sha256'] or m['template_sha256'] != template_hash:
                raise ValueError('FAISS shard input changed')
        else:
            partial = target.parent / ('.' + target.name + '.partial')
            recover_partial(output, partial); partial.mkdir(parents=True)
            index = faiss.read_index(str(trained / 'template.faiss'))
            raw = np.memmap(Path(row['path']) / 'usrcat.f32.bin', dtype='<f4', mode='r').reshape(-1, 60)
            if len(raw) != row['conformers']:
                raise ValueError('USRCAT row count mismatch')
            for offset in range(0, len(raw), 50000):
                batch = np.asarray(raw[offset:offset+50000], dtype=np.float32).copy()
                if not np.isfinite(batch).all():
                    raise ValueError('Non-finite USRCAT batch')
                batch -= mean; batch /= std
                ids = np.arange(row['global_id_start']+offset,
                                row['global_id_start']+offset+len(batch), dtype=np.int64)
                index.add_with_ids(batch, ids)
            if index.ntotal != row['conformers']:
                raise ValueError('FAISS shard count mismatch')
            fp = partial / 'index.faiss'; faiss.write_index(index, str(fp))
            write(partial / 'manifest.json', dict(source_manifest_sha256=row['manifest_sha256'],
                  template_sha256=template_hash, ntotal=int(index.ntotal),
                  files={'index.faiss': dict(bytes=fp.stat().st_size, sha256=digest(fp))}))
            partial.rename(target); del index, raw
        records.append(dict(row, index_path=str(target / 'index.faiss')))
        log('FAISS shard ready: ' + row['name'])
    # Each source is encoded once. Restarting an interrupted merge only rereads
    # compact indices; no molecular parsing or vector encoding is repeated.
    final = directory / 'manifest.json'
    if final.exists():
        existing = read(final)
        if (existing['catalog_sha256'] != stamp or digest(directory / 'index.faiss') != existing['index_sha256']
                or digest(directory / 'transform.npz') != digest(trained / 'transform.npz')):
            raise ValueError('Final FAISS integrity failure')
        return existing
    merged = faiss.read_index(str(trained / 'template.faiss'))
    for row in records:
        other = faiss.read_index(row['index_path'])
        faiss.merge_into(merged, other, False)
        del other
    if int(merged.ntotal) != total:
        raise ValueError('Merged FAISS count mismatch')
    merged.nprobe = min(256, merged.nlist)
    temp = directory / 'index.faiss.partial'; faiss.write_index(merged, str(temp))
    os.replace(temp, directory / 'index.faiss')
    shutil.copyfile(trained / 'transform.npz', directory / 'transform.npz')
    result = dict(format='aidd-incremental-faiss', version=1,
                  catalog=str(output / 'artifacts' / 'catalog.json'), catalog_sha256=stamp,
                  final_ntotal=total, index_sha256=digest(directory / 'index.faiss'),
                  additions=[dict(index_path=str(directory / 'index.faiss'), ntotal=total)],
                  parameters=dict(dimension=60, nlist=int(merged.nlist), m=20, nbits=8, seed=20260911),
                  shard_indices=records, recall_status='pending independent exact-vs-approximate calibration')
    write(final, result)
    return result


WORKERS = 20


def execute(args):
    import numpy as np
    from aidd_agent.conformer_artifacts import build_shard_parallel, META_DTYPE
    from aidd_agent.mol2 import discover_mol2_files
    global WORKERS
    WORKERS = args.workers
    output = args.output.resolve(); source = args.source.resolve()
    if output == source or source.is_relative_to(output) or output.is_relative_to(source):
        raise ValueError('Keep outputs separate from source tree')
    if args.output.exists() and not (output / 'batch.json').exists():
        raise ValueError('Output exists but is not an owned batch; choose a new directory')
    files = [p.resolve() for p in discover_mol2_files(source)]
    if not files or len(set(files)) != len(files) or any(not p.is_relative_to(source) for p in files):
        raise ValueError('Empty source, aliased source, or source symlink escapes root')
    stems = [p.stem for p in files]
    if len(stems) != len(set(stems)):
        raise ValueError('Repeated file stems across folders: resolve explicitly before building shards')
    plan = dict(version=VERSION, source=str(source), files=[snapshot(p) for p in files],
                old_artifacts=str(args.old_artifacts.resolve()), old_artifacts_sha256=digest(args.old_artifacts),
                old_chemical=str(args.old_chemical.resolve()), old_chemical_sha256=digest(args.old_chemical))
    if not args.run:
        print(json.dumps(dict(files=len(files), source_gib=sum(p.stat().st_size for p in files)/2**30,
                              free_gib=shutil.disk_usage(source).free/2**30,
                              workers=WORKERS, output=str(output), note='Plan only; use --run'), indent=2))
        return
    # Fail before any expensive source scan if the workstation environment is
    # incomplete. The scientific dependencies are never installed implicitly.
    import rdkit
    import faiss
    log(f'Runtime Python={sys.version.split()[0]} RDKit={rdkit.__version__} FAISS={faiss.__version__}')
    if hasattr(os, 'sched_getaffinity') and WORKERS > len(os.sched_getaffinity(0)):
        raise ValueError('workers exceeds CPUs assigned to this process')
    output.mkdir(parents=True, exist_ok=True)
    if (output / 'batch.json').exists():
        if read(output / 'batch.json') != plan:
            raise ValueError('Frozen input inventory changed; do not mix batches')
    else:
        write(output / 'batch.json', plan)
    # Advisory POSIX lock releases on termination, without persistent stale locks.
    import fcntl
    with (output / 'run.lock').open('a') as lock:
        fcntl.flock(lock.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        old = resolve_shards(args.old_artifacts)
        oldchem = read(args.old_chemical)
        if oldchem['library_id'] != old['library_id']:
            raise ValueError('Old chemistry belongs to a different library')
        chemmap = {r['name']: r for r in oldchem['shards']}
        for r in chemmap.values():
            if not Path(r['path']).is_dir():
                r['path'] = str(args.old_chemical.resolve().parent / r['name'])
        bystem = {p.stem: p for p in files}
        oldsources = {}
        for row in old['shards']:
            p = bystem.get(row['name'])
            if p is None:
                raise ValueError('Original MOL2 required to restore stable identity: ' + row['name'])
            if digest(p) != read(Path(row['path'])/'manifest.json')['source_sha256']:
                raise ValueError('Original MOL2 changed: ' + str(p))
            oldsources[str(p)] = row
        dbpath = output / 'registry.sqlite3'
        private_registry(dbpath, old['library_id'])
        for path, row in oldsources.items():
            log('Restore old identity without recomputing descriptors: ' + row['name'])
            register_source(dbpath, old['library_id'], Path(path),
                            read(Path(row['path'])/'manifest.json')['source_sha256'], Path(row['path']))
        rows = list(old['shards']); completed_new = 0
        for item in plan['files']:
            check_source(item)
            p = Path(item['path'])
            if str(p) in oldsources:
                continue
            if args.max_new_files and completed_new >= args.max_new_files:
                log('Pilot complete. Repeat without --max-new-files to process the frozen full batch.')
                return
            if shutil.disk_usage(output).free < 32 * 2**30:
                raise ValueError('Stopped safely: less than 32 GiB free')
            log('Start ' + p.name)
            t = time.monotonic(); sha = digest(p)
            report = register_source(dbpath, old['library_id'], p, sha)
            check_source(item)
            if report['inserted']:
                if report['inserted'] > 1000000:
                    raise ValueError('Source exceeds the 1M-conformer memory guard: '+str(p)+
                                     '; prepare explicit bounded shards before launching this file')
                start = sum(int(r['conformers']) for r in rows)
                target = output / 'artifacts' / p.stem
                if target.exists():
                    m = verify(target)
                    if (m['source_sha256'], m['global_id_start'], m['conformers'], m['library_id']) != (
                            sha, start, report['inserted'], old['library_id']):
                        raise ValueError('Existing v1 does not match source/range')
                else:
                    recover_partial(output, target.parent / ('.'+target.name+'.partial'))
                    build_shard_parallel(dbpath, old['library_id'], p, output/'artifacts', start, workers=WORKERS)
                row = artifact_row(target)
                m = np.memmap(target/'meta.bin', dtype=META_DTYPE, mode='r')
                if len(m) != row['conformers'] or not np.array_equal(m['global_id'], np.arange(start, start+len(m))):
                    raise ValueError('Artifact global IDs are not contiguous')
                del m
                auxiliary(output, row, old['library_id'], dbpath)
                rows.append(row)
            completed_new += 1
            report['elapsed_this_run_seconds'] = time.monotonic()-t
            write(output/'progress'/f'{p.stem}.json', report)
            log(f'Ready {p.name}: inserted={report["inserted"]:,}, duplicates={report["duplicates"]:,}, seconds={report["elapsed_this_run_seconds"]:.1f}')
        # Only publish expanded catalogs after EVERY frozen source is complete.
        catalog = dict(format='aidd-conformer-artifact-catalog', version=1, library_id=old['library_id'],
                       conformers=sum(int(r['conformers']) for r in rows),
                       heavy_atoms=sum(int(r['heavy_atoms']) for r in rows),
                       features=sum(int(r['features']) for r in rows), shards=rows)
        write(output/'artifacts'/'catalog.json', catalog)
        chemical, pharma = [], []
        for row in rows:
            c, p = auxiliary(output, row, old['library_id'], dbpath,
                             old_chem=chemmap.get(row['name']) if row in old['shards'] else None)
            chemical.append(c); pharma.append(p)
        catpath = output/'artifacts'/'catalog.json'
        write(output/'chemical'/'catalog.json', dict(format='aidd-chemical-companion-catalog', version=1,
              library_id=old['library_id'], artifact_v1_catalog=str(catpath),
              artifact_v1_catalog_sha256=digest(catpath), shards=chemical))
        from aidd_agent.pharmacophore_index import INDEX_FORMAT
        write(output/'pharmacophore'/'catalog.json', dict(format=INDEX_FORMAT+'-catalog', version=1,
              library_id=old['library_id'], artifact_catalog=str(catpath), artifact_catalog_sha256=digest(catpath),
              bin_width_angstrom=0.5, max_distance_angstrom=20.0, conformers=catalog['conformers'], shards=pharma))
        make_indices(output, catalog, WORKERS)
        for item in plan['files']:
            check_source(item)
        if digest(args.old_artifacts) != plan['old_artifacts_sha256'] or digest(args.old_chemical) != plan['old_chemical_sha256']:
            raise ValueError('Original catalog changed during build')
        write(output/'COMPLETE.json', dict(status='precomputed_recall_calibration_pending',
              conformers=catalog['conformers'], shards=len(rows), batch_sha256=digest(output/'batch.json')))
        log(f'COMPLETE: {catalog["conformers"]:,} conformers. Recall calibration remains separate.')


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--old-artifacts', type=Path, required=True)
    p.add_argument('--old-chemical', type=Path, required=True)
    p.add_argument('--workers', type=int, default=20)
    p.add_argument('--max-new-files', type=int, default=0)
    p.add_argument('--run', action='store_true')
    args = p.parse_args()
    if args.workers < 1 or args.max_new_files < 0:
        p.error('workers must be positive; max-new-files must be nonnegative')
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'src'))
    for name in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
        os.environ[name] = '1'
    execute(args)


if __name__ == '__main__':
    main()
