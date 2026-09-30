"""Standalone NumPy-only export of every stored conformer of selected molecules.

Original coordinates and record names are preserved. No chemistry reconstruction,
search, docking, or pose transformation is performed. Resume uses sealed parts.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import sqlite3

import numpy as np


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda: f.read(8 * 1024 * 1024), b''):
            h.update(b)
    return h.hexdigest()


def read(path):
    return json.loads(Path(path).read_text(encoding='utf-8'))


def save(path, value):
    temp = path.with_suffix(path.suffix + '.partial')
    temp.write_text(json.dumps(value, indent=2), encoding='utf-8')
    temp.replace(path)


def require(ok, message):
    if not ok:
        raise ValueError(message)


def blocks(path):
    # Same text normalization as the original registry ingestion.
    current = []
    index = 0
    with path.open(encoding='utf-8', errors='replace', newline=None) as f:
        for line in f:
            if line.strip().upper() == '@<TRIPOS>MOLECULE':
                if current:
                    yield index, ''.join(current)
                    index += 1
                current = [line]
            elif current:
                current.append(line)
    if current:
        yield index, ''.join(current)


def run(batch, search, output, resume=False):
    batch, search, output = [Path(p).resolve() for p in (batch, search, output)]
    for source in (batch, search):
        require(not output.is_relative_to(source) and not source.is_relative_to(output),
                'Output must be separate from inputs')
    report = read(search / 'report.json')
    require(report['status'] == 'complete', 'Search is not complete')
    paths = [search / 'selected-molecules.npy', search / 'selected-conformers.npy',
             search / 'ranking.sqlite']
    expected = {**report['retrieval']['outputs'], **report['outputs']}
    hashes = {}
    for p in paths:
        print('Verifying ' + str(p), flush=True)
        hashes[str(p)] = sha(p)
        require(hashes[str(p)] == expected.get(str(p)), 'Input hash mismatch: ' + str(p))
    catalog_path = batch / 'artifacts/catalog.json'
    catalog = read(catalog_path)
    registry_stat = (batch / 'registry.sqlite3').stat()
    signature = dict(inputs=hashes, catalog_sha256=sha(catalog_path),
                     registry_size=registry_stat.st_size, registry_mtime_ns=registry_stat.st_mtime_ns,
                     script_sha256=sha(__file__), batch=str(batch), search=str(search))
    if output.exists():
        require(resume and (output / 'signature.json').exists(), 'Use fresh output or --resume')
        require(read(output / 'signature.json') == signature, 'Resume inputs or code changed')
    else:
        output.mkdir(parents=True)
        save(output / 'signature.json', signature)
    mids = np.load(paths[0], allow_pickle=False)
    gids = np.load(paths[1], allow_pickle=False)
    require(mids.ndim == gids.ndim == 1 and len(mids) > 0 and len(gids) > 0, 'Empty or invalid selection')
    require(mids.dtype.kind in 'SU' and gids.dtype.kind in 'iu', 'Invalid ID types')
    require(len(np.unique(mids)) == len(mids) and len(np.unique(gids)) == len(gids), 'Duplicate selected IDs')
    gids = np.sort(gids)
    require(gids[0] >= 0, 'Negative global ID')
    db = sqlite3.connect((output / 'manifest.sqlite').as_uri(), uri=True)
    db.execute('PRAGMA temp_store=FILE')
    db.execute('ATTACH DATABASE ? AS registry', ((batch / 'registry.sqlite3').as_uri() + '?mode=ro',))
    db.execute('ATTACH DATABASE ? AS ranks', (paths[2].as_uri() + '?mode=ro',))
    # Rebuild only an incomplete metadata plan. Completed export parts are separate.
    if not (output / 'plan.json').exists():
        db.executescript('DROP TABLE IF EXISTS selected; DROP TABLE IF EXISTS wanted; DROP TABLE IF EXISTS member;'
                         'CREATE TABLE member(mid TEXT PRIMARY KEY);'
                         'CREATE TABLE selected(gid INTEGER PRIMARY KEY,cid TEXT UNIQUE,mid TEXT);')
        decode = lambda x: bytes(x).decode() if isinstance(x, (bytes, np.bytes_)) else str(x)
        db.executemany('INSERT INTO member VALUES (?)', ((decode(m),) for m in mids))
        offset = 0
        for row in catalog['shards']:
            start, n = int(row['global_id_start']), int(row['conformers'])
            require(start == offset, 'Noncontiguous catalog')
            offset += n
            folder = Path(row['path'])
            if not folder.is_dir():
                folder = catalog_path.parent / row['name']
            manifest = read(folder / 'manifest.json')
            for name in ('conformer_ids.bin', 'molecule_ids.bin'):
                p = folder / name
                require(p.stat().st_size == n * 16, 'ID file size mismatch')
                require(sha(p) == manifest['files'][name]['sha256'], 'ID file hash mismatch')
            c = np.memmap(folder / 'conformer_ids.bin', dtype='S16', mode='r')
            m = np.memmap(folder / 'molecule_ids.bin', dtype='S16', mode='r')
            ids = gids[np.searchsorted(gids, start):np.searchsorted(gids, start + n)]
            db.executemany('INSERT INTO selected VALUES (?,?,?)',
                           ((int(g), decode(c[g-start]), decode(m[g-start])) for g in ids))
            del c, m
        require(gids[-1] < offset, 'Global ID outside catalog')
        require(db.execute('SELECT count(*) FROM selected').fetchone()[0] == len(gids), 'Missing selected IDs')
        require(db.execute('SELECT count(*) FROM selected s LEFT JOIN member m ON s.mid=m.mid WHERE m.mid IS NULL').fetchone()[0] == 0,
                'Unexpected selected molecule')
        require(db.execute('SELECT count(DISTINCT mid) FROM selected').fetchone()[0] == len(mids), 'Missing selected molecule')
        # Count all registry conformers of selected molecules, not just ranking representatives.
        require(db.execute('SELECT count(*) FROM registry.conformer c JOIN member m ON c.molecule_id=m.mid').fetchone()[0] == len(gids),
                'Selection does not cover all registered conformers of selected molecules')
        db.execute('CREATE TABLE wanted AS SELECT s.gid,s.cid,s.mid,r.rank,m.source_name AS molecule_name,'
                   'c.source_path,c.source_record_index,c.source_record_name,c.content_sha256 '
                   'FROM selected s JOIN registry.conformer c ON c.id=s.cid AND c.molecule_id=s.mid '
                   'JOIN registry.molecule m ON m.id=s.mid JOIN ranks.ranking r ON r.mid=s.mid')
        require(db.execute('SELECT count(*) FROM wanted').fetchone()[0] == len(gids), 'Registry/ranking identity mismatch')
        db.execute('CREATE UNIQUE INDEX wanted_source ON wanted(source_path,source_record_index)')
        db.execute('CREATE UNIQUE INDEX wanted_gid ON wanted(gid)')
        db.commit()
        save(output / 'plan.json', dict(molecules=len(mids), conformers=len(gids),
                                       manifest_sha256=sha(output / 'manifest.sqlite')))
    plan = read(output / 'plan.json')
    require(plan['molecules'] == len(mids) and plan['conformers'] == len(gids), 'Plan count mismatch')
    require(plan['manifest_sha256'] == sha(output / 'manifest.sqlite'), 'Plan database changed')
    names = output / 'molecules.csv'
    with names.open('w', encoding='utf-8', newline='') as f:
        w = csv.writer(f)
        w.writerow(['molecule_id', 'rank', 'original_name', 'conformers'])
        w.writerows(db.execute('SELECT mid,rank,molecule_name,count(*) FROM wanted GROUP BY mid ORDER BY rank'))
    sources = [r[0] for r in db.execute('SELECT DISTINCT source_path FROM wanted ORDER BY source_path')]
    receipts = []
    total = 0
    for number, source in enumerate(sources, 1):
        stem = f'part-{number:05d}'
        receipt = output / (stem + '.json')
        count = db.execute('SELECT count(*) FROM wanted WHERE source_path=?', (source,)).fetchone()[0]
        if receipt.exists():
            part = read(receipt)
            require(part['source'] == source and part['conformers'] == count, 'Part identity mismatch')
            for name, digest in part['outputs'].items():
                require(sha(output / name) == digest, 'Completed part hash mismatch')
        else:
            target, mapping = output / (stem + '.mol2'), output / (stem + '.csv')
            raw_temp, map_temp = Path(str(target) + '.partial'), Path(str(mapping) + '.partial')
            cursor = db.execute('SELECT * FROM wanted WHERE source_path=? ORDER BY source_record_index', (source,))
            row = next(cursor, None)
            written = 0
            with raw_temp.open('wb') as raw_out, map_temp.open('w', encoding='utf-8', newline='') as meta:
                writer = csv.writer(meta)
                writer.writerow(['global_id', 'conformer_id', 'molecule_id', 'rank', 'molecule_name',
                                 'source_path', 'source_record_index', 'record_name', 'content_sha256',
                                 'output_record_index', 'byte_offset', 'byte_length'])
                for idx, text in blocks(Path(source)):
                    if row is None:
                        break
                    if idx != row[6]:
                        continue
                    payload = text.encode('utf-8')
                    require(hashlib.sha256(payload).hexdigest() == row[8], f'Source record changed: {source}:{idx}')
                    writer.writerow([*row, written, raw_out.tell(), len(payload)])
                    raw_out.write(payload)
                    # Separate a final record lacking newline when concatenating.
                    if not payload.endswith(b'\n'):
                        raw_out.write(b'\n')
                    written += 1
                    row = next(cursor, None)
            require(row is None and written == count, 'Missing source records: ' + source)
            raw_temp.replace(target)
            map_temp.replace(mapping)
            part = dict(source=source, conformers=count,
                        outputs={p.name: sha(p) for p in (target, mapping)})
            save(receipt, part)
        total += count
        receipts.append(dict(receipt=receipt.name, sha256=sha(receipt)))
        save(output / 'progress.json', dict(completed_sources=number, sources=len(sources),
                                            exported_conformers=total, expected_conformers=len(gids)))
        print(f'Exported {total:,}/{len(gids):,} conformers; sources {number}/{len(sources)}', flush=True)
    require(total == len(gids), 'Final count mismatch')
    db.close()
    result = dict(status='complete', molecules=len(mids), conformers=total, parts=len(sources),
                  scope='All selected stored conformers; original source coordinates, not docked poses',
                  validation='Selected IDs, registry identity, all-conformer counts and every source record hash checked',
                  outputs={p.name: sha(p) for p in (names, output / 'manifest.sqlite')}, parts_manifest=receipts)
    save(output / 'report.json', result)
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    for key in ('batch', 'search', 'output'):
        parser.add_argument('--' + key, required=True, type=Path)
    parser.add_argument('--resume', action='store_true')
    args = parser.parse_args()
    run(args.batch, args.search, args.output, args.resume)
