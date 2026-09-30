"""Reconstruct conformer-level RRF from sealed pose chunks and export original MOL2.

Requires NumPy and export_selected_conformers.py beside this script. No rescoring.
"""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import sqlite3

import numpy as np

from export_selected_conformers import read, require, save, sha


def verify(path, expected):
    require(sha(path) == expected, 'Hash mismatch: ' + str(path))


def rank_template(db, template, quota, k):
    db.execute('DROP TABLE IF EXISTS trank')
    db.execute('CREATE TABLE trank AS SELECT gid,contact,gaussian,'
               'ROW_NUMBER() OVER(ORDER BY contact DESC,gaussian DESC,gid) AS rank FROM scores')
    db.execute('''INSERT INTO aggregate
        SELECT gid,1.0/(?+rank),1,contact>0,rank<=?,rank,?,contact,gaussian
        FROM trank WHERE true
        ON CONFLICT(gid) DO UPDATE SET
        rrf=aggregate.rrf+excluded.rrf,support=aggregate.support+1,
        has_contact=MAX(aggregate.has_contact,excluded.has_contact),
        tier=MAX(aggregate.tier,excluded.tier),
        best_template=CASE WHEN excluded.best_rank<aggregate.best_rank THEN excluded.best_template ELSE aggregate.best_template END,
        contact=CASE WHEN excluded.best_rank<aggregate.best_rank THEN excluded.contact ELSE aggregate.contact END,
        gaussian=CASE WHEN excluded.best_rank<aggregate.best_rank THEN excluded.gaussian ELSE aggregate.gaussian END,
        best_rank=MIN(aggregate.best_rank,excluded.best_rank)''', (k, quota, template))
    db.execute('DROP TABLE trank')
    db.execute('DELETE FROM scores')


def run(search, exported, output, count=1000000, resume=False, part_size=10000):
    search, exported, output = [Path(p).resolve() for p in (search, exported, output)]
    require(count > 0 and part_size > 0, 'Positive count and part size required')
    for source in (search, exported):
        require(not output.is_relative_to(source) and not source.is_relative_to(output), 'Output overlaps input')
    sr, er = read(search / 'report.json'), read(exported / 'report.json')
    require(sr['status'] == er['status'] == 'complete', 'Completed search and export required')
    signature = dict(search=str(search), exported=str(exported), search_report=sha(search / 'report.json'),
                     export_report=sha(exported / 'report.json'), count=count, part_size=part_size,
                     implementation=sha(__file__), protocol=sha(search / 'protocol.json'),
                     design=sha(search / 'adopted-design/report.json'))
    if output.exists():
        require(resume and read(output / 'signature.json') == signature, 'Use fresh output or matching --resume')
    else:
        output.mkdir(parents=True)
        save(output / 'signature.json', signature)
    gidpath = search / 'selected-conformers.npy'
    verify(gidpath, sr['retrieval']['outputs'][str(gidpath)])
    gids = np.sort(np.load(gidpath, allow_pickle=False))
    require(len(gids) == er['conformers'] and len(np.unique(gids)) == len(gids), 'Selection/export count mismatch')
    verify(exported / 'manifest.sqlite', er['outputs']['manifest.sqlite'])
    design = read(search / 'adopted-design/report.json')
    protocol = read(search / 'protocol.json')
    quota, k = int(protocol['template_quota']), float(protocol['rrf_k'])
    require(quota > 0 and math.isfinite(k) and k > 0, 'Invalid saved ranking policy')
    templates = [t['query_id'] for t in design['templates']]
    require(len(templates) == len(set(templates)) and len(templates) > 0, 'Invalid template IDs')
    require(templates == [t['template'] for t in sr['retrieval']['templates']], 'Template order mismatch')
    # Export signature links the source selection to the completed MOL2 export.
    es = read(exported / 'signature.json')
    require(es['inputs'][str(gidpath)] == sr['retrieval']['outputs'][str(gidpath)], 'Export selection lineage mismatch')
    db = sqlite3.connect((output / 'selection.sqlite').as_uri(), uri=True)
    db.execute('PRAGMA temp_store=FILE')
    db.execute('ATTACH DATABASE ? AS source', ((exported / 'manifest.sqlite').as_uri() + '?mode=ro',))
    db.executescript('''CREATE TABLE IF NOT EXISTS aggregate(
        gid INTEGER PRIMARY KEY,rrf REAL,support INT,has_contact INT,tier INT,
        best_rank INT,best_template TEXT,contact REAL,gaussian REAL);
        CREATE TABLE IF NOT EXISTS scores(gid INTEGER PRIMARY KEY,contact REAL,gaussian REAL);
        CREATE TABLE IF NOT EXISTS completed(template TEXT PRIMARY KEY,receipt_digest TEXT,scored INT);''')
    # Compact exact identity map for all selected global IDs.
    mapped_gids = np.empty(len(gids), dtype=np.int64)
    mids = np.empty(len(gids), dtype='S16')
    cids = np.empty(len(gids), dtype='S16')
    n = 0
    for gid, mid, cid in db.execute('SELECT gid,mid,cid FROM source.wanted ORDER BY gid'):
        require(n < len(gids), 'Export has extra conformers')
        mapped_gids[n], mids[n] = gid, mid.encode()
        cids[n] = cid.encode()
        n += 1
    require(n == len(gids) and np.array_equal(mapped_gids, gids), 'Export global IDs differ')
    for ti, template in enumerate(templates):
        root = search / 'chunks' / f'{ti:02d}'
        files = sorted(root.glob('*.npz'))
        require(files, 'Missing per-conformer score chunks: ' + str(root))
        digest = hashlib.sha256()
        for p in files:
            receipt = p.with_suffix('.receipt.json')
            require(receipt.is_file(), 'Missing chunk receipt: ' + str(receipt))
            digest.update(p.name.encode())
            digest.update(sha(receipt).encode())
        receipt_digest = digest.hexdigest()
        done = db.execute('SELECT receipt_digest FROM completed WHERE template=?', (template,)).fetchone()
        if done:
            require(done[0] == receipt_digest, 'Completed template receipts changed')
            print('Reusing completed template ' + template, flush=True)
            continue
        print(f'Reading template {ti+1}/{len(templates)}: {template}; {len(files):,} chunks', flush=True)
        with db:
            db.execute('DELETE FROM scores')
            offset = 0
            for fi, p in enumerate(files):
                receipt = read(p.with_suffix('.receipt.json'))
                poses = p.with_suffix('.poses.jsonl')
                for artifact in (p, poses):
                    require(str(artifact) in receipt['files'], 'Chunk artifact missing from receipt')
                    verify(artifact, receipt['files'][str(artifact)])
                with np.load(p, allow_pickle=False) as z:
                    ids = z['global_ids']
                    require(len(ids) > 0 and np.array_equal(ids, gids[offset:offset+len(ids)]), 'Missing/duplicate/out-of-order chunk IDs')
                    chunk_mids = np.asarray(z['molecule_ids'], dtype='S16')
                    require(np.array_equal(chunk_mids, mids[offset:offset+len(ids)]), 'Chunk molecule identity mismatch')
                    stages = z['stage_levels']
                    require(len(stages) == len(ids) and np.isin(stages, np.arange(8)).all(), 'Chunk stage coverage mismatch')
                rows = []
                seen = set()
                with poses.open(encoding='utf-8') as f:
                    for line in f:
                        row = json.loads(line)
                        gid = row['global_id']
                        pos = int(np.searchsorted(ids, gid))
                        require(pos < len(ids) and int(ids[pos]) == gid, 'Pose outside chunk')
                        require(row['molecule_id'].encode() == chunk_mids[pos], 'Pose identity mismatch')
                        require(row['conformer_id'].encode() == cids[offset+pos], 'Pose conformer identity mismatch')
                        require(stages[pos] == 7, 'Pose stage mismatch')
                        require(gid not in seen, 'Multiple poses per conformer/template; unsupported historical schema')
                        seen.add(gid)
                        contact, gaussian = float(row['optional_score']), float(row['gaussian_same_pose'])
                        require(math.isfinite(contact) and math.isfinite(gaussian), 'Nonfinite ranking terms')
                        rows.append((gid, contact, gaussian))
                require(len(seen) == int(np.count_nonzero(stages == 7)), 'Missing scored poses from chunk')
                db.executemany('INSERT INTO scores VALUES (?,?,?)', rows)
                offset += len(ids)
                if (fi+1) % 2000 == 0:
                    print(f'  {offset:,}/{len(gids):,} conformers checked', flush=True)
            require(offset == len(gids), 'Incomplete template computation coverage')
            scored = db.execute('SELECT count(*) FROM scores').fetchone()[0]
            rank_template(db, template, quota, k)
            db.execute('INSERT INTO completed VALUES (?,?,?)', (template, receipt_digest, scored))
        print(f'Template complete; {scored:,} scored conformers', flush=True)
    scored = db.execute('SELECT count(*) FROM aggregate').fetchone()[0]
    require(scored >= count, f'Only {scored} conformers have saved scores; requested {count}. No unscored padding.')
    if not (output / 'ranking.json').exists():
        with db:
            db.execute('DROP TABLE IF EXISTS main.selected')
            db.execute('CREATE TABLE main.selected AS SELECT ROW_NUMBER() OVER(ORDER BY has_contact DESC,tier DESC,rrf DESC,gid) AS rank,* '
                       'FROM aggregate ORDER BY has_contact DESC,tier DESC,rrf DESC,gid LIMIT ?', (count,))
            db.execute('CREATE UNIQUE INDEX selected_gid ON selected(gid)')
        save(output / 'ranking.json', dict(scored=scored, selected=count, database_sha256=sha(output / 'selection.sqlite')))
    verify(output / 'selection.sqlite', read(output / 'ranking.json')['database_sha256'])
    # Build byte locators using hash-verified export manifests; read original coordinates only.
    wanted = {r[0] for r in db.execute('SELECT gid FROM selected')}
    locators = {}
    for item in er['parts_manifest']:
        rp = exported / item['receipt']
        verify(rp, item['sha256'])
        part = read(rp)
        csv_names = [name for name in part['outputs'] if name.endswith('.csv')]
        mol_names = [name for name in part['outputs'] if name.endswith('.mol2')]
        require(len(csv_names) == len(mol_names) == 1, 'Ambiguous export part')
        mapping, mol2 = exported / csv_names[0], exported / mol_names[0]
        require(mapping.resolve().parent == exported and mol2.resolve().parent == exported, 'Unsafe part path')
        verify(mapping, part['outputs'][csv_names[0]])
        used = False
        with mapping.open(encoding='utf-8', newline='') as f:
            for row in csv.DictReader(f):
                gid = int(row['global_id'])
                if gid in wanted:
                    require(gid not in locators, 'Duplicate exported conformer')
                    locators[gid] = (mol2, row)
                    used = True
        if used:
            verify(mol2, part['outputs'][mol_names[0]])
    require(len(locators) == count, 'Selected conformers missing from export')
    print(f'Exporting {count:,} conformers in rank order', flush=True)
    receipts = []
    unique_molecules = set()
    cursor = db.execute('SELECT s.rank,s.gid,w.cid,w.mid,w.molecule_name,s.rrf,s.support,s.has_contact,s.tier,'
                        's.best_template,s.contact,s.gaussian,w.content_sha256 '
                        'FROM selected s JOIN source.wanted w ON s.gid=w.gid ORDER BY s.rank')
    fields = ['rank','global_id','conformer_id','molecule_id','name','rrf','template_support','has_contact','quota_tier',
              'representative_template','contact','gaussian','source_sha256','output_record_index','byte_offset','byte_length']
    written = 0
    for part_index in range((count+part_size-1)//part_size):
        rows = cursor.fetchmany(part_size)
        require(rows, 'Missing selected rank rows')
        unique_molecules.update(r[3] for r in rows)
        stem = f'top-{part_index+1:05d}'
        receipt_path = output / (stem + '.json')
        if receipt_path.exists():
            part = read(receipt_path)
            require(part['start_rank'] == rows[0][0] and part['end_rank'] == rows[-1][0], 'Output rank range changed')
            for name, h in part['outputs'].items():
                verify(output / name, h)
        else:
            molpath, csvpath = output / (stem + '.mol2'), output / (stem + '.csv')
            mt, ct = Path(str(molpath)+'.partial'), Path(str(csvpath)+'.partial')
            with mt.open('wb') as m, ct.open('w', encoding='utf-8', newline='') as c:
                writer = csv.writer(c)
                writer.writerow(fields)
                # Read in source order within each output part, then write rank order.
                payloads = {}
                ordered = sorted(rows, key=lambda r: (str(locators[r[1]][0]), int(locators[r[1]][1]['byte_offset'])))
                current, handle = None, None
                try:
                    for r in ordered:
                        path, locator = locators[r[1]]
                        require(locator['conformer_id'] == r[2] and locator['molecule_id'] == r[3], 'Locator identity mismatch')
                        if current != path:
                            if handle:
                                handle.close()
                            handle, current = path.open('rb'), path
                        handle.seek(int(locator['byte_offset']))
                        payload = handle.read(int(locator['byte_length']))
                        require(hashlib.sha256(payload).hexdigest() == r[12], 'Selected MOL2 record hash mismatch')
                        payloads[r[1]] = payload
                finally:
                    if handle:
                        handle.close()
                for i, r in enumerate(rows):
                    payload = payloads.pop(r[1])
                    writer.writerow([*r, i, m.tell(), len(payload)])
                    m.write(payload)
                    if not payload.endswith(b'\n'):
                        m.write(b'\n')
            mt.replace(molpath)
            ct.replace(csvpath)
            part = dict(start_rank=rows[0][0], end_rank=rows[-1][0], conformers=len(rows),
                        outputs={p.name: sha(p) for p in (molpath, csvpath)})
            save(receipt_path, part)
        receipts.append(dict(receipt=receipt_path.name, sha256=sha(receipt_path)))
        written += len(rows)
        print(f'Exported {written:,}/{count:,}', flush=True)
    require(written == count, 'Final export count mismatch')
    summary = output / 'selected-conformers.csv'
    temp_summary = output / 'selected-conformers.csv.partial'
    with temp_summary.open('w', encoding='utf-8', newline='') as f:
        writer = csv.writer(f)
        writer.writerow([*fields, 'mol2_file'])
        summary_count = 0
        for item in receipts:
            stem = Path(item['receipt']).stem
            with (output / (stem + '.csv')).open(encoding='utf-8', newline='') as part_csv:
                reader = csv.reader(part_csv)
                require(next(reader) == fields, 'Output CSV schema mismatch')
                for row in reader:
                    summary_count += 1
                    require(int(row[0]) == summary_count, 'Output ranks are not contiguous')
                    writer.writerow([*row, stem + '.mol2'])
    require(summary_count == count, 'Combined manifest count mismatch')
    temp_summary.replace(summary)
    db.close()
    result = dict(status='complete', input_conformers=len(gids), scored_conformers=scored,
                  unscored_conformers=len(gids)-scored, exported_conformers=written,
                  distinct_molecules=len(unique_molecules), templates=templates,
                  policy=dict(unit='conformer', within_template='contact DESC, gaussian DESC, global_id ASC',
                              across_templates='any-positive-contact tier, quota tier, RRF DESC, global_id ASC',
                              template_quota=quota, rrf_k=k),
                  scope='New conformer-level ranking over saved scored ANN candidates; original coordinates, no docking',
                  parts_manifest=receipts, outputs={'selection.sqlite': sha(output / 'selection.sqlite'),
                                                   summary.name: sha(summary)})
    save(output / 'report.json', result)
    return result


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    for key in ('search', 'exported', 'output'):
        p.add_argument('--'+key, type=Path, required=True)
    p.add_argument('--count', type=int, default=1000000)
    p.add_argument('--part-size', type=int, default=10000)
    p.add_argument('--resume', action='store_true')
    a = p.parse_args()
    run(a.search, a.exported, a.output, a.count, a.resume, a.part_size)
