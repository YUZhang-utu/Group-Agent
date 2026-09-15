"""E033: read-only library acceptance and bounded-memory USRCAT calibration."""
from __future__ import annotations

import argparse
from contextlib import contextmanager
import csv
import hashlib
import json
import os
from pathlib import Path
import platform
import sqlite3
import sys
import time

import numpy as np

from .conformer_artifacts import META_DTYPE
from .chemical_companion import CHEM_META_DTYPE


def log(message):
    print(time.strftime('%Y-%m-%d %H:%M:%S'), message, flush=True)


def read(path):
    return json.loads(Path(path).read_text(encoding='utf-8'))


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def write(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False), encoding='utf-8')


def require(condition, message):
    if not condition:
        raise ValueError(message)


@contextmanager
def mapped(path, dtype, shape=None):
    array = np.memmap(path, dtype=dtype, mode='r', shape=shape)
    try:
        yield array
    finally:
        array._mmap.close()


def shard_path(row, catalog_path):
    p = Path(row['path'])
    if not p.is_dir():
        p = catalog_path.parent / row['name']
    return p.resolve()


def verify_manifest(path, full):
    manifest = read(path / 'manifest.json')
    payloads = list(manifest.get('files', {}).items())
    payloads += [(r['file'], dict(bytes=int(r['count']) * 8, sha256=r['sha256']))
                 for r in manifest.get('keys', [])]
    for name, record in payloads:
        target = (path / name).resolve()
        require(target.is_relative_to(path.resolve()), 'Manifest path escapes shard')
        require(target.stat().st_size == record['bytes'], f'Byte count mismatch: {target}')
        if full:
            require(sha(target) == record['sha256'], f'Checksum mismatch: {target}')
    return manifest


def accept_library(batch, full=True):
    """Verify immutable manifests/counts/identity, without parsing original MOL2."""
    started = time.perf_counter()
    complete_path = batch / 'COMPLETE.json'
    require(complete_path.is_file(),
            'E032 is not complete: COMPLETE.json missing. Resume precompute first; registration alone is insufficient.')
    complete, inventory = read(complete_path), read(batch / 'batch.json')
    require(complete['batch_sha256'] == sha(batch / 'batch.json'), 'COMPLETE/inventory hash mismatch')
    catpaths = {k: batch / k / 'catalog.json' for k in ('artifacts', 'chemical', 'pharmacophore')}
    catalogs = {k: read(p) for k, p in catpaths.items()}
    catalog = catalogs['artifacts']; total = int(catalog['conformers'])
    require(total > 1, 'Need at least two conformers')
    require(complete['conformers'] == total and complete['shards'] == len(catalog['shards']),
            'COMPLETE/catalog count mismatch')
    maps = {}
    for kind, cat in catalogs.items():
        require(cat['library_id'] == catalog['library_id'], f'{kind}: cross-library catalog')
        maps[kind] = {r['name']: r for r in cat['shards']}
        require(len(maps[kind]) == len(cat['shards']), f'{kind}: duplicate shard names')
        require(set(maps[kind]) == {r['name'] for r in catalog['shards']}, f'{kind}: shard coverage mismatch')
        if kind != 'artifacts':
            require(cat['artifact_catalog_sha256'] == sha(catpaths['artifacts'])
                    if kind == 'pharmacophore' else
                    cat['artifact_v1_catalog_sha256'] == sha(catpaths['artifacts']),
                    f'{kind}: artifact catalog lineage mismatch')
            require(sum(int(r['conformers']) for r in cat['shards']) == total, f'{kind}: total mismatch')
    with sqlite3.connect((batch / 'registry.sqlite3').resolve().as_uri() + '?mode=ro', uri=True) as db:
        db.row_factory = sqlite3.Row
        require(db.execute('PRAGMA quick_check').fetchone()[0] == 'ok', 'SQLite quick_check failed')
        sources = {r['path']: dict(r) for r in db.execute('SELECT * FROM batch_source')}
        require(set(sources) == {r['path'] for r in inventory['files']}, 'Not all frozen sources registered')
        by_stem = {Path(p).stem: p for p in sources}
        require(len(by_stem) == len(sources), 'Repeated source stems')
        require(all(r['records'] == r['inserted'] + r['duplicates'] for r in sources.values()),
                'Source accounting mismatch')
        count = db.execute('SELECT COUNT(*) FROM conformer c JOIN molecule m ON m.id=c.molecule_id '
                           'WHERE m.library_id=?', (catalog['library_id'],)).fetchone()[0]
        molecules = db.execute('SELECT COUNT(*) FROM molecule WHERE library_id=?',
                               (catalog['library_id'],)).fetchone()[0]
        conflicts = db.execute('SELECT COUNT(*) FROM conformer c JOIN molecule m ON m.id=c.molecule_id '
                               'WHERE m.library_id=? AND c.conformer_index<0',
                               (catalog['library_id'],)).fetchone()[0]
        require(count == total == sum(r['inserted'] for r in sources.values()), 'Registry/artifact total mismatch')
    expected = 0
    for position, row in enumerate(catalog['shards'], 1):
        log(f'Accept shard {position}/{len(catalog["shards"])}: {row["name"]}')
        n = int(row['conformers']); require(n > 0 and row['global_id_start'] == expected, 'Noncontiguous catalog')
        paths, manifests = {}, {}
        for kind in catalogs:
            peer = maps[kind][row['name']]; p = shard_path(peer, catpaths[kind])
            require(sha(p / 'manifest.json') == peer['manifest_sha256'], f'{kind}: manifest hash mismatch')
            m = verify_manifest(p, full)
            require(m['library_id'] == catalog['library_id'] and m['conformers'] == n
                    and m['global_id_start'] == expected, f'{kind}: range/library mismatch')
            paths[kind], manifests[kind] = p, m
        # E032 restores old IDs from relocated MOL2 by unique stem + frozen
        # content hash; immutable old manifests may still contain the old path.
        source = by_stem.get(row['name'])
        require(source in sources and sources[source]['inserted'] == n, 'Source/shard registration mismatch')
        require(manifests['artifacts']['source_sha256'] == sources[source]['sha'], 'Source hash lineage mismatch')
        require(manifests['chemical']['source_sha256'] == sources[source]['sha'], 'Chemical source hash mismatch')
        for kind, key in [('chemical', 'artifact_v1_manifest_sha256'), ('pharmacophore', 'source_manifest_sha256')]:
            require(manifests[kind][key] == row['manifest_sha256'], f'{kind}: source manifest lineage mismatch')
        a, c = paths['artifacts'], paths['chemical']
        for filename, dtype in [('molecule_ids.bin', 'S16'), ('conformer_ids.bin', 'S16')]:
            require((a / filename).stat().st_size == n * 16, 'ID array size mismatch')
            require((c / filename).stat().st_size == n * 16, 'Chemical ID array size mismatch')
            with mapped(a / filename, dtype) as left, mapped(c / filename, dtype) as right:
                for start in range(0, n, 65536):
                    require(np.array_equal(left[start:start+65536], right[start:start+65536]), 'Cross-stage ID mismatch')
        with sqlite3.connect((batch / 'registry.sqlite3').resolve().as_uri() + '?mode=ro', uri=True) as db:
            cursor = db.execute('SELECT c.id,c.molecule_id FROM conformer c JOIN molecule m ON m.id=c.molecule_id '
                                'WHERE m.library_id=? AND c.source_path=? ORDER BY c.source_record_index',
                                (catalog['library_id'], source))
            with mapped(a / 'conformer_ids.bin', 'S16') as confs, mapped(a / 'molecule_ids.bin', 'S16') as mols:
                offset = 0
                while True:
                    registered = cursor.fetchmany(65536)
                    if not registered:
                        break
                    ids = np.asarray(registered, dtype='S16')
                    require(np.array_equal(confs[offset:offset+len(ids)], ids[:, 0]) and
                            np.array_equal(mols[offset:offset+len(ids)], ids[:, 1]), 'Registry/artifact identity mismatch')
                    offset += len(ids)
                require(offset == n, 'Registry source count mismatch')
        require((a / 'meta.bin').stat().st_size == n * META_DTYPE.itemsize, 'Artifact metadata size mismatch')
        require((c / 'chem-meta.bin').stat().st_size == n * CHEM_META_DTYPE.itemsize, 'Chemical metadata size mismatch')
        with mapped(a / 'meta.bin', META_DTYPE) as left, mapped(c / 'chem-meta.bin', CHEM_META_DTYPE) as right:
            for start in range(0, n, 65536):
                truth = np.arange(expected+start, expected+min(n, start+65536))
                require(np.array_equal(left['global_id'][start:start+65536], truth)
                        and np.array_equal(right['global_id'][start:start+65536], truth), 'Global ID mapping mismatch')
        require((a / 'usrcat.f32.bin').stat().st_size == n * 60 * 4, 'USRCAT size mismatch')
        row['path'] = str(a); expected += n
    require(expected == total, 'Catalog total mismatch')
    fm = read(batch / 'faiss' / 'manifest.json')
    require(fm['final_ntotal'] == total and fm['catalog_sha256'] == sha(catpaths['artifacts']), 'FAISS lineage mismatch')
    if full:
        require(sha(batch / 'faiss' / 'index.faiss') == fm['index_sha256'], 'Final FAISS checksum mismatch')
    require(sha(batch / 'faiss' / 'transform.npz') == sha(batch / 'faiss' / 'trained' / 'transform.npz'),
            'Transform differs from trained transform')
    training = verify_manifest(batch / 'faiss' / 'trained', full)
    require(training['catalog_sha256'] == sha(catpaths['artifacts']), 'Training catalog lineage mismatch')
    return catalog, dict(status='passed' if full else 'metadata_passed_hashes_not_checked',
                         library_conformers=total, library_molecules=molecules,
                         registered_sources=len(sources), shards=len(catalog['shards']),
                         source_records=sum(r['records'] for r in sources.values()),
                         exact_duplicates=sum(r['duplicates'] for r in sources.values()),
                         preserved_index_conflicts=conflicts, wall_seconds=time.perf_counter()-started,
                         hashes='full_generated_payloads' if full else 'manifests_only',
                         original_mol2_bytes_rehashed=False, chemical_identity_qc='not_assessed',
                         sqlite_check='quick_check', faiss_id_check='count_and_returned_id_range; not_full_ID_permutation')


class Corpus:
    def __init__(self, catalog):
        self.rows = catalog['shards']
        self.total = int(catalog['conformers'])
        self.starts = np.asarray([r['global_id_start'] for r in self.rows], dtype=np.int64)

    def fetch(self, ids, vectors=False):
        ids = np.asarray(ids, dtype=np.int64)
        require(ids.ndim == 1 and np.all((ids >= 0) & (ids < self.total)), 'Invalid retrieved IDs')
        groups = np.searchsorted(self.starts, ids, side='right') - 1
        mols = np.empty(len(ids), dtype='S16')
        raw = np.empty((len(ids), 60), dtype=np.float32) if vectors else None
        for group in np.unique(groups):
            slots = np.flatnonzero(groups == group); row = self.rows[int(group)]
            local = ids[slots] - row['global_id_start']; p = Path(row['path'])
            with mapped(p / 'molecule_ids.bin', 'S16') as data:
                mols[slots] = data[local]
            if vectors:
                with mapped(p / 'usrcat.f32.bin', '<f4', (row['conformers'], 60)) as data:
                    raw[slots] = data[local]
        return mols, raw


def topk(distances, ids, k):
    """Deterministic distance/ID Top-K, including ties at the partition edge."""
    finite = np.isfinite(distances); distances, ids = distances[finite], ids[finite]
    if len(ids) > k:
        edge = np.partition(distances, k-1)[k-1]
        low = np.flatnonzero(distances < edge)
        equal = np.flatnonzero(distances == edge)
        equal = equal[np.argsort(ids[equal], kind='stable')[:k-len(low)]]
        chosen = np.concatenate((low, equal)); distances, ids = distances[chosen], ids[chosen]
    order = np.lexsort((ids, distances))
    return distances[order], ids[order]


def exact_truth(corpus, queries, excluded, mean, std, k, chunk_size, progress=None):
    best = [(np.empty(0, dtype=np.float32), np.empty(0, dtype=np.int64)) for _ in queries]
    excluded_counts = [0] * len(queries); compute = [0.] * len(queries)
    started = time.perf_counter()
    for row in corpus.rows:
        p, n = Path(row['path']), int(row['conformers'])
        log(f'Exact descriptor truth: {row["name"]}')
        with mapped(p / 'usrcat.f32.bin', '<f4', (n, 60)) as vectors, mapped(p / 'molecule_ids.bin', 'S16') as mols:
            for start in range(0, n, chunk_size):
                raw = np.asarray(vectors[start:start+chunk_size], dtype=np.float32)
                require(np.all(np.isfinite(raw)), f'Nonfinite descriptors: {p}')
                data = (raw - mean) / std
                require(np.all(np.isfinite(data)), 'Nonfinite transformed descriptors')
                ids = np.arange(row['global_id_start']+start, row['global_id_start']+start+len(data), dtype=np.int64)
                for qi, query in enumerate(queries):
                    t = time.perf_counter(); delta = data-query
                    dist = np.einsum('ij,ij->i', delta, delta)
                    if excluded[qi]:
                        mask = mols[start:start+len(data)] == excluded[qi]
                        excluded_counts[qi] += int(mask.sum()); dist[mask] = np.inf
                    d, ix = topk(dist, ids, k)
                    best[qi] = topk(np.concatenate((best[qi][0], d)), np.concatenate((best[qi][1], ix)), k)
                    compute[qi] += time.perf_counter()-t
        if progress:
            progress(row['name'])
    return best, excluded_counts, dict(wall_seconds=time.perf_counter()-started,
                                      per_query_distance_topk_seconds=compute,
                                      timing_note='Shared single full-corpus streaming pass; per-query compute excludes shared I/O/normalization.')


def unique_representatives(ids, molecules, cap):
    _, positions = np.unique(molecules, return_index=True)
    positions = np.sort(positions)[:cap]
    return ids[positions], molecules[positions]


def retained_counts(ids, mols, total, molecule_total):
    unique = len(np.unique(mols))
    return dict(retained_conformers=len(ids), retained_molecules=unique,
                conformer_reduction_fraction=1-len(ids)/total,
                molecule_reduction_fraction=1-unique/molecule_total)


def evaluate_setting(index, corpus, query, exclude, excluded_count, truth, mean, std,
                     nprobe, budget, repeats, truth_ks, caps, molecule_total):
    index.nprobe = nprobe
    search_k = min(corpus.total, budget + excluded_count)
    t = time.perf_counter(); index.search(query[None], search_k)
    first_call = time.perf_counter()-t
    timings = []
    for _ in range(repeats):
        t = time.perf_counter(); _, found = index.search(query[None], search_k)
        timings.append(time.perf_counter()-t)
    ids = found[0]; ids = ids[ids >= 0]
    require(len(np.unique(ids)) == len(ids), 'FAISS returned duplicate global IDs')
    t = time.perf_counter(); mols, raw = corpus.fetch(ids, vectors=True)
    eligible = mols != exclude if exclude else np.ones(len(ids), dtype=bool)
    ids, mols, raw = ids[eligible][:budget], mols[eligible][:budget], raw[eligible][:budget]
    delta = (raw-mean)/std-query
    distances = np.einsum('ij,ij->i', delta, delta)
    require(np.all(np.isfinite(distances)), 'Nonfinite candidate distances')
    order = np.lexsort((ids, distances)); ids, mols, distances = ids[order], mols[order], distances[order]
    post_seconds = time.perf_counter()-t
    truth_dist, truth_ids = truth
    truth_mols, _ = corpus.fetch(truth_ids)
    recall = {}
    for k in truth_ks:
        actual = min(k, len(truth_ids))
        require(actual > 0, 'No eligible exact neighbors after query-molecule exclusion')
        strict = len(np.intersect1d(ids, truth_ids[:actual])) / actual
        ties = min(actual, int(np.sum(distances <= truth_dist[actual-1]))) / actual
        recall[f'top{k}_strict'] = strict
        recall[f'top{k}_boundary_tie_aware'] = ties
    cap_results, cap_arrays = [], {}
    reference_molecules = np.unique(truth_mols)
    for cap in caps:
        chosen, chosen_mols = unique_representatives(ids, mols, cap)
        cap_results.append(dict(molecule_cap=cap, **retained_counts(chosen, chosen_mols, corpus.total, molecule_total),
                                exact_truth_molecule_coverage=len(np.intersect1d(chosen_mols, reference_molecules))/len(reference_molecules)))
        cap_arrays[f'cap_{cap}_global_ids'] = chosen
        cap_arrays[f'cap_{cap}_molecule_ids'] = chosen_mols
    result = dict(nprobe=nprobe, requested_candidate_budget=budget, faiss_requested_k=search_k,
                  first_call_seconds=first_call, warm_search_seconds=timings,
                  warm_search_p50_seconds=float(np.median(timings)),
                  warm_search_p95_seconds=float(np.percentile(timings, 95)),
                  candidate_fetch_and_exact_descriptor_rerank_seconds=post_seconds,
                  warm_search_plus_single_fetch_rerank_seconds=float(np.median(timings))+post_seconds,
                  end_to_end_timing_note='Sum of warm search median and one fetch/rerank; not a repeated end-to-end percentile.',
                  **retained_counts(ids, mols, corpus.total, molecule_total),
                  recall=recall, cap_scenarios=cap_results,
                  eligible_library_conformers=corpus.total-excluded_count)
    return result, dict(global_ids=ids, molecule_ids=mols, squared_l2=distances, **cap_arrays)


def summarize(rows, gate, gate_k):
    groups = {}
    for row in rows:
        groups.setdefault((row['requested_candidate_budget'], row['nprobe']), []).append(row)
    summaries = []
    for (budget, nprobe), group in sorted(groups.items()):
        recalls = [r['recall'][f'top{gate_k}_strict'] for r in group]
        times = [t for r in group for t in r['warm_search_seconds']]
        summaries.append(dict(candidate_budget=budget, nprobe=nprobe, query_count=len(group),
                              recall_min=min(recalls), recall_mean=float(np.mean(recalls)),
                              gate_passed=min(recalls) >= gate,
                              search_p50_seconds=float(np.median(times)), search_p95_seconds=float(np.percentile(times, 95)),
                              search_plus_fetch_rerank_median_seconds=float(np.median([
                                  r['warm_search_plus_single_fetch_rerank_seconds'] for r in group])),
                              retained_molecules_min=min(r['retained_molecules'] for r in group),
                              retained_molecules_max=max(r['retained_molecules'] for r in group),
                              molecule_reduction_min=min(r['molecule_reduction_fraction'] for r in group),
                              molecule_reduction_max=max(r['molecule_reduction_fraction'] for r in group)))
    passing = [r for r in summaries if r['gate_passed']]
    selected = min(passing, key=lambda r: (r['candidate_budget'], r['search_plus_fetch_rerank_median_seconds'])) if passing else None
    return summaries, selected


def render_report(report, output):
    a = report['acceptance']
    lines = ['# E033 library acceptance and 3D coarse-search calibration', '',
             f"Acceptance: {a['status']}; conformers: {a['library_conformers']:,}; source-grouped molecules: {a['library_molecules']:,}.",
             '', 'This measures standardized-USRCAT retrieval, not activity enrichment or final pose/docking quality.',
             'Reduction is caused by explicit candidate budgets; it is not a measured chemical rejection rate.', '',
             '| Candidate budget | nprobe | Worst query recall | Search p50 / p95 (s) | Retained molecules min–max | Molecule reduction min–max | Gate |',
             '|---:|---:|---:|---:|---:|---:|---|']
    for r in report['summary']:
        lines.append(f"| {r['candidate_budget']} | {r['nprobe']} | {r['recall_min']:.4f} | "
                     f"{r['search_p50_seconds']:.4f} / {r['search_p95_seconds']:.4f} | "
                     f"{r['retained_molecules_min']}–{r['retained_molecules_max']} | "
                     f"{r['molecule_reduction_min']:.2%}–{r['molecule_reduction_max']:.2%} | {r['gate_passed']} |")
    lines += ['', '## Provisional calibration selection', '',
              json.dumps(report['provisional_setting'], indent=2) if report['provisional_setting'] else
              'No tested configuration met the recall gate. Do not accept a low-recall setting solely because it filters more.',
              '', '## Interpretation and next step', '',
              'Candidate NPZ files preserve IDs and exact descriptor distances. Caps keep one descriptor-ranked representative per molecule.',
              'Run target-specific Gaussian pose refinement / E031 analysis on calibrated candidates next; measure its retention and cost separately.',
              'Held-out active/decoy labels are required for biological enrichment claims. None were inferred by this benchmark.',
              'Exact reference scanning and integrity checking are one-time evaluation costs, not online query latency.',
              'First-call timing is not guaranteed cold-cache timing. Warm p95 uses a small descriptive panel, not a service SLA.',
              'See report.json and metrics.csv for every query, recall denominator, cap scenario and timing.', '']
    (output / 'report.md').write_text('\n'.join(lines), encoding='utf-8')
    flat = []
    for r in report['results']:
        for c in r['cap_scenarios']:
            flat.append(dict(query=r['query'], nprobe=r['nprobe'], candidate_budget=r['requested_candidate_budget'],
                             candidate_conformers=r['retained_conformers'], candidate_molecules=r['retained_molecules'],
                             query_search_p50_seconds=r['warm_search_p50_seconds'],
                             fetch_rerank_seconds=r['candidate_fetch_and_exact_descriptor_rerank_seconds'],
                             **r['recall'], **{'cap_'+k: v for k, v in c.items()}))
    with (output / 'metrics.csv').open('w', newline='', encoding='utf-8') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(flat[0])); writer.writeheader(); writer.writerows(flat)


def run(args):
    import faiss
    batch, output = args.batch.resolve(), args.output.resolve()
    require(not output.is_relative_to(batch) and not batch.is_relative_to(output), 'Keep evaluation output separate from batch')
    require(not output.exists(), 'Output already exists: choose a new --output directory to preserve prior results')
    # No file creation in the library, including its existing advisory lock.
    import fcntl
    with (batch / 'run.lock').open('r') as lock:
        fcntl.flock(lock.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        output.mkdir(parents=True)
        try:
            _run_locked(args, batch, output, faiss)
        except Exception as exc:
            write(output / 'FAILED.json', dict(error=str(exc), type=type(exc).__name__))
            raise


def _run_locked(args, batch, output, faiss):
    started = time.perf_counter()
    catalog, acceptance = accept_library(batch, not args.metadata_only)
    write(output / 'acceptance.json', acceptance)
    corpus = Corpus(catalog)
    faiss.omp_set_num_threads(args.threads)
    with np.load(batch / 'faiss' / 'transform.npz', allow_pickle=False) as tr:
        mean, std = np.asarray(tr['mean'], dtype=np.float32), np.asarray(tr['std'], dtype=np.float32)
    require(mean.shape == std.shape == (60,) and np.all(np.isfinite(mean))
            and np.all(np.isfinite(std)) and np.all(std > 0), 'Invalid frozen USRCAT transform')
    t = time.perf_counter(); index = faiss.read_index(str(batch / 'faiss' / 'index.faiss'))
    load_seconds = time.perf_counter()-t
    require(index.ntotal == corpus.total and index.d == 60, 'FAISS count/dimension mismatch')
    query_rows, raw_queries, excluded = [], [], []
    if args.queries:
        with np.load(args.queries, allow_pickle=False) as q:
            raw_queries = np.asarray(q['vectors'], dtype=np.float32)
            names = q['names'].astype(str)
            excluded = list(q['exclude_molecule_ids'].astype('S16')) if 'exclude_molecule_ids' in q else [b''] * len(names)
        require(len(names) == len(raw_queries) == len(excluded), 'Query metadata length mismatch')
        query_rows = [dict(name=str(n), exclude_molecule_id=e.decode(), kind='external_raw_descriptor') for n, e in zip(names, excluded)]
    else:
        rng = np.random.default_rng(args.seed); used = set()
        for _ in range(max(1000, args.query_count*100)):
            gid = int(rng.integers(corpus.total)); mols, raw = corpus.fetch([gid], vectors=True)
            mid = bytes(mols[0])
            if mid in used:
                continue
            used.add(mid); raw_queries.append(raw[0]); excluded.append(mid)
            query_rows.append(dict(name=f'library_gid_{gid}', global_id=gid, exclude_molecule_id=mid.decode(),
                                   kind='library_calibration_all_same_molecule_conformers_excluded'))
            if len(query_rows) == args.query_count:
                break
        require(len(query_rows) == args.query_count, 'Unable to sample requested distinct molecules; reduce --query-count')
        raw_queries = np.asarray(raw_queries, dtype=np.float32)
    require(raw_queries.ndim == 2 and raw_queries.shape[1] == 60 and len(raw_queries) > 0
            and np.all(np.isfinite(raw_queries)), 'Queries must be finite raw shape (n,60)')
    queries = np.ascontiguousarray((raw_queries-mean)/std, dtype=np.float32)
    np.savez(output / 'queries.npz', vectors=raw_queries, names=np.asarray([r['name'] for r in query_rows]),
             exclude_molecule_ids=np.asarray(excluded, dtype='S16'))
    params = dict(nprobes=args.nprobes, budgets=args.budgets, truth_ks=args.truth_ks, caps=args.caps,
                  repeats=args.repeats, seed=args.seed, chunk_size=args.chunk_size, threads=args.threads,
                  recall_gate=args.recall_gate, metadata_only=args.metadata_only)
    write(output / 'protocol.json', dict(parameters=params, queries=query_rows,
          artifact_catalog_sha256=sha(batch / 'artifacts' / 'catalog.json'),
          faiss_manifest_sha256=sha(batch / 'faiss' / 'manifest.json'), transform_sha256=sha(batch / 'faiss' / 'transform.npz'),
          query_sha256=sha(output / 'queries.npz'), python=sys.version, numpy=np.__version__, faiss=faiss.__version__))
    truth, excluded_counts, exact_timing = exact_truth(corpus, queries, excluded, mean, std,
                                                       max(args.truth_ks), args.chunk_size)
    for qi, (d, ids) in enumerate(truth):
        require(len(ids) > 0, 'No eligible exact neighbors')
        np.savez(output / f'truth_{qi:03d}.npz', global_ids=ids, squared_l2=d)
        query_rows[qi]['excluded_conformers'] = excluded_counts[qi]
        query_rows[qi]['truth_denominators'] = {str(k): min(k, len(ids)) for k in args.truth_ks}
    results = []
    actual_nprobes = sorted({min(n, int(index.nlist)) for n in args.nprobes})
    for nprobe in actual_nprobes:
        for budget in args.budgets:
            for qi, q in enumerate(queries):
                log(f'Query {qi+1}/{len(queries)}: nprobe={nprobe}, candidate budget={budget}')
                result, arrays = evaluate_setting(index, corpus, q, excluded[qi], excluded_counts[qi], truth[qi], mean, std,
                    nprobe, budget, args.repeats, args.truth_ks, args.caps, acceptance['library_molecules'])
                result['query'] = query_rows[qi]['name']; result['query_index'] = qi
                filename = f'candidates_q{qi:03d}_np{nprobe}_k{budget}.npz'
                np.savez(output / filename, **arrays); result['candidate_file'] = filename
                results.append(result)
                write(output / 'partial-results.json', results)
    summary, selected = summarize(results, args.recall_gate, max(args.truth_ks))
    peak = None
    try:
        import resource
        peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * (1 if sys.platform == 'darwin' else 1024)
    except ImportError:
        pass
    report = dict(format='aidd-e033-library-acceptance', version=1, acceptance=acceptance, parameters=params,
                  queries=query_rows, exact_reference_timing=exact_timing, index_load_seconds=load_seconds,
                  results=results, summary=summary, provisional_setting=selected,
                  calibration_status='gate_passed_on_panel' if selected else 'recall_gate_not_met',
                  biological_quality='not_evaluated_no_active_decoy_labels',
                  pose_refinement_filtering='not_evaluated_USRCAT_coarse_stage_only',
                  query_panel_status='calibration_not_independent_holdout',
                  peak_process_rss_bytes=peak, wall_seconds=time.perf_counter()-started,
                  environment=dict(platform=platform.platform(), cpu_count=os.cpu_count(), faiss=faiss.__version__))
    write(output / 'report.json', report); render_report(report, output)
    write(output / 'EVALUATION_COMPLETE.json', dict(report_sha256=sha(output / 'report.json'),
          calibration_status=report['calibration_status'], acceptance_status=acceptance['status']))
    log(f'Evaluation finished: {output / "report.md"}; {report["calibration_status"]}')


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--batch', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--metadata-only', action='store_true')
    p.add_argument('--queries', type=Path, help='Optional NPZ vectors[n,60], names[n], exclude_molecule_ids[n]')
    p.add_argument('--query-count', type=int, default=8)
    p.add_argument('--seed', type=int, default=20260915)
    p.add_argument('--threads', type=int, default=20)
    p.add_argument('--chunk-size', type=int, default=32768)
    p.add_argument('--nprobes', nargs='+', type=int, default=[64, 128, 256])
    p.add_argument('--budgets', nargs='+', type=int, default=[1000, 10000, 100000])
    p.add_argument('--truth-ks', nargs='+', type=int, default=[100, 1000])
    p.add_argument('--caps', nargs='+', type=int, default=[100, 1000, 10000])
    p.add_argument('--repeats', type=int, default=3)
    p.add_argument('--recall-gate', type=float, default=.95)
    args = p.parse_args()
    require(all(x > 0 for x in [args.query_count, args.threads, args.chunk_size, args.repeats,
                               *args.nprobes, *args.budgets, *args.truth_ks, *args.caps]), 'All numeric counts must be positive')
    require(0 < args.recall_gate <= 1, 'Recall gate must be in (0,1]')
    if hasattr(os, 'sched_getaffinity'):
        require(args.threads <= len(os.sched_getaffinity(0)), 'threads exceeds available CPU affinity')
    run(args)


if __name__ == '__main__':
    main()
