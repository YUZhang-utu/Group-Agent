"""Uncapped per-conformer pose/feature evaluation with bounded, resumable chunks."""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
import copy
import json
import os
from pathlib import Path
import sqlite3
import time

import numpy as np

from . import library_acceptance as ev
from .chunk_execution import bounded_results
from .classified_features import write_tables
from .expanded_wee1 import fingerprint, ensure_file_descriptor_limit
from .gaussian_batch import ArtifactCatalogReader, _load_query, _score_ids, _atomic_json, _atomic_savez
from .interaction_fast import score_batched
from .interaction_review import OBJECTIVE
from .screening_selection import archive, check_hashes

THRESHOLDS = (.25, .5, .75)
POSE_PARAMETERS = dict(sigma=1., cutoff=4.5, pair_tolerance=2., axial_samples=6,
                       max_pair_seeds=512, bounded_pair_seeds=True)
_STATE = None


def catalog_count(catalog):
    end = 0
    for row in sorted(catalog['shards'], key=lambda r: r['global_id_start']):
        ev.require(row['global_id_start'] == end and row['conformers'] > 0,
                   'Catalog ranges have gaps or overlap')
        end += row['conformers']
    ev.require(end == catalog['conformers'] and end > 0, 'Invalid catalog total')
    return end


def initialize(query):
    global _STATE
    catalog = Path(query['artifact_catalog'])
    if os.name == 'posix': ensure_file_descriptor_limit(len(ev.read(catalog)['shards']))
    _, original = _load_query(Path(query['query_npz']))
    indices = sorted({a['feature_index'] for a in query['anchors']})
    expanded = dict(original, anchor_feature_indices=np.array(indices, dtype=np.int64))
    _STATE = (ArtifactCatalogReader(catalog), original, expanded, query)


def compute_chunk(task):
    start, stop, target = task
    reader, original, expanded, query = _STATE
    ids = np.arange(start, stop, dtype=np.int64)
    started = time.perf_counter()
    # No descriptor Top-K and no Gaussian Top-N: every supplied ID is evaluated.
    rigid = _score_ids(reader, original, ids, **POSE_PARAMETERS)
    _, assignments, scores = score_batched(Path(query['artifact_catalog']),
        Path(query['chemical_companion']), expanded, ids, rigid['molecule_ids'],
        rigid['conformer_ids'], [OBJECTIVE], {OBJECTIVE: rigid[OBJECTIVE+'__transform']},
        sigma=1., cutoff=4.5, angular_power=2.)
    arrays = {k: rigid[k] for k in ('global_ids', 'molecule_ids', 'conformer_ids',
                                    OBJECTIVE+'__objective', OBJECTIVE+'__transform')}
    arrays.update({OBJECTIVE+'__anchor_scores': scores[OBJECTIVE],
                   OBJECTIVE+'__anchor_assignments': assignments[OBJECTIVE],
                   'query_anchor_feature_indices': expanded['anchor_feature_indices']})
    _atomic_savez(Path(target), arrays)
    receipt = dict(start=start, stop=stop, sha256=ev.sha(target),
                   wall_seconds=time.perf_counter()-started)
    _atomic_json(Path(target).with_suffix('.receipt.json'), receipt)
    return receipt


def load_chunk(path, start, stop, indices):
    receipt = ev.read(path.with_suffix('.receipt.json'))
    ev.require(receipt['start'] == start and receipt['stop'] == stop and
               ev.sha(path) == receipt['sha256'], 'Changed full-library chunk receipt')
    data = archive(path)
    ev.require(np.array_equal(data['global_ids'], np.arange(start, stop)) and
               np.array_equal(data['query_anchor_feature_indices'], indices), 'Incomplete chunk IDs/features')
    scores, assignments = data[OBJECTIVE+'__anchor_scores'], data[OBJECTIVE+'__anchor_assignments']
    ev.require(scores.shape == assignments.shape == (stop-start, len(indices)) and
               np.isfinite(scores).all(), 'Invalid feature scores')
    for key in ('molecule_ids', 'conformer_ids', OBJECTIVE+'__objective', OBJECTIVE+'__transform'):
        ev.require(len(data[key]) == stop-start, 'Incomplete chunk arrays')
    ev.require(np.isfinite(data[OBJECTIVE+'__objective']).all() and
               np.isfinite(data[OBJECTIVE+'__transform']).all(), 'Invalid pose arrays')
    return data


class Counts:
    """Disk-backed molecule union; never use molecule masks for same-pose ALL."""
    def __init__(self, path, width):
        ev.require(0 < width <= 63, 'Diagnostic mask supports 1..63 distinct query features')
        self.db = sqlite3.connect(path)
        # This derived index is rebuilt from verified chunks on every invocation.
        self.db.execute('DROP TABLE IF EXISTS molecules')
        self.db.execute('CREATE TABLE molecules (id TEXT PRIMARY KEY, a INTEGER, b INTEGER, c INTEGER)')
        self.molecules = np.zeros((3, width), dtype=np.int64)
        self.conformers = np.zeros((3, width), dtype=np.int64)
        self.processed = 0

    def add(self, data):
        scores = data[OBJECTIVE+'__anchor_scores']
        assigned = data[OBJECTIVE+'__anchor_assignments'] >= 0
        flags = np.asarray([assigned & (scores >= t) for t in THRESHOLDS])
        self.conformers += flags.sum(axis=1)
        weights = np.left_shift(np.int64(1), np.arange(scores.shape[1], dtype=np.int64))
        masks = flags.astype(np.int64) @ weights
        grouped = {}
        for i, mid in enumerate(data['molecule_ids']):
            row = grouped.setdefault(str(mid), [0, 0, 0])
            for j in range(3): row[j] |= int(masks[j, i])
        with self.db:
            for mid, row in grouped.items():
                old = self.db.execute('SELECT a,b,c FROM molecules WHERE id=?', (mid,)).fetchone() or (0,0,0)
                merged = [a | b for a, b in zip(old, row)]
                for j, (a,b) in enumerate(zip(old, merged)):
                    new = a ^ b
                    while new:
                        bit = new & -new
                        self.molecules[j, bit.bit_length()-1] += 1
                        new ^= bit
                self.db.execute('INSERT OR REPLACE INTO molecules VALUES (?,?,?,?)', (mid, *merged))
        self.processed += len(scores)

    def total_molecules(self):
        return self.db.execute('SELECT COUNT(*) FROM molecules').fetchone()[0]


def preview(source, output, policy):
    """Stream same-pose conditions, then keep one representative per molecule."""
    from .screening_selection import validate_selection, select_rows
    validate_selection(policy)
    report = ev.read(source)
    ev.require(report.get('status') == 'complete' and report.get('kind') == 'full_library_conditions',
               'Selection requires completed full-library condition evidence')
    check_hashes(report['sources'])
    available = {a['anchor_id']: (q, a) for q in report['queries'] for a in q['anchors']}
    ev.require(set(policy['required_anchors']) <= available.keys(), 'Unknown anchor ID')
    ev.require(len({available[a][0]['query_id'] for a in policy['required_anchors']}) == 1,
               'Select anchors from one query')
    q = available[policy['required_anchors'][0]][0]
    columns = sorted({available[a][1]['score_column'] for a in policy['required_anchors']})
    output = Path(output); output.mkdir(parents=True, exist_ok=True)
    db = sqlite3.connect(output/'selection.sqlite')
    conformers = 0
    try:
        db.execute('DROP TABLE IF EXISTS selected')
        db.execute('CREATE TABLE selected (mid TEXT PRIMARY KEY, score REAL, gid INTEGER, payload TEXT)')
        end = 0
        for chunk in q['chunks']:
            ev.require(chunk['start'] == end, 'Noncontiguous full-library evidence')
            data = load_chunk(Path(chunk['path']), chunk['start'], chunk['stop'],
                              sorted({a['feature_index'] for a in q['anchors']}))
            rows, matched = select_rows(data, data, columns, policy, q['query_id'])
            conformers += matched
            with db:
                for row in rows:
                    row['rigid_chunk'] = chunk['path']
                    db.execute('INSERT INTO selected VALUES (?,?,?,?) ON CONFLICT(mid) DO UPDATE SET '
                        'score=excluded.score,gid=excluded.gid,payload=excluded.payload '
                        'WHERE excluded.score>selected.score OR (excluded.score=selected.score AND excluded.gid<selected.gid)',
                        (row['molecule_id'], row['gaussian_score'], row['global_id'], json.dumps(row, sort_keys=True)))
            end = chunk['stop']
        ev.require(end == q['counts']['total_conformers'] and q['full_coverage'], 'Incomplete coverage')
        matched = db.execute('SELECT COUNT(*) FROM selected').fetchone()[0]
        limit = policy.get('max_molecules', matched)
        representatives = output/'representatives.jsonl'
        with representatives.open('w', encoding='utf-8', newline='\n') as handle:
            for (payload,) in db.execute('SELECT payload FROM (SELECT payload,gid FROM selected ORDER BY score DESC,gid ASC LIMIT ?) ORDER BY gid', (limit,)):
                handle.write(payload+'\n')
    finally:
        db.close()
    result = dict(status='complete', kind='selection_preview', query=q, policy=policy,
        counts=dict(**q['counts'], matching_conformers=conformers, matching_molecules=matched,
                    selected_molecules=min(matched,limit), selected_representatives=min(matched,limit)),
        representatives_jsonl=str(representatives), evidence_report=str(Path(source).resolve()),
        sources={**report['sources'], **fingerprint([source, representatives])},
        e031_changes_ranking=False, biological_quality='not_evaluated',
        approval='Preview only; separately request export of this selection')
    _atomic_json(output/'report.json', result)
    return result


def export(source, output):
    import tempfile
    from .chemical_companion import ChemicalCompanionReader
    from .gaussian_overlay import apply_transform
    from .predocking_qc import _chemical_sdf
    selected = ev.read(source)
    check_hashes(selected['sources'])
    ev.require(str(Path(selected['evidence_report']).resolve()) in selected['sources'], 'Unsealed evidence')
    with tempfile.TemporaryDirectory() as tmp:
        expected = preview(selected['evidence_report'], tmp, selected['policy'])
        ev.require(expected['counts'] == selected['counts'] and
            ev.sha(expected['representatives_jsonl']) == ev.sha(selected['representatives_jsonl']), 'Selection changed')
    q = selected['query']; output = Path(output); output.mkdir(parents=True, exist_ok=True)
    reader = chemistry = None
    if selected['counts']['selected_molecules']:
        if os.name == 'posix': ensure_file_descriptor_limit(len(ev.read(q['artifact_catalog'])['shards']))
        reader, chemistry = ArtifactCatalogReader(Path(q['artifact_catalog'])), ChemicalCompanionReader(Path(q['chemical_companion']))
    last, arrays = None, None
    with open(selected['representatives_jsonl'], encoding='utf-8') as rows, \
         (output/'selected-poses.sdf').open('w', encoding='utf-8') as sdf, \
         (output/'selected-ids.json').open('w', encoding='utf-8') as ids:
        ids.write('['); first = True
        for line in rows:
            row = json.loads(line)
            ev.require(row['rigid_chunk'] in selected['sources'], 'Unsealed rigid chunk')
            if last != row['rigid_chunk']:
                last = row['rigid_chunk']; arrays = archive(last)
            i, gid = row['row_index'], row['global_id']
            candidate, chem = reader.get(gid), chemistry.get(gid)
            ev.require(int(arrays['global_ids'][i]) == gid and all(x.molecule_id == row['molecule_id'] and
                x.conformer_id == row['conformer_id'] for x in (candidate,chem)), 'Export identity mismatch')
            points = apply_transform(candidate.shape_points, arrays[OBJECTIVE+'__transform'][i])
            sdf.write(_chemical_sdf(points,chem,q['query_id']).replace('$$$$',
                '>  <AIDD_MOLECULE_ID>\n'+row['molecule_id']+'\n\n$$$$'))
            if not first: ids.write(',\n')
            ids.write(json.dumps(row)); first = False
        ids.write(']\n')
    check_hashes(selected['sources'])
    result = dict(status='complete', kind='docking_handoff', counts=selected['counts'], policy=selected['policy'], query=q,
        sdf=str(output/'selected-poses.sdf'), ids=str(output/'selected-ids.json'),
        sources={**selected['sources'], **fingerprint([source])},
        outputs=fingerprint([output/'selected-poses.sdf',output/'selected-ids.json']), e031_changes_ranking=False,
        docking_status='not_run_requires_engine_specific_preparation',
        preparation='Heavy-atom rigid poses only; engine-specific preparation required')
    _atomic_json(output/'report.json', result)
    return result


def run(source, output, *, workers=8, chunk_size=2048, max_chunks=None):
    from .prompt_workflow import file_lock
    source, output = Path(source).resolve(), Path(output).resolve()
    ev.require(workers > 0 and chunk_size > 0 and (max_chunks is None or max_chunks > 0), 'Invalid scheduling')
    evidence = ev.read(source)
    ev.require(evidence.get('status') == 'complete' and evidence.get('classification') and evidence.get('kind') == 'screening_evidence',
               'Use a completed classified evidence report')
    query_ids = [q['query_id'] for q in evidence.get('queries', [])]
    ev.require(query_ids and len(set(query_ids)) == len(query_ids), 'Expected distinct nonempty query definitions')
    for q in evidence['queries']:
        for protected in (source.parent, Path(q['artifact_catalog']).resolve().parent,
                          Path(q['chemical_companion']).resolve().parent):
            ev.require(not output.is_relative_to(protected) and not protected.is_relative_to(output),
                       'Output overlaps protected inputs')
    output.mkdir(parents=True, exist_ok=True)
    with file_lock(output/'run.lock'):
        check_hashes(evidence['sources'])
        sources = {**evidence['sources'], **fingerprint([source])}
        seen_catalogs = set()
        for q in evidence['queries']:
            for name in ('artifact_catalog', 'chemical_companion'):
                catalog_path = Path(q[name]).resolve()
                if catalog_path in seen_catalogs: continue
                seen_catalogs.add(catalog_path)
                ev.log('Full-library integrity preflight: '+str(catalog_path))
                for row in ev.read(catalog_path)['shards']:
                    directory = Path(row['path'])
                    if not directory.is_dir(): directory = catalog_path.parent/row['name']
                    manifest = directory/'manifest.json'
                    ev.require(ev.sha(manifest) == row['manifest_sha256'], 'Shard manifest changed')
                    manifest_data = ev.verify_manifest(directory, full=True)
                    for filename, record in manifest_data.get('files', {}).items():
                        sources[str((directory/filename).resolve())] = record['sha256']
                    sources.update(fingerprint([manifest]))
        protocol = dict(source_hashes=sources, code=fingerprint(sorted(Path(__file__).parent.glob('*.py'))),
                        chunk_size=chunk_size, pose_parameters=POSE_PARAMETERS, thresholds=THRESHOLDS,
                        numpy=np.__version__, candidate_top_k=None, gaussian_top_n=None)
        # JSON normalizes tuple values, including thresholds.
        protocol = json.loads(json.dumps(protocol))
        path = output/'protocol.json'
        if path.exists(): ev.require(ev.read(path) == protocol, 'Changed protocol; use a fresh output')
        _atomic_json(path, protocol)
        started = time.perf_counter()
        result = dict(kind='full_library_conditions', status='running', queries=[], sources=sources,
            library=evidence['library'], classification=evidence['classification'], classes=evidence['classes'],
            unsupported_classes=evidence['unsupported_classes'], e031_changes_ranking=False,
            policy='All catalog conformers; one heuristic anchored-Gaussian pose per conformer, not all possible poses. No Top-K/Top-N membership truncation.',
            biological_quality='not_evaluated', candidate_top_k=None, gaussian_top_n=None)
        try:
            for qi, original in enumerate(evidence['queries']):
                q = copy.deepcopy(original)
                total = catalog_count(ev.read(q['artifact_catalog']))
                ev.require(total == evidence['library']['library_conformers'], 'Accepted library size differs')
                indices = sorted({a['feature_index'] for a in q['anchors']})
                for a in q['anchors']: a['score_column'] = indices.index(a['feature_index'])
                root = output/f'query-{qi}'; root.mkdir(exist_ok=True)
                count = Counts(root/'molecule-counts.sqlite', len(indices))
                chunk_total = (total+chunk_size-1)//chunk_size
                selected = min(chunk_total, max_chunks) if max_chunks else chunk_total
                def tasks():
                    for ci in range(selected):
                        start, stop = ci*chunk_size, min(total, (ci+1)*chunk_size)
                        target = root/f'chunk-{ci:08d}.npz'
                        if target.with_suffix('.receipt.json').exists():
                            load_chunk(target, start, stop, indices)
                        else: yield start, stop, str(target)
                q['chunks'] = []
                q['counts'] = dict(evaluated_conformers=0, evaluated_molecules=0, total_conformers=total)
                q['full_coverage'] = False
                result['queries'].append(q)
                _atomic_json(output/'report.json', result)
                try:
                    # Persistent workers reuse the mapped library and query across chunks.
                    with ProcessPoolExecutor(max_workers=workers, initializer=initialize, initargs=(q,)) as pool:
                        for _, receipt in bounded_results(pool, compute_chunk, tasks(), workers*2):
                            ev.log(f"{q['query_id']}: evaluated IDs {receipt['start']}..{receipt['stop']-1} / {total}")
                            _atomic_json(output/'progress.json', dict(query_id=q['query_id'], latest_completed_chunk=receipt,
                                total_conformers=total, status='running', note='Chunks can finish out of order; this is not cumulative coverage'))
                    for ci in range(selected):
                        start, stop = ci*chunk_size, min(total, (ci+1)*chunk_size)
                        target = root/f'chunk-{ci:08d}.npz'
                        count.add(load_chunk(target, start, stop, indices))
                        q['chunks'].append(dict(path=str(target), start=start, stop=stop))
                        sources.update(fingerprint([target, target.with_suffix('.receipt.json')]))
                    q['counts'] = dict(evaluated_conformers=count.processed, evaluated_molecules=count.total_molecules(),
                                       total_conformers=total)
                    q['full_coverage'] = count.processed == total
                    if q['full_coverage']:
                        ev.require(count.total_molecules() == evidence['library']['library_molecules'],
                                   'Full-library molecule count differs from acceptance')
                    for a in q['anchors']:
                        column = a['score_column']
                        a['diagnostic_counts'] = [dict(minimum_score=t, conformers=int(count.conformers[j,column]),
                            molecules=int(count.molecules[j,column])) for j,t in enumerate(THRESHOLDS)]
                finally:
                    count.db.close()
                # Old budgeted candidate artifacts are evidence only, never the new result set.
                for key in ('rigid','sidecar','retrieval'): q.pop(key, None)
                _atomic_json(output/'report.json', result)
            result['status'] = 'complete' if all(q['full_coverage'] for q in result['queries']) else 'partial'
            result['timing_scope'] = 'Includes scoring, chunk persistence, count aggregation and final integrity checks; excludes integrity preflight; resume is not fresh timing'
            check_hashes(sources)
            result['wall_seconds_this_invocation'] = time.perf_counter()-started
            _atomic_json(output/'report.json', result)
            write_tables(result, output)
            _atomic_json(output/'RUN_STATUS.json', dict(status=result['status'], report_sha256=ev.sha(output/'report.json')))
        except Exception as exc:
            result.update(status='failed', error=str(exc))
            _atomic_json(output/'report.json', result)
            _atomic_json(output/'RUN_STATUS.json', dict(status='failed'))
            raise
        return result


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--classification', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--workers', type=int, default=8)
    p.add_argument('--chunk-size', type=int, default=2048)
    p.add_argument('--max-chunks', type=int, help='Explicit pilot limit per query; incomplete coverage is labeled partial')
    a = p.parse_args()
    result = run(a.classification, a.output, workers=a.workers, chunk_size=a.chunk_size, max_chunks=a.max_chunks)
    print(json.dumps(dict(status=result['status'], report=str(a.output/'report.json'))))


if __name__ == '__main__': main()
