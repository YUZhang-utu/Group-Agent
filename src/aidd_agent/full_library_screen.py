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
_FILTER_READER = None
_BOUND = None
_JOINT = None


def catalog_count(catalog):
    end = 0
    for row in sorted(catalog['shards'], key=lambda r: r['global_id_start']):
        ev.require(row['global_id_start'] == end and row['conformers'] > 0,
                   'Catalog ranges have gaps or overlap')
        end += row['conformers']
    ev.require(end == catalog['conformers'] and end > 0, 'Invalid catalog total')
    return end


def initialize(query):
    global _STATE, _FILTER_READER, _BOUND, _JOINT
    catalog = Path(query['artifact_catalog'])
    if os.name == 'posix': ensure_file_descriptor_limit(len(ev.read(catalog)['shards']))
    _, original = _load_query(Path(query['query_npz']))
    indices = sorted({a['feature_index'] for a in query['anchors']})
    pocket = _load_query(Path(query['consensus_npz']))[1] if query.get('consensus_npz') else original
    expanded = dict(pocket, anchor_feature_indices=np.array(indices, dtype=np.int64))
    _STATE = (ArtifactCatalogReader(catalog), original, expanded, query)
    _FILTER_READER = _BOUND = _JOINT = None
    if query.get('condition_policy'):
        from .interaction_fast import InteractionFeatureReader
        from .necessary_conditions import NecessaryConditions
        policy = query['condition_policy']
        required = [a['feature_index'] for a in query['anchors'] if a['anchor_id'] in policy['required_anchors']]
        _FILTER_READER = InteractionFeatureReader(catalog, Path(query['chemical_companion']))
        _BOUND = NecessaryConditions(expanded,required,policy['match_mode'],policy['minimum_score'])
        if 'coarse_constraints' in policy:
            from .joint_coarse import JointCoarse
            _JOINT = JointCoarse(original, policy['coarse_constraints'])


def coarse_decisions(ids, records):
    """Joint eligibility precedes anchor bounds and any seed generation."""
    joint = np.ones(len(ids), dtype=bool)
    decisions = []
    for i, (gid, record) in enumerate(zip(ids, records)):
        if _JOINT is not None:
            candidate = _STATE[0].get(int(gid))
            ev.require(candidate.molecule_id == record.molecule_id and
                       candidate.conformer_id == record.conformer_id, 'Joint coarse identity mismatch')
            passed, reason = _JOINT.check(candidate)
            joint[i] = passed
            if not passed:
                decisions.append((False, reason))
                continue
        decisions.append(_BOUND.check(record))
    return decisions, joint


def compute_chunk(task):
    start, stop, target = task
    reader, original, expanded, query = _STATE
    ids = np.arange(start, stop, dtype=np.int64)
    started = time.perf_counter()
    if _BOUND is not None:
        return compute_funnel_chunk(task, ids, started)
    # No descriptor Top-K and no Gaussian Top-N: every supplied ID is evaluated.
    rigid = _score_ids(reader, original, ids, backend=query.get('pose_backend','reference'), **POSE_PARAMETERS)
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


def compute_funnel_chunk(task, ids, started):
    from collections import Counter
    from .necessary_conditions import pose_mask
    start, stop, target = task
    reader, original, expanded, query = _STATE
    policy = query['condition_policy']
    columns = sorted({a['score_column'] for a in query['anchors'] if a['anchor_id'] in policy['required_anchors']})
    records = [_FILTER_READER.get(int(gid)) for gid in ids]
    decisions, joint = coarse_decisions(ids, records)
    keep = np.array([x[0] for x in decisions],dtype=bool)
    reasons = Counter(x[1] for x in decisions if not x[0])
    prefilter_seconds = time.perf_counter()-started
    n, width = len(ids), len(expanded['anchor_feature_indices'])
    arrays = dict(global_ids=ids, molecule_ids=np.array([c.molecule_id for c in records]),
        conformer_ids=np.array([c.conformer_id for c in records]),
        query_anchor_feature_indices=expanded['anchor_feature_indices'], pose_evaluated=keep,
        prefilter_passed=keep, **{OBJECTIVE+'__objective':np.zeros(n),
        OBJECTIVE+'__transform':np.tile(np.eye(4).reshape(1,16),(n,1)),
        OBJECTIVE+'__anchor_scores':np.zeros((n,width)),
        OBJECTIVE+'__anchor_assignments':np.full((n,width),-2,dtype=np.int32)})
    if 'coarse_constraints' in policy:
        arrays['joint_eligible'] = joint
    pose_seconds = annotation_seconds = 0.
    reused = False
    if query.get('reuse_directory'):
        ev.require('coarse_constraints' not in policy, 'Joint rules require fresh chunks; legacy chunks lack joint eligibility')
        legacy = Path(query['reuse_directory'])/f"chunk-{start//query['reuse_chunk_size']:08d}.npz"
        if legacy.with_suffix('.receipt.json').exists():
            prior = load_chunk(legacy,start,stop,expanded['anchor_feature_indices'])
            ev.require('condition_passed' not in prior, 'Only unfiltered E044 chunks can be reused')
            ev.require(np.array_equal(prior['molecule_ids'],arrays['molecule_ids']) and
                       np.array_equal(prior['conformer_ids'],arrays['conformer_ids']), 'Reused chunk identity mismatch')
            prior_matches = pose_mask(prior,columns,policy,OBJECTIVE)
            ev.require(not np.any(prior_matches & ~keep), 'Necessary-condition bound rejected an existing passing pose')
            for key in (OBJECTIVE+'__objective',OBJECTIVE+'__anchor_scores',OBJECTIVE+'__anchor_assignments'):
                arrays[key] = prior[key]
            arrays[OBJECTIVE+'__transform'] = prior[OBJECTIVE+'__transform'].reshape(-1,16)
            arrays['pose_evaluated'] = np.ones(n,dtype=bool)
            reused = True
    prepared = None
    feasibility = {}
    invariant_keep = keep.copy()
    if query.get('pose_feasibility') and not reused:
        from .pose_feasibility import prepare_survivors
        keep,prepared,feasibility = prepare_survivors(reader,original,_BOUND,records,ids,keep,POSE_PARAMETERS)
        arrays['pose_evaluated'] = keep.copy()
    arrays['pose_feasibility_passed'] = keep.copy()
    if keep.any() and not reused:
        stage_started = time.perf_counter()
        rigid = _score_ids(reader,original,ids[keep],backend=query.get('pose_backend','reference'),prepared=prepared,**POSE_PARAMETERS)
        ev.require(np.array_equal(rigid['molecule_ids'],arrays['molecule_ids'][keep]) and
                   np.array_equal(rigid['conformer_ids'],arrays['conformer_ids'][keep]), 'Prefilter/pose identity mismatch')
        pose_seconds = time.perf_counter()-stage_started
        stage_started = time.perf_counter()
        _, assignments, scores = score_batched(Path(query['artifact_catalog']),Path(query['chemical_companion']),
            expanded,ids[keep],rigid['molecule_ids'],rigid['conformer_ids'],[OBJECTIVE],
            {OBJECTIVE:rigid[OBJECTIVE+'__transform']},sigma=1.,cutoff=4.5,angular_power=2.)
        annotation_seconds = time.perf_counter()-stage_started
        arrays[OBJECTIVE+'__objective'][keep] = rigid[OBJECTIVE+'__objective']
        arrays[OBJECTIVE+'__transform'][keep] = rigid[OBJECTIVE+'__transform'].reshape(-1,16)
        arrays[OBJECTIVE+'__anchor_scores'][keep] = scores[OBJECTIVE]
        arrays[OBJECTIVE+'__anchor_assignments'][keep] = assignments[OBJECTIVE]
    arrays['condition_passed'] = pose_mask(arrays,columns,policy,OBJECTIVE)
    _atomic_savez(Path(target),arrays)
    receipt = dict(start=start,stop=stop,sha256=ev.sha(target),wall_seconds=time.perf_counter()-started,
        prefilter_seconds=prefilter_seconds,pose_seconds=pose_seconds,annotation_seconds=annotation_seconds,pose_backend=query.get('pose_backend','reference'),prefilter_passed=int(invariant_keep.sum()),
        prefilter_rejected=int((~invariant_keep).sum()),rejection_reasons=dict(reasons),
        matching_conformers=int(arrays['condition_passed'].sum()))
    receipt.update(feasibility)
    receipt['pose_candidates'] = int(keep.sum())
    receipt['legacy_chunk_reused'] = reused
    receipt['new_pose_evaluations'] = 0 if reused else int(keep.sum())
    if reused:
        receipt['legacy_source'] = fingerprint([legacy,legacy.with_suffix('.receipt.json')])
    _atomic_json(Path(target).with_suffix('.receipt.json'),receipt)
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
    if 'condition_passed' in data:
        if 'joint_eligible' in data:
            mask = data['joint_eligible']
            ev.require(mask.shape == (stop-start,) and mask.dtype == bool, 'Invalid joint eligibility mask')
            ev.require(not np.any(data['condition_passed'] & ~mask), 'A match cannot fail joint eligibility')
        for key in ('condition_passed','pose_evaluated','prefilter_passed'):
            ev.require(data[key].shape == (stop-start,) and data[key].dtype == bool, 'Invalid funnel masks')
        ev.require(not np.any(data['condition_passed'] & ~data['pose_evaluated']), 'A match needs an evaluated pose')
        if 'pose_feasibility_passed' in data:
            mask=data['pose_feasibility_passed']
            ev.require(mask.shape==(stop-start,) and mask.dtype==bool,'Invalid pose feasibility mask')
            ev.require(not np.any(data['condition_passed'] & ~mask),'A match cannot fail pose feasibility')
    return data


class Counts:
    """Disk-backed molecule union; never use molecule masks for same-pose ALL."""
    def __init__(self, path, width):
        ev.require(0 < width <= 63, 'Diagnostic mask supports 1..63 distinct query features')
        self.db = sqlite3.connect(path)
        # This derived index is rebuilt from verified chunks on every invocation.
        self.db.execute('DROP TABLE IF EXISTS molecules')
        self.db.execute('CREATE TABLE molecules (id TEXT PRIMARY KEY, a INTEGER, b INTEGER, c INTEGER, matched INTEGER)')
        self.molecules = np.zeros((3, width), dtype=np.int64)
        self.conformers = np.zeros((3, width), dtype=np.int64)
        self.processed = 0
        self.prefilter_passed = self.pose_evaluated = self.matching_conformers = self.pose_candidates = 0

    def add(self, data):
        scores = data[OBJECTIVE+'__anchor_scores']
        assigned = data[OBJECTIVE+'__anchor_assignments'] >= 0
        flags = np.asarray([assigned & (scores >= t) for t in THRESHOLDS])
        condition = data.get('condition_passed',np.zeros(len(scores),dtype=bool))
        if 'condition_passed' in data: flags &= condition[None,:,None]
        self.prefilter_passed += int(data.get('prefilter_passed',np.ones(len(scores),dtype=bool)).sum())
        self.pose_candidates += int(data.get('pose_feasibility_passed',data.get('prefilter_passed',np.ones(len(scores),dtype=bool))).sum())
        self.pose_evaluated += int(data.get('pose_evaluated',np.ones(len(scores),dtype=bool)).sum())
        self.matching_conformers += int(condition.sum())
        self.conformers += flags.sum(axis=1)
        weights = np.left_shift(np.int64(1), np.arange(scores.shape[1], dtype=np.int64))
        masks = flags.astype(np.int64) @ weights
        grouped = {}
        for i, mid in enumerate(data['molecule_ids']):
            row = grouped.setdefault(str(mid), [0, 0, 0, 0])
            for j in range(3): row[j] |= int(masks[j, i])
            row[3] |= int(condition[i])
        with self.db:
            for mid, row in grouped.items():
                old = self.db.execute('SELECT a,b,c,matched FROM molecules WHERE id=?', (mid,)).fetchone() or (0,0,0,0)
                merged = [a | b for a, b in zip(old, row)]
                for j, (a,b) in enumerate(zip(old[:3], merged[:3])):
                    new = a ^ b
                    while new:
                        bit = new & -new
                        self.molecules[j, bit.bit_length()-1] += 1
                        new ^= bit
                self.db.execute('INSERT OR REPLACE INTO molecules VALUES (?,?,?,?,?)', (mid, *merged))
        self.processed += len(scores)

    def total_molecules(self):
        return self.db.execute('SELECT COUNT(*) FROM molecules').fetchone()[0]


def preview(source, output, policy):
    """Stream same-pose conditions, then keep one representative per molecule."""
    from .screening_selection import validate_selection, select_rows
    validate_selection(policy)
    report = ev.read(source)
    ev.require(report.get('status') == 'complete' and report.get('kind') in ('full_library_conditions','condition_funnel'),
               'Selection requires completed full-library condition evidence')
    if report.get('kind') == 'condition_funnel':
        ev.require(normalize_policy(policy) == report['condition_policy'],
                   'Changed rule requires a new funnel run; rejected poses were not scored for other conditions')
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
        counts=dict(q['counts'], matching_conformers=conformers, matching_molecules=matched,
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


def normalize_policy(policy):
    from .screening_selection import validate_selection
    validate_selection(policy)
    result = dict(required_anchors=sorted(policy['required_anchors']), match_mode=policy['match_mode'],
                  minimum_score=policy['minimum_score'])
    if 'coarse_constraints' in policy:
        result['coarse_constraints'] = copy.deepcopy(policy['coarse_constraints'])
    return result


def run_funnel(selection, output, **kwargs):
    selected = ev.read(selection)
    ev.require(selected.get('kind') == 'selection_preview', 'Use an explicit condition selection preview')
    check_hashes(selected['sources'])
    source = Path(selected['evidence_report']).resolve()
    ev.require(str(source) in selected['sources'], 'Unsealed condition evidence')
    return run(source,output,condition=normalize_policy(selected['policy']),selection_source=selection,**kwargs)


def run(source, output, *, workers=None, chunk_size=256, max_chunks=None, condition=None, selection_source=None, reuse=None, backend=None, pose_feasibility=None):
    from .prompt_workflow import file_lock
    from .pose_acceleration import hardware_options, array_module
    backend, workers = hardware_options(backend,workers)
    if pose_feasibility is None:
        setting=os.environ.get('AIDD_POSE_FEASIBILITY','1')
        ev.require(setting in ('0','1'),'AIDD_POSE_FEASIBILITY must be 0 or 1')
        pose_feasibility=setting=='1'
    if backend == 'cupy': array_module(backend)
    ev.log(f'Pose backend: {backend}; workers: {workers}; chunk size: {chunk_size}')
    source, output = Path(source).resolve(), Path(output).resolve()
    ev.require(workers > 0 and chunk_size > 0 and (max_chunks is None or max_chunks > 0), 'Invalid scheduling')
    evidence = ev.read(source)
    ev.require(evidence.get('status') == 'complete' and evidence.get('classification') and evidence.get('kind') == 'screening_evidence',
               'Use a completed classified evidence report')
    query_ids = [q['query_id'] for q in evidence.get('queries', [])]
    ev.require(query_ids and len(set(query_ids)) == len(query_ids), 'Expected distinct nonempty query definitions')
    if condition:
        condition = normalize_policy(condition)
        available = {a['anchor_id']:q['query_id'] for q in evidence['queries'] for a in q['anchors']}
        ev.require(set(condition['required_anchors']) <= available.keys(), 'Unknown condition feature')
        selected_queries = {available[a] for a in condition['required_anchors']}
        ev.require(len(selected_queries) == 1, 'A funnel rule must use one crystal query')
        evidence['queries'] = [q for q in evidence['queries'] if q['query_id'] in selected_queries]
        for q in evidence['queries']: q['condition_policy'] = condition
    for q in evidence['queries']:
        q['pose_backend'] = backend
        q['pose_feasibility'] = bool(condition and pose_feasibility)
    for q in evidence['queries']:
        for protected in (source.parent, Path(q['artifact_catalog']).resolve().parent,
                          Path(q['chemical_companion']).resolve().parent):
            ev.require(not output.is_relative_to(protected) and not protected.is_relative_to(output),
                       'Output overlaps protected inputs')
    output.mkdir(parents=True, exist_ok=True)
    with file_lock(output/'run.lock'):
        check_hashes(evidence['sources'])
        sources = {**evidence['sources'], **fingerprint([source])}
        if selection_source: sources.update(fingerprint([selection_source]))
        if reuse:
            ev.require(condition is not None, 'Legacy reuse is only supported for a condition funnel')
            legacy = ev.read(Path(reuse)/'protocol.json')
            ev.require(legacy.get('condition_policy') is None and legacy['chunk_size'] == chunk_size and
                       legacy['pose_parameters'] == POSE_PARAMETERS and legacy['numpy'] == np.__version__,
                       'Legacy pose protocol differs')
            original_queries = ev.read(source)['queries']
            scientific_files = ['gaussian_batch.py','gaussian_overlay.py','interaction_fast.py',
                'interaction_matching.py','chemical_geometry.py','conformer_artifacts.py','chemical_companion.py']
            for name in scientific_files:
                code = Path(__file__).parent/name
                ev.require(legacy['code'].get(str(code.resolve())) == ev.sha(code), 'Legacy scoring code differs: '+name)
            for q in evidence['queries']:
                for key in ('artifact_catalog','chemical_companion','query_npz'):
                    ev.require(legacy['source_hashes'].get(str(Path(q[key]).resolve())) == ev.sha(q[key]),
                               'Legacy query/library input differs')
                ev.require(legacy['source_hashes'].get(str(source)) == ev.sha(source), 'Legacy classification differs')
                qi = next(i for i, old in enumerate(original_queries) if old['query_id'] == q['query_id'])
                q['reuse_directory'] = str(Path(reuse).resolve()/f'query-{qi}')
                q['reuse_chunk_size'] = chunk_size
            sources.update(fingerprint([Path(reuse)/'protocol.json']))
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
                        chunk_size=chunk_size, pose_parameters=POSE_PARAMETERS, thresholds=THRESHOLDS, pose_backend=backend, pose_feasibility=pose_feasibility,
                        numpy=np.__version__, candidate_top_k=None, gaussian_top_n=None,condition_policy=condition,
                        reuse=str(Path(reuse).resolve()) if reuse else None)
        # JSON normalizes tuple values, including thresholds.
        protocol = json.loads(json.dumps(protocol))
        path = output/'protocol.json'
        if path.exists(): ev.require(ev.read(path) == protocol, 'Changed protocol; use a fresh output')
        _atomic_json(path, protocol)
        started = time.perf_counter()
        result = dict(kind='condition_funnel' if condition else 'full_library_conditions', condition_policy=condition,
            status='running', queries=[], sources=sources, hardware=dict(pose_backend=backend,workers=workers,chunk_size=chunk_size),
            library=evidence['library'], classification=evidence['classification'], classes=evidence['classes'],
            unsupported_classes=evidence['unsupported_classes'], e031_changes_ranking=False,
            policy='All catalog conformers; one heuristic anchored-Gaussian pose per conformer, not all possible poses. No Top-K/Top-N membership truncation.',
            biological_quality='not_evaluated', candidate_top_k=None, gaussian_top_n=None)
        if condition:
            result['coarse_eligibility'] = dict(
                constraints=condition.get('coarse_constraints'),
                scope='Explicit whole-ligand eligibility AND requested anchor necessary conditions before seeds; not receptor clash evaluation',
                rejection_target=.9, rejection_target_status='Not guaranteed; measure actual counts')
            result['policy'] = 'All catalog conformers checked by necessary conditions; existing seed poses receive optimistic rule checks; only possible conformers receive unchanged Gaussian winner selection over ALL original seeds. No Top-K/Top-N. Feature counts are conditional on the requested rule.'
            if not pose_feasibility:
                result['policy'] = 'Invariant necessary conditions only; survivors receive unchanged Gaussian winner selection. No Top-K/Top-N. Conditional feature counts.'
            result['unscored_rows'] = 'pose_evaluated=false, assignments=-2; zero scores and identity transforms are placeholders, not calculated results'
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
                progress = dict(query_id=q['query_id'], total_conformers=total,
                    checked_conformers=0, checked_this_invocation=0, reused_completed_rows=0, new_pose_evaluations=0,
                    prefilter_passed=0, prefilter_rejected=0, pose_candidates=0, pose_feasibility_rejected=0, matching_conformers=0, status='running')
                def update_progress(receipt, reused):
                    n = receipt['stop']-receipt['start']
                    progress['checked_conformers'] += n
                    progress['reused_completed_rows' if reused else 'checked_this_invocation'] += n
                    if not reused: progress['new_pose_evaluations'] += receipt.get('new_pose_evaluations',n)
                    progress['prefilter_passed'] += receipt.get('prefilter_passed',n)
                    progress['prefilter_rejected'] += receipt.get('prefilter_rejected',0)
                    progress['pose_candidates'] += receipt.get('pose_candidates',receipt.get('prefilter_passed',n))
                    progress['pose_feasibility_rejected'] += receipt.get('pose_feasibility_rejected',0)
                    progress['matching_conformers'] += receipt.get('matching_conformers',0)
                    _atomic_json(output/'progress.json',progress)
                def tasks():
                    for ci in range(selected):
                        start, stop = ci*chunk_size, min(total, (ci+1)*chunk_size)
                        target = root/f'chunk-{ci:08d}.npz'
                        if target.with_suffix('.receipt.json').exists():
                            load_chunk(target, start, stop, indices)
                            update_progress(ev.read(target.with_suffix('.receipt.json')),True)
                        else: yield start, stop, str(target)
                q['chunks'] = []
                q['counts'] = dict(evaluated_conformers=0, evaluated_molecules=0, total_conformers=total)
                q['full_coverage'] = False
                result['queries'].append(q)
                _atomic_json(output/'report.json', result)
                try:
                    # One CUDA owner; CPU execution uses persistent mapped workers.
                    from contextlib import nullcontext
                    gpu = backend == 'cupy'
                    if gpu: initialize(q)
                    context = nullcontext(None) if gpu else ProcessPoolExecutor(max_workers=workers, initializer=initialize, initargs=(q,))
                    with context as pool:
                        results = ((task,compute_chunk(task)) for task in tasks()) if gpu else bounded_results(pool, compute_chunk, tasks(), workers*2)
                        for _, receipt in results:
                            update_progress(receipt,False)
                            elapsed = time.perf_counter()-started
                            progress['elapsed_seconds'] = elapsed
                            progress['checked_per_second'] = progress['checked_this_invocation']/max(elapsed,1e-9)
                            progress['prefilter_survival_fraction'] = progress['prefilter_passed']/max(progress['checked_conformers'],1)
                            _atomic_json(output/'progress.json',progress)
                            if condition and progress['checked_conformers'] >= chunk_size*workers and progress['prefilter_rejected'] == 0:
                                if progress['checked_conformers'] == chunk_size*workers:
                                    ev.log('Invariant prefilter has rejected zero rows; inspect Gaussian candidate counts to assess the additional pose-feasibility stage. Final matches are separate.')
                            detail = (f"; prefilter kept {receipt['prefilter_passed']}/{receipt['stop']-receipt['start']}; "
                                      f"Gaussian candidates {receipt.get('pose_candidates',receipt['prefilter_passed'])}; matches {receipt['matching_conformers']}" if condition else '')
                            ev.log(f"{q['query_id']}: checked {progress['checked_conformers']} / {total}; "
                                   f"finished IDs {receipt['start']}..{receipt['stop']-1}{detail}")
                    progress['stage'] = 'aggregating_verified_chunks'
                    _atomic_json(output/'progress.json',progress)
                    for ci in range(selected):
                        start, stop = ci*chunk_size, min(total, (ci+1)*chunk_size)
                        target = root/f'chunk-{ci:08d}.npz'
                        count.add(load_chunk(target, start, stop, indices))
                        q['chunks'].append(dict(path=str(target), start=start, stop=stop))
                        sources.update(fingerprint([target, target.with_suffix('.receipt.json')]))
                    q['counts'] = dict(evaluated_conformers=count.processed, evaluated_molecules=count.total_molecules(),
                                       total_conformers=total)
                    if condition:
                        q['counts'].update(prefilter_passed=count.prefilter_passed,
                            prefilter_rejected=count.processed-count.prefilter_passed,
                            pose_candidates=count.pose_candidates,
                            pose_feasibility_rejected=count.prefilter_passed-count.pose_candidates,
                            pose_evaluated_conformers=count.pose_evaluated,
                            matching_conformers=count.matching_conformers,
                            matching_molecules=count.db.execute('SELECT COUNT(*) FROM molecules WHERE matched=1').fetchone()[0])
                        q['diagnostic_scope'] = 'Feature counts only within the molecules/conformers passing the selected rule'
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
            progress.update(status=result['status'],stage='finished')
            _atomic_json(output/'progress.json',progress)
        except Exception as exc:
            result.update(status='failed', error=str(exc))
            _atomic_json(output/'report.json', result)
            _atomic_json(output/'RUN_STATUS.json', dict(status='failed'))
            raise
        return result


def main():
    p = argparse.ArgumentParser(description=__doc__)
    source = p.add_mutually_exclusive_group(required=True)
    source.add_argument('--classification', type=Path)
    source.add_argument('--selection', type=Path, help='Explicit rule preview for conservative full-library funnel')
    p.add_argument('--reuse', type=Path, help='Optional stopped E044 output with compatible verified chunks')
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--workers', type=int, help='CPU default: affinity minus two, capped at 24; GPU: one')
    p.add_argument('--backend', choices=['reference','numpy','cupy'], help='Default: AIDD_POSE_BACKEND or numpy')
    p.add_argument('--chunk-size', type=int, default=256)
    p.add_argument('--no-pose-feasibility',action='store_true',help='Reference comparison: disable rule-aware seed rejection')
    p.add_argument('--max-chunks', type=int, help='Explicit pilot limit per query; incomplete coverage is labeled partial')
    a = p.parse_args()
    options = dict(workers=a.workers,chunk_size=a.chunk_size,max_chunks=a.max_chunks,reuse=a.reuse,backend=a.backend,pose_feasibility=False if a.no_pose_feasibility else None)
    result = run_funnel(a.selection,a.output,**options) if a.selection else run(a.classification,a.output,**options)
    print(json.dumps(dict(status=result['status'], report=str(a.output/'report.json'))))


if __name__ == '__main__': main()
