"""E118 finite-pose soft ranking inside the completed KRAS direct-contact pool."""
import argparse
from concurrent.futures import ProcessPoolExecutor, wait, FIRST_COMPLETED
import hashlib
import json
import math
import os
from pathlib import Path
import sqlite3
import sys
import time

for variable in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[variable] = '1'
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'src'))
import numpy as np
import run_kras_weekend as run
from aidd_agent.prompt_workflow import file_lock

POLICY = dict(version='E118-soft-v1', minimum_distance=2.0, retained_conformers=3,
              score='mean(shape, residue_feature_match, mean_three_region_coverage)',
              color='diagnostic_only', score_fraction=.8, diversity='graph_simhash_12bit')
COUNTS = None


def readonly(path):
    return sqlite3.connect(Path(path).resolve().as_uri() + '?mode=ro', uri=True)


def code_hash():
    return hashlib.sha256(Path(__file__).read_bytes().replace(b'\r\n', b'\n')).hexdigest()


def diversity_bucket(chem):
    """Coarse graph similarity bucket, independent of atom order; not a scaffold ID."""
    labels = [f'{int(z)}:{int(q)}:{int(f)}' for z, q, f in
              zip(chem.atomic_numbers, chem.formal_charges, chem.atom_flags)]
    adjacency = [[] for _ in labels]
    for bond in chem.bonds:
        a, b, order = int(bond['begin']), int(bond['end']), int(bond['order'])
        adjacency[a].append((b, order)); adjacency[b].append((a, order))
    tokens = list(labels)
    for _ in range(2):
        labels = [hashlib.sha256((labels[i] + '|' + '|'.join(sorted(
            f'{order}:{labels[j]}' for j, order in neighbors))).encode()).hexdigest()
                  for i, neighbors in enumerate(adjacency)]
        tokens.extend(labels)
    votes = [0] * 12
    for token in tokens:
        value = int.from_bytes(hashlib.sha256(token.encode()).digest()[:2], 'big')
        for bit in range(12):
            votes[bit] += 1 if value & (1 << bit) else -1
    return sum((1 << i) for i, value in enumerate(votes) if value >= 0)


def score_candidate(candidate, features, chemistry, state):
    kw = run.kw
    sets = kw.chemical_sets(chemistry, candidate.feature_types)
    stats = dict(seeds=0, direct_poses=0, geometry_poses=0)
    if not all(len(x) for x in sets.values()):
        return None, stats
    query = state['query']
    seeds, _ = kw.prepare_seeds(candidate, query, **kw.POSE_PARAMETERS)
    matrices = np.array([s.transform_matrix for s in seeds]).reshape(-1, 4, 4)
    stats['seeds'] = len(matrices)
    matrices = matrices[kw.direct_contact_masks(candidate.shape_points, sets, state, matrices)['all']]
    stats['direct_poses'] = len(matrices)
    if not len(matrices):
        return None, stats
    moved = np.einsum('bij,nj->bni', matrices[:, :3, :3], candidate.shape_points) + matrices[:, None, :3, 3]
    nearest = state['tree'].query(moved.reshape(-1, 3))[0].reshape(len(matrices), -1)
    minimum = nearest.min(axis=1)
    keep = minimum >= POLICY['minimum_distance']
    matrices, moved, minimum = matrices[keep], moved[keep], minimum[keep]
    stats['geometry_poses'] = len(matrices)
    if not len(matrices):
        return None, stats
    indices = query['anchor_feature_indices']
    q = tuple(query[k][indices] for k in ('feature_points', 'feature_types', 'feature_directions', 'feature_direction_kinds'))
    types = features.feature_types.copy()
    for i, kind in enumerate(types):
        if kind in (1, 2) and not any(chemistry.atomic_numbers[j] in (7, 8) for j in chemistry.feature_members[i]):
            types[i] = 0
    n = len(matrices)
    _, assignments, matches = kw.match_batch((*q, np.ones(len(indices))),
        np.repeat(features.feature_points[None], n, axis=0), np.repeat(types[None], n, axis=0),
        np.repeat(features.feature_directions[None], n, axis=0),
        np.repeat(features.feature_kinds[None], n, axis=0), matrices)
    matches = np.where(assignments >= 0, matches, 0)
    lookup = {int(v): i for i, v in enumerate(indices)}
    groups = {}
    for anchor in state['record']['anchors']:
        groups.setdefault(anchor['target_residue'], set()).add(lookup[anchor['feature_index']])
    feature_score = np.mean([matches[:, sorted(group)].max(axis=1) for group in groups.values()], axis=0)
    best = None
    for j, matrix in enumerate(matrices):
        tree = kw.cKDTree(moved[j])
        geometry = state['geometry']
        regions = [float(np.mean(tree.query(geometry[f'region{i}'])[0] <= 1.5)) for i in range(3)]
        envelope = float(np.mean(tree.query(geometry['contact'])[0] <= 1.5))
        overlay = kw.score_overlay(query['shape_points'], candidate.shape_points,
            query['feature_points'], query['feature_types'], candidate.feature_points,
            candidate.feature_types, matrix.ravel(), cutoff=None)
        shape, color = float(overlay['shape']['tanimoto']), float(overlay['color']['tanimoto'])
        score = (shape + float(feature_score[j]) + float(np.mean(regions))) / 3
        if not all(math.isfinite(x) for x in (score, shape, color, minimum[j])):
            raise ValueError('Non-finite pose score')
        pose = dict(score=score, shape=shape, color=color, feature_match=float(feature_score[j]),
                    envelope=envelope, regions=regions, minimum_distance=float(minimum[j]),
                    transform=matrix.ravel().tolist(), seed_order=j)
        if best is None or pose['score'] > best['score']:
            best = pose
    best.update(gid=int(candidate.global_id), mid=str(candidate.molecule_id), cid=str(candidate.conformer_id))
    return best, stats


def initialize(batch, prepared, counts):
    global COUNTS
    run.initialize(batch, prepared)
    COUNTS = readonly(counts)


def process(task):
    index, ids = task
    reader, feature_reader, states, chemistry_reader = run.WORK
    output = {t: dict(pool_conformers=0, seeds=0, direct_poses=0, geometry_poses=0,
                      geometry_conformers=0, rows=[]) for t in states}
    for gid in ids:
        candidate = reader.get(gid)
        targets = [t for t in states if COUNTS.execute(
            'SELECT 1 FROM molecules WHERE target=? AND mid=?', (t, candidate.molecule_id)).fetchone()]
        if not targets:
            continue
        chemistry = chemistry_reader.get(gid); features = feature_reader.get(gid)
        identity = (candidate.molecule_id, candidate.conformer_id)
        if any((r.molecule_id, r.conformer_id) != identity for r in (chemistry, features)):
            raise ValueError('Candidate/chemical/feature identity mismatch')
        bucket = None
        for target in targets:
            row = output[target]; row['pool_conformers'] += 1
            pose, stats = score_candidate(candidate, features, chemistry, states[target])
            for key, value in stats.items():
                row[key] += value
            if pose is not None:
                row['geometry_conformers'] += 1
                if bucket is None:
                    bucket = diversity_bucket(chemistry)
                pose['bucket'] = bucket
                row['rows'].append(pose)
    return index, len(ids), output


def create_tables(db):
    db.executescript('''CREATE TABLE IF NOT EXISTS chunks(id INTEGER PRIMARY KEY, processed INTEGER, counts TEXT);
    CREATE TABLE IF NOT EXISTS molecules(target TEXT, mid TEXT, score REAL, bucket INTEGER,
        payload TEXT, PRIMARY KEY(target,mid));''')


def commit_chunk(db, result):
    index, processed, output = result
    with db:
        if db.execute('SELECT 1 FROM chunks WHERE id=?', (index,)).fetchone():
            raise ValueError('Chunk already committed')
        for target, rows in output.items():
            for pose in rows['rows']:
                old = db.execute('SELECT payload FROM molecules WHERE target=? AND mid=?', (target, pose['mid'])).fetchone()
                poses = {p['cid']: p for p in json.loads(old[0])} if old else {}
                if pose['cid'] not in poses or pose['score'] > poses[pose['cid']]['score']:
                    poses[pose['cid']] = pose
                poses = sorted(poses.values(), key=lambda p: (-p['score'], p['gid'], p['cid']))[:POLICY['retained_conformers']]
                best = poses[0]
                db.execute('INSERT OR REPLACE INTO molecules VALUES (?,?,?,?,?)',
                           (target, pose['mid'], best['score'], best['bucket'], json.dumps(poses)))
        counts = {t: {k: v for k, v in rows.items() if k != 'rows'} for t, rows in output.items()}
        db.execute('INSERT INTO chunks VALUES (?,?,?)', (index, processed, json.dumps(counts)))


def summary(db, scheduled):
    processed = 0; totals = {t: {} for t in run.kw.REFERENCES}
    for n, raw in db.execute('SELECT processed,counts FROM chunks'):
        processed += n
        for target, counts in json.loads(raw).items():
            for key, value in counts.items():
                totals[target][key] = totals[target].get(key, 0) + value
    for target in totals:
        totals[target]['eligible_molecules'] = db.execute('SELECT COUNT(*) FROM molecules WHERE target=?', (target,)).fetchone()[0]
        totals[target]['score_ranges'] = {}
        for metric in ('score', 'shape', 'color', 'feature_match', 'envelope', 'minimum_distance'):
            # SQLite JSON support is standard in supported workstation Python builds.
            totals[target]['score_ranges'][metric] = db.execute(
                "SELECT MIN(json_extract(payload,?)),MAX(json_extract(payload,?)) FROM molecules WHERE target=?",
                (f'$[0].{metric}', f'$[0].{metric}', target)).fetchone()
    return dict(status='complete' if processed == scheduled else 'partial', processed=processed,
                scheduled=scheduled, targets=totals)


def controls(prepared):
    from types import SimpleNamespace
    output = {}
    for target in run.kw.REFERENCES:
        state = run.kw.load_state(prepared, target); q = state['query']
        rotation = np.array([[0, -1, 0], [0, 0, -1], [1, 0, 0]])
        move = lambda x: x @ rotation.T + np.array([17, -11, 6])
        candidate = SimpleNamespace(global_id=-1, molecule_id='control', conformer_id='control',
            shape_points=move(q['shape_points']), feature_points=move(q['feature_points']), feature_types=q['feature_types'])
        features = SimpleNamespace(feature_points=candidate.feature_points, feature_types=q['feature_types'],
            feature_directions=q['feature_directions'] @ rotation.T, feature_kinds=q['feature_direction_kinds'])
        chemistry = SimpleNamespace(**state['record']['chemistry'])
        pose, stats = score_candidate(candidate, features, chemistry, state)
        if pose is None or pose['shape'] < .99:
            raise ValueError('Native control failed: ' + target)
        output[target] = dict(pose=pose, counts=stats)
    return output


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--direct', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--sample', type=int, default=10000)
    p.add_argument('--workers', type=int, default=20)
    p.add_argument('--chunk-size', type=int, default=128)
    p.add_argument('--hours', type=float, default=0)
    p.add_argument('--pilot', type=Path)
    args = p.parse_args()
    if min(args.sample, args.hours) < 0 or min(args.workers, args.chunk_size) < 1:
        p.error('Invalid budget')
    source = args.direct.resolve(); out = args.output.resolve()
    report = run.kw.read(source / 'report.json'); protocol = run.kw.read(source / 'protocol.json')
    if report.get('status') != 'complete' or not report.get('full_catalog') or report.get('scope') != 'direct_contacts_only_not_final_screening_hits':
        p.error('Complete full direct-contact count required')
    original = protocol['original']; batch = Path(original['batch']); prepared = Path(original['prepared'])
    for protected in (source, batch.resolve(), prepared.resolve()):
        if out.is_relative_to(protected) or protected.is_relative_to(out):
            p.error('Output overlaps protected input')
    if run.seal_inputs(batch, prepared, original['chunk_size'], original['sample_size']) != original:
        p.error('Original input/code lineage changed; do not bypass this check')
    counts_path = source / 'counts.sqlite'
    with readonly(counts_path) as counts:
        for target, expected in report['targets'].items():
            actual = counts.execute('SELECT COUNT(*) FROM molecules WHERE target=?', (target,)).fetchone()[0]
            if actual != expected['distinct_molecules']:
                raise ValueError('Direct-pool database count mismatch')
    common = dict(original=original, policy=POLICY, code_hash=code_hash(),
                  direct_protocol_sha256=run.kw.sha(source / 'protocol.json'),
                  direct_report_sha256=run.kw.sha(source / 'report.json'),
                  counts_sha256=run.kw.sha(counts_path))
    if not args.sample:
        if not args.pilot:
            p.error('Full scoring requires --pilot PATH/report.json')
        pilot = run.kw.read(args.pilot)
        pilot_protocol = run.kw.read(args.pilot.parent / 'protocol.json')
        if pilot.get('status') != 'complete' or pilot.get('scope') != 'spread_pilot' or pilot_protocol['common'] != common:
            p.error('Pilot missing, incomplete or lineage mismatch')
        if any(pilot.get('targets', {}).get(t, {}).get('eligible_molecules', 0) == 0 for t in run.kw.REFERENCES):
            p.error('Pilot has a target with zero geometry-qualified molecules; inspect before a full run')
    total = run.catalog_count(run.kw.read(batch / 'artifacts/catalog.json'))
    scheduled = min(total, args.sample) if args.sample else total
    signature = dict(common=common, sample=args.sample, chunk_size=args.chunk_size)
    out.mkdir(parents=True, exist_ok=True)
    with file_lock(out / 'run.lock'):
        if (out / 'protocol.json').exists() and run.kw.read(out / 'protocol.json') != signature:
            raise ValueError('Resume mismatch; use a new output directory')
        run.kw.save(out / 'protocol.json', signature)
        native = controls(prepared); run.kw.save(out / 'controls.json', native)
        db = sqlite3.connect(out / 'ranking.sqlite'); create_tables(db)
        done = {r[0] for r in db.execute('SELECT id FROM chunks')}
        pending = (t for t in run.tasks(total, args.chunk_size, args.sample) if t[0] not in done)
        started = time.monotonic(); last = started; added = 0; stopped = False
        run.kw.save(out / 'report.json', dict(status='running', scope='spread_pilot' if args.sample else 'full_direct_pool_reranking'))
        with ProcessPoolExecutor(max_workers=args.workers, initializer=initialize,
                                 initargs=(batch, prepared, counts_path)) as pool:
            active = set(); exhausted = False
            while active or not exhausted:
                while len(active) < args.workers and not exhausted:
                    if args.hours and time.monotonic() - started >= args.hours * 3600:
                        exhausted = True; stopped = True; break
                    task = next(pending, None)
                    if task is None:
                        exhausted = True; break
                    active.add(pool.submit(process, task))
                if not active:
                    break
                finished, active = wait(active, return_when=FIRST_COMPLETED)
                for future in finished:
                    result = future.result(); commit_chunk(db, result); added += result[1]
                if time.monotonic() - last >= 30:
                    progress = dict(invocation_new_conformers=added, invocation_seconds=time.monotonic()-started,
                                    committed_chunks=db.execute('SELECT COUNT(*) FROM chunks').fetchone()[0])
                    run.kw.save(out / 'progress.json', progress); print(json.dumps(progress), flush=True); last=time.monotonic()
        result = summary(db, scheduled); db.close()
        result.update(scope='spread_pilot' if args.sample else 'full_direct_pool_reranking',
                      invocation_seconds=time.monotonic()-started, invocation_new_conformers=added,
                      full_catalog=not args.sample, stopped_on_budget=stopped, native_controls_passed=True,
                      limitations=['Finite rigid seed recovery for pool molecules; not cached-pose rescoring.',
                                   'Exploratory soft ranking; no activity or allele selectivity validation.',
                                   'Distance floor is a gross-clash proxy; coverage is reference-envelope coverage.',
                                   'Three conformers retained per molecule; exported representative is not an energy minimum.'])
        result['ranking_sha256'] = run.kw.sha(out / 'ranking.sqlite')
        if not args.sample and result['status'] == 'complete':
            result['direct_pose_reconciliation'] = {t: result['targets'][t].get('direct_poses', 0) == report['targets'][t]['direct_poses']
                                                   for t in run.kw.REFERENCES}
            if not all(result['direct_pose_reconciliation'].values()):
                result['status'] = 'failed_direct_pose_reconciliation'
        if args.sample and added:
            result['full_hours_linear_estimate'] = result['invocation_seconds'] / added * total / 3600
        run.kw.save(out / 'report.json', result); print(json.dumps(result, indent=2), flush=True)
        if result['status'].startswith('failed'):
            raise SystemExit(2)


if __name__ == '__main__':
    main()
