"""Full-library ANY-anchor funnel with pose combinations and scaffold groups.

No pilot prerequisite, Top-K, or candidate cap. Gaussian winner scores rank
conformers; they do not erase matching alternative original-seed poses.
"""
from __future__ import annotations

import argparse
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
import copy
import csv
import hashlib
import json
import multiprocessing
from pathlib import Path
import sqlite3
import time

import numpy as np

from . import full_library_screen as full
from .chunk_execution import bounded_results
from .gaussian_batch import prepare_seeds, _score_ids
from .interaction_fast import match_batch
from .joint_coarse import extents
from .pose_feasibility import possible_seed_mask
from .prompt_workflow import file_lock

_GUIDED_CACHE = None


def record_json(path, document):
    full._atomic_json(Path(path), document)


def coarse_stages(candidate, features, original, rule, bound, joint=None, timings=None):
    """Independent eligibility rules, in the recorded execution order."""
    timings = Counter() if timings is None else timings
    t = time.perf_counter()
    low, high = rule['heavy_atom_ratio']
    size_ok = low <= len(candidate.shape_points) / len(original['shape_points']) <= high
    timings['size_filter'] += time.perf_counter()-t
    if not size_ok:
        return 0
    t = time.perf_counter()
    types, counts = (np.unique(original['feature_types'], return_counts=True) if joint is None
                     else (joint.types, joint.counts))
    ct, cc = np.unique(candidate.feature_types, return_counts=True)
    table = dict(zip(ct.tolist(), cc.tolist()))
    coverage = sum(min(int(n), table.get(int(t), 0)) for t, n in zip(types, counts)) / counts.sum()
    timings['typed_coverage_filter'] += time.perf_counter()-t
    if coverage < rule['minimum_feature_coverage']:
        return 1
    t = time.perf_counter()
    qe = extents(original['shape_points']) if joint is None else joint.extent
    distance = np.linalg.norm(extents(candidate.shape_points) - qe) / np.linalg.norm(qe)
    timings['extent_filter'] += time.perf_counter()-t
    if distance > rule['maximum_extent_distance'] + 1e-7:
        return 2
    return 4 if bound is not None and bound.check(features)[0] else 3


def pose_representatives(features, seeds, possible, expanded, columns, threshold, design=None, anchor_order=None, gaussian_scores=None, shape_points=None, budget_mode=False):
    """Keep one actual pose per exact required-anchor bitmask, never a union pose."""
    indices = expanded['anchor_feature_indices']
    query = tuple(expanded[k][indices] for k in (
        'feature_points', 'feature_types', 'feature_directions',
        'feature_direction_kinds', 'anchored_weights'))
    selected = np.flatnonzero(possible)
    best = {}
    matching_seeds = 0
    for start in range(0, len(selected), 32):
        numbers = selected[start:start+32]
        matrices = np.asarray([seeds[i].transform_matrix for i in numbers]).reshape(-1, 4, 4)
        arrays = [np.repeat(getattr(features, name)[None], len(numbers), axis=0)
                  for name in ('feature_points', 'feature_types', 'feature_directions', 'feature_kinds')]
        _, assignments, scores = match_batch(query, *arrays, matrices)
        for i, number in enumerate(numbers):
            values, assigned = scores[i, columns], assignments[i, columns]
            hits = (assigned >= 0) & (values >= threshold)
            mask = sum(1 << j for j, passed in enumerate(hits) if passed)
            if not mask and not budget_mode:
                continue
            if design is not None and not budget_mode:
                from .guided_filters import rule_passes
                if not rule_passes(mask,anchor_order,design):continue
            quality = float(values[hits].min()) if hits.any() else 0.
            row = dict(mask=mask, matched_anchor_count=int(hits.sum()),
                       min_matched_score=quality, seed_index=int(number),
                       seed_id=str(seeds[number].seed_id), scores=values.tolist(),
                       assignments=assigned.tolist(), transform=matrices[i].reshape(-1).tolist())
            if gaussian_scores is not None:
                from .guided_filters import pose_rank
                row.update(pose_rank(values,assigned,anchor_order,gaussian_scores[number],shape_points,matrices[i],design))
                if budget_mode:row['legacy_score_threshold_passed'] = row['composite_score']>=design['minimum_pose_score']
                elif row['composite_score']<design['minimum_pose_score']:continue
            matching_seeds += 1
            rank=(row.get('optional_score',quality),row.get('gaussian_same_pose',0.)) if budget_mode else row.get('composite_score',quality)
            signature=json.dumps(sorted(row.get('occupied_spatial_groups',[])),separators=(',',':'))
            key=(0,'') if budget_mode else (mask,signature)
            old_rank=((best[key].get('optional_score',best[key]['min_matched_score']),best[key].get('gaussian_same_pose',0.))
                      if budget_mode and key in best else best[key].get('composite_score',best[key]['min_matched_score']) if key in best else None)
            if key not in best or rank > old_rank:
                best[key] = row
    return [best[k] for k in sorted(best)], matching_seeds


def compute(task):
    global _GUIDED_CACHE
    start, stop, target = task[:3]
    target = Path(target)
    reader, original, expanded, q = full._STATE
    policy = q['condition_policy']
    mode=q.get('selection_mode','threshold')
    if mode not in ('threshold','budget'):raise ValueError('Unknown selection mode')
    budget_mode=mode=='budget'
    guided=None
    if q.get('guided_design'):
        from .guided_filters import GuidedFilter
        key=json.dumps(q['guided_design'],sort_keys=True)
        if _GUIDED_CACHE is None or _GUIDED_CACHE[0]!=key:
            _GUIDED_CACHE=(key,GuidedFilter(q,expanded,q['guided_design']))
        guided=_GUIDED_CACHE[1]
    columns = [next(a['score_column'] for a in q['anchors'] if a['anchor_id'] == aid)
               for aid in policy['required_anchors']]
    counts = Counter({key:0 for key in ('size_passed','coverage_passed','extent_passed',
        'anchor_bound_passed','original_seeds','possible_seeds','pose_feasibility_passed',
        'gaussian_evaluated_conformers','exact_seed_annotations','matching_seeds',
        'matching_conformers','pose_combination_records')})
    counts['input_conformers'] = stop-start
    seconds = Counter()
    tick = time.perf_counter()
    ids = np.arange(start, stop, dtype=np.int64) if len(task)==3 else np.asarray(task[3],dtype=np.int64)
    if ids.ndim!=1 or len(ids)!=stop-start or len(np.unique(ids))!=len(ids) or np.any(ids<0):
        raise ValueError('Invalid explicit pilot conformer IDs')
    mids = []
    levels = []
    records = []
    for gid in ids:
        t = time.perf_counter()
        candidate = reader.get(int(gid))
        mids.append(candidate.molecule_id)
        seconds['artifact_read'] += time.perf_counter()-t
        if budget_mode:
            stage=4
        elif guided is not None and getattr(guided,'design',{}).get('adaptive_coarse'):
            from .guided_filters import adaptive_stage
            stage=adaptive_stage(candidate,guided.design['adaptive_coarse'])
        else:
            stage = coarse_stages(candidate, None, original, policy['coarse_constraints'], None,
                                  joint=full._JOINT, timings=seconds)
        if stage >= 3:
            t = time.perf_counter()
            features = full._FILTER_READER.get(int(gid))
            if (candidate.molecule_id != features.molecule_id or
                candidate.conformer_id != features.conformer_id or
                not np.array_equal(candidate.feature_types, features.feature_types) or
                not np.array_equal(candidate.feature_points, features.feature_points)):
                raise ValueError(f'Artifact/chemical mismatch: {gid}')
            stage = 4 if budget_mode or full._BOUND.check(features)[0] else 3
            seconds['chemical_read_and_anchor_bound'] += time.perf_counter()-t
        levels.append(stage)
        for j, name in enumerate(('size_passed', 'coverage_passed', 'extent_passed', 'anchor_bound_passed')):
            counts[name] += int(stage > j)
        if stage != 4:
            continue
        if guided is not None and not budget_mode:
            t=time.perf_counter()
            geometry_ok=guided.geometry_possible(features)
            seconds['mandatory_feature_geometry']+=time.perf_counter()-t
            if not geometry_ok:
                counts['guided_geometry_rejected']+=1
                continue
            counts['guided_geometry_passed']+=1
        t = time.perf_counter()
        seeds, pairs = prepare_seeds(candidate, original, **full.POSE_PARAMETERS)
        seconds['seed_generation'] += time.perf_counter()-t
        counts['original_seeds'] += len(seeds)
        t = time.perf_counter()
        possible = np.ones(len(seeds),dtype=bool) if budget_mode else possible_seed_mask(full._BOUND, features, seeds)
        if guided is not None and not budget_mode:possible=guided.geometry_seeds(features,seeds,possible)
        seconds['seed_feasibility'] += time.perf_counter()-t
        counts['possible_seeds'] += int(possible.sum())
        if not possible.any():
            continue
        if guided is not None:
            t=time.perf_counter()
            positions=np.flatnonzero(possible)
            transforms=np.asarray([seeds[i].transform_matrix for i in positions]).reshape(-1,4,4)
            pocket=guided.pocket_mask(candidate.shape_points,transforms)
            counts['pocket_tested_seeds']+=len(positions)
            counts['pocket_rejected_seeds']+=int((~pocket).sum())
            possible[positions]=pocket
            seconds['pocket_exclusion']+=time.perf_counter()-t
            if not possible.any():
                counts['pocket_rejected_conformers']+=1
                continue
            counts['pocket_passed_conformers']+=1
        counts['pose_feasibility_passed'] += 1
        levels[-1] = 5
        t = time.perf_counter()
        scoring_seeds = seeds if guided is None else [seeds[i] for i in np.flatnonzero(possible)]
        gaussian_scores=None
        if q.get('consensus_npz'):
            from .guided_filters import same_pose_gaussian
            gaussian_scores=np.full(len(seeds),np.nan)
            gaussian_scores[np.flatnonzero(possible)]=same_pose_gaussian(original,candidate,scoring_seeds)
            gaussian_best=float(np.nanmax(gaussian_scores))
        else:
            rigid = _score_ids(reader, original, np.array([gid]), backend='numpy',
                               prepared={int(gid): (candidate, scoring_seeds, pairs)}, **full.POSE_PARAMETERS)
            gaussian_best=float(rigid[full.OBJECTIVE+'__objective'][0])
        counts['gaussian_evaluated_seeds'] += len(scoring_seeds)
        seconds['gaussian'] += time.perf_counter()-t
        counts['gaussian_evaluated_conformers'] += 1
        levels[-1] = 6
        t = time.perf_counter()
        poses, matching_seeds = pose_representatives(features, seeds, possible, expanded, columns,
                                                     policy['minimum_score'],q.get('guided_design'),policy['required_anchors'],gaussian_scores,candidate.shape_points,budget_mode=budget_mode)
        seconds['anchor_assignment'] += time.perf_counter()-t
        counts['exact_seed_annotations'] += int(possible.sum())
        counts['matching_seeds'] += matching_seeds
        counts['matching_conformers'] += int(bool(poses))
        if poses:
            levels[-1] = 7
        counts['pose_combination_records'] += len(poses)
        for pose in poses:
            pose.update(global_id=int(gid), molecule_id=candidate.molecule_id,
                        conformer_id=candidate.conformer_id,
                        matched_anchors=[aid for j, aid in enumerate(policy['required_anchors'])
                                         if pose['mask'] & (1 << j)],
                        conformer_gaussian_best=gaussian_best)
            records.append(pose)
    t = time.perf_counter()
    full._atomic_savez(target, dict(global_ids=ids, molecule_ids=np.asarray(mids),stage_levels=np.asarray(levels,dtype=np.uint8)))
    poses_path = target.with_suffix('.poses.jsonl')
    temporary = poses_path.with_suffix('.partial')
    with temporary.open('w', encoding='utf-8') as stream:
        for row in records:
            stream.write(json.dumps(row)+'\n')
    temporary.replace(poses_path)
    seconds['persistence'] += time.perf_counter()-t
    receipt = dict(start=start, stop=stop, counts=dict(counts), worker_seconds=dict(seconds),
                   wall_seconds=time.perf_counter()-tick,
                   files=full.fingerprint([target, poses_path]))
    record_json(target.with_suffix('.receipt.json'), receipt)
    return receipt


def verify_chunk(target, start, stop):
    receipt = full.ev.read(target.with_suffix('.receipt.json'))
    if (receipt['start'], receipt['stop']) != (start, stop):
        raise ValueError('Chunk range mismatch')
    expected = {str(target.resolve()), str(target.with_suffix('.poses.jsonl').resolve())}
    if set(receipt['files']) != expected:
        raise ValueError('Chunk file identity mismatch')
    full.check_hashes(receipt['files'])
    with np.load(target, allow_pickle=False) as z:
        if (not np.array_equal(z['global_ids'], np.arange(start, stop)) or len(z['molecule_ids']) != stop-start
                or z['stage_levels'].shape != (stop-start,) or np.any(z['stage_levels'] > 7)):
            raise ValueError('Chunk coverage mismatch')
    return receipt


def schema(db):
    db.executescript('''
        CREATE TABLE IF NOT EXISTS library_molecules(mid TEXT PRIMARY KEY, stage INTEGER);
        CREATE TABLE IF NOT EXISTS poses(mid TEXT, mask INTEGER, quality REAL, gid INTEGER,
                                        payload TEXT, spatial TEXT NOT NULL, PRIMARY KEY(mid,mask,spatial));
        CREATE TABLE IF NOT EXISTS members(mid TEXT PRIMARY KEY, cluster TEXT, scaffold TEXT,
                                          error TEXT, union_mask INTEGER, pose_count INTEGER);
    ''')
    if 'spatial' not in {r[1] for r in db.execute('PRAGMA table_info(poses)')}:
        raise ValueError('Legacy pose database requires a fresh run after grouped-scoring code changes')


def merge_pose(db, row):
    """SQL union is only by molecule; each payload remains one real pose."""
    signature=json.dumps(sorted(row.get('occupied_spatial_groups',[])),separators=(',',':'))
    old = db.execute('SELECT quality,gid FROM poses WHERE mid=? AND mask=? AND spatial=?',
                     (row['molecule_id'], row['mask'],signature)).fetchone()
    quality, gid = row.get('composite_score',row['min_matched_score']), row['global_id']
    if old is None or quality > old[0] or (quality == old[0] and gid < old[1]):
        db.execute('INSERT OR REPLACE INTO poses VALUES (?,?,?,?,?,?)',
                   (row['molecule_id'], row['mask'], quality, gid, json.dumps(row),signature))


def scaffold_key(chemistry):
    """Exact non-stereochemical Murcko scaffold grouping, not similarity clustering."""
    from rdkit import Chem
    from rdkit.Chem.Scaffolds import MurckoScaffold
    mol = Chem.RWMol()
    for index,(number, charge) in enumerate(zip(chemistry.atomic_numbers, chemistry.formal_charges)):
        atom = Chem.Atom(int(number)); atom.SetFormalCharge(int(charge))
        hydrogens=getattr(chemistry,'total_hydrogens',None)
        if hydrogens is not None:
            atom.SetNumExplicitHs(int(hydrogens[index]));atom.SetNoImplicit(True)
        mol.AddAtom(atom)
    orders = {1: Chem.BondType.SINGLE, 2: Chem.BondType.DOUBLE,
              3: Chem.BondType.TRIPLE, 4: Chem.BondType.AROMATIC}
    for bond in chemistry.bonds:
        left, right, order = int(bond['begin']), int(bond['end']), int(bond['order'])
        mol.AddBond(left, right, orders[order])
        if order == 4:
            mol.GetAtomWithIdx(left).SetIsAromatic(True)
            mol.GetAtomWithIdx(right).SetIsAromatic(True)
    mol = mol.GetMol(); Chem.SanitizeMol(mol)
    scaffold = MurckoScaffold.GetScaffoldForMol(mol)
    # Do not collapse all acyclic molecules into one empty-scaffold group.
    smiles = Chem.MolToSmiles(scaffold if scaffold.GetNumAtoms() else mol, isomericSmiles=False)
    return 'scaffold-'+hashlib.sha256(smiles.encode()).hexdigest(), smiles


def aggregate(root, q, total, chunk_size, expected_molecules):
    from .chemical_companion import ChemicalCompanionReader
    db = sqlite3.connect(root/'candidates.sqlite')
    schema(db)
    # Derived outputs are rebuilt; sealed scan chunks remain untouched.
    with db:
        for table in ('library_molecules', 'poses', 'members'):
            db.execute('DELETE FROM '+table)
    try:
        for start in range(0, total, chunk_size):
            target = root/'chunks'/f'chunk-{start//chunk_size:08d}.npz'
            verify_chunk(target, start, min(start+chunk_size, total))
            with np.load(target, allow_pickle=False) as z, db:
                db.executemany('INSERT INTO library_molecules VALUES (?,?) ON CONFLICT(mid) DO UPDATE SET stage=MAX(stage,excluded.stage)',
                               ((str(mid),int(stage)) for mid,stage in zip(z['molecule_ids'],z['stage_levels'])))
                with target.with_suffix('.poses.jsonl').open(encoding='utf-8') as stream:
                    for line in stream:
                        row = json.loads(line)
                        if not start <= row['global_id'] < min(start+chunk_size, total):
                            raise ValueError('Pose global ID outside chunk')
                        merge_pose(db, row)
        n = db.execute('SELECT COUNT(*) FROM library_molecules').fetchone()[0]
        if n != expected_molecules:
            raise ValueError(f'Molecule coverage mismatch: {n} != {expected_molecules}')
        coarse_end = time.perf_counter()
        chemistry = ChemicalCompanionReader(Path(q['chemical_companion']))
        failures = 0
        with (root/'candidate-poses.jsonl').open('w', encoding='utf-8') as stream, db:
            for mid, gid in db.execute('SELECT mid,MIN(gid) FROM poses GROUP BY mid ORDER BY mid'):
                chem = chemistry.get(gid)
                if chem.molecule_id != mid:
                    raise ValueError('Scaffold source molecule identity mismatch')
                error = ''
                try:
                    key, smiles = scaffold_key(chem)
                except Exception as exc:
                    # Preserve every failed member as its own explicit group.
                    failures += 1; error = str(exc); smiles = ''
                    key = 'unresolved-'+hashlib.sha256(mid.encode()).hexdigest()
                union = 0; count = 0
                for mask, payload in db.execute('SELECT mask,payload FROM poses WHERE mid=? ORDER BY mask', (mid,)):
                    union |= mask; count += 1
                    row = json.loads(payload); row['cluster_id'] = key
                    stream.write(json.dumps(row)+'\n')
                db.execute('INSERT INTO members VALUES (?,?,?,?,?,?)', (mid,key,smiles,error,union,count))
        with (root/'cluster-members.csv').open('w', newline='', encoding='utf-8') as stream:
            writer = csv.writer(stream)
            writer.writerow(['molecule_id','cluster_id','scaffold','grouping_error','union_across_poses_not_simultaneous','pose_combination_count'])
            writer.writerows(db.execute('SELECT * FROM members ORDER BY cluster,mid'))
        with (root/'clusters.csv').open('w', newline='', encoding='utf-8') as stream:
            writer = csv.writer(stream); writer.writerow(['cluster_id','scaffold','molecules','display_member'])
            writer.writerows(db.execute('SELECT cluster,MIN(scaffold),COUNT(*),MIN(mid) FROM members GROUP BY cluster ORDER BY COUNT(*) DESC,cluster'))
        return dict(evaluated_molecules=n,
                    molecules_by_stage=[db.execute('SELECT COUNT(*) FROM library_molecules WHERE stage>=?',(i,)).fetchone()[0] for i in range(8)],
                    matching_molecules=db.execute('SELECT COUNT(*) FROM members').fetchone()[0],
                    pose_combination_records=db.execute('SELECT COUNT(*) FROM poses').fetchone()[0],
                    clusters=db.execute('SELECT COUNT(DISTINCT cluster) FROM members').fetchone()[0],
                    scaffold_failures_preserved_as_singletons=failures,
                    scaffold_and_export_seconds=time.perf_counter()-coarse_end,
                    combinations=[dict(mask=mask, molecules=count) for mask,count in
                                  db.execute('SELECT mask,COUNT(*) FROM poses GROUP BY mask ORDER BY mask')])
    finally:
        db.close()


def run(selection, output, workers=22, chunk_size=2048):
    from rdkit import rdBase
    selection, output = Path(selection).resolve(), Path(output).resolve()
    selected = full.ev.read(selection)
    if selected.get('kind') != 'selection_preview' or workers < 1 or chunk_size < 1:
        raise ValueError('Need selection_preview and positive workers/chunk size')
    full.check_hashes(selected['sources'])
    source = Path(selected['evidence_report']).resolve()
    if str(source) not in selected['sources']:
        raise ValueError('Classification is not sealed by the selection')
    evidence = full.ev.read(source)
    if evidence.get('status') != 'complete' or evidence.get('kind') != 'screening_evidence':
        raise ValueError('Need completed classified evidence')
    policy = full.normalize_policy(dict(selected['policy'], match_mode='any'))
    if 'coarse_constraints' not in policy or len(policy['required_anchors']) > 20:
        raise ValueError('Need explicit joint constraints and at most 20 specified anchors')
    q = copy.deepcopy(next(q for q in evidence['queries'] if
        set(policy['required_anchors']) <= {a['anchor_id'] for a in q['anchors']}))
    q['condition_policy'] = policy
    if selected.get('guided_design'):
        q['guided_design']=selected['guided_design']
        if str(Path(q['guided_design']['receptor_npz']).resolve()) not in selected['sources']:
            raise ValueError('Guided receptor is not sealed by the adopted design')
    indices = sorted({a['feature_index'] for a in q['anchors']})
    for a in q['anchors']:
        a['score_column'] = indices.index(a['feature_index'])
    for p in (selection.parent,source.parent,Path(q['artifact_catalog']).resolve().parent,
              Path(q['chemical_companion']).resolve().parent):
        if output.is_relative_to(p) or p.is_relative_to(output):
            raise ValueError('Output overlaps protected inputs')
    output.mkdir(parents=True, exist_ok=True)
    with file_lock(output.parent/'preselection-family.lock'), file_lock(output/'run.lock'):
        started = time.perf_counter()
        sources = {**evidence['sources'], **full.fingerprint([selection, source,
            Path(q['query_npz']), Path(q['artifact_catalog']), Path(q['chemical_companion'])])}
        full.check_hashes(sources)
        protocol = dict(version=1, sources=sources, code=full.fingerprint(sorted(Path(__file__).parent.glob('*.py'))),
                        workers=workers, chunk_size=chunk_size, policy=policy,
                        seed_parameters=full.POSE_PARAMETERS, numpy=np.__version__, rdkit=rdBase.rdkitVersion,
                        pose_retention='One actual pose per exact anchor mask per molecule; all members retained',
                        clustering='Exact non-stereochemical Murcko scaffold groups; acyclic full-connectivity groups',
                        gaussian='Unchanged best original-seed conformer score, ranking only; no cutoff',
                        reference_acceptance='Not an unpruned full-library equivalence benchmark')
        if q.get('guided_design'):
            protocol['guided_design']=q['guided_design']
            protocol['gaussian']='Best geometry/pocket-surviving seed score, ranking only; no cutoff; not legacy ranking equivalence'
        if q.get('consensus_npz'):
            if str(Path(q['consensus_npz']).resolve()) not in sources:raise ValueError('Consensus coordinates are not sealed')
            protocol['gaussian']='Same-pose mean shape/color Tanimoto (sigma=1, no distance cutoff), weighted optional contacts minus explicit soft exclusions; optional composite cutoff'
        protocol = json.loads(json.dumps(protocol))
        protocol_path = output/'protocol.json'
        if protocol_path.exists() and full.ev.read(protocol_path) != protocol:
            raise ValueError('Changed protocol/code: use a new output directory')
        record_json(protocol_path, protocol)
        report = dict(status='running', stage='integrity_preflight', full_library=True,
                      policy=policy, anchor_order=policy['required_anchors'], counts={},
                      candidate_top_k=None, gaussian_top_n=None, full_coverage=False,
                      scientific_scope='Original bounded seeds; explicit coarse eligibility; no docking or affinity validation')
        record_json(output/'suite-report.json', report)
        try:
            manifests = {}
            for name in ('artifact_catalog','chemical_companion'):
                catalog = Path(q[name])
                for row in full.ev.read(catalog)['shards']:
                    directory = Path(row['path'])
                    if not directory.is_dir(): directory = catalog.parent/row['name']
                    manifest = directory/'manifest.json'
                    if full.ev.sha(manifest) != row['manifest_sha256']:
                        raise ValueError('Shard manifest changed')
                    full.ev.verify_manifest(directory, full=True)
                    manifests.update(full.fingerprint([manifest]))
            report['integrity_preflight_seconds'] = time.perf_counter()-started
            total = full.catalog_count(full.ev.read(q['artifact_catalog']))
            if total != evidence['library']['library_conformers']:
                raise ValueError('Library conformer count differs')
            (output/'chunks').mkdir(exist_ok=True)
            counts, worker_seconds = Counter(), Counter()
            reused = 0
            report.update(stage='full_library_scan', total_conformers=total)
            record_json(output/'suite-report.json', report)
            def update(receipt, cached):
                nonlocal reused
                counts.update(receipt['counts']); worker_seconds.update(receipt['worker_seconds'])
                if cached: reused += receipt['stop']-receipt['start']
                report.update(counts=dict(counts), summed_worker_seconds=dict(worker_seconds),
                              reused_conformers=reused, elapsed_seconds_this_invocation=time.perf_counter()-started)
                record_json(output/'progress.json', report)
            def tasks():
                for start in range(0,total,chunk_size):
                    stop = min(total,start+chunk_size)
                    target = output/'chunks'/f'chunk-{start//chunk_size:08d}.npz'
                    if target.with_suffix('.receipt.json').exists():
                        update(verify_chunk(target,start,stop), True)
                    else:
                        yield start,stop,str(target)
            scan_start = time.perf_counter()
            with ProcessPoolExecutor(max_workers=workers, initializer=full.initialize, initargs=(q,),
                    mp_context=multiprocessing.get_context('spawn')) as pool:
                for _, receipt in bounded_results(pool,compute,tasks(),workers*2):
                    update(receipt,False)
                    full.ev.log(f"Checked {counts['input_conformers']}/{total}; coarse {counts['anchor_bound_passed']}; Gaussian {counts['gaussian_evaluated_conformers']}; matching {counts['matching_conformers']}")
                    if q.get('guided_design'):
                        full.ev.log(f"Guided geometry {counts['guided_geometry_passed']}; pocket conformers {counts['pocket_passed_conformers']}; seeds generated {counts['original_seeds']}, geometry possible {counts['possible_seeds']}, pocket rejected {counts['pocket_rejected_seeds']}, Gaussian evaluated {counts['gaussian_evaluated_seeds']}")
            report['scan_wall_seconds_this_invocation'] = time.perf_counter()-scan_start
            if counts['input_conformers'] != total:
                raise ValueError('Incomplete conformer coverage')
            report.update(stage='molecule_aggregation_and_scaffold_groups')
            record_json(output/'progress.json', report)
            record_json(output/'suite-report.json', report)
            t = time.perf_counter()
            report['aggregation'] = aggregate(output,q,total,chunk_size,evidence['library']['library_molecules'])
            report['aggregation_wall_seconds'] = time.perf_counter()-t
            full.check_hashes(sources); full.check_hashes(manifests); full.check_hashes(protocol['code'])
            previous = total; stages = []
            for level,key in enumerate(('size_passed','coverage_passed','extent_passed','anchor_bound_passed',
                        'pose_feasibility_passed','gaussian_evaluated_conformers','matching_conformers'),1):
                current = counts[key]
                stages.append(dict(stage=key,input_conformers=previous,output_conformers=current,
                                   output_molecules=report['aggregation']['molecules_by_stage'][level],
                                   rejected=previous-current,rejection_fraction=(previous-current)/previous if previous else None))
                previous = current
            report.update(status='complete',stage='finished',full_coverage=True,stage_counts=stages,
                          wall_seconds_this_invocation=time.perf_counter()-started,
                          timing_scope='Includes integrity, scan, persistence, aggregation and grouping; worker sums are not wall time; resumed runs are not fresh timings',
                          outputs={name:str(output/name) for name in ('candidate-poses.jsonl','cluster-members.csv','clusters.csv','candidates.sqlite')})
            report['output_hashes'] = full.fingerprint([output/name for name in report['outputs']])
            record_json(output/'suite-report.json', report); record_json(output/'progress.json', report)
            return report
        except Exception as exc:
            report.update(status='failed',error=str(exc))
            record_json(output/'suite-report.json', report)
            raise


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--selection',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--workers',type=int,default=22)
    parser.add_argument('--chunk-size',type=int,default=2048)
    args = parser.parse_args()
    run(args.selection,args.output,args.workers,args.chunk_size)


if __name__ == '__main__':
    main()
