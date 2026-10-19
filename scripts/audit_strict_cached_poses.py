"""E111: bounded, read-only strict review of saved E109 pose representatives."""
import argparse
import csv
from contextlib import closing
import hashlib
import json
from pathlib import Path
import random
import sqlite3
import sys
import time

import numpy as np
from scipy.spatial import cKDTree

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'src'))
ILE61 = 'pocket:74a6061918ece571'
GATES = ['contacts', 'shape', 'color', 'steric', 'envelope', 'regions']
POLICY = dict(anchor_threshold=0.5, minimum_residues=3, shape_minimum=0.50,
              color_minimum=0.30, minimum_distance=2.0, envelope_fraction=0.65,
              region_fraction=0.50, coverage_distance=1.5, receptor_contact_distance=4.5,
              region_count=3, minimum_region_atoms=3, minimum_seed_separation=2.0)


def read(path):
    return json.loads(Path(path).read_text(encoding='utf-8'))


def save(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False), encoding='utf-8')


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for part in iter(lambda: stream.read(1024*1024), b''):
            h.update(part)
    return h.hexdigest()


def require(condition, message):
    if not condition:
        raise ValueError(message)


def verify(path, digest):
    require(sha(path) == digest, 'Hash mismatch: ' + str(path))


def ro(path):
    return closing(sqlite3.connect(Path(path).resolve().as_uri() + '?mode=ro', uri=True))


def points(value):
    a = np.asarray(value, dtype=float)
    require(a.ndim == 2 and a.shape[1] == 3 and len(a) > 0 and np.isfinite(a).all(),
            'Expected finite, nonempty Nx3 coordinates')
    return a


def contact_gate(scores, assignments, order, anchors):
    require(len(scores) == len(assignments) == len(order), 'Anchor vector length mismatch')
    require(np.isfinite(scores).all(), 'Nonfinite anchor score')
    lookup = {a['anchor_id']: a for a in anchors}
    hits = [key for key, score, assignment in zip(order, scores, assignments)
            if assignment >= 0 and score >= POLICY['anchor_threshold']]
    residues = sorted({int(lookup[key]['target_residue']) for key in hits})
    return ILE61 in hits and len(residues) >= 3, residues


def make_regions(reference, receptor_tree):
    reference = points(reference)
    contact = reference[receptor_tree.query(reference)[0] <= POLICY['receptor_contact_distance']]
    if len(contact) < 9:
        return contact, None
    center = contact.mean(axis=0)
    seeds = [int(np.argmax(np.sum((contact-center)**2, axis=1)))]
    for _ in range(2):
        nearest = np.min(np.sum((contact[:, None]-contact[seeds][None, :])**2, axis=2), axis=1)
        if nearest.max() < POLICY['minimum_seed_separation']**2:
            return contact, None
        seeds.append(int(np.argmax(nearest)))
    labels = np.argmin(np.sum((contact[:, None]-contact[seeds][None, :])**2, axis=2), axis=1)
    regions = [contact[labels == i] for i in range(3)]
    return contact, regions if all(len(r) >= 3 for r in regions) else None


def geometry_metrics(moved, receptor_tree, contact, regions):
    moved = points(moved)
    minimum = float(receptor_tree.query(moved)[0].min())
    tree = cKDTree(moved)
    coverage = lambda p: float(np.mean(tree.query(p)[0] <= POLICY['coverage_distance'])) if len(p) else 0.0
    fractions = [coverage(r) for r in regions] if regions is not None else []
    return dict(minimum_distance=minimum, envelope_coverage=coverage(contact),
                region_coverages=fractions, regions_supported=regions is not None)


def gates(contact_ok, shape, color, metrics):
    require(np.isfinite([shape, color]).all(), 'Nonfinite shape/color score')
    return dict(contacts=bool(contact_ok), shape=bool(shape >= POLICY['shape_minimum']),
                color=bool(color >= POLICY['color_minimum']),
                steric=bool(metrics['minimum_distance'] >= POLICY['minimum_distance']),
                envelope=bool(metrics['envelope_coverage'] >= POLICY['envelope_fraction']),
                regions=bool(metrics['regions_supported'] and
                             min(metrics['region_coverages']) >= POLICY['region_fraction']))


def accepted(flags):
    return all(flags[key] for key in GATES)


def write_csv(path, rows):
    if not rows:
        return
    with Path(path).open('w', newline='', encoding='utf-8') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def audit(args):
    start = time.monotonic()
    prepared, reference, output = args.prepared.resolve(), args.reference.resolve(), args.output.resolve()
    require(args.sample > 0 and args.seconds > 0, 'Positive sample size and time budget required')
    require(all(not output.is_relative_to(p) and not p.is_relative_to(output)
                for p in [prepared, reference]), 'Output overlaps inputs')
    output.mkdir(parents=True, exist_ok=False)
    save(output/'report.json', dict(status='incomplete'))
    prep, ref = read(prepared/'report.json'), read(reference/'report.json')
    require(prep['status'] == ref['status'] == 'complete', 'Incomplete input')
    require(ref['prepared_report_sha256'] == sha(prepared/'report.json'), 'Reference/prepared mismatch')
    verify(prepared/'query.json', prep['query_sha256'])
    require(ref['query_sha256'] == prep['query_sha256'], 'Query mismatch')
    require(ref['mode'] == 'full' and not ref['include_training'], 'Need held-out full reference')
    for name in ['schedule.sqlite', 'scores.sqlite']:
        verify(reference/name, ref['output_hashes'][name])
    query = read(prepared/'query.json')
    from aidd_agent.contact_groups import selected_anchors
    from aidd_agent.gaussian_batch import ArtifactCatalogReader, _load_query
    from aidd_agent.gaussian_overlay import score_overlay, apply_transform
    from aidd_agent.block_sampling import snapshot
    batch = Path(prep['batch'])
    require(snapshot(batch/'registry.sqlite3') == prep['registry_snapshot'], 'Registry changed')
    for name in ['gaussian_overlay.py', 'gaussian_batch.py', 'contact_groups.py']:
        raw = (Path(__file__).resolve().parents[1]/'src'/'aidd_agent'/name).read_bytes()
        normalized = raw.replace(b'\r\n', b'\n')
        require(prep['code'][name] in [hashlib.sha256(raw).hexdigest(), hashlib.sha256(normalized).hexdigest()],
                'Pinned scientific helper changed: ' + name)
    order = query.get('anchor_order', selected_anchors(query['design']))
    require(order == selected_anchors(query['design']), 'Unexpected anchor order')
    require(ILE61 in order, 'Missing specified Ile61 anchor')
    resources = [Path(query['design']['receptor_npz'])] + [Path(t['query_npz']) for t in query['templates']]
    for resource in resources:
        digest = query['sources'].get(str(resource)) or query['sources'].get(str(resource.resolve()))
        require(digest is not None, 'Resource lacks source seal: ' + str(resource))
        verify(resource, digest)
    with np.load(resources[0], allow_pickle=False) as z:
        receptor_tree = cKDTree(points(z['points']))
    templates = {}
    region_report = []
    for t in query['templates']:
        data = _load_query(Path(t['query_npz']))[1]
        contact, regions = make_regions(data['shape_points'], receptor_tree)
        templates[t['query_id']] = data, contact, regions
        region_report.append(dict(template=t['query_id'], contact_atoms=len(contact),
                                  region_atoms=[len(r) for r in regions] if regions else [],
                                  region_points=[r.tolist() for r in regions] if regions else [],
                                  supported=regions is not None))
    save(output/'regions.json', region_report)
    settings = dict(policy=POLICY, sample=args.sample, seed=args.seed, seconds=args.seconds,
                    prepared=str(prepared), reference=str(reference),
                    query_sha256=prep['query_sha256'], script_sha256=sha(__file__),
                    prepared_report_sha256=sha(prepared/'report.json'),
                    reference_report_sha256=sha(reference/'report.json'))
    save(output/'protocol.json', settings)
    controls = []
    for control in query['crystal_self_controls']:
        if control['query_id'] not in templates:
            continue
        data, contact, regions = templates[control['query_id']]
        ok, residues = contact_gate(control['scores'], [0]*len(order), order, query['anchors'])
        metric = geometry_metrics(data['shape_points'], receptor_tree, contact, regions)
        flags = gates(ok, 1.0, 1.0, metric)
        controls.append(dict(query_id=control['query_id'], residues=residues, **metric, **flags,
                             accepted=accepted(flags), scope='identity pose; self shape/color equal one by definition'))
    save(output/'reference_controls.json', controls)
    setup_seconds = time.monotonic()-start
    tick = time.monotonic()
    rng = random.Random(args.seed)
    selected = []
    with ro(reference/'schedule.sqlite') as schedule:
        for i, row in enumerate(schedule.execute('SELECT gid,cid,mid FROM selected ORDER BY gid')):
            if i < args.sample:
                selected.append(row)
            else:
                j = rng.randrange(i+1)
                if j < args.sample:
                    selected[j] = row
            if i % 10000 == 0 and time.monotonic()-tick > args.seconds:
                raise TimeoutError('Sampling exceeded budget; no pass rate is claimed')
    rng.shuffle(selected)
    write_csv(output/'sample.csv', [dict(gid=g, cid=c, mid=m) for g, c, m in selected])
    reader = ArtifactCatalogReader(batch/'artifacts/catalog.json')
    results, pose_rows = [], []
    with ro(reference/'scores.sqlite') as db:
        for gid, cid, mid in selected:
            if time.monotonic()-tick > args.seconds:
                break
            saved = list(db.execute('SELECT template,payload FROM poses WHERE gid=? ORDER BY template', (gid,)))
            incremental = {k: False for k in GATES}
            standalone = {k: False for k in GATES}
            if saved:
                candidate = reader.get(gid)
                require(str(candidate.conformer_id) == str(cid), 'Conformer identity mismatch')
                for template, payload in saved:
                    pose = json.loads(payload)
                    require(pose['global_id'] == gid and str(pose['conformer_id']) == str(cid), 'Pose identity mismatch')
                    require(str(pose['molecule_id']) == str(candidate.molecule_id), 'Pose molecule identity mismatch')
                    transform = np.asarray(pose['transform'], dtype=float).reshape(4, 4)
                    require(np.isfinite(transform).all() and np.allclose(transform[3], [0,0,0,1]) and
                            np.allclose(transform[:3,:3] @ transform[:3,:3].T, np.eye(3), atol=1e-4) and
                            np.linalg.det(transform[:3,:3]) > 0, 'Invalid rigid transform')
                    data, contact, regions = templates[template]
                    ok, residues = contact_gate(pose['scores'], pose['assignments'], order, query['anchors'])
                    score = score_overlay(data['shape_points'], candidate.shape_points,
                                          data['feature_points'], data['feature_types'],
                                          candidate.feature_points, candidate.feature_types,
                                          transform.ravel(), cutoff=None)
                    shape, color = score['shape']['tanimoto'], score['color']['tanimoto']
                    require(abs(0.5*(shape+color)-pose['gaussian_same_pose']) < 1e-4,
                            'Gaussian decomposition differs from stored pose')
                    metric = geometry_metrics(apply_transform(candidate.shape_points, transform.ravel()),
                                              receptor_tree, contact, regions)
                    flags = gates(ok, shape, color, metric)
                    prefix = True
                    for key in GATES:
                        standalone[key] |= flags[key]
                        prefix = prefix and flags[key]
                        incremental[key] |= prefix
                    pose_rows.append(dict(gid=gid, mid=mid, template=template, residues=json.dumps(residues),
                                          shape_tanimoto=shape, color_tanimoto=color, minimum_distance=metric['minimum_distance'],
                                          envelope_coverage=metric['envelope_coverage'],
                                          region_coverages=json.dumps(metric['region_coverages']),
                                          **flags, accepted=accepted(flags)))
            results.append(dict(gid=gid, mid=mid, old_pass=bool(saved), **incremental,
                                **{'standalone_'+k:v for k,v in standalone.items()}))
            if len(results) % 100 == 0:
                print('Audited ' + str(len(results)) + '/' + str(len(selected)), flush=True)
    if hasattr(reader, 'close'):
        reader.close()
    write_csv(output/'conformers.csv', results)
    write_csv(output/'poses.csv', pose_rows)
    complete = len(results) == len(selected)
    count = len(results)
    by_molecule = {}
    for row in results:
        by_molecule.setdefault(row['mid'], []).append(row)
    report = dict(status='complete' if complete else 'time_budget_exhausted', kind='strict_cached_pose_audit',
                  requested_conformers=args.sample, selected_conformers=len(selected), processed_conformers=count,
                  old_pass_conformers=sum(r['old_pass'] for r in results),
                  cumulative_conformer_counts={k:sum(r[k] for r in results) for k in GATES},
                  standalone_conformer_counts={k:sum(r['standalone_'+k] for r in results) for k in GATES},
                  strict_pass_fraction=sum(r['regions'] for r in results)/count if complete and count else None,
                  sampled_molecules=len(by_molecule),
                  sampled_molecules_with_surviving_sampled_conformer=sum(any(r['regions'] for r in rs) for rs in by_molecule.values()),
                  selected_template_controls=len(controls), template_controls_passed=sum(c['accepted'] for c in controls),
                  setup_seconds=setup_seconds, audit_seconds=time.monotonic()-tick, total_seconds=time.monotonic()-start,
                  policy=POLICY, limitations=[
                      'Exploratory heuristic thresholds; no activity or production acceptance claim.',
                      'Saved best-per-template poses only; alternative valid poses may have been discarded.',
                      'Conformers sampled uniformly; molecule output is not full molecule-level pass rate.',
                      'Reference contact-envelope groups are geometric proxies, not validated canonical subpockets.',
                      'Steric rule uses untyped distance, not atom-specific VDW or energy.',
                      'No new alignment, docking, or full-library run.',
                      'Partial time-limited output has no overall pass-rate estimate.'])
    report['output_hashes'] = {p.name:sha(p) for p in output.iterdir() if p.name != 'report.json' and p.is_file()}
    save(output/'report.json', report)
    print(json.dumps(report, indent=2))
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--prepared', type=Path, required=True)
    parser.add_argument('--reference', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--sample', type=int, default=2000)
    parser.add_argument('--seed', type=int, default=20261019)
    parser.add_argument('--seconds', type=float, default=600)
    args = parser.parse_args()
    audit(args)


if __name__ == '__main__':
    main()
