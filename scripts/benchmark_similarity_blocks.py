"""E110 similarity-aligned sampling and held-out yield replay; E109 engine unchanged."""
import argparse
from collections import defaultdict
import csv
import json
from pathlib import Path
import random
import sqlite3
import sys
import time

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'src'))
from aidd_agent import matched_block_benchmark as b


def own_hash():
    return b.sha(Path(__file__))


def fresh(path, sources):
    path = path.resolve()
    b.require(all(not path.is_relative_to(s.resolve()) and not s.resolve().is_relative_to(path)
                  for s in sources), 'Output overlaps an input directory')
    path.mkdir(parents=True, exist_ok=False)
    b.save(path / 'report.json', dict(status='incomplete'))
    return path


def inputs(prepared, reference):
    prep, ref = b.checked(prepared), b.checked(reference)
    b.require(prep['code'] == b.code_identity(), 'Scientific code changed; use the E109 pinned engine')
    b.require(ref['prepared_report_sha256'] == b.sha(prepared / 'report.json'), 'Prepared receipt mismatch')
    b.require(ref['query_sha256'] == prep['query_sha256'], 'Query mismatch')
    b.require(ref['mode'] == 'full' and not ref['include_training'], 'Need a held-out full pilot reference')
    protocol = b.read(reference / 'protocol.json')
    b.require(protocol['code'] == prep['code'], 'Reference evaluator code mismatch')
    return prep, ref


def sample_blocks(db, schemes, size, seed):
    db.execute('CREATE TABLE draws(scheme TEXT,block TEXT,mid TEXT,PRIMARY KEY(scheme,block,mid))')
    db.create_function('sample_hash', 2, lambda block, mid: b.molecule_key(block + '/' + mid, seed))
    for scheme in schemes:
        column = 'block_id' if scheme == 'E094' else scheme
        print('Sampling ' + scheme, flush=True)
        db.execute('CREATE TEMP TABLE pool AS SELECT ' + column + ' AS block,canonical_mid AS mid '
                   'FROM population.member WHERE excluded=0 GROUP BY ' + column + ',canonical_mid')
        counts = list(db.execute('SELECT block,count(*) FROM pool GROUP BY block'))
        b.require(counts and min(n for _, n in counts) >= size,
                  'Every block needs at least the requested distinct eligible molecules')
        db.execute('INSERT INTO draws SELECT ?,block,mid FROM (SELECT block,mid,ROW_NUMBER() OVER '
                   '(PARTITION BY block ORDER BY sample_hash(block,mid),mid) AS position FROM pool) '
                   'WHERE position<=?', (scheme, size))
        db.execute('DROP TABLE pool')
    db.execute('CREATE TABLE fit_mid AS SELECT DISTINCT mid FROM draws')
    db.execute('CREATE UNIQUE INDEX fit_mid_pk ON fit_mid(mid)')
    db.execute('CREATE TABLE selected AS SELECT gid,artifact_cid AS cid,canonical_mid AS mid,'
               'block_id AS E094,E095,E096 FROM population.member '
               'WHERE excluded=0 AND canonical_mid IN (SELECT mid FROM fit_mid)')
    db.execute('CREATE UNIQUE INDEX selected_gid ON selected(gid)')
    db.execute('CREATE INDEX selected_mid ON selected(mid)')
    db.commit()


def plan(args):
    started = time.perf_counter()
    prepared, reference = args.prepared.resolve(), args.reference.resolve()
    prep, ref = inputs(prepared, reference)
    out = fresh(args.output, [prepared, reference])
    with sqlite3.connect((out / 'schedule.sqlite').as_uri(), uri=True) as db:
        b.attach(db, prepared / 'universe.sqlite', 'population')
        b.attach(db, reference / 'schedule.sqlite', 'pilot')
        b.attach(db, reference / 'scores.sqlite', 'scores')
        sample_blocks(db, args.schemes, args.sample_size, args.seed)
        b.require(db.execute('SELECT count(*) FROM pilot.selected s LEFT JOIN population.member p USING(gid) '
                             'WHERE p.gid IS NULL OR s.mid!=p.canonical_mid OR p.excluded!=0').fetchone()[0] == 0,
                  'Pilot identities do not match prepared population')
        db.execute('CREATE TABLE missing AS SELECT s.gid FROM selected s '
                   'WHERE NOT EXISTS(SELECT 1 FROM pilot.selected p WHERE p.gid=s.gid)')
        db.execute('CREATE TABLE positive(gid INTEGER PRIMARY KEY)')
        db.execute('INSERT INTO positive SELECT DISTINCT p.gid FROM scores.poses p JOIN selected s USING(gid)')
        db.execute('CREATE TABLE validation AS SELECT s.*,EXISTS(SELECT 1 FROM scores.poses p WHERE p.gid=s.gid) '
                   'AS positive FROM pilot.selected s WHERE s.mid NOT IN (SELECT mid FROM fit_mid)')
        db.execute('CREATE UNIQUE INDEX validation_gid ON validation(gid)')
        db.commit()
        n = db.execute('SELECT count(*) FROM selected').fetchone()[0]
        missing = db.execute('SELECT count(*) FROM missing').fetchone()[0]
        validation = db.execute('SELECT count(*) FROM validation').fetchone()[0]
        b.require(validation > 0, 'No validation molecules remain after excluding all fitting molecules')
        counts = [dict(scheme=s, blocks=k, sampled_molecule_memberships=k*args.sample_size)
                  for s, k in db.execute('SELECT scheme,count(DISTINCT block) FROM draws GROUP BY scheme')]
        mids = db.execute('SELECT count(*) FROM fit_mid').fetchone()[0]
    settings = dict(prepared=str(prepared), reference=str(reference),
                    prepared_report_sha256=b.sha(prepared / 'report.json'),
                    reference_report_sha256=b.sha(reference / 'report.json'), script_sha256=own_hash(),
                    scientific_code=prep['code'], schemes=args.schemes, sample_size=args.sample_size, seed=args.seed)
    b.save(out / 'protocol.json', settings)
    return b.seal(out, dict(status='complete', kind='similarity_sampling_plan', **settings,
                 blocks=counts, fit_molecules=mids, fit_conformers=n, reused_conformers=n-missing,
                 new_conformers=missing, validation_conformers=validation,
                 estimated_new_refinement_hours=missing/ref['input_conformers']*ref['refinement_seconds']/3600,
                 plan_seconds=time.perf_counter()-started,
                 limitations=['Runtime estimate assumes pilot throughput; not measured.',
                              'Same-pilot follow-up is exploratory; validation molecules are excluded from fitting.']),
                 ['protocol.json', 'schedule.sqlite'])


def fit(args):
    started = time.perf_counter()
    plan_dir = args.plan.resolve()
    planned = b.checked(plan_dir)
    b.require(planned['script_sha256'] == own_hash(), 'Benchmark script changed since planning')
    prepared, reference = Path(planned['prepared']), Path(planned['reference'])
    prep, ref = inputs(prepared, reference)
    b.require(planned['prepared_report_sha256'] == b.sha(prepared / 'report.json') and
              planned['reference_report_sha256'] == b.sha(reference / 'report.json'), 'Plan inputs changed')
    settings = dict(plan_report_sha256=b.sha(plan_dir / 'report.json'), workers=args.workers,
                    chunk_size=args.chunk_size, hardware=args.hardware, script_sha256=own_hash())
    out = args.output.resolve()
    resumed = out.exists()
    if resumed:
        b.require(args.resume and b.read(out / 'protocol.json') == settings, 'Use fresh output or identical --resume')
        if b.read(out / 'report.json').get('status') == 'complete':
            return b.checked(out)
    else:
        fresh(out, [plan_dir, prepared, reference])
        b.save(out / 'protocol.json', settings)
    from aidd_agent.block_sampling import snapshot
    from aidd_agent.screening_selection import check_hashes
    from aidd_agent.contact_policy import require_protein_contacts
    from aidd_agent.budget_screen import run_refinement
    import numpy as np
    b.require(snapshot(Path(prep['batch']) / 'registry.sqlite3') == prep['registry_snapshot'], 'Library registry changed')
    query = b.read(prepared / 'query.json')
    check_hashes(query['sources'])
    require_protein_contacts(query['design'], query['anchors'])
    b.require(query.get('readiness') == 'ready_for_consensus_funnel', 'Query is not ready')
    query = dict(query, block_search_mode='threshold', assignment_backend='numpy')
    query.pop('_block_cache', None)
    with b.readonly(plan_dir / 'schedule.sqlite') as db:
        ids = np.fromiter((r[0] for r in db.execute('SELECT gid FROM missing ORDER BY gid')), dtype=np.int64)
    b.save(out / 'report.json', dict(status='running', stage='missing_sample_refinement', conformers=len(ids)))
    print('New fitting conformers: ' + str(len(ids)), flush=True)
    tick = time.perf_counter()
    if len(ids):
        run_refinement(Path(prep['batch']), out, query, ids, args.workers, args.chunk_size)
        b.collect_scores(out, ids, query, args.chunk_size)
    else:
        with sqlite3.connect(out / 'scores.sqlite') as db:
            db.execute('CREATE TABLE IF NOT EXISTS poses(gid INTEGER,template TEXT,contact REAL,gaussian REAL,payload TEXT)')
    refinement_seconds = time.perf_counter()-tick
    rows = []
    with sqlite3.connect((out / 'fit.sqlite').as_uri(), uri=True) as db:
        b.attach(db, plan_dir / 'schedule.sqlite', 'plan')
        b.attach(db, out / 'scores.sqlite', 'scores')
        db.execute('DROP TABLE IF EXISTS main.positive')
        db.execute('CREATE TABLE positive(gid INTEGER PRIMARY KEY)')
        db.execute('INSERT INTO positive SELECT gid FROM plan.positive')
        db.execute('INSERT OR IGNORE INTO positive SELECT DISTINCT gid FROM scores.poses')
        for scheme in planned['schemes']:
            counts = db.execute('SELECT d.block,count(DISTINCT d.mid),count(DISTINCT CASE WHEN p.gid IS NOT NULL '
                'THEN d.mid END) FROM plan.draws d JOIN plan.selected s ON s.mid=d.mid AND s.' + scheme +
                '=d.block LEFT JOIN positive p USING(gid) WHERE d.scheme=? GROUP BY d.block', (scheme,)).fetchall()
            b.require(all(n == planned['sample_size'] for _, n, _ in counts), 'Incomplete fitting membership')
            for rank, (block, n, positives) in enumerate(sorted(counts, key=lambda r: (-r[2], r[0])), 1):
                rows.append(dict(scheme=scheme, rank=rank, block_id=block, sampled_molecules=n,
                                 qualified_molecules=positives, qualified_fraction=positives/n))
        db.commit()
    b.table(out / 'priorities.csv', rows)
    return b.seal(out, dict(status='complete', kind='similarity_fit', **settings,
                 fit_conformers=planned['fit_conformers'], new_conformers=len(ids),
                 reused_conformers=planned['reused_conformers'], resumed_execution=resumed,
                 refinement_and_collection_seconds=refinement_seconds, wall_seconds=time.perf_counter()-started,
                 timing_eligible=False, limitations=['Cache-assisted fitting; not a fresh latency comparison.',
                 'Qualification is stage 7 on any template, counted once per in-block molecule.']),
                 ['protocol.json', 'scores.sqlite', 'fit.sqlite', 'priorities.csv'])


def curves(order, populations, positives):
    searched, found = 0, set()
    result = [(0, 0, 0)]
    for index, block in enumerate(order, 1):
        searched += populations.get(block, 0)
        found.update(positives.get(block, ()))
        result.append((index, searched, len(found)))
    return result


def replay(args):
    started = time.perf_counter()
    plan_dir, fit_dir = args.plan.resolve(), args.fit.resolve()
    planned, fitted = b.checked(plan_dir), b.checked(fit_dir)
    b.require(planned['script_sha256'] == own_hash() == fitted['script_sha256'], 'Benchmark script changed')
    b.require(fitted['plan_report_sha256'] == b.sha(plan_dir / 'report.json'), 'Fit/plan mismatch')
    prepared = Path(planned['prepared'])
    b.require(b.sha(prepared / 'report.json') == planned['prepared_report_sha256'], 'Preparation changed')
    b.verified(prepared / 'priorities.csv', b.read(prepared / 'report.json')['output_hashes']['priorities.csv'])
    old = b.load_priorities(prepared)
    with (fit_dir / 'priorities.csv').open(newline='') as stream:
        fitting = list(csv.DictReader(stream))
    out = fresh(args.output, [plan_dir, fit_dir, prepared, Path(planned['reference'])])
    all_curves, budgets, targets = [], [], []
    with b.readonly(plan_dir / 'schedule.sqlite') as db:
        total = db.execute('SELECT count(*) FROM validation').fetchone()[0]
        total_hits = db.execute('SELECT count(DISTINCT mid) FROM validation WHERE positive=1').fetchone()[0]
        for scheme in planned['schemes']:
            pops = dict(db.execute('SELECT ' + scheme + ',count(*) FROM validation GROUP BY ' + scheme))
            hits = defaultdict(set)
            for block, mid in db.execute('SELECT ' + scheme + ',mid FROM validation WHERE positive=1'):
                hits[block].add(mid)
            ordered = [r['block_id'] for r in sorted((r for r in fitting if r['scheme']==scheme), key=lambda r:int(r['rank']))]
            b.require(set(pops) <= set(ordered), 'Validation block has no fitted priority')
            methods = [('similarity', -1, ordered)]
            for scorer in b.SCORERS:
                order = [block for _, block in old[scheme, scorer]]
                b.require(set(order) == set(ordered), 'Historical and similarity block sets differ')
                methods.append((scorer, -1, order))
            for repeat in range(args.random_repeats):
                order = sorted(ordered)
                random.Random(args.seed + repeat).shuffle(order)
                methods.append(('random', args.seed + repeat, order))
            for method, seed, order in methods:
                route = curves(order, pops, hits)
                overhead = planned['fit_conformers'] if method == 'similarity' else 0
                common = dict(scheme=scheme, method=method, seed=seed)
                for blocks, work, count in route:
                    all_curves.append(dict(**common, blocks=blocks, conformers=work,
                                           qualified_molecules=count, work_fraction=work/total,
                                           cold_conformer_work=work+overhead))
                for fraction in args.work_fractions:
                    blocks, work, count = [r for r in route if r[1] <= int(total*fraction)][-1]
                    budgets.append(dict(**common, budget_fraction=fraction, blocks=blocks, conformers=work,
                                        qualified_molecules=count, recall_all_qualified=count/total_hits if total_hits else 0,
                                        cold_conformer_work=work+overhead))
                for target in args.targets:
                    match = next((r for r in route if r[2] >= target), None)
                    targets.append(dict(**common, target=target, reached=match is not None,
                                        blocks=match[0] if match else '', conformers=match[1] if match else '',
                                        qualified_molecules=match[2] if match else '',
                                        cold_conformer_work=match[1]+overhead if match else ''))
    b.table(out / 'yield_curves.csv', all_curves)
    b.table(out / 'yield_at_work.csv', budgets)
    b.table(out / 'work_to_candidates.csv', targets)
    report = dict(status='complete', kind='similarity_yield_replay', script_sha256=own_hash(),
                  plan_report_sha256=b.sha(plan_dir / 'report.json'), fit_report_sha256=b.sha(fit_dir / 'report.json'),
                  validation_conformers=total, validation_qualified_molecules=total_hits,
                  fit_conformers=planned['fit_conformers'], fit_new_conformers=planned['new_conformers'],
                  fit_reused_conformers=planned['reused_conformers'], fit_wall_seconds=fitted['wall_seconds'],
                  random_repeats=args.random_repeats, seed=args.seed, work_fractions=args.work_fractions,
                  targets=args.targets, replay_seconds=time.perf_counter()-started,
                  limitations=['Exploratory pilot validation; no full-library result or measured search speedup.',
                      'Cold work charges all logical fit conformers, even reused ones; it is not wall time.',
                      'Multi-scheme runs charge the shared fitting union to each similarity arm conservatively.',
                      'Historical docking priority construction costs are excluded; random has no priority fitting cost.',
                      'Targets count validation positives only; sample positives do not enter validation.',
                      'Whole blocks are atomic; unreached targets remain explicitly marked.'])
    lines = ['# E110 similarity-aligned yield replay', '',
             'Exploratory held-out pilot result. Work is conformer evaluations, not measured latency.', '',
             '| Scheme | Method | Yield at 20% work |', '|---|---|---:|']
    for scheme in planned['schemes']:
        for method in ('similarity', 'chemplp', 'equiscore', 'random'):
            values = [r['qualified_molecules'] for r in budgets if r['scheme']==scheme and r['method']==method and r['budget_fraction']==.2]
            if values:
                lines.append('| ' + scheme + ' | ' + method + ' | ' + str(round(sum(values)/len(values), 2)) + ' |')
    lines += ['', 'Logical fitting conformers: ' + str(planned['fit_conformers']) + '.',
              'New fitting conformers: ' + str(planned['new_conformers']) + '.', '',
              'See work_to_candidates.csv for reached/shortfall status and cold work including fitting.',
              'Random rows are 20 orders by default; the table uses their mean. No speedup claim is made.']
    (out / 'summary.md').write_text('\n'.join(lines)+'\n', encoding='utf-8')
    return b.seal(out, report, ['yield_curves.csv', 'yield_at_work.csv', 'work_to_candidates.csv', 'summary.md'])


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='stage', required=True)
    cmd = sub.add_parser('plan')
    cmd.add_argument('--prepared', type=Path, required=True)
    cmd.add_argument('--reference', type=Path, required=True)
    cmd.add_argument('--schemes', nargs='+', choices=b.SCHEMES, default=['E094'])
    cmd.add_argument('--sample-size', type=int, default=500)
    cmd.add_argument('--seed', type=int, default=20261015)
    cmd.set_defaults(func=plan)
    cmd = sub.add_parser('fit')
    cmd.add_argument('--plan', type=Path, required=True)
    cmd.add_argument('--workers', type=int, default=20)
    cmd.add_argument('--chunk-size', type=int, default=4048)
    cmd.add_argument('--hardware', required=True)
    cmd.add_argument('--resume', action='store_true')
    cmd.set_defaults(func=fit)
    cmd = sub.add_parser('replay')
    cmd.add_argument('--plan', type=Path, required=True)
    cmd.add_argument('--fit', type=Path, required=True)
    cmd.add_argument('--random-repeats', type=int, default=20)
    cmd.add_argument('--seed', type=int, default=20261015)
    cmd.add_argument('--work-fractions', nargs='+', type=float, default=[.1,.2,.5,1.])
    cmd.add_argument('--targets', nargs='+', type=int, default=[1000,10000,50000])
    cmd.set_defaults(func=replay)
    for cmd in sub.choices.values():
        cmd.add_argument('--output', type=Path, required=True)
    args = parser.parse_args(argv)
    for key in ('sample_size', 'workers', 'chunk_size', 'random_repeats'):
        if hasattr(args, key):
            b.require(getattr(args, key) > 0, key + ' must be positive')
    if hasattr(args, 'schemes'):
        b.require(len(args.schemes)==len(set(args.schemes)), 'Duplicate schemes')
    if hasattr(args, 'work_fractions'):
        b.require(all(0 < x <= 1 for x in args.work_fractions), 'Invalid work fractions')
        b.require(all(x > 0 for x in args.targets), 'Targets must be positive')
    print(json.dumps(args.func(args), indent=2))


if __name__ == '__main__':
    main()
