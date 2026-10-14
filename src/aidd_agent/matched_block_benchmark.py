"""Matched threshold reference, block-prefix replay and fresh measured executions."""
import argparse
from collections import defaultdict
import csv
import hashlib
import json
import math
from pathlib import Path
import platform
import random
import sqlite3
import statistics
import time

from .routing_benchmark import read, require, save, sha, table, verified


SCHEMES = ('E094', 'E095', 'E096')
SCORERS = ('chemplp', 'equiscore')
DEFAULT_CAMPAIGN = '/mnt/local/hand/yuzhang/aidd/e097-chat-workspace/users/workstation/projects/prj-adf8a1b9f4f8-prompt-aidd/runs/PROMPT-5f7c28e4bd704616/execution/campaign/blocks'


def readonly(path):
    return sqlite3.connect(Path(path).resolve().as_uri() + '?mode=ro', uri=True)


def attach(connection, path, alias):
    require(alias.isidentifier(), 'Invalid database alias')
    connection.execute(f'ATTACH DATABASE ? AS {alias}', (Path(path).resolve().as_uri() + '?mode=ro',))


def code_identity():
    return {path.name: sha(path) for path in sorted(Path(__file__).parent.glob('*.py'))}


def seal(folder, report, names):
    report['output_hashes'] = {name: sha(folder / name) for name in names}
    save(folder / 'report.json', report)
    return report


def checked(folder):
    report = read(folder / 'report.json')
    require(report.get('status') == 'complete', f'Incomplete input: {folder}')
    for name, digest in report['output_hashes'].items():
        path = (folder / name).resolve()
        require(path.is_relative_to(folder.resolve()), 'Output seal escapes source folder')
        verified(path, digest)
    return report


def molecule_key(mid, seed):
    return int.from_bytes(hashlib.sha256(f'{seed}/{mid}'.encode()).digest()[:8], 'big') >> 1


def fraction_rule(fraction, alias='member'):
    require(math.isfinite(fraction) and 0 < fraction <= 1, 'Fraction must be in (0,1]')
    return ('1', ()) if fraction == 1 else (f'{alias}.sample_key < ?', (int(fraction * (1 << 63)),))


def prepare(args):
    from .block_sampling import snapshot
    from .block_conformer_search import members, map_artifacts, profile_paths
    started = time.perf_counter()
    campaign = args.campaign.resolve()
    receipt = read(campaign / 'report.json')
    require(receipt.get('status') == 'complete' and receipt.get('search_engine') == 'consensus-threshold-v1',
            'Choose an existing completed threshold campaign')
    for name in ('protocol.json', 'query/report.json'):
        verified(campaign / name, receipt['output_hashes'][name])
    protocol = read(campaign / 'protocol.json')
    for name, digest in protocol['inputs'].items():
        verified(Path(name), digest)
    analysis = Path(protocol['request']['analysis'])
    full = read(Path(read(analysis)['source_report']))
    sample_report = Path(full['sample_report'])
    profiles = []
    for name in protocol['inputs']:
        path = Path(name)
        if path.suffix == '.json':
            value = read(path)
            if isinstance(value, dict) and {'profiles', 'descriptors', *SCHEMES} <= value.keys():
                profiles.append(path)
    require(len(profiles) == 1, 'Cannot uniquely locate the sealed sampling profile')
    paths = profile_paths(profiles[0])
    batch = Path(protocol['library_registry']['path']).parent
    require(snapshot(batch / 'registry.sqlite3') == protocol['library_registry'], 'Source registry changed')
    signature = read(sample_report.parent / 'signature.json')
    recorded = {Path(entry['path']).resolve() for entry in signature['database_snapshots']}
    require({(paths[scheme] / 'property_blocks.sqlite').resolve() for scheme in ('E095', 'E096')} <= recorded,
            'Both child partitions must have frozen sampling snapshots')
    for entry in signature['database_snapshots']:
        require(snapshot(entry['path']) == entry, 'Source partition/profile snapshot changed')
    with readonly(sample_report.parent / 'samples.sqlite') as samples:
        populations = {scheme: dict(samples.execute('SELECT block_id,n FROM population WHERE scheme=?', (scheme,)))
                       for scheme in SCHEMES}
    require(all(populations.values()), 'All three regular partitions are required')
    require(len({sum(values.values()) for values in populations.values()}) == 1, 'Partition totals differ')
    priorities = []
    for scorer, location in [('equiscore', analysis.parent), ('chemplp', analysis.parent / 'chemplp_baseline')]:
        rank_report = read(location / 'report.json')
        verified(location / 'block_ranking.sqlite', rank_report['output_hashes']['block_ranking.sqlite'])
        with readonly(location / 'block_ranking.sqlite') as ranks:
            for scheme in SCHEMES:
                rows = list(ranks.execute('SELECT block_id,rank,population FROM rankings WHERE scheme=? AND receptor=? '
                                         'AND top_n=10 AND status=? ORDER BY rank,block_id',
                                         (scheme, args.receptor, 'eligible')))
                require({block: count for block, _, count in rows} == populations[scheme],
                        f'Missing/ineligible priority blocks or wrong population: {scheme}/{scorer}')
                require([rank for _, rank, _ in rows] == list(range(1, len(rows)+1)), 'Noncontiguous priority ranks')
                priorities.extend(dict(scheme=scheme, scorer=scorer, rank=rank, block_id=block)
                                  for block, rank, _ in rows)
    output = args.output.resolve()
    protected = [batch, campaign, analysis.parent, sample_report.parent, *paths.values()]
    require(all(not output.is_relative_to(path.resolve()) and not path.resolve().is_relative_to(output)
                for path in protected), 'Output overlaps a scientific source directory')
    output.mkdir(parents=True, exist_ok=False)
    save(output / 'report.json', dict(status='incomplete', stage='mapping_common_population'))
    database = output / 'universe.sqlite'
    members(profiles[0], dict(scheme='E094', block_ids=sorted(populations['E094'])), sample_report, database)
    map_artifacts(batch, database, paths['descriptors'])
    with sqlite3.connect(database.as_uri(), uri=True) as connection:
        attach(connection, paths['E095'] / 'property_blocks.sqlite', 'child95')
        attach(connection, paths['E096'] / 'property_blocks.sqlite', 'child96')
        attach(connection, batch / 'registry.sqlite3', 'registry')
        attach(connection, sample_report.parent / 'samples.sqlite', 'training')
        for definition in ('E095 TEXT', 'E096 TEXT', 'canonical_mid TEXT', 'sample_key INTEGER', 'excluded INTEGER'):
            connection.execute('ALTER TABLE member ADD COLUMN ' + definition)
        connection.execute('UPDATE member SET E095=(SELECT block_id FROM child95.membership WHERE cid=member.cid), '
                           'E096=(SELECT block_id FROM child96.membership WHERE cid=member.cid), '
                           'canonical_mid=(SELECT molecule_id FROM registry.conformer WHERE id=member.artifact_cid)')
        require(connection.execute('SELECT count(*) FROM member WHERE E095 IS NULL OR E096 IS NULL '
                                   'OR canonical_mid IS NULL').fetchone()[0] == 0, 'Incomplete common membership/identity')
        connection.execute('CREATE INDEX canonical_mid_index ON member(canonical_mid)')
        connection.execute('CREATE TABLE molecules(mid TEXT PRIMARY KEY,sample_key INTEGER,excluded INTEGER)')
        connection.execute('INSERT INTO molecules(mid,excluded) SELECT canonical_mid,0 FROM member GROUP BY canonical_mid')
        connection.execute('CREATE TABLE training_mid AS SELECT DISTINCT mid FROM training.selected')
        connection.execute('CREATE UNIQUE INDEX training_mid_index ON training_mid(mid)')
        connection.execute('UPDATE molecules SET excluded=1 WHERE mid IN '
                           '(SELECT canonical_mid FROM member JOIN training_mid USING(mid))')
        connection.executemany('UPDATE molecules SET sample_key=? WHERE mid=?',
                               ((molecule_key(mid, args.seed), mid) for mid, in connection.execute('SELECT mid FROM molecules')))
        connection.execute('UPDATE member SET sample_key=(SELECT sample_key FROM molecules WHERE mid=canonical_mid), '
                           'excluded=(SELECT excluded FROM molecules WHERE mid=canonical_mid)')
        for column in ('block_id', 'E095', 'E096', 'sample_key'):
            connection.execute(f'CREATE INDEX benchmark_{column} ON member({column})')
        for scheme, column in [('E094', 'block_id'), ('E095', 'E095'), ('E096', 'E096')]:
            actual = dict(connection.execute(f'SELECT {column},count(*) FROM member GROUP BY {column}'))
            require(actual == populations[scheme], f'Full membership conservation failed: {scheme}')
        total, excluded = connection.execute('SELECT count(*),sum(excluded) FROM member').fetchone()
        molecule_count, training_molecules = connection.execute('SELECT count(*),sum(excluded) FROM molecules').fetchone()
        connection.commit()
    query = read(campaign / 'query/report.json')
    require(query.get('readiness') == 'ready_for_consensus_funnel' and not query.get('failed_template_self_controls'),
            'Query is not ready for the frozen threshold evaluator')
    save(output / 'query.json', query)
    table(output / 'priorities.csv', priorities)
    for entry in signature['database_snapshots']:
        require(snapshot(entry['path']) == entry, 'Source changed during preparation')
    return seal(output, dict(status='complete', kind='matched_benchmark_preparation', batch=str(batch),
                query_sha256=sha(output / 'query.json'), source_campaign=str(campaign),
                source_protocol_sha256=sha(campaign / 'protocol.json'), source_inputs=protocol['inputs'],
                source_snapshots=signature['database_snapshots'], registry_snapshot=protocol['library_registry'],
                seed=args.seed, priority_receptor=args.receptor, population_scope='common_regular_partitions_only',
                regular_conformers=total, regular_molecules=molecule_count, priority_fit_conformers=excluded,
                priority_fit_molecules=training_molecules, heldout_conformers=total-excluded,
                code=code_identity(), preparation_seconds=time.perf_counter()-started,
                limitations=['Special and chemistry-review populations excluded.',
                             'Priorities were fitted previously; query-level holdout is not established.',
                             'Preparation time excludes earlier partition construction and target scoring.']),
                ['universe.sqlite', 'query.json', 'priorities.csv'])


def load_priorities(prepared):
    result = defaultdict(list)
    with (prepared / 'priorities.csv').open(newline='', encoding='utf-8') as stream:
        for row in csv.DictReader(stream):
            result[row['scheme'], row['scorer']].append((int(row['rank']), row['block_id']))
    require(set(result) == {(scheme, scorer) for scheme in SCHEMES for scorer in SCORERS}, 'Incomplete priority matrix')
    for rows in result.values():
        rows.sort()
        require([rank for rank, _ in rows] == list(range(1, len(rows)+1)) and
                len({block for _, block in rows}) == len(rows), 'Invalid priority list')
    return result


def select_members(connection, prepared, fraction, include_training, scheme=None, blocks=None):
    attach(connection, prepared / 'universe.sqlite', 'population')
    condition, values = fraction_rule(fraction)
    connection.execute('CREATE TABLE selected AS SELECT gid,artifact_cid AS cid,canonical_mid AS mid, '
                       'block_id AS E094,E095,E096 FROM population.member WHERE ' + condition +
                       ('' if include_training else ' AND excluded=0'), values)
    if scheme:
        require(scheme in SCHEMES, 'Unknown partition')
        connection.execute('CREATE TEMP TABLE allowed(block TEXT PRIMARY KEY)')
        connection.executemany('INSERT INTO allowed VALUES(?)', ((block,) for block in blocks))
        connection.execute(f'DELETE FROM selected WHERE {scheme} NOT IN (SELECT block FROM allowed)')
    connection.execute('CREATE UNIQUE INDEX selected_gid ON selected(gid)')
    connection.execute('CREATE UNIQUE INDEX selected_cid ON selected(cid)')
    connection.commit()


def aggregate(connection):
    connection.execute('DROP TABLE IF EXISTS main.template_ranks')
    connection.execute('DROP TABLE IF EXISTS main.ranking')
    connection.execute('CREATE TABLE template_ranks AS SELECT gid,template,ROW_NUMBER() OVER '
                       '(PARTITION BY template ORDER BY contact DESC,gaussian DESC,gid) AS rank FROM poses')
    connection.execute('CREATE TABLE ranking AS SELECT ROW_NUMBER() OVER(ORDER BY positive DESC,rrf DESC,gid) AS rank,'
                       'gid,rrf,support FROM (SELECT t.gid,MAX(p.contact>0) AS positive,'
                       'SUM(1.0/(60+t.rank)) AS rrf,COUNT(*) AS support FROM template_ranks t '
                       'JOIN poses p USING(gid,template) GROUP BY t.gid)')
    connection.execute('CREATE UNIQUE INDEX rank_gid ON ranking(gid)')
    connection.commit()


def collect_scores(output, ids, query, chunk_size):
    import numpy as np
    with sqlite3.connect(output / 'scores.sqlite') as connection:
        connection.executescript('DROP TABLE IF EXISTS poses; CREATE TABLE poses('
                                'gid INTEGER,template TEXT,contact REAL,gaussian REAL,payload TEXT,PRIMARY KEY(gid,template));')
        for ordinal, template in enumerate(query['templates']):
            for start in range(0, len(ids), chunk_size):
                path = output / 'chunks' / f'{ordinal:02d}' / f'{start:010d}.npz'
                receipt = read(path.with_suffix('.receipt.json'))
                for name, digest in receipt['files'].items():
                    verified(Path(name), digest)
                expected = ids[start:start+chunk_size]
                with np.load(path, allow_pickle=False) as arrays:
                    require(np.array_equal(arrays['global_ids'], expected), 'Chunk schedule mismatch')
                    passed = set(map(int, arrays['global_ids'][arrays['stage_levels'] == 7]))
                seen = set()
                with path.with_suffix('.poses.jsonl').open(encoding='utf-8') as stream:
                    for line in stream:
                        pose = json.loads(line)
                        gid = pose['global_id']
                        require(type(gid) is int and gid in passed, 'Pose is outside passed schedule')
                        require(all(math.isfinite(pose[key]) for key in ('optional_score', 'gaussian_same_pose', 'composite_score')),
                                'Nonfinite pose scores')
                        seen.add(gid)
                        connection.execute('INSERT INTO poses VALUES(?,?,?,?,?) ON CONFLICT(gid,template) DO UPDATE SET '
                                           'contact=excluded.contact,gaussian=excluded.gaussian,payload=excluded.payload '
                                           'WHERE (excluded.contact,excluded.gaussian)>(poses.contact,poses.gaussian)',
                                           (gid, template['query_id'], pose['optional_score'], pose['gaussian_same_pose'],
                                            json.dumps(pose, sort_keys=True, separators=(',', ':'))))
                require(seen == passed, 'Qualified schedule has missing pose evidence')
                connection.commit()
        aggregate(connection)


def scientific_digest(database):
    digest = hashlib.sha256()
    with readonly(database) as connection:
        for row in connection.execute('SELECT gid,template,contact,gaussian,payload FROM poses ORDER BY gid,template'):
            digest.update(json.dumps(row, separators=(',', ':'), allow_nan=False).encode() + b'\n')
    return digest.hexdigest()


def evaluate(args):
    import numpy as np
    from .budget_screen import run_refinement
    from .screening_selection import check_hashes
    from .block_sampling import snapshot
    from .contact_policy import require_protein_contacts
    started = time.perf_counter()
    prepared = args.prepared.resolve()
    preparation = checked(prepared)
    require(preparation['code'] == code_identity(), 'Code changed after preparation; use the pinned checkout')
    require(snapshot(Path(preparation['batch']) / 'registry.sqlite3') == preparation['registry_snapshot'], 'Library registry changed')
    query = read(prepared / 'query.json')
    check_hashes(query['sources'])
    require_protein_contacts(query['design'], query['anchors'])
    query = dict(query, block_search_mode='threshold', assignment_backend='numpy')
    query.pop('_block_cache', None)
    require(query.get('readiness') == 'ready_for_consensus_funnel', 'Threshold query is not ready')
    priorities = load_priorities(prepared)
    blocks = None
    if args.mode == 'blocks':
        require(args.scheme in SCHEMES and args.scorer in SCORERS, 'Choose a valid partition and scorer')
        ordering = [block for _, block in priorities[args.scheme, args.scorer]]
        count = len(ordering) if args.top_blocks == 'all' else int(args.top_blocks)
        require(0 < count <= len(ordering), 'Invalid block-prefix length')
        blocks = ordering[:count]
    output = args.output.resolve()
    require(not output.is_relative_to(prepared) and not prepared.is_relative_to(output), 'Execution output overlaps preparation')
    settings = dict(prepared=str(prepared), prepared_report_sha256=sha(prepared / 'report.json'),
                    query_sha256=preparation['query_sha256'], code=code_identity(), fraction=args.fraction,
                    include_training=args.include_training, mode=args.mode, scheme=args.scheme, scorer=args.scorer,
                    blocks=blocks, workers=args.workers, chunk_size=args.chunk_size, output_budget=args.output_budget,
                    hardware=args.hardware, repeat=args.repeat)
    resumed = output.exists()
    if resumed:
        require(args.resume and read(output / 'protocol.json') == settings, 'Use fresh output or identical --resume')
        if (output / 'report.json').exists() and read(output / 'report.json').get('status') == 'complete':
            return checked(output)
    else:
        output.mkdir(parents=True)
        save(output / 'protocol.json', settings)
        with sqlite3.connect((output / 'schedule.sqlite').as_uri(), uri=True) as connection:
            select_members(connection, prepared, args.fraction, args.include_training,
                           args.scheme if blocks is not None else None, blocks)
        save(output / 'schedule-seal.json', {'sha256': sha(output / 'schedule.sqlite')})
    verified(output / 'schedule.sqlite', read(output / 'schedule-seal.json')['sha256'])
    with readonly(output / 'schedule.sqlite') as connection:
        ids = np.fromiter((row[0] for row in connection.execute('SELECT gid FROM selected ORDER BY gid')), dtype=np.int64)
        molecule_count = connection.execute('SELECT count(DISTINCT mid) FROM selected').fetchone()[0]
    require(len(ids) > 0, 'No evaluation members at this scale/prefix')
    identity_digest = hashlib.sha256(ids.astype('<i8', copy=False).tobytes()).hexdigest()
    save(output / 'report.json', dict(status='running', stage='threshold_evaluation', conformers=len(ids)))
    prepared_seconds = time.perf_counter()-started
    refinement_started = time.perf_counter()
    run_refinement(Path(preparation['batch']), output, query, ids, args.workers, args.chunk_size)
    refinement_seconds = time.perf_counter()-refinement_started
    aggregation_started = time.perf_counter()
    collect_scores(output, ids, query, args.chunk_size)
    with sqlite3.connect((output / 'scores.sqlite').as_uri(), uri=True) as connection:
        attach(connection, output / 'schedule.sqlite', 'schedule')
        qualified = connection.execute('SELECT count(*) FROM ranking').fetchone()[0]
        with (output / 'candidates.csv').open('w', newline='', encoding='utf-8') as stream:
            writer = csv.writer(stream)
            writer.writerow(['rank', 'global_id', 'conformer_id', 'molecule_id', 'rrf', 'template_support'])
            writer.writerows(connection.execute('SELECT rank,gid,cid,mid,rrf,support FROM ranking '
                                                 'JOIN schedule.selected USING(gid) WHERE rank<=? ORDER BY rank', (args.output_budget,)))
    semantic_hash = scientific_digest(output / 'scores.sqlite')
    aggregation_seconds = time.perf_counter()-aggregation_started
    report = dict(status='complete', kind='matched_threshold_execution', method='full' if blocks is None else
                  f'{args.scheme}/{args.scorer}/Top-{args.top_blocks}', fraction=args.fraction,
                  include_training=args.include_training, input_conformers=len(ids), input_molecules=molecule_count,
                  qualified_conformers=qualified, exported_conformers=min(qualified, args.output_budget),
                  shortfall=max(0, args.output_budget-qualified), schedule_sha256=identity_digest,
                  semantic_pose_sha256=semantic_hash, mode=args.mode, top_blocks=args.top_blocks if blocks is not None else None,
                  prepared_report_sha256=settings['prepared_report_sha256'], query_sha256=preparation['query_sha256'],
                  cache_reused=False, resumed_execution=resumed, timing_eligible=not resumed,
                  preflight_seconds=prepared_seconds, refinement_seconds=refinement_seconds,
                  aggregation_seconds=aggregation_seconds, wall_seconds_before_sealing=time.perf_counter()-started,
                  hardware=args.hardware, machine=platform.platform(), workers=args.workers, chunk_size=args.chunk_size,
                  repeat=args.repeat, output_budget=args.output_budget, policy='current_threshold_contact_first_rrf60',
                  limitations=['Scientific evaluation plus ranking/CSV only; no MOL2 export or docking.',
                               'Prepared block identities/priorities reused; setup and target-scoring costs reported separately.',
                               'Fixed prefix dispatched together; no oracle/adaptive recall stop or result cache.'])
    seal(output, report, ['protocol.json', 'schedule.sqlite', 'scores.sqlite', 'candidates.csv'])
    report['total_wall_seconds'] = time.perf_counter()-started
    save(output / 'report.json', report)
    return report


def subset_ranking(connection, source, gids=None, select_sql=None):
    connection.execute('DROP TABLE IF EXISTS main.poses')
    connection.execute('DROP TABLE IF EXISTS main.keep')
    connection.execute('CREATE TABLE keep(gid INTEGER PRIMARY KEY)')
    if select_sql is not None:
        connection.execute('INSERT INTO keep ' + select_sql)
    else:
        connection.executemany('INSERT INTO keep VALUES(?)', ((int(gid),) for gid in gids))
    connection.execute(f'CREATE TABLE poses AS SELECT p.gid,p.template,p.contact,p.gaussian '
                       f'FROM {source}.poses p JOIN keep USING(gid)')
    aggregate(connection)


def reference_sets(connection, cutoffs):
    ordered = [row[0] for row in connection.execute('SELECT gid FROM ranking ORDER BY rank LIMIT ?', (max(cutoffs),))]
    molecules = [row[0] for row in connection.execute('SELECT mid,min(rank) first_rank FROM ranking '
                 'JOIN scale.selected USING(gid) GROUP BY mid ORDER BY first_rank LIMIT ?', (max(cutoffs),))]
    return {(unit, cutoff): set(values[:cutoff]) for unit, values in [('conformer', ordered), ('molecule', molecules)]
            for cutoff in cutoffs if cutoff <= len(values)}


def route_statistics(connection, references):
    populations = {scheme: dict(connection.execute(f'SELECT {scheme},count(*) FROM scale.selected GROUP BY {scheme}'))
                   for scheme in SCHEMES}
    hits = {scheme: defaultdict(lambda: defaultdict(set)) for scheme in SCHEMES}
    for gid, mid, first, second, third in connection.execute('SELECT gid,mid,E094,E095,E096 FROM ranking JOIN scale.selected USING(gid)'):
        for key, reference in references.items():
            identity = gid if key[0] == 'conformer' else mid
            if identity in reference:
                for scheme, block in zip(SCHEMES, (first, second, third)):
                    hits[scheme][block][key].add(identity)
    return populations, hits


def replay(args):
    prepared, reference = args.prepared.resolve(), args.reference.resolve()
    preparation, receipt = checked(prepared), checked(reference)
    require(receipt['kind'] == 'matched_threshold_execution' and receipt['mode'] == 'full', 'Replay requires a complete matched full reference')
    require(receipt['prepared_report_sha256'] == sha(prepared / 'report.json'), 'Preparation/reference mismatch')
    require(all(fraction <= receipt['fraction'] for fraction in args.fractions), 'Scale exceeds evaluated reference coverage')
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    save(output / 'report.json', dict(status='incomplete'))
    priorities = load_priorities(prepared)
    checkpoints, targets, export_rows, scales, budgets = [], [], [], [], []
    curve_fields = ['fraction', 'scheme', 'scorer', 'seed', 'blocks', 'unit', 'reference_k',
                    'searched', 'work_fraction', 'hits', 'recall']
    with (output / 'routing_curves.csv').open('w', newline='', encoding='utf-8') as stream:
        csv.DictWriter(stream, fieldnames=curve_fields).writeheader()
    for fraction in sorted(set(args.fractions)):
        scale_dir = output / ('scale-' + str(fraction))
        scale_dir.mkdir()
        with sqlite3.connect((scale_dir / 'members.sqlite').as_uri(), uri=True) as database:
            select_members(database, prepared, fraction, receipt['include_training'])
            total, molecule_count = database.execute('SELECT count(*),count(DISTINCT mid) FROM selected').fetchone()
            attach(database, reference / 'schedule.sqlite', 'coverage')
            require(database.execute('SELECT count(*) FROM selected LEFT JOIN coverage.selected c USING(gid) '
                                     'WHERE c.gid IS NULL').fetchone()[0] == 0, 'Reference misses scale members')
        require(total, f'Empty scale: {fraction}')
        with sqlite3.connect((scale_dir / 'reference.sqlite').as_uri(), uri=True) as ranking:
            attach(ranking, reference / 'scores.sqlite', 'source')
            attach(ranking, scale_dir / 'members.sqlite', 'scale')
            subset_ranking(ranking, 'source', select_sql='SELECT gid FROM scale.selected ORDER BY gid')
            references = reference_sets(ranking, args.top_k)
            require(references, f'No requested K fits qualified reference at fraction {fraction}')
            qualified_count = ranking.execute('SELECT count(*) FROM ranking').fetchone()[0]
            populations, block_hits = route_statistics(ranking, references)
            scales.append(dict(fraction=fraction, conformers=total, molecules=molecule_count,
                               qualified=qualified_count, measured_full_runtime_available=fraction == receipt['fraction'] and receipt['timing_eligible']))
        schedules = []
        for scheme in SCHEMES:
            for scorer in SCORERS:
                schedules.append((scheme, scorer, -1, [block for _, block in priorities[scheme, scorer]]))
            for seed in range(args.seed, args.seed+args.random_repeats):
                blocks = [block for _, block in priorities[scheme, 'chemplp']]
                random.Random(f'{seed}/{scheme}').shuffle(blocks)
                schedules.append((scheme, 'random', seed, blocks))
        with sqlite3.connect((scale_dir / 'prefix-work.sqlite').as_uri(), uri=True) as scratch, \
                (output / 'routing_curves.csv').open('a', newline='', encoding='utf-8') as stream:
            attach(scratch, scale_dir / 'reference.sqlite', 'reference')
            attach(scratch, scale_dir / 'members.sqlite', 'scale')
            scratch.execute('CREATE TABLE allowed_blocks(block TEXT PRIMARY KEY)')
            curve_writer = csv.DictWriter(stream, fieldnames=curve_fields)
            for scheme, scorer, seed, ordering in schedules:
                require(set(populations[scheme]) <= set(ordering), 'Present block lacks priority')
                searched, found, curve = 0, {key: set() for key in references}, []
                scratch.execute('DELETE FROM allowed_blocks')
                prefixes = set(args.top_blocks) | {len(ordering)}
                for step, block in enumerate(ordering, 1):
                    searched += populations[scheme].get(block, 0)
                    scratch.execute('INSERT INTO allowed_blocks VALUES(?)', (block,))
                    for key in references:
                        found[key].update(block_hits[scheme].get(block, {}).get(key, set()))
                    for (unit, cutoff), expected in references.items():
                        row = dict(fraction=fraction, scheme=scheme, scorer=scorer, seed=seed, blocks=step,
                                   unit=unit, reference_k=cutoff, searched=searched, work_fraction=searched/total,
                                   hits=len(found[unit, cutoff]), recall=len(found[unit, cutoff])/len(expected))
                        curve_writer.writerow(row)
                        curve.append(row)
                        if step in prefixes:
                            checkpoints.append(row)
                    if scorer != 'random' and step in prefixes and args.rerank_exports:
                        if searched == total:
                            rank_table = 'reference.ranking'
                        else:
                            subset_ranking(scratch, 'reference', select_sql=
                                f'SELECT gid FROM scale.selected WHERE {scheme} IN (SELECT block FROM allowed_blocks) ORDER BY gid')
                            rank_table = 'ranking'
                        delivered_rows = list(scratch.execute(f'SELECT gid,mid FROM {rank_table} JOIN scale.selected USING(gid) '
                                                             'WHERE rank<=? ORDER BY rank', (args.output_budget,)))
                        delivered = {row[0] for row in delivered_rows}
                        for (unit, cutoff), expected in references.items():
                            selected = delivered if unit == 'conformer' else {row[1] for row in delivered_rows}
                            export_rows.append(dict(fraction=fraction, scheme=scheme, scorer=scorer, blocks=step, unit=unit,
                                                    reference_k=cutoff, output_budget=args.output_budget, exported=len(delivered),
                                                    hits=len(selected & expected), recall=len(selected & expected)/len(expected)))
                require(searched == total, 'All-block schedule did not cover the common universe')
                for (unit, cutoff) in references:
                    points = [row for row in curve if row['unit'] == unit and row['reference_k'] == cutoff]
                    for target in args.target_recall:
                        reached = next(row for row in points if row['recall'] >= target)
                        targets.append(dict(fraction=fraction, scheme=scheme, scorer=scorer, seed=seed, unit=unit,
                                            reference_k=cutoff, target_recall=target, blocks=reached['blocks'],
                                            searched=reached['searched'], work_fraction=reached['work_fraction'],
                                            work_reduction_factor=total/reached['searched'], achieved_recall=reached['recall']))
                    for budget_fraction in (0.01, 0.05, 0.1, 0.2, 0.5, 1.0):
                        eligible = [row for row in points if row['searched'] <= math.floor(budget_fraction*total)]
                        last = eligible[-1] if eligible else dict(searched=0, recall=0)
                        budgets.append(dict(fraction=fraction, scheme=scheme, scorer=scorer, seed=seed, unit=unit,
                                            reference_k=cutoff, work_budget_fraction=budget_fraction,
                                            searched=last['searched'], recall=last['recall']))
    table(output / 'scale_summary.csv', scales)
    table(output / 'top_block_checkpoints.csv', checkpoints)
    table(output / 'work_to_recall.csv', targets)
    table(output / 'recall_at_work.csv', budgets)
    if export_rows:
        table(output / 'delivered_recall.csv', export_rows)
    return seal(output, dict(status='complete', kind='matched_routing_replay', reference=str(reference),
                reference_report_sha256=sha(reference / 'report.json'), prepared_report_sha256=sha(prepared / 'report.json'),
                fractions=args.fractions, top_k=args.top_k, top_blocks=args.top_blocks,
                target_recall=args.target_recall, seed=args.seed, random_repeats=args.random_repeats,
                rerank_exports=args.rerank_exports, output_budget=args.output_budget,
                limitations=['Replay measures candidate discovery/work, not runtime.',
                             'Scales are nested molecule-hash subsets of frozen full-library partitions, not independently rebuilt indexes.',
                             'Priorities are fixed; molecules used to fit them are excluded unless explicitly included.',
                             'Work-to-recall prefixes are retrospective operating points, not deployable oracle stopping rules.',
                             'Common original block rank is retained even when a small scale has empty blocks.']),
                ['scale_summary.csv', 'routing_curves.csv', 'top_block_checkpoints.csv', 'work_to_recall.csv', 'recall_at_work.csv'] +
                (['delivered_recall.csv'] if export_rows else []))


def verify_equivalence(args):
    first, second = checked(args.full), checked(args.blocks)
    require(first['mode'] == 'full' and second['mode'] == 'blocks' and second['top_blocks'] == 'all', 'Require full and all-block runs')
    for key in ('prepared_report_sha256', 'query_sha256', 'schedule_sha256', 'semantic_pose_sha256', 'output_budget'):
        require(first[key] == second[key], f'Equivalence failed: {key}')
    require(sha(args.full / 'candidates.csv') == sha(args.blocks / 'candidates.csv'), 'Ordered candidate exports differ')
    args.output.mkdir(parents=True, exist_ok=False)
    result = dict(status='complete', equivalence='passed', full=str(args.full), blocks=str(args.blocks),
                  full_report_sha256=sha(args.full / 'report.json'), blocks_report_sha256=sha(args.blocks / 'report.json'),
                  scope='Same evaluated IDs, best-pose payloads and ordered candidate export; no early stopping')
    save(args.output / 'report.json', result)
    return result


def timing(args):
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    save(output / 'report.json', dict(status='incomplete'))
    records = [(folder, checked(folder)) for folder in args.runs]
    groups = defaultdict(list)
    for folder, report in records:
        require(report.get('kind') == 'matched_threshold_execution' and report['timing_eligible'] and
                not report['cache_reused'] and not report['resumed_execution'], 'Timing contains reused/resumed or invalid execution')
        protocol = read(folder / 'protocol.json')
        require(protocol['code'] == read(records[0][0] / 'protocol.json')['code'], 'Timing code versions differ')
        key = (report['fraction'], report['prepared_report_sha256'], report['query_sha256'], report['include_training'],
               report['hardware'], report['workers'], report['chunk_size'], report['output_budget'])
        groups[key].append((folder, report))
    rows = []
    for key, members in groups.items():
        baselines = {row['repeat']: (folder, row) for folder, row in members if row['mode'] == 'full'}
        require(len(baselines) >= 3 and len(baselines) == sum(row['mode'] == 'full' for _, row in members),
                'Need at least three unique full repeats per condition/scale')
        for repeat, (folder, baseline) in baselines.items():
            require(baseline['semantic_pose_sha256'] == next(iter(baselines.values()))[1]['semantic_pose_sha256'],
                    'Repeated baselines differ scientifically')
        methods = defaultdict(dict)
        for folder, row in members:
            if row['mode'] == 'blocks':
                require(row['repeat'] not in methods[row['method']], 'Duplicate method/repeat')
                methods[row['method']][row['repeat']] = (folder, row)
        for method, repeats in methods.items():
            require(set(repeats) == set(baselines), 'Need matched repeats for each prefix')
            ratios, recalls, output_matches = [], [], []
            for repeat, (folder, row) in repeats.items():
                base_folder, baseline = baselines[repeat]
                with readonly(base_folder / 'scores.sqlite') as database:
                    reference = {gid for gid, in database.execute('SELECT gid FROM ranking WHERE rank<=?', (args.reference_k,))}
                require(len(reference) == args.reference_k, 'Reference K exceeds baseline qualified count')
                with readonly(folder / 'scores.sqlite') as database:
                    delivered = {gid for gid, in database.execute('SELECT gid FROM ranking WHERE rank<=?', (row['output_budget'],))}
                recall = len(reference & delivered)/len(reference)
                recalls.append(recall)
                output_matches.append(row['exported_conformers'] == baseline['exported_conformers'])
                ratios.append(baseline['total_wall_seconds']/row['total_wall_seconds'])
            rows.append(dict(fraction=key[0], method=method, repeats=len(repeats), reference_k=args.reference_k,
                             target_recall=args.target_recall, min_delivered_recall=min(recalls),
                             meets_target=all(value >= args.target_recall for value in recalls),
                             meets_output_budget=all(output_matches),
                             median_time_ratio=statistics.median(ratios),
                             accepted_speedup=statistics.median(ratios) if min(recalls) >= args.target_recall and all(output_matches) else None,
                             median_full_seconds=statistics.median(row['total_wall_seconds'] for _, row in baselines.values()),
                             median_prefix_seconds=statistics.median(row['total_wall_seconds'] for _, row in repeats.values())))
    require(rows, 'No prefix timing methods supplied')
    table(output / 'timing_comparison.csv', rows)
    return seal(output, dict(status='complete', kind='matched_measured_timing_summary',
                inputs={str(folder): sha(folder / 'report.json') for folder, _ in records},
                limitations=['Includes preflight/evaluation/aggregation/sealing, excludes preparation and upstream priority fitting.',
                             'Hardware label is user-declared; OS caches and competing workloads require operator control.',
                             'Accepted speedup requires the declared delivered reference recall, not biological activity.']),
                ['timing_comparison.csv'])


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest='command', required=True)
    prepare_command = commands.add_parser('prepare')
    prepare_command.add_argument('--campaign', type=Path, default=Path(DEFAULT_CAMPAIGN))
    prepare_command.add_argument('--seed', type=int, default=20261014)
    prepare_command.add_argument('--receptor', default='7NA2_A')
    prepare_command.set_defaults(handler=prepare)
    run = commands.add_parser('evaluate')
    run.add_argument('--prepared', type=Path, required=True)
    run.add_argument('--fraction', type=float, default=0.01)
    run.add_argument('--mode', choices=['full', 'blocks'], required=True)
    run.add_argument('--scheme', choices=SCHEMES)
    run.add_argument('--scorer', choices=SCORERS)
    run.add_argument('--top-blocks', default='20')
    run.add_argument('--workers', type=int, default=20)
    run.add_argument('--chunk-size', type=int, default=4048)
    run.add_argument('--output-budget', type=int, default=100000)
    run.add_argument('--hardware', required=True)
    run.add_argument('--repeat', type=int, default=1)
    run.add_argument('--include-training', action='store_true')
    run.add_argument('--resume', action='store_true')
    run.set_defaults(handler=evaluate)
    replay_command = commands.add_parser('replay')
    replay_command.add_argument('--prepared', type=Path, required=True)
    replay_command.add_argument('--reference', type=Path, required=True)
    replay_command.add_argument('--fractions', nargs='+', type=float, default=[0.01])
    replay_command.add_argument('--top-k', nargs='+', type=int, default=[100, 1000, 10000])
    replay_command.add_argument('--top-blocks', nargs='+', type=int, default=[1, 5, 10, 20, 50])
    replay_command.add_argument('--target-recall', nargs='+', type=float, default=[0.9, 0.95, 0.99])
    replay_command.add_argument('--random-repeats', type=int, default=20)
    replay_command.add_argument('--seed', type=int, default=20261014)
    replay_command.add_argument('--rerank-exports', action='store_true')
    replay_command.add_argument('--output-budget', type=int, default=100000)
    replay_command.set_defaults(handler=replay)
    verify = commands.add_parser('verify-equivalence')
    verify.add_argument('--full', type=Path, required=True)
    verify.add_argument('--blocks', type=Path, required=True)
    verify.set_defaults(handler=verify_equivalence)
    timings = commands.add_parser('timing-summary')
    timings.add_argument('--runs', nargs='+', type=Path, required=True)
    timings.add_argument('--reference-k', type=int, default=10000)
    timings.add_argument('--target-recall', type=float, default=0.95)
    timings.set_defaults(handler=timing)
    for command in (prepare_command, run, replay_command, verify, timings):
        command.add_argument('--output', type=Path, required=True)
    args = parser.parse_args(argv)
    if args.command == 'evaluate':
        fraction_rule(args.fraction)
        require(args.workers > 0 and args.chunk_size > 0 and args.output_budget > 0 and args.repeat > 0, 'Invalid execution size')
    if args.command == 'replay':
        for fraction in args.fractions:
            fraction_rule(fraction)
        require(all(value > 0 for value in args.top_k + args.top_blocks) and args.output_budget > 0 and
                args.random_repeats > 0 and all(0 < value <= 1 for value in args.target_recall), 'Invalid replay options')
    if args.command == 'timing-summary':
        require(args.reference_k > 0 and 0 < args.target_recall <= 1, 'Invalid timing endpoint')
    args.handler(args)
    print('Saved', args.command, 'to', args.output)


if __name__ == '__main__':
    main()
