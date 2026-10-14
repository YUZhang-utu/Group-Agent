"""Read-only historical recovery, complete-universe routing replay and timing receipts."""
import argparse
import csv
import hashlib
import html
import itertools
import json
import math
from pathlib import Path
import random
import sqlite3
import statistics
import subprocess
import time


RUNS = ('5f7c28e4bd704616 3a90516120f74140 4f0fad82c9714b1e '
        '34bf4cc4a192449e 7a7fddb92b744a38 9bd13fe58ad14df0 '
        'dcc64417cf534000 317795cf58ad4f04 803b293dd57c4aec '
        '82f400cb6c1c46df c71019343cb04b62 bfa50ed8eb9a489e').split()
SCHEMES = ('E094', 'E095', 'E096')
METHODS = ('chemplp', 'equiscore')


def require(condition, message):
    if not condition:
        raise ValueError(message)


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(chunk)
    return digest.hexdigest()


def read(path):
    return json.loads(Path(path).read_text(encoding='utf-8'))


def save(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False), encoding='utf-8')


def table(path, rows, fields=None):
    with Path(path).open('w', encoding='utf-8', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=fields or list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def csv_rows(path):
    with Path(path).open(encoding='utf-8-sig', newline='') as stream:
        yield from csv.DictReader(stream)


def verified(path, digest):
    require(isinstance(digest, str) and sha(path) == digest, f'Missing hash or changed file: {path}')


def unique(values):
    return list(dict.fromkeys(values))


def new_output(path):
    path = Path(path).resolve()
    path.mkdir(parents=True, exist_ok=False)
    save(path / 'audit.json', {'status': 'incomplete'})
    return path


class Registry:
    def __init__(self, path):
        self.path = Path(path).resolve()
        self.connection = sqlite3.connect(self.path.as_uri() + '?mode=ro', uri=True)
        self.connection.execute('PRAGMA query_only=ON')
        self.connection.execute('BEGIN')
        self.cache = {}

    def resolve(self, cid, digest, molecule=None):
        require(isinstance(cid, str) and cid and digest, 'Missing artifact identity/hash')
        if cid not in self.cache:
            record = self.connection.execute(
                'SELECT molecule_id,content_sha256 FROM conformer WHERE id=?', (cid,)).fetchone()
            require(record is not None, f'Artifact CID absent from registry: {cid}')
            self.cache[cid] = record
        mid, expected = self.cache[cid]
        require(expected == digest, f'Source hash mismatch: {cid}')
        require(molecule is None or molecule == mid, f'Canonical molecule mismatch: {cid}')
        return mid

    def close(self):
        self.connection.close()


def load_reference(folder, registry):
    report = read(folder / 'report.json')
    require(report.get('status') == 'complete', 'Historical export is incomplete')
    path = folder / 'selected-conformers.csv'
    verified(path, report.get('outputs', {}).get(path.name))
    conformers, molecules, identities = [], [], {}
    for expected_rank, row in enumerate(csv_rows(path), 1):
        require(int(row['rank']) == expected_rank, 'Historical ranks must be contiguous and ordered')
        cid = row['conformer_id']
        require(cid not in identities, 'Duplicate historical conformer')
        mid = registry.resolve(cid, row['source_sha256'], row['molecule_id'])
        identities[cid] = (int(row['global_id']), mid)
        conformers.append(cid)
        molecules.append(mid)
    require(len(conformers) == report['exported_conformers'] and conformers,
            'Historical export count mismatch or empty reference')
    molecules = unique(molecules)
    require(len(molecules) == report['distinct_molecules'], 'Historical molecule count mismatch')
    require(len({record[0] for record in identities.values()}) == len(identities),
            'Duplicate historical global ID')
    return report, conformers, molecules, identities


def load_arm(folder, registry, reference_ids):
    report = read(folder / 'report.json')
    require(report.get('status') == 'complete' and report.get('kind') == 'block_search_dock',
            f'Incomplete or unexpected campaign: {folder}')
    require(report.get('search_engine') == 'consensus-threshold-v1'
            and report.get('docking_requested') is False, 'Expected threshold search-only arm')
    for name in ('protocol.json', 'selected/report.json', 'query/report.json'):
        verified(folder / name, report.get('output_hashes', {}).get(name))
    protocol = read(folder / 'protocol.json')
    export = read(folder / 'selected/report.json')
    require(export.get('status') == 'complete', 'Incomplete candidate export')
    path = folder / 'selected/candidates.jsonl'
    verified(path, export.get('output_hashes', {}).get('candidates.jsonl'))
    candidates, seen, structural_mids, gids = [], set(), set(), set()
    with path.open(encoding='utf-8') as stream:
        for rank, line in enumerate(stream, 1):
            row = json.loads(line)
            require(row['rank'] == rank, 'Candidate order is not sealed rank order')
            cid = row['artifact_cid']
            require(cid not in seen and row['gid'] not in gids, 'Duplicate candidate artifact/global ID')
            mid = registry.resolve(cid, row['source_sha256'])
            if cid in reference_ids:
                require(reference_ids[cid] == (row['gid'], mid), 'Historical/current global ID mismatch')
            seen.add(cid)
            gids.add(row['gid'])
            structural_mids.add(row['molecule_id'])
            candidates.append((cid, mid))
    require(len(candidates) == report['exported_conformers'] == export['unique_conformers'],
            'Candidate count mismatch')
    require(len(structural_mids) == report['unique_molecules'] == export['unique_molecules'],
            'Reported structural molecule count mismatch')
    require(report['requested_conformers'] - len(candidates) == report['shortfall'] >= 0,
            'Inconsistent shortfall')
    require(0 <= len(candidates) <= report['ranked_conformers'] <= report['searched_conformers'],
            'Invalid search/rank/export counts')
    selection = report['selection']
    name = f"{selection['scheme']}/{selection['method']}/Top-{selection['blocks_per_method']}"
    blocks = report.get('searched_blocks', [])
    require(sum(block['searched_conformers'] for block in blocks) == report['searched_conformers'],
            'Block history work count mismatch')
    return dict(name=name, report=report, protocol=protocol, candidates=candidates,
                path=str(folder), query_sha256=sha(folder / 'query/report.json'),
                report_sha256=sha(folder / 'report.json'), candidate_sha256=sha(path))


def stage_coverage(arm, reference_ids):
    folder = Path(arm['path'])
    report = arm['report']
    for name in ('members.sqlite', 'ranking.sqlite'):
        verified(folder / name, report.get('output_hashes', {}).get(name))
    visited = {row['block_id'] for row in report['searched_blocks']}
    require(len(visited) == len(report['searched_blocks']), 'Repeated searched block')
    selected, searched, qualified = set(), set(), set()
    with sqlite3.connect((folder / 'ranking.sqlite').resolve().as_uri() + '?mode=ro', uri=True) as ranking:
        ranked_gids = {row[0] for row in ranking.execute('SELECT gid FROM ranking')}
    require(len(ranked_gids) == report['ranked_conformers'], 'Ranked database count mismatch')
    examined_count, member_gids, visited_gids = 0, set(), set()
    with sqlite3.connect((folder / 'members.sqlite').resolve().as_uri() + '?mode=ro', uri=True) as members:
        for cid, gid, block in members.execute('SELECT artifact_cid,gid,block_id FROM member'):
            require(gid not in member_gids, 'Duplicate member global ID')
            member_gids.add(gid)
            if block in visited:
                examined_count += 1
                visited_gids.add(gid)
            if cid in reference_ids:
                require(reference_ids[cid][0] == gid, 'Membership/reference global ID mismatch')
                selected.add(cid)
                if block in visited:
                    searched.add(cid)
                if gid in ranked_gids:
                    qualified.add(cid)
    require(examined_count == report['searched_conformers'], 'Visited member count mismatch')
    require(ranked_gids <= visited_gids, 'Ranked candidate outside visited blocks')
    exported = {row[0] for row in arm['candidates']} & reference_ids.keys()
    require(exported <= qualified <= searched <= selected, 'Reference stage sets are not nested')
    return dict(selected_blocks=selected, searched_blocks=searched, qualified=qualified,
                exported=exported)


def overlap(first, second):
    common = len(first & second)
    return dict(intersection=common, jaccard=common / len(first | second) if first | second else None,
                identical=first == second)


def recovery_rows(arm, reference_conformers, reference_molecules, cutoffs, budgets):
    candidates = arm['candidates']
    rows = []
    for budget in sorted(set(budgets) | {len(candidates)}):
        if budget > len(candidates):
            continue
        chosen = candidates[:budget]
        for unit, reference, selected in (
                ('conformer', reference_conformers, {row[0] for row in chosen}),
                ('molecule', reference_molecules, {row[1] for row in chosen})):
            for cutoff in cutoffs:
                if cutoff > len(reference):
                    continue
                hits = len(selected & set(reference[:cutoff]))
                rows.append(dict(arm=arm['name'], unit=unit, reference_k=cutoff,
                                 exported_prefix=budget, full_output=budget == len(candidates),
                                 selected_units=len(selected), reference_hits=hits,
                                 reference_recall=hits / cutoff,
                                 selected_reference_fraction=hits / len(selected) if selected else None,
                                 searched_conformers=arm['report']['searched_conformers']))
    return rows


def scatter(path, rows):
    points = [row for row in rows if row['full_output'] and row['unit'] == 'conformer']
    cutoffs = sorted({row['reference_k'] for row in points})
    width, height = 980, max(1, len(cutoffs)) * 300
    parts = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}">',
             '<rect width="100%" height="100%" fill="white"/>']
    colors = ('#2563eb', '#dc2626', '#059669', '#9333ea', '#d97706', '#0891b2')
    for panel, cutoff in enumerate(cutoffs):
        selected = [row for row in points if row['reference_k'] == cutoff]
        maximum = max(row['searched_conformers'] for row in selected) or 1
        top = panel * 300
        parts.append(f'<text x="60" y="{top+25}">Historical reference Top-{cutoff}: recovery vs searched work (not time)</text>')
        parts.append(f'<path d="M60 {top+50} V{top+250} H600" fill="none" stroke="black"/>')
        parts.append(f'<text x="8" y="{top+60}">100%</text><text x="20" y="{top+250}">0%</text>')
        parts.append(f'<text x="60" y="{top+275}">0</text><text x="475" y="{top+275}">{maximum} conformers</text>')
        for index, row in enumerate(selected):
            horizontal = 60 + 520 * row['searched_conformers'] / maximum
            vertical = top + 250 - 190 * row['reference_recall']
            color = colors[index // 2 % len(colors)]
            label = html.escape(row['arm'])
            parts.append(f'<circle cx="{horizontal:.2f}" cy="{vertical:.2f}" r="{4+index%2*2}" fill="{color}" fill-opacity="0.6"><title>{label}: {row["reference_recall"]:.4f}</title></circle>')
            parts.append(f'<text x="620" y="{top+45+index*18}" fill="{color}" font-size="12">{label}: {row["reference_recall"]:.4f}</text>')
    path.write_text('\n'.join(parts + ['</svg>']), encoding='utf-8')


def historical(args):
    output = new_output(args.output)
    run_ids = args.run_ids or RUNS
    folders = [args.runs_root / ('PROMPT-' + run) / 'execution/campaign/blocks' for run in run_ids]
    registry_path = args.registry or Path(read(folders[0] / 'protocol.json')['library_registry']['path'])
    registry = Registry(registry_path)
    try:
        reference_report, ref_cids, ref_mids, ref_ids = load_reference(args.reference, registry)
        arms = [load_arm(folder, registry, ref_ids) for folder in folders]
    finally:
        registry.close()
    expected = {f'{scheme}/{method}/Top-{limit}' for scheme in SCHEMES for method in METHODS for limit in (5, 10)}
    require(len(arms) == 12 and {arm['name'] for arm in arms} == expected, 'Expected exact twelve-arm matrix')
    arms.sort(key=lambda arm: arm['name'])
    require(len({arm['query_sha256'] for arm in arms}) == 1, 'Current arms use different queries')
    require(len({json.dumps(arm['protocol']['library_registry'], sort_keys=True) for arm in arms}) == 1,
            'Current arms use different library registry snapshots')
    global_ids = {}
    for arm in arms:
        for cid, mid in arm['candidates']:
            require(global_ids.setdefault(cid, mid) == mid, 'Inconsistent canonical identity')
    common_budget = min(len(arm['candidates']) for arm in arms)
    summaries, recovery, pairs, gains, stages = [], [], [], [], []
    for arm in arms:
        report = arm['report']
        cids = {record[0] for record in arm['candidates']}
        summaries.append(dict(arm=arm['name'], searched=report['searched_conformers'],
                              qualified=report['ranked_conformers'], exported=len(cids),
                              canonical_molecules=len({record[1] for record in arm['candidates']}),
                              shortfall=report['shortfall'], blocks=len(report['searched_blocks']),
                              cached_blocks=sum(bool(block.get('cache_reused')) for block in report['searched_blocks']),
                              recorded_seconds=report.get('elapsed_seconds'),
                              timing_scope=report.get('timing_scope', 'unverified'),
                              historical_reference_overlap=len(cids & ref_ids.keys()),
                              outside_historical_reference=len(cids - ref_ids.keys())))
        recovery.extend(recovery_rows(arm, ref_cids, ref_mids, args.top_k, args.budgets + [common_budget]))
        if not args.skip_stages:
            coverage = stage_coverage(arm, ref_ids)
            for cutoff in args.top_k:
                if cutoff > len(ref_cids):
                    continue
                reference = set(ref_cids[:cutoff])
                counts = {stage: len(reference & candidates) for stage, candidates in coverage.items()}
                require(counts['exported'] <= counts['qualified'] <= counts['searched_blocks'] <= counts['selected_blocks'],
                        'Invalid reference stage containment')
                stages.append(dict(arm=arm['name'], reference_k=cutoff, **counts,
                                   outside_selected_blocks=cutoff-counts['selected_blocks'],
                                   unvisited_after_stop=counts['selected_blocks']-counts['searched_blocks'],
                                   not_currently_qualified=counts['searched_blocks']-counts['qualified'],
                                   dropped_by_output_ranking=counts['qualified']-counts['exported']))
    require(recovery, 'No reference cutoffs fit the available reference')
    for first, second in itertools.combinations(arms, 2):
        for index, unit in enumerate(('conformer', 'molecule')):
            left = {row[index] for row in first['candidates']}
            right = {row[index] for row in second['candidates']}
            pairs.append(dict(arm_a=first['name'], arm_b=second['name'], unit=unit, **overlap(left, right)))
    by_name = {arm['name']: arm for arm in arms}
    for scheme, method in itertools.product(SCHEMES, METHODS):
        first, second = [by_name[f'{scheme}/{method}/Top-{limit}'] for limit in (5, 10)]
        for index, unit in enumerate(('conformer', 'molecule')):
            left, right = [{row[index] for row in arm['candidates']} for arm in (first, second)]
            gains.append(dict(scheme=scheme, scorer=method, unit=unit, added=len(right-left),
                              lost=len(left-right), identical=left == right,
                              extra_searched=second['report']['searched_conformers']-first['report']['searched_conformers']))
    table(output / 'arm_summary.csv', summaries)
    table(output / 'historical_recovery.csv', recovery)
    table(output / 'pairwise_overlap.csv', pairs)
    table(output / 'top5_top10_gains.csv', gains)
    if stages:
        table(output / 'reference_stage_coverage.csv', stages)
    scatter(output / 'recovery_work.svg', recovery)
    warnings = [
        'Historical ANN-selected reference is not an exhaustive same-policy library ground truth.',
        'Old quota-tier RRF and current threshold-funnel/RRF policies differ.',
        'Prefix results share the full arm search cost; they are not partial-runtime measurements.',
        'Molecule reference rank is first appearance in historical conformer rank, not the old molecule rank.',
        'Outside-reference candidates are unknown, not false positives.',
        'Existing cached wall times cannot establish cold speedup. No biological enrichment was evaluated.',
        'Registry content identities were checked; whole-library/query equivalence to the old search remains unverified.']
    audit = dict(status='complete', mode='historical_reference_recovery', common_export_prefix=common_budget,
                 reference_report=reference_report, reference_report_sha256=sha(args.reference / 'report.json'),
                 reference_csv_sha256=sha(args.reference / 'selected-conformers.csv'),
                 registry=str(registry.path), arms=[{key: arm[key] for key in ('name', 'path', 'query_sha256', 'report_sha256', 'candidate_sha256')} for arm in arms],
                 arguments=dict(top_k=args.top_k, budgets=args.budgets, skip_stages=args.skip_stages), warnings=warnings,
                 missing_stages=['complete_universe_random_control_replay', 'matched_actual_runtime'],
                 source_sha256=sha(Path(__file__)), biological_validation=False)
    save(output / 'audit.json', audit)
    (output / 'report.md').write_text('# Historical routing benchmark\n\n'
        'All twelve arms were audited against canonical registry identities.\n'
        f'Common exported prefix: {common_budget}. See historical_recovery.csv and recovery_work.svg.\n\n'
        + '\n'.join('- ' + warning for warning in warnings) + '\n', encoding='utf-8')
    return audit


def load_universe(path):
    rows, seen, ranks = [], set(), []
    for row in csv_rows(path):
        cid, mid = row['conformer_id'], row['molecule_id']
        require(cid and mid and cid not in seen, 'Empty or duplicate universe identity')
        require(row['qualified'] in ('0', '1'), 'qualified must be 0 or 1, including rejected candidates')
        qualified = row['qualified'] == '1'
        rank = int(row['reference_rank']) if row['reference_rank'] else None
        require((qualified and rank is not None and rank > 0) or (not qualified and rank is None),
                'Every qualified candidate needs a rank; rejected candidates must not have one')
        require(all(row.get(scheme) for scheme in SCHEMES), 'Missing regular block membership')
        seen.add(cid)
        if qualified:
            ranks.append(rank)
        rows.append(dict(cid=cid, mid=mid, qualified=qualified, rank=rank,
                         **{scheme: row[scheme] for scheme in SCHEMES}))
    require(rows and sorted(ranks) == list(range(1, len(ranks)+1)) and ranks,
            'Reference ranks must cover all qualified rows exactly once')
    return rows


def route_curve(groups, reference, total, method, seed, unit, cutoff):
    found, searched, qualified = set(), 0, 0
    curves = []
    for step, members in enumerate(groups, 1):
        searched += len(members)
        for row in members:
            qualified += row['qualified']
            if row['qualified']:
                identity = row['cid'] if unit == 'conformer' else row['mid']
                if identity in reference:
                    found.add(identity)
        curves.append(dict(method=method, seed=seed, unit=unit, reference_k=cutoff, step=step,
                           searched=searched, work_fraction=searched/total, qualified=qualified,
                           reference_hits=len(found), reference_recall=len(found)/len(reference)))
    return curves


def replay_schedules(populations, priorities, rows, args):
    for (scheme, scorer), ordering in sorted(priorities.items()):
        yield (f'{scheme}/{scorer}', -1,
               [populations[scheme][block] for _, block in ordering if block in populations[scheme]])
    for seed in range(args.seed, args.seed + args.repeats):
        for scheme in SCHEMES:
            blocks = sorted(populations[scheme])
            random.Random(f'{seed}/{scheme}').shuffle(blocks)
            yield f'{scheme}/random_blocks', seed, [populations[scheme][block] for block in blocks]
        shuffled = list(rows)
        random.Random(seed).shuffle(shuffled)
        yield ('uniform_random_conformers', seed,
               [shuffled[start:start+args.uniform_chunk] for start in range(0, len(rows), args.uniform_chunk)])


def replay(args):
    output = new_output(args.output)
    manifest = read(args.manifest)
    require(manifest.get('scope') in ('complete_regular_library', 'complete_fixed_sublibrary'),
            'Declare a complete evaluated regular universe, not the selected union')
    for field in ('library_identity', 'query_identity', 'evaluation_policy', 'priority_fit_exclusion'):
        require(bool(manifest.get(field)), f'Missing universe provenance: {field}')
    verified(args.universe, manifest.get('universe_sha256'))
    verified(args.priorities, manifest.get('priorities_sha256'))
    rows = load_universe(args.universe)
    require(len(rows) == manifest['conformers'], 'Universe population count mismatch')
    priorities = {}
    for row in csv_rows(args.priorities):
        key = (row['scheme'], row['scorer'])
        require(key in set(itertools.product(SCHEMES, METHODS)), 'Unknown priority scheme/scorer')
        priorities.setdefault(key, []).append((int(row['rank']), row['block_id']))
    require(set(priorities) == set(itertools.product(SCHEMES, METHODS)), 'Need all six frozen priority lists')
    populations = {scheme: {} for scheme in SCHEMES}
    for row in rows:
        for scheme in SCHEMES:
            populations[scheme].setdefault(row[scheme], []).append(row)
    for (scheme, scorer), ordering in sorted(priorities.items()):
        ordering.sort()
        require([rank for rank, _ in ordering] == list(range(1, len(ordering)+1)), 'Priority ranks must be contiguous')
        require(len({block for _, block in ordering}) == len(ordering), 'Duplicate priority block')
        require(set(populations[scheme]) <= {block for _, block in ordering}, 'Unranked universe block')
    ranked = sorted((row for row in rows if row['qualified']), key=lambda row: row['rank'])
    summaries, budgets = [], []
    writer = None
    with (output / 'routing_curves.csv').open('w', encoding='utf-8', newline='') as stream:
        for method, seed, groups in replay_schedules(populations, priorities, rows, args):
            for unit, ordered in [('conformer', [row['cid'] for row in ranked]),
                                  ('molecule', unique(row['mid'] for row in ranked))]:
                for cutoff in args.top_k:
                    if cutoff > len(ordered):
                        continue
                    curve = route_curve(groups, set(ordered[:cutoff]), len(rows), method, seed, unit, cutoff)
                    if writer is None:
                        writer = csv.DictWriter(stream, fieldnames=list(curve[0]))
                        writer.writeheader()
                    writer.writerows(curve)
                    reached = next(row for row in curve if row['reference_recall'] >= args.target_recall)
                    summaries.append(dict(method=method, seed=seed, unit=unit, reference_k=cutoff,
                                          target_recall=args.target_recall, work_to_target=reached['searched'],
                                          work_fraction_to_target=reached['work_fraction'],
                                          achieved_recall=reached['reference_recall']))
                    for fraction in (0.01, 0.05, 0.1, 0.2, 0.5, 1.0):
                        completed = [row for row in curve if row['searched'] <= math.floor(fraction * len(rows))]
                        last = completed[-1] if completed else dict(searched=0, reference_recall=0)
                        budgets.append(dict(method=method, seed=seed, unit=unit, reference_k=cutoff,
                                            work_budget_fraction=fraction, actual_work=last['searched'],
                                            reference_recall=last['reference_recall']))
    require(summaries, 'No cutoffs fit qualified reference')
    table(output / 'recall_at_work.csv', budgets)
    table(output / 'work_to_recall.csv', summaries)
    aggregates = []
    for key in sorted({(row['method'], row['unit'], row['reference_k']) for row in summaries}):
        values = [row['work_to_target'] for row in summaries if (row['method'], row['unit'], row['reference_k']) == key]
        aggregates.append(dict(method=key[0], unit=key[1], reference_k=key[2], repetitions=len(values),
                               median_work=statistics.median(values), min_work=min(values), max_work=max(values)))
    table(output / 'replay_summary.csv', aggregates)
    audit = dict(status='complete', mode='routing_discovery_replay_not_runtime_or_final_export',
                 manifest=manifest, manifest_sha256=sha(args.manifest), seed=args.seed, repeats=args.repeats,
                 top_k=args.top_k, target_recall=args.target_recall, uniform_chunk=args.uniform_chunk,
                 source_sha256=sha(Path(__file__)), actual_search_launched=False,
                 limitations=['Completeness and fit/test exclusion are manifest declarations, not independently proven.',
                              'Whole-block boundary costs, no interpolation; uniform baseline is batched separately.',
                              'Qualified discovery is before RRF/export truncation; no measured wall-time speedup.',
                              'Random repetitions quantify random-order variation, not uncertainty across targets.',
                              'Special/excluded pools and setup costs are outside this regular-universe replay.'])
    save(output / 'audit.json', audit)
    return audit


def measure(args):
    output = new_output(args.output)
    spec = read(args.spec)
    require(spec.get('phase') in ('setup', 'online'), 'Timing phase must be setup or online')
    require(spec.get('cache_policy') in ('fresh_application_cache', 'warm_application_cache'), 'Declare cache policy')
    for field in ('method', 'comparison_id', 'hardware', 'workers', 'protocol_identity', 'repeat', 'command', 'cwd', 'result_receipt'):
        require(field in spec and spec[field] not in ('', None), f'Missing timing specification: {field}')
    require(isinstance(spec['command'], list) and all(isinstance(item, str) for item in spec['command']) and spec['command'],
            'command must be an argv list, never a shell string')
    receipt = Path(spec['result_receipt'])
    if not receipt.is_absolute():
        receipt = Path(spec['cwd']) / receipt
    require(not receipt.exists(), 'Result receipt already exists; do not time a resumed/completed run')
    started = time.perf_counter()
    with (output / 'stdout.log').open('wb') as stdout, (output / 'stderr.log').open('wb') as stderr:
        process = subprocess.run(spec['command'], cwd=spec['cwd'], stdout=stdout, stderr=stderr, shell=False)
    elapsed = time.perf_counter() - started
    result = read(receipt) if receipt.exists() else {}
    valid = (process.returncode == 0 and result.get('status') == 'complete'
             and result.get('cache_reused') is False and result.get('resumed_execution') is False)
    audit = dict(status='complete' if valid else 'ineligible', mode='actual_command_timing',
                 spec=spec, spec_sha256=sha(args.spec), wall_seconds=elapsed, returncode=process.returncode,
                 receipt_path=str(receipt), receipt_sha256=sha(receipt) if receipt.exists() else None,
                 result_cache_reused=result.get('cache_reused'), result_resumed=result.get('resumed_execution'),
                 eligible_fresh_execution=valid, source_sha256=sha(Path(__file__)),
                 recovery_and_budget_equivalence_verified=False,
                 limitations=['OS page cache is not cleared or measured.',
                              'Hardware and comparison conditions are user declarations.',
                              'Analyze achieved recovery and output budget before interpreting timing ratios.'])
    save(output / 'audit.json', audit)
    require(valid, 'Measurement retained but ineligible: failure, cache reuse, resume, or missing explicit receipt flags')
    return audit


def timing_summary(args):
    output = new_output(args.output)
    records = [read(path) for path in args.receipts]
    require(records and all(row.get('eligible_fresh_execution') is True and row.get('status') == 'complete' for row in records),
            'Only successful fresh-execution timing receipts may be summarized')
    conditions = ('comparison_id', 'hardware', 'workers', 'protocol_identity', 'phase', 'cache_policy')
    require(len({json.dumps({key: row['spec'][key] for key in conditions}, sort_keys=True) for row in records}) == 1,
            'Timing conditions differ; do not pool these receipts')
    groups = {}
    for row in records:
        spec = row['spec']
        require(math.isfinite(row['wall_seconds']) and row['wall_seconds'] > 0, 'Invalid wall time')
        trials = groups.setdefault(spec['method'], {})
        require(spec['repeat'] not in trials, 'Duplicate method/repeat')
        trials[spec['repeat']] = row['wall_seconds']
    require(args.baseline in groups, 'Missing baseline method')
    trials = set(groups[args.baseline])
    require(len(trials) >= 3 and all(set(group) == trials for group in groups.values()),
            'Require at least three matched repeats per method')
    rows = []
    for method, values in sorted(groups.items()):
        ratios = [groups[args.baseline][trial] / values[trial] for trial in sorted(trials)]
        rows.append(dict(method=method, repeats=len(values), median_seconds=statistics.median(values.values()),
                         min_seconds=min(values.values()), max_seconds=max(values.values()),
                         median_paired_time_ratio=statistics.median(ratios)))
    table(output / 'timing_summary.csv', rows)
    audit = dict(status='complete', mode='declared_matched_timing_only', baseline=args.baseline,
                 receipts={str(path): sha(path) for path in args.receipts},
                 recovery_and_budget_equivalence_verified=False,
                 limitation='Time ratios are not accepted screening speedups until recovery and budgets are matched.')
    save(output / 'audit.json', audit)
    return audit


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest='mode', required=True)
    historic = commands.add_parser('historical', help='Audit all twelve arms against the existing historical export')
    historic.add_argument('--reference', type=Path, required=True)
    historic.add_argument('--runs-root', type=Path, required=True)
    historic.add_argument('--registry', type=Path)
    historic.add_argument('--run-ids', nargs='+')
    historic.add_argument('--skip-stages', action='store_true', help='Skip membership/ranking DB audit if those files were not transferred')
    historic.add_argument('--top-k', nargs='+', type=int, default=[100, 1000, 10000, 100000, 1000000])
    historic.add_argument('--budgets', nargs='+', type=int, default=[10000, 25000, 50000, 100000])
    historic.set_defaults(handler=historical)
    simulation = commands.add_parser('replay', help='Replay frozen routing over a complete evaluated universe')
    for name in ('manifest', 'universe', 'priorities'):
        simulation.add_argument('--' + name, type=Path, required=True)
    simulation.add_argument('--top-k', nargs='+', type=int, default=[100, 1000, 10000])
    simulation.add_argument('--repeats', type=int, default=20)
    simulation.add_argument('--seed', type=int, default=20261014)
    simulation.add_argument('--uniform-chunk', type=int, default=1000)
    simulation.add_argument('--target-recall', type=float, default=0.9)
    simulation.set_defaults(handler=replay)
    timing = commands.add_parser('measure', help='Explicitly execute and time one user-specified command')
    timing.add_argument('--spec', type=Path, required=True)
    timing.set_defaults(handler=measure)
    summary = commands.add_parser('timing-summary', help='Summarize at least three declared-matched repeats')
    summary.add_argument('--receipts', nargs='+', type=Path, required=True)
    summary.add_argument('--baseline', required=True)
    summary.set_defaults(handler=timing_summary)
    for command in (historic, simulation, timing, summary):
        command.add_argument('--output', type=Path, required=True)
    args = parser.parse_args(argv)
    for name in ('top_k', 'budgets'):
        require(all(value > 0 for value in getattr(args, name, [])), 'Cutoffs and budgets must be positive')
    if args.mode == 'replay':
        require(args.repeats > 0 and args.uniform_chunk > 0 and 0 < args.target_recall <= 1, 'Invalid replay parameters')
    args.handler(args)
    print('Saved benchmark to', args.output)


if __name__ == '__main__':
    main()
