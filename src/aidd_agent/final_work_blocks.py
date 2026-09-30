"""Validated E091 overlays and size-controlled, enumerable execution blocks."""
import argparse
from collections import Counter, defaultdict
import csv
import hashlib
import json
import math
from pathlib import Path
import sqlite3
import time

from .boundary_pair_review import readonly

VERSION = 'final-work-blocks-v1'


def sha(path):
    h = hashlib.sha256()
    processed = 0
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda: stream.read(8 * 1024 * 1024), b''):
            h.update(chunk)
            processed += len(chunk)
            if processed % (4 * 1024 ** 3) == 0:
                print(f'Hash progress {Path(path).name}: {processed // 1024 ** 3} GiB', flush=True)
    return h.hexdigest()


def read(path):
    return json.loads(Path(path).read_text(encoding='utf-8'))


def sealed(directory, names):
    directory = Path(directory).resolve()
    report = read(directory / 'report.json')
    if report.get('status') != 'complete':
        raise ValueError('Completed input required: ' + str(directory))
    for name in names:
        print('Verifying ' + str(directory / name), flush=True)
        if sha(directory / name) != report['output_hashes'][name]:
            raise ValueError('Input hash changed: ' + name)
    return report


def identifier(value):
    return 'work-' + hashlib.sha256(json.dumps(value, sort_keys=True).encode()).hexdigest()[:24]


def allocate(classes, minimum):
    """Pool tails deterministically without erasing original class identities."""
    if type(minimum) is not int or minimum < 1:
        raise ValueError('Minimum size must be a positive integer')
    pools, small = [], defaultdict(list)
    for row in sorted(classes, key=lambda r: r['hard_group']):
        if row['conformers'] <= 0:
            raise ValueError('Empty class')
        if row['conformers'] >= minimum:
            pools.append([row])
        else:
            small[(row['unit_count'], row['omega_pattern'])].append(row)
    tails = []
    for _, rows in sorted(small.items()):
        (pools if sum(r['conformers'] for r in rows) >= minimum else tails).append(rows)
    if not pools and tails:
        pools.append([r for rows in tails for r in rows]); tails = []
    for rows in tails:
        first = rows[0]
        def preference(pool):
            exact = any(r['omega_pattern'] == first['omega_pattern'] for r in pool)
            same_length = any(r['unit_count'] == first['unit_count'] for r in pool)
            return (not exact, not same_length, sum(r['conformers'] for r in pool),
                    min(r['hard_group'] for r in pool))
        min(pools, key=preference).extend(rows)
    result = []
    for rows in pools:
        groups = sorted(r['hard_group'] for r in rows)
        patterns = sorted(set(r['omega_pattern'] for r in rows))
        lengths = sorted(set(r['unit_count'] for r in rows))
        n = sum(r['conformers'] for r in rows)
        result.append(dict(block_id=identifier([VERSION, minimum, groups]),
            kind='single_class' if len(rows) == 1 else 'same_pattern_pool' if len(patterns) == 1 else 'mixed_execution_pool',
            conformers=n, logical_classes=len(rows), unit_counts=lengths,
            omega_patterns=patterns, groups=groups, below_minimum=n < minimum))
    if len(result) > 1 and any(r['below_minimum'] for r in result):
        raise AssertionError('Tiny regular block escaped pooling')
    return sorted(result, key=lambda r: r['block_id'])


def write_csv(path, rows, columns):
    with Path(path).open('w', newline='', encoding='utf-8') as stream:
        writer = csv.DictWriter(stream, fieldnames=columns); writer.writeheader()
        writer.writerows({k: json.dumps(v) if isinstance(v, (list, dict)) else v
                         for k, v in row.items() if k in columns} for row in rows)


def build_blocks(model, diagnostic, candidates, routing, output, minimum=5000):
    started = time.monotonic()
    if type(minimum) is not int or minimum < 1:
        raise ValueError('Positive minimum size required')
    model, diagnostic, candidates, routing, output = [Path(p).resolve() for p in (model, diagnostic, candidates, routing, output)]
    if any(output == p or p in output.parents or output in p.parents for p in (model, diagnostic, candidates, routing)):
        raise ValueError('Output must be a separate fresh directory')
    if output.exists():
        raise FileExistsError(output)
    mr = sealed(model, ['blocks.sqlite'])
    dr = sealed(diagnostic, ['hard_groups.csv'])
    cr = sealed(candidates, ['candidates.sqlite'])
    rr = sealed(routing, ['boundary_assignments.sqlite', 'routing_models.json'])
    if dr['provenance']['model_report_sha256'] != sha(model / 'report.json'):
        raise ValueError('Diagnostic/model mismatch')
    if rr['sources'].get(str(candidates / 'report.json')) != sha(candidates / 'report.json'):
        raise ValueError('Routing/candidate mismatch')
    for receipt in (cr, rr):
        for path, digest in receipt['sources'].items():
            if sha(Path(path)) != digest:
                raise ValueError('Provenance changed: ' + path)
    with (diagnostic / 'hard_groups.csv').open(newline='', encoding='utf-8') as stream:
        groups = {}
        for row in csv.DictReader(stream):
            g = row['hard_group']; states = row['omega_pattern'].split(';')
            if g in groups or set(states) - {'cis', 'trans', 'boundary'}:
                raise ValueError('Invalid or duplicate hard group')
            groups[g] = dict(hard_group=g, conformers=int(row['conformers']),
                omega_pattern=row['omega_pattern'], unit_count=len(states), boundary='boundary' in states)
    models = read(routing / 'routing_models.json')
    output.mkdir()
    try:
        with readonly(model / 'blocks.sqlite') as original, sqlite3.connect(output / 'work_blocks.sqlite') as dest:
            dest.execute('PRAGMA temp_store=FILE')
            dest.executescript('''CREATE TABLE class_map(group_id TEXT PRIMARY KEY,block_id TEXT NOT NULL,
                original_n INTEGER NOT NULL,assigned_n INTEGER NOT NULL,omega_pattern TEXT NOT NULL);
                CREATE TABLE override(cid TEXT PRIMARY KEY,original_group TEXT NOT NULL,target_group TEXT,
                state TEXT NOT NULL,rotation INTEGER,score REAL,block_id TEXT);
                CREATE TABLE block(block_id TEXT PRIMARY KEY,kind TEXT NOT NULL,n INTEGER NOT NULL,metadata TEXT NOT NULL);''')
            populations = dict(original.execute('SELECT group_id,count(*) FROM point GROUP BY group_id'))
            if populations != {g: r['conformers'] for g, r in groups.items()}:
                raise ValueError('Actual source populations differ from diagnostic')
            if sum(populations.values()) != rr['total_conformers'] or sum(populations.values()) != mr['counts']['conformers']:
                raise ValueError('Source population mismatch')
            original.execute('ATTACH DATABASE ? AS route', ((routing / 'boundary_assignments.sqlite').as_uri() + '?mode=ro',))
            original.execute('ATTACH DATABASE ? AS cand', ((candidates / 'candidates.sqlite').as_uri() + '?mode=ro',))
            counts, by_original, assigned = Counter(), Counter(), Counter()
            sql = '''SELECT a.cid,a.original_group,a.state,a.target_group,a.rotation,a.score,a.evidence,
                p.group_id,c.original_group,c.maximum_deviation,c.candidates
                FROM route.assignment a LEFT JOIN point p ON p.cid=a.cid
                LEFT JOIN cand.candidate c ON c.cid=a.cid'''
            for row in original.execute(sql):
                cid, g, state, target, rotation, score, encoded, source_group, candidate_group, deviation, matches_json = row
                if source_group != g or candidate_group != g or g not in groups or not groups[g]['boundary']:
                    raise ValueError('Boundary CID missing or original class changed: ' + str(cid))
                if state == 'assigned_provisional':
                    if target not in groups or groups[target]['boundary']:
                        raise ValueError('Invalid target class')
                    if groups[g]['unit_count'] != groups[target]['unit_count'] or type(rotation) is not int or not 0 <= rotation < groups[g]['unit_count']:
                        raise ValueError('Invalid directed rotation')
                    if deviation is None or deviation > 45 or score is None or not math.isfinite(score) or not 0 <= score <= 1:
                        raise ValueError('Invalid routing gate or score')
                    allowed = any(m['hard_group'] == target and m['rotation'] == rotation for m in json.loads(matches_json))
                    evidence = json.loads(encoded)
                    accepted = [e for e in evidence.get('accepted', []) if e['hard_group'] == target and e['rotation'] == rotation and abs(e['score'] - score) < 1e-10]
                    if not allowed or not accepted or models.get(target, {}).get('status') != 'routable':
                        raise ValueError('Assignment lacks eligible calibrated target evidence')
                    if any(e.get('prototype_id') not in models[target]['prototype_ids'] for e in accepted):
                        raise ValueError('Unknown supporting prototype')
                    assigned[target] += 1
                elif state == 'special':
                    if target is not None or rotation is not None or score is not None:
                        raise ValueError('Special member has an adopted target')
                else:
                    raise ValueError('Unknown routing state')
                dest.execute('INSERT INTO override VALUES(?,?,?,?,?,?,NULL)', (cid, g, target, state, rotation, score))
                counts[state] += 1; by_original[g] += 1
                if sum(counts.values()) % 50000 == 0:
                    dest.commit(); print('Validated boundary rows: ' + str(sum(counts.values())), flush=True)
            if dict(by_original) != {g: n for g, n in populations.items() if groups[g]['boundary']}:
                raise ValueError('Boundary rows are missing or duplicated')
            if {k: v for k, v in rr['boundary_counts'].items() if v} != dict(counts):
                raise ValueError('Routing count mismatch')
            regular = [dict(r, conformers=r['conformers'] + assigned[g]) for g, r in groups.items() if not r['boundary']]
            blocks = allocate(regular, minimum)
            class_to_block = {g: b['block_id'] for b in blocks for g in b['groups']}
            for row in regular:
                g = row['hard_group']
                dest.execute('INSERT INTO class_map VALUES(?,?,?,?,?)',
                    (g, class_to_block[g], populations[g], assigned[g], row['omega_pattern']))
            for block in blocks:
                dest.execute('INSERT INTO block VALUES(?,?,?,?)',
                    (block['block_id'], block['kind'], block['conformers'], json.dumps(block)))
            if counts['special']:
                dest.execute('INSERT INTO block VALUES(?,?,?,?)', ('special-exhaustive', 'special_exhaustive', counts['special'], '{}'))
            dest.execute("UPDATE override SET block_id=CASE WHEN state='special' THEN 'special-exhaustive' ELSE (SELECT block_id FROM class_map WHERE group_id=target_group) END")
            if dest.execute('SELECT count(*) FROM override WHERE block_id IS NULL').fetchone()[0]:
                raise ValueError('Unmapped boundary members')
            dest.execute('CREATE INDEX override_block ON override(block_id)')
            dest.execute('CREATE INDEX class_block ON class_map(block_id)'); dest.commit()
            regular_n = sum(b['conformers'] for b in blocks)
            if regular_n != rr['provisional_regular_conformers'] or regular_n + counts['special'] != rr['total_conformers']:
                raise ValueError('Final population conservation failed')
        write_csv(output / 'blocks.csv', blocks, ['block_id', 'kind', 'conformers', 'logical_classes', 'unit_counts', 'omega_patterns', 'below_minimum'])
        class_rows = [dict(r, block_id=class_to_block[r['hard_group']], assigned_boundary=assigned[r['hard_group']]) for r in regular]
        write_csv(output / 'classes.csv', class_rows, ['hard_group', 'block_id', 'conformers', 'assigned_boundary', 'omega_pattern', 'unit_count'])
        result = dict(status='complete', version=VERSION, readiness='offline_enumerable_blocks_not_validated_search_dispatch',
            model=str(model), diagnostic=str(diagnostic), candidates=str(candidates), routing=str(routing),
            minimum_regular_block_size=minimum, logical_capacity_limit=None,
            source_conformers=rr['total_conformers'], regular_conformers=regular_n, special_conformers=counts['special'],
            special_fraction=counts['special']/rr['total_conformers'], original_blocks=mr['counts']['blocks'],
            regular_work_blocks=len(blocks), special_work_blocks=int(bool(counts['special'])),
            smallest_regular_block=min((b['conformers'] for b in blocks), default=0),
            largest_regular_block=max((b['conformers'] for b in blocks), default=0),
            blocks_below_minimum=sum(b['below_minimum'] for b in blocks),
            kind_counts=dict(Counter(b['kind'] for b in blocks)),
            class_identity_preserved=True, complete_boundary_cid_validation=True,
            original_inputs_modified=False, wall_seconds=time.monotonic()-started,
            sources={str(p / 'report.json'): sha(p / 'report.json') for p in (model, diagnostic, candidates, routing)},
            output_hashes={n: sha(output / n) for n in ['work_blocks.sqlite', 'blocks.csv', 'classes.csv']},
            implementation_sha256=sha(Path(__file__)),
            limitations=['Mixed execution pools do not imply chemical equivalence; original class and chirality hashes remain available',
                'No maximum population cap; runtime batching is separate',
                'Special members remain exhaustive; no timing or recall claim',
                'Source-only CIDs require registry joining; old radii must not reject expanded memberships',
                'Property refinement requires full source-verified profile coverage',
                'The separate source chemistry-review population is outside this artifact'])
        (output / 'report.json').write_text(json.dumps(result, indent=2), encoding='utf-8')
        return result
    except Exception:
        (output / 'FAILED.txt').write_text('Build incomplete. Use a fresh output for retry.', encoding='utf-8')
        raise


def members(directory, block_id=None):
    """Stream exact identities without rewriting the full original library."""
    directory = Path(directory).resolve(); report = read(directory / 'report.json')
    if report['status'] != 'complete':
        raise ValueError('Completed work blocks required')
    with readonly(directory / 'work_blocks.sqlite') as overlay, readonly(Path(report['model']) / 'blocks.sqlite') as source:
        if block_id and not overlay.execute('SELECT 1 FROM block WHERE block_id=?', (block_id,)).fetchone():
            raise ValueError('Unknown work block')
        clauses = ' WHERE block_id=?' if block_id else ''; args = (block_id,) if block_id else ()
        for group, block in overlay.execute('SELECT group_id,block_id FROM class_map' + clauses, args):
            for cid, mid, node in source.execute('SELECT cid,mid,node FROM point WHERE group_id=?', (group,)):
                yield dict(cid=cid, mid=mid, original_group=group, logical_class=group, block_id=block,
                           state='unchanged_definite', rotation=0, old_leaf=node)
        source.execute('ATTACH DATABASE ? AS overlay', ((directory / 'work_blocks.sqlite').as_uri() + '?mode=ro',))
        where = ' WHERE o.block_id=?' if block_id else ''
        for cid, mid, node, original, target, block, state, rotation in source.execute('''SELECT o.cid,p.mid,p.node,
                o.original_group,o.target_group,o.block_id,o.state,o.rotation FROM overlay.override o
                JOIN point p ON p.cid=o.cid''' + where, args):
            yield dict(cid=cid, mid=mid, original_group=original, logical_class=target,
                       block_id=block, state=state, rotation=rotation, old_leaf=node)


def validate(directory, property_blocks=None):
    """Independently enumerate all reported memberships, without coordinate recomputation."""
    directory = Path(directory).resolve()
    report = sealed(directory, ['work_blocks.sqlite', 'blocks.csv', 'classes.csv'])
    mr = sealed(report['model'], ['blocks.sqlite'])
    for path, digest in report['sources'].items():
        if sha(path) != digest:
            raise ValueError('Upstream receipt changed')
    if report['source_conformers'] != mr['counts']['conformers']:
        raise ValueError('Source population changed')
    actual = Counter(); n = 0
    for row in members(directory):
        actual[row['block_id']] += 1; n += 1
        if n % 1000000 == 0:
            print(f'Validated work memberships: {n}', flush=True)
    with readonly(directory / 'work_blocks.sqlite') as db:
        expected = dict(db.execute('SELECT block_id,n FROM block'))
    if dict(actual) != expected or n != report['source_conformers']:
        raise ValueError('Enumerated membership counts differ from report')
    result = dict(status='complete', structural_gate='passed', conformers=n, work_blocks=len(actual),
                  source_coordinates_recomputed=False, scientific_recall='not_validated')
    if property_blocks:
        property_blocks = Path(property_blocks).resolve()
        pr = sealed(property_blocks, ['property_blocks.sqlite', 'blocks.csv'])
        if Path(pr['parent_blocks']).resolve() != directory:
            raise ValueError('Wrong parent work blocks')
        for path, digest in pr['sources'].items():
            if sha(path) != digest:
                raise ValueError('Property source receipt changed')
        with readonly(property_blocks / 'property_blocks.sqlite') as db:
            db.execute('ATTACH DATABASE ? AS source', ((Path(report['model']) / 'blocks.sqlite').as_uri() + '?mode=ro',))
            db.execute('ATTACH DATABASE ? AS work', ((directory / 'work_blocks.sqlite').as_uri() + '?mode=ro',))
            counts = Counter()
            for cid, child, parent, truth, child_parent in db.execute('''SELECT m.cid,m.block_id,m.parent_block,
                COALESCE(o.block_id,c.block_id),b.parent_block FROM membership m
                LEFT JOIN source.point p ON p.cid=m.cid LEFT JOIN work.override o ON o.cid=m.cid
                LEFT JOIN work.class_map c ON c.group_id=p.group_id LEFT JOIN block b ON b.block_id=m.block_id'''):
                if truth != parent or parent == 'special-exhaustive' or child_parent != parent:
                    raise ValueError('Property member has a wrong parent or unknown CID')
                counts[child] += 1
            if dict(counts) != dict(db.execute('SELECT block_id,n FROM block')) or sum(counts.values()) != report['regular_conformers']:
                raise ValueError('Property memberships do not conserve the regular library')
            if len(counts) > 1 and min(counts.values()) < report['minimum_regular_block_size']:
                raise ValueError('Property refinement created a tiny block')
        result.update(property_gate='passed', property_blocks=len(counts), regular_conformers=sum(counts.values()))
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ['model', 'diagnostic', 'candidates', 'routing', 'output']:
        parser.add_argument('--' + name, type=Path, required=True)
    parser.add_argument('--minimum-block-size', type=int, default=5000)
    a = parser.parse_args()
    print(json.dumps(build_blocks(a.model, a.diagnostic, a.candidates, a.routing, a.output, a.minimum_block_size), indent=2))


if __name__ == '__main__':
    main()
