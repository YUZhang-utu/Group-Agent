"""Top-N block prioritization over sealed saved poses; no docking execution."""
from collections import defaultdict
from contextlib import closing
import csv
import json
import math
from pathlib import Path
import sqlite3

from .block_sampling import readonly
from .final_work_blocks import sha
from .joint_spatial_profiles import save


def options(top_n=10, ranking_unit='molecule'):
    if type(top_n) is not int or not 1 <= top_n <= 100:
        raise ValueError('top_n must be an integer from 1 to 100')
    if ranking_unit not in {'molecule', 'conformer'}:
        raise ValueError('ranking_unit must be molecule or conformer')
    return dict(top_n=top_n, ranking_unit=ranking_unit)


def add_rankings(source, output, *, top_n=10, ranking_unit='molecule'):
    """Extend an already verified block analysis with bounded candidate tables.

    Keep existing engine/source code seals unchanged. Independent molecule mode
    uses the best scored in-block conformer per molecule and receptor. Incomplete
    blocks and blocks with fewer than N scored units have no primary rank.
    """
    config = options(top_n, ranking_unit)
    output, source = Path(output), Path(source)
    report = json.loads((output / 'report.json').read_text(encoding='utf-8'))
    run = json.loads(source.read_text(encoding='utf-8'))
    signature = json.loads((output / 'signature.json').read_text(encoding='utf-8'))
    if signature.get('sha256') != sha(source):
        raise ValueError('Analysis does not match source report')
    if run.get('score') in {'chemplp', 'ChemPLP', 'ChemPLP TOTAL_SCORE (lower is better)'}:
        score_name, better, direction = 'ChemPLP', 'lower', 1
    elif run.get('kind') == 'block_equiscore_run' and run.get('score') == 'EquiScore' and run.get('better') == 'higher':
        score_name, better, direction = 'EquiScore', 'higher', -1
    else:
        raise ValueError('Unsupported score or scoring direction')
    for name in ('block_scores.csv', 'block_summary.csv'):
        if sha(output / name) != report['output_hashes'][name]:
            raise ValueError('Analysis score table changed')
    sample = Path(run['sample_report']).parent / 'samples.sqlite'
    with readonly(sample) as db:
        molecules = dict(db.execute('SELECT cid,mid FROM selected'))
    with (output / 'block_summary.csv').open(encoding='utf-8', newline='') as stream:
        summaries = {(r['scheme'], r['block_id'], r['receptor']): r for r in csv.DictReader(stream)}
    groups = defaultdict(lambda: dict(seen=set(), units=set(), scored=set(), best={}))
    retained = max(10, top_n)
    def order(row):
        return direction * row['score'], row['molecule_id'], row['cid']
    with (output / 'block_scores.csv').open(encoding='utf-8', newline='') as stream:
        for row in csv.DictReader(stream):
            key = row['scheme'], row['block_id'], row['receptor']
            if key not in summaries or row['cid'] not in molecules:
                raise ValueError('Unknown block or conformer')
            group = groups[key]
            if row['cid'] in group['seen']:
                raise ValueError('Duplicate block conformer/receptor')
            group['seen'].add(row['cid'])
            mid = molecules[row['cid']]
            if not mid:
                raise ValueError('Missing molecule identity')
            unit = mid if ranking_unit == 'molecule' else row['cid']
            group['units'].add(unit)
            if row['status'] != 'ok':
                continue
            row['score'] = float(row['score'])
            if not math.isfinite(row['score']):
                raise ValueError('Nonfinite successful score')
            if score_name == 'EquiScore' and not 0 <= row['score'] <= 1:
                raise ValueError('EquiScore probability is outside [0,1]')
            row['molecule_id'] = mid
            group['scored'].add(unit)
            best = group['best']
            if unit not in best or order(row) < order(best[unit]):
                best[unit] = row
                if len(best) > retained:
                    del best[max(best, key=lambda k: order(best[k]))]
    ranking_rows, candidates = [], []
    for key, summary in sorted(summaries.items()):
        group = groups[key]
        if len(group['seen']) != int(summary['sampled']):
            raise ValueError('Sample count differs from analysis')
        best = sorted(group['best'].values(), key=order)
        for rank, row in enumerate(best, 1):
            candidates.append(dict(scheme=key[0], block_id=key[1], receptor=key[2],
                candidate_rank=rank, molecule_id=row['molecule_id'], cid=row['cid'],
                score=row['score'], pose_name=row['pose_name'], pose_file=row['pose_file'],
                pose_index=row['pose_index']))
        for n in sorted({5, 10, top_n}):
            enough = len(best) >= n
            complete = int(summary['missing']) == 0 and run['status'] == 'complete'
            ranking_rows.append(dict(scheme=key[0], block_id=key[1], receptor=key[2],
                top_n=n, rank=None, top_n_mean=math.fsum(r['score'] for r in best[:n])/n if enough else None,
                status='eligible' if enough and complete else 'insufficient_candidates' if not enough else 'incomplete_panel',
                sampled_conformers=int(summary['sampled']), scored_conformers=int(summary['scored']),
                sampled_units=len(group['units']), scored_units=len(group['scored']),
                missing=int(summary['missing']), population=int(summary['population'])))
    partitions = defaultdict(list)
    for row in ranking_rows:
        if row['status'] == 'eligible':
            partitions[row['scheme'], row['receptor'], row['top_n']].append(row)
    for rows in partitions.values():
        for rank, row in enumerate(sorted(rows, key=lambda r: (direction * r['top_n_mean'], r['block_id'])), 1):
            row['rank'] = rank
    dbpath = output / 'block_ranking.sqlite'
    with closing(sqlite3.connect(dbpath)) as db:
        db.executescript('''DROP TABLE IF EXISTS rankings; DROP TABLE IF EXISTS candidates;
            CREATE TABLE rankings(scheme TEXT, block_id TEXT, receptor TEXT, top_n INTEGER,
                rank INTEGER, top_n_mean REAL, status TEXT, sampled_conformers INTEGER,
                scored_conformers INTEGER, sampled_units INTEGER, scored_units INTEGER,
                missing INTEGER, population INTEGER, PRIMARY KEY(scheme,block_id,receptor,top_n));
            CREATE TABLE candidates(scheme TEXT, block_id TEXT, receptor TEXT, candidate_rank INTEGER,
                molecule_id TEXT, cid TEXT, score REAL, pose_name TEXT, pose_file TEXT, pose_index TEXT,
                PRIMARY KEY(scheme,block_id,receptor,candidate_rank));
            CREATE INDEX block_rank ON rankings(scheme,receptor,top_n,rank);''')
        db.executemany('INSERT INTO rankings VALUES(' + ','.join('?' * 13) + ')', [tuple(r.values()) for r in ranking_rows])
        db.executemany('INSERT INTO candidates VALUES(' + ','.join('?' * 10) + ')', [tuple(r.values()) for r in candidates])
        db.commit()
    for name, rows in [('block_rankings.csv', ranking_rows), ('block_top_candidates.csv', candidates)]:
        with (output / name).open('w', encoding='utf-8', newline='') as stream:
            if rows:
                writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
                writer.writeheader(); writer.writerows(rows)
    report['ranking'] = dict(**config, score=score_name, better=better,
        available_top_n=sorted({5, 10, top_n}), retained_candidates=retained,
        group_by=['scheme', 'receptor'], source_report_sha256=sha(source),
        implementation_sha256=sha(Path(__file__)), source_root=str(Path(run.get('source_report', source)).resolve().parent),
        ranked_blocks=sum(r['rank'] is not None and r['top_n'] == top_n for r in ranking_rows),
        leaders=[r for r in ranking_rows if r['rank'] == 1 and r['top_n'] == top_n],
        model_executed=bool(run.get('model_executed', False)))
    review = [f'# Top-{top_n} block prioritization', '',
        f'Primary score: mean of the best {top_n} {ranking_unit} scores per block and receptor; {better} {score_name} is better.',
        'The whole-block mean is diagnostic only and does not determine this ranking.',
        'Molecule mode selects the best in-block conformer per molecule. Ties use molecule/CID and block ID.',
        'Top-5 and Top-10 summaries are available. Incomplete panels and insufficient candidates are unranked.',
        'Rankings are within each scheme/receptor. Different sampling budgets preclude declaring a winning partition scheme.',
        'These are sampled-candidate priorities, not experimental affinity, enrichment or held-out recall.', '',
        '## Leading blocks', '']
    review += [f'- {r["scheme"]} / {r["receptor"]}: {r["block_id"]}, Top-{top_n} mean {r["top_n_mean"]:.6g}.' for r in report['ranking']['leaders']]
    (output / 'ranking_review.md').write_text('\n'.join(review) + '\n', encoding='utf-8')
    names = ['block_ranking.sqlite', 'block_rankings.csv', 'block_top_candidates.csv', 'ranking_review.md']
    report['output_hashes'].update({name: sha(output / name) for name in names})
    report.setdefault('outputs', {}).update({name: str((output / name).resolve()) for name in names})
    save(output / 'report.json', report)
    return report


def inspect(report_path, *, operation='ranks', scheme=None, receptor=None,
            block_id=None, top_n=None, limit=5):
    """Read bounded ranked blocks or their saved top candidates, with no inference."""
    if operation not in {'ranks', 'top'}:
        raise ValueError('Choose ranks or top')
    if type(limit) is not int or not 1 <= limit <= 100:
        raise ValueError('limit must be 1..100')
    path = Path(report_path)
    report = json.loads(path.read_text(encoding='utf-8'))
    ranking = report.get('ranking')
    if not ranking:
        raise ValueError('This older analysis has no Top-N ranking; queue a new analysis')
    n = ranking['top_n'] if top_n is None else top_n
    options(n, ranking['ranking_unit'])
    if n not in ranking['available_top_n']:
        raise ValueError('Queue analysis with this top_n first')
    from .project_context import ensure_within
    dbpath = ensure_within(path.parent / 'block_ranking.sqlite', path.parent)
    if sha(dbpath) != report['output_hashes']['block_ranking.sqlite']:
        raise ValueError('Ranking artifact changed')
    with readonly(dbpath) as db:
        db.row_factory = sqlite3.Row
        groups = [dict(r) for r in db.execute('SELECT DISTINCT scheme,receptor FROM rankings ORDER BY scheme,receptor')]
        if scheme is None or receptor is None:
            return dict(status='select_group', groups=groups, ranking=ranking,
                        message='Choose a scheme and receptor; scores are not pooled across receptors.')
        if dict(scheme=scheme, receptor=receptor) not in groups:
            raise ValueError('Unknown scheme/receptor')
        if operation == 'ranks':
            rows = [dict(r) for r in db.execute('SELECT * FROM rankings WHERE scheme=? AND receptor=? AND top_n=? AND rank IS NOT NULL ORDER BY rank LIMIT ?', (scheme, receptor, n, limit))]
            excluded = [dict(r) for r in db.execute('SELECT status,count(*) AS blocks FROM rankings WHERE scheme=? AND receptor=? AND top_n=? AND rank IS NULL GROUP BY status', (scheme, receptor, n))]
            return dict(status='complete', score=ranking['score'], better=ranking['better'], ranking_unit=ranking['ranking_unit'], top_n=n, blocks=rows, excluded=excluded)
        if not isinstance(block_id, str) or not block_id:
            raise ValueError('Choose a block_id from ranks')
        block = db.execute('SELECT * FROM rankings WHERE scheme=? AND receptor=? AND top_n=? AND block_id=?', (scheme, receptor, n, block_id)).fetchone()
        if block is None:
            raise ValueError('Unknown block')
        rows = [dict(r) for r in db.execute('SELECT * FROM candidates WHERE scheme=? AND receptor=? AND block_id=? ORDER BY candidate_rank LIMIT ?', (scheme, receptor, block_id, min(limit, ranking['retained_candidates'])))]
    return dict(status='complete', score=ranking['score'], better=ranking['better'], ranking_unit=ranking['ranking_unit'], block=dict(block), candidates=rows,
                source_root=ranking['source_root'], retained_candidates=ranking['retained_candidates'],
                message='Saved pose references only; no new docking or affinity inference.')


def main():
    import argparse
    from .block_plants import analyze
    from .prompt_workflow import file_lock
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', required=True, help='Completed PLANTS report.json')
    parser.add_argument('--output', required=True, help='Fresh analysis directory')
    parser.add_argument('--top-n', type=int, default=10)
    parser.add_argument('--ranking-unit', choices=['molecule', 'conformer'], default='molecule')
    args = parser.parse_args()
    config = options(args.top_n, args.ranking_unit)
    source, output = Path(args.source).resolve(), Path(args.output).resolve()
    if output.exists():
        raise ValueError('Use a fresh ranking output directory')
    with file_lock(source.parent.parent / ('.' + source.parent.name + '.e097.lock')):
        digest = sha(source)
        analyze(source, output)
        result = add_rankings(source, output, **config)
        if sha(source) != digest:
            raise ValueError('Docking report changed during analysis')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
