"""Export a bounded, source-sealed review of existing EquiScore rankings. No inference."""
import argparse
from collections import defaultdict, Counter
import csv
import hashlib
import json
import math
from pathlib import Path
import zipfile

DEFAULT = '/mnt/local/hand/yuzhang/aidd/e097-chat-workspace/users/workstation/projects/prj-adf8a1b9f4f8-prompt-aidd/runs/PROMPT-1d4f82d50d294925/execution/equiscore-v3/equiscore-analysis'


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for part in iter(lambda: stream.read(1024 * 1024), b''): digest.update(part)
    return digest.hexdigest()


def rows(path):
    with path.open(encoding='utf-8', newline='') as stream: return list(csv.DictReader(stream))


def write_csv(path, data):
    if not data: return
    with path.open('w', encoding='utf-8', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(data[0]))
        writer.writeheader(); writer.writerows(data)


def review(analysis, output, limit=10):
    if limit not in (5, 10): raise ValueError('Choose 5 or 10 leading blocks per group')
    analysis, output = Path(analysis).resolve(), Path(output).resolve()
    if output.exists() or output.with_suffix('.zip').exists(): raise ValueError('Use a new review output directory')
    if output == analysis or analysis in output.parents: raise ValueError('Keep review outputs outside the sealed analysis')
    report_path = analysis / 'report.json'
    report_hash = sha(report_path)
    report = json.loads(report_path.read_text(encoding='utf-8'))
    ranking = report.get('ranking', {})
    if report.get('status') != 'complete' or ranking.get('score') != 'EquiScore' or ranking.get('better') != 'higher' or ranking.get('ranking_unit') != 'molecule' or ranking.get('top_n') != 10:
        raise ValueError('Requires a completed EquiScore Top-10 distinct-molecule analysis')
    source = Path(report['source_report'])
    if sha(source) != ranking['source_report_sha256']: raise ValueError('Full scoring report changed')
    full = json.loads(source.read_text(encoding='utf-8'))
    if full.get('scope') != 'full' or full.get('status') != 'complete' or full.get('failed_pairs') != 0 or full.get('pairs', 0) <= 0 or full.get('pairs') != full.get('scored') or full.get('worker_exit_code') != 0:
        raise ValueError('Requires a successful full scoring receipt')
    names = ('block_rankings.csv', 'block_top_candidates.csv', 'rank_comparison.csv')
    seals = {}
    for name in names:
        seals[name] = sha(analysis / name)
        if seals[name] != report['output_hashes'].get(name): raise ValueError('Missing or changed sealed artifact: ' + name)
    rankings, candidates, comparison = [rows(analysis / name) for name in names]
    key = lambda r: (r['scheme'], r['receptor'], r['block_id'])
    comparisons = {key(r): r for r in comparison}
    if len(comparisons) != len(comparison): raise ValueError('Duplicate comparison record')
    by_block = defaultdict(list)
    for row in candidates: by_block[key(row)].append(row)
    top5 = {key(r): r for r in rankings if int(r['top_n']) == 5}
    groups = defaultdict(list); excluded = defaultdict(Counter)
    primary = [r for r in rankings if int(r['top_n']) == 10]
    if len({key(r) for r in primary}) != len(primary): raise ValueError('Duplicate primary block record')
    for row in primary:
        group = row['scheme'], row['receptor']
        if row['status'] == 'eligible' and row['rank']: groups[group].append(row)
        else: excluded[group][row['status']] += 1
    if sum(map(len, groups.values())) != ranking['ranked_blocks']: raise ValueError('Ranked block count differs from receipt')
    leaders, selected_candidates, summaries = [], [], []
    markdown = ['# EquiScore leading-block review', '',
        'Primary ranking: best 10 distinct molecules per block, each represented by its best in-block conformer. Higher EquiScore is better.',
        'Top-5/Top-10 leading blocks are display limits; they do not change the primary Top-10 mean. Top-5-mean ranks are shown as a sensitivity check.',
        'Groups remain separate. No experimental enrichment, affinity, pose-quality or partition-superiority claim is made.', '']
    for group in sorted(set(groups) | set(excluded)):
        ordered = sorted(groups[group], key=lambda r: int(r['rank']))
        if [int(r['rank']) for r in ordered] != list(range(1, len(ordered)+1)): raise ValueError('Non-contiguous ranks')
        means = [float(r['top_n_mean']) for r in ordered]
        if any(not math.isfinite(v) or not 0 <= v <= 1 for v in means) or means != sorted(means, reverse=True): raise ValueError('Invalid score ordering')
        lead = ordered[:limit]
        leading5 = {r['block_id'] for r in ordered[:5]}
        secondary5 = {r['block_id'] for r in top5.values() if (r['scheme'], r['receptor']) == group and r['rank'] and int(r['rank']) <= 5}
        summary = dict(scheme=group[0], receptor=group[1], ranked_blocks=len(ordered), excluded=dict(excluded[group]),
            top5_block_overlap_between_top10_and_top5_means=len(leading5 & secondary5),
            top1_minus_top2=means[0]-means[1] if len(means)>1 else None)
        summaries.append(summary)
        markdown += ['## ' + ' / '.join(group), '',
            '| Rank | Block | Top-10 mean | Top-5-mean rank | ChemPLP rank | Shared top-10 molecules |',
            '|---:|---|---:|---:|---:|---:|']
        for row in lead:
            block_key = key(row)
            values = sorted(by_block[block_key], key=lambda r: int(r['candidate_rank']))[:10]
            scores = [float(r['score']) for r in values]
            if len(values) != 10 or len({r['molecule_id'] for r in values}) != 10 or [int(r['candidate_rank']) for r in values] != list(range(1,11)):
                raise ValueError('Leading block does not contain 10 distinct ranked molecules')
            if any(not math.isfinite(v) or not 0 <= v <= 1 for v in scores) or scores != sorted(scores, reverse=True): raise ValueError('Invalid candidate ordering')
            if not math.isclose(math.fsum(scores)/10, float(row['top_n_mean']), rel_tol=0, abs_tol=1e-12): raise ValueError('Candidate mean does not reproduce block mean')
            comp = comparisons[block_key]
            if comp['equiscore_rank'] != row['rank'] or not math.isclose(float(comp['equiscore_top10_mean']), float(row['top_n_mean']), rel_tol=0, abs_tol=1e-12):
                raise ValueError('Comparison differs from ranking')
            second = top5[block_key]
            item = dict(row, top5_mean_rank=second['rank'], top5_mean=second['top_n_mean'],
                chemplp_rank=comp['chemplp_rank'], chemplp_top10_mean=comp['chemplp_top10_mean'],
                rank_improvement=comp['rank_improvement'], shared_top10_molecules=comp['shared_top10_molecules'],
                exact_equal_top10_scores=len(set(scores)) == 1)
            leaders.append(item)
            selected_candidates.extend(dict(r, block_rank=row['rank'], inspect_first=int(r['candidate_rank']) <= 5) for r in values)
            markdown.append('| {rank} | {block_id} | {top_n_mean} | {top5_mean_rank} | {chemplp_rank} | {shared_top10_molecules} |'.format(**item))
        markdown += ['', 'Top-5 block overlap between the two within-block summaries: '+str(summary['top5_block_overlap_between_top10_and_top5_means'])+'.',
            'Excluded blocks: '+json.dumps(summary['excluded'],sort_keys=True)+'.', '']
    if sha(report_path) != report_hash or sha(source) != ranking['source_report_sha256'] or any(sha(analysis / n) != s for n,s in seals.items()):
        raise ValueError('Inputs changed during review')
    output.mkdir(parents=True)
    write_csv(output / 'leading_blocks.csv', leaders)
    write_csv(output / 'leading_block_candidates.csv', selected_candidates)
    (output / 'review.md').write_text('\n'.join(markdown),encoding='utf-8')
    result = dict(status='complete', operation='read_only_leader_review', analysis_report=str(report_path),
        analysis_sha256=report_hash, input_hashes=seals, blocks_per_group=limit, primary_top_n=10,
        groups=summaries, selected_blocks=len(leaders), selected_candidate_rows=len(selected_candidates),
        model_executed=False, docking_executed=False, pose_geometry_inspected=False,
        limitations=['Saved numeric priorities only; original pose references retained, no pose extraction or geometry assessment.',
                     'Agreement and rank shifts are descriptive, not experimental enrichment or predictive validation.'])
    result['output_hashes']={p.name:sha(p) for p in output.iterdir() if p.is_file()}
    (output / 'report.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    with zipfile.ZipFile(output.with_suffix('.zip'),'x',compression=zipfile.ZIP_DEFLATED) as archive:
        for path in sorted(output.iterdir()): archive.write(path,path.name)
    return result


if __name__ == '__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--analysis',default=DEFAULT)
    parser.add_argument('--output',required=True)
    parser.add_argument('--blocks',type=int,choices=(5,10),default=10)
    args=parser.parse_args()
    result=review(args.analysis,args.output,args.blocks)
    print(json.dumps({k:result[k] for k in ('status','selected_blocks','selected_candidate_rows','groups')},indent=2))
    print('Review bundle: '+str(Path(args.output).resolve().with_suffix('.zip')))
