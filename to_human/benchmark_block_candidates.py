"""Read-only E107 membership audit; writes only a new user-selected output folder.

Usage: python benchmark_block_candidates.py --runs-root PATH --output NEW_DIRECTORY
No third-party dependencies, model calls, docking, or search execution.
"""
import argparse
import csv
import hashlib
import itertools
import json
from pathlib import Path

RUNS = ('5f7c28e4bd704616 3a90516120f74140 4f0fad82c9714b1e '
        '34bf4cc4a192449e 7a7fddb92b744a38 9bd13fe58ad14df0 '
        'dcc64417cf534000 317795cf58ad4f04 803b293dd57c4aec '
        '82f400cb6c1c46df c71019343cb04b62 bfa50ed8eb9a489e').split()


def read(path):
    return json.loads(path.read_text(encoding='utf-8'))


def load_arm(folder):
    report = read(folder / 'report.json')
    if report.get('kind') != 'block_search_dock' or report.get('status') != 'complete':
        raise ValueError(f'Not a completed search: {folder}')
    if report.get('docking_requested') is not False:
        raise ValueError(f'Not an explicitly search-only arm: {folder}')
    if report.get('search_engine') != 'consensus-threshold-v1':
        raise ValueError(f'Unexpected search engine: {folder}')
    selected = folder / 'selected'
    selection_report = read(selected / 'report.json')
    raw = (selected / 'candidates.jsonl').read_bytes()
    digest = hashlib.sha256(raw).hexdigest()
    if digest != selection_report.get('output_hashes', {}).get('candidates.jsonl'):
        raise ValueError(f'Candidate hash mismatch: {folder}')
    mapping = {}
    for line in raw.decode('utf-8').splitlines():
        row = json.loads(line)
        cid, mid = row['cid'], row['molecule_id']
        if not isinstance(cid, str) or not cid or not isinstance(mid, str) or not mid:
            raise ValueError('Missing or invalid candidate identity')
        if cid in mapping:
            raise ValueError(f'Duplicate conformer: {cid}')
        mapping[cid] = mid
    molecules = set(mapping.values())
    if len(mapping) != report['exported_conformers'] or len(molecules) != report['unique_molecules']:
        raise ValueError(f'Candidate counts disagree with report: {folder}')
    s = report['selection']
    name = f"{s['scheme']}/{s['method']}/Top-{s['blocks_per_method']}"
    blocks = report.get('searched_blocks', [])
    summary = dict(arm=name, searched=report['searched_conformers'],
                   ranked=report['ranked_conformers'], exported=len(mapping),
                   molecules=len(molecules), shortfall=report['shortfall'],
                   blocks_searched=len(blocks), seconds=report.get('elapsed_seconds'),
                   cached_blocks=sum(b.get('cache_reused') is True for b in blocks),
                   survival_fraction=report['ranked_conformers']/report['searched_conformers']
                   if report['searched_conformers'] else None,
                   conformers_per_molecule=len(mapping)/len(molecules) if molecules else None,
                   candidate_sha256=digest, report_path=str(folder/'report.json'))
    return summary, mapping, molecules


def overlap(a, b):
    intersection = len(a & b)
    return dict(intersection=intersection, union=len(a | b),
                jaccard=intersection/len(a | b) if a | b else None,
                fraction_of_a=intersection/len(a) if a else None,
                fraction_of_b=intersection/len(b) if b else None,
                identical=a == b)


def write_csv(path, rows):
    with path.open('w', encoding='utf-8', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--runs-root', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError('Use a new output directory; existing outputs are preserved')
    arms = [load_arm(args.runs_root / ('PROMPT-'+run) / 'execution/campaign/blocks') for run in RUNS]
    if len({a[0]['arm'] for a in arms}) != 12:
        raise ValueError('Expected twelve different arms')
    global_mapping = {}
    for _, mapping, _ in arms:
        for cid, mid in mapping.items():
            if global_mapping.setdefault(cid, mid) != mid:
                raise ValueError(f'Inconsistent cross-arm molecule identity: {cid}')
    pairs = []
    for a, b in itertools.combinations(arms, 2):
        for unit, sa, sb in [('conformer', set(a[1]), set(b[1])), ('molecule', a[2], b[2])]:
            pairs.append(dict(arm_a=a[0]['arm'], arm_b=b[0]['arm'], unit=unit, **overlap(sa, sb)))
    args.output.mkdir(parents=True, exist_ok=False)
    write_csv(args.output/'arm_summary.csv', [a[0] for a in arms])
    write_csv(args.output/'pairwise_overlap.csv', pairs)
    (args.output/'audit.json').write_text(json.dumps(dict(
        status='complete', scope='candidate_identity_and_membership_only', arms=12,
        candidate_hashes_verified=True, full_query_and_library_identity_audit=False,
        biological_enrichment_evaluated=False, scientific_computation_launched=False,
        limitations=['Set overlap is not activity or structural diversity.',
                     'Cache timings are not independent cold-run benchmarks.',
                     'No verification of all original MOL2 file contents in this audit.']), indent=2), encoding='utf-8')
    print('Saved arm_summary.csv, pairwise_overlap.csv and audit.json to', args.output)


if __name__ == '__main__':
    main()
