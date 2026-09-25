"""Replay adaptive block sampling on a fully evaluated molecular panel."""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import random


def prepare(rows):
    members = {}
    for row in rows:
        mid, block = row['molecule_id'], row['block_id']
        if not mid or not block or mid in members:
            raise ValueError('Require one row per unique molecule and a nonempty block ID')
        raw = row['score'].strip()
        score = float(raw) if raw else None
        if score is not None and not math.isfinite(score):
            raise ValueError('Nonfinite score')
        if score is None and row.get('status') != 'no_surviving_pose':
            raise ValueError('Blank scores require explicit no_surviving_pose status')
        rank = int(row['reference_rank']) if row.get('reference_rank', '').strip() else None
        if rank is not None and (rank <= 0 or score is None):
            raise ValueError('Invalid reference rank')
        members[mid] = dict(block=block, score=score, rank=rank)
    if not members:
        raise ValueError('Empty panel')
    valid = [mid for mid, row in members.items() if row['score'] is not None]
    ranked = [mid for mid in valid if members[mid]['rank'] is not None]
    if ranked and len(ranked) != len(valid):
        raise ValueError('Provide reference ranks for all scored molecules or none')
    if len({members[mid]['rank'] for mid in ranked}) != len(ranked):
        raise ValueError('Duplicate reference ranks')
    reference = sorted(valid, key=(lambda mid: members[mid]['rank']) if ranked else
                       (lambda mid: (-members[mid]['score'], mid)))
    return members, reference, 'supplied_full_run_rank' if ranked else 'maximum_contact_score'


def replay(rows, initial=32, top_k=5, batch=32, fraction=.25,
           exploration=.2, top=100, seed=20260924):
    if min(initial, top_k, batch, top) <= 0 or not 0 < fraction <= 1 or not 0 <= exploration <= 1:
        raise ValueError('Invalid sampling budgets')
    members, reference, objective = prepare(rows)
    blocks = {}
    for mid in sorted(members):
        blocks.setdefault(members[mid]['block'], []).append(mid)
    rng = random.Random(seed)
    for block in sorted(blocks):
        rng.shuffle(blocks[block])
    budget = math.ceil(len(members) * fraction)
    initial_cost = sum(min(initial, len(ids)) for ids in blocks.values())
    if initial_cost > budget:
        raise ValueError('Initial per-block samples exceed total budget; lower initial or raise fraction')
    sampled = {block: ids[:initial] for block, ids in blocks.items()}
    pending = {block: ids[initial:] for block, ids in blocks.items()}
    used = initial_cost
    while used < budget:
        available = sorted(block for block, ids in pending.items() if ids)
        if rng.random() < exploration:
            chosen = rng.choice(available)
        else:
            def priority(block):
                scores = sorted((members[mid]['score'] if members[mid]['score'] is not None else -math.inf
                                 for mid in sampled[block]), reverse=True)[:top_k]
                return sum(scores) / len(scores)
            best = max(priority(block) for block in available)
            chosen = rng.choice([block for block in available if priority(block) == best])
        count = min(batch, budget-used, len(pending[chosen]))
        sampled[chosen].extend(pending[chosen][:count])
        del pending[chosen][:count]
        used += count
    found = {mid for ids in sampled.values() for mid in ids}
    uniform = set(random.Random(seed).sample(sorted(members), budget))
    target = set(reference[:top])
    def metrics(selected):
        return dict(evaluated_molecules=len(selected), recovered_top=len(target & selected),
                    reference_top_recall=len(target & selected)/len(target) if target else None,
                    missed_top_ids=sorted(target-selected))
    return dict(seed=seed, panel_molecules=len(members), blocks=len(blocks),
                reference_objective=objective, reference_top_size=len(target),
                initial_evaluations=initial_cost, adaptive=metrics(found), uniform=metrics(uniform),
                block_allocations={block:len(ids) for block, ids in sorted(sampled.items())})


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--csv', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--initial', type=int, default=32)
    parser.add_argument('--top-k', type=int, default=5)
    parser.add_argument('--batch', type=int, default=32)
    parser.add_argument('--fraction', type=float, default=.25)
    parser.add_argument('--exploration', type=float, default=.2)
    parser.add_argument('--top', type=int, default=100)
    parser.add_argument('--seed', type=int, default=20260924)
    parser.add_argument('--repeats', type=int, default=5)
    args = parser.parse_args()
    if args.repeats < 1:
        parser.error('Positive repeat count required')
    raw = args.csv.read_bytes()
    rows = list(csv.DictReader(raw.decode('utf-8-sig').splitlines()))
    settings = {key:getattr(args,key) for key in
                ('initial','top_k','batch','fraction','exploration','top')}
    results = [replay(rows, **settings, seed=args.seed+i) for i in range(args.repeats)]
    report = dict(kind='retrospective_block_exploration', settings=settings, results=results,
                  source_sha256=hashlib.sha256(raw).hexdigest(),
                  implementation_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  scope='Fully evaluated input panel only; no live scheduling, ANN recall, activity or wall-time claim')
    with args.output.open('x', encoding='utf-8') as stream:
        json.dump(report, stream, indent=2, allow_nan=False)
    print(args.output)


if __name__ == '__main__':
    main()
