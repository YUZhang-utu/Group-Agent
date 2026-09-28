"""Bounded retrospective exploration of conformer blocks with overlapping molecules."""
import argparse
import csv
import json
import math
from pathlib import Path
import random

from .conformer_block_store import open_model
from .macrocycle_identity_audit import sha


def replay(model, scores, initial=100, fraction=.25, tail_fraction=.05,
           exploration=.2, top=100000, seed=20260928, max_conformers=1000000):
    if initial < 1 or top < 1 or not 0 < fraction <= 1 or not 0 < tail_fraction <= 1 or not 0 <= exploration <= 1:
        raise ValueError('Invalid replay settings')
    db, model_report = open_model(model)
    try:
        count = db.execute('SELECT count(*) FROM point').fetchone()[0]
        if not 0 < count <= max_conformers:
            raise ValueError('Replay panel is empty or exceeds explicit memory limit')
        members = {cid: (mid, block) for cid, mid, block in db.execute('SELECT cid,mid,node FROM point')}
    finally:
        db.close()
    values = {}
    with Path(scores).open(encoding='utf-8-sig', newline='') as stream:
        for row in csv.DictReader(stream):
            cid = row['conformer_id']
            if cid not in members or cid in values:
                raise ValueError('Duplicate or out-of-model score identity')
            raw = row['score'].strip()
            score = float(raw) if raw else None
            if score is not None and not math.isfinite(score):
                raise ValueError('Nonfinite reference score')
            if score is None and row.get('status') != 'no_surviving_pose':
                raise ValueError('Missing scores require no_surviving_pose status')
            values[cid] = score
    if len(values) != len(members):
        raise ValueError('Reference must cover every model conformer')
    blocks = {}
    for cid, (mid, block) in sorted(members.items()):
        blocks.setdefault(block, {}).setdefault(mid, []).append(cid)
    rng = random.Random(seed)
    pending = {}; selected = {}; found = set()
    for block in sorted(blocks):
        names = sorted(blocks[block]); rng.shuffle(names)
        selected[block] = names[:initial]; pending[block] = names[initial:]
        found.update(cid for mid in selected[block] for cid in blocks[block][mid])
    budget = math.ceil(len(members)*fraction)
    if len(found) > budget:
        raise ValueError('Initial unique-molecule block samples exceed conformer budget')
    initial_cost = len(found)

    def tails(block):
        # Each molecule contributes its best evaluated in-block conformer score.
        maxima = [max((values[cid] for cid in blocks[block][mid] if values[cid] is not None),
                      default=-math.inf) for mid in selected[block]]
        k = max(1, math.ceil(len(maxima)*tail_fraction))
        upper = sorted(maxima, reverse=True)[:k]
        return sum(upper)/len(upper), k, len(maxima)

    while len(found) < budget:
        remaining = budget-len(found)
        available = {block: [mid for mid in names if len(blocks[block][mid]) <= remaining]
                     for block, names in pending.items()}
        available = {block: names for block, names in available.items() if names}
        if not available:
            break
        options = sorted(available)
        if rng.random() < exploration:
            block = rng.choice(options)
        else:
            best = max(tails(block)[0] for block in options)
            block = rng.choice([block for block in options if tails(block)[0] == best])
        mid = available[block][0]
        pending[block].remove(mid); selected[block].append(mid)
        found.update(blocks[block][mid])

    def ranked(conformers):
        best = {}
        for cid in conformers:
            score = values[cid]
            if score is not None:
                mid = members[cid][0]; best[mid] = max(best.get(mid, -math.inf), score)
        return sorted(best, key=lambda mid: (-best[mid], mid))

    target = set(ranked(members)[:top])
    def metrics(ids):
        ranked_ids = set(ranked(ids)[:top])
        seen = {members[cid][0] for cid in ids}
        return dict(scored_conformers=len(ids), evaluated_unique_molecules=len(seen),
                    molecule_discovery_recall=len(target & seen)/len(target) if target else None,
                    molecule_ranking_recall=len(target & ranked_ids)/len(target) if target else None)
    uniform = set(random.Random(seed).sample(sorted(members), len(found)))
    tail_rows = []
    for block in sorted(blocks):
        mean, k, n = tails(block)
        tail_rows.append(dict(block_id=block, sampled_molecules=n, effective_k=k,
                              upper_tail_mean=mean if math.isfinite(mean) else None,
                              tail_contains_no_surviving_pose=not math.isfinite(mean)))
    return dict(model_id=model_report['model_id'], seed=seed, initial=initial,
        fraction=fraction, tail_fraction=tail_fraction, exploration_probability=exploration,
        reference_score_sha256=sha(Path(scores)), reference_top_molecules=len(target),
        panel_conformers=len(members), conformer_budget=budget, initial_conformer_cost=initial_cost,
        unused_budget=budget-len(found), adaptive=metrics(found), uniform_conformer=metrics(uniform),
        block_tails=tail_rows, scope='Retrospective complete supplied panel; larger scores are better',
        limitations=['No activity or out-of-panel recall claim',
                    'Exploration probability is not a guaranteed minimum per-block allocation',
                    'No uncertainty interval from one seed; repeat seeds and use held-out molecules',
                    'Explicit bounded in-memory replay; full source extraction and fitting are disk-backed'])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--model', type=Path, required=True)
    parser.add_argument('--scores', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--initial', type=int, choices=(100, 500), default=100)
    parser.add_argument('--fraction', type=float, default=.25)
    parser.add_argument('--tail-fraction', type=float, default=.05)
    parser.add_argument('--exploration', type=float, default=.2)
    parser.add_argument('--top', type=int, default=100000)
    parser.add_argument('--seed', type=int, default=20260928)
    parser.add_argument('--max-conformers', type=int, default=1000000)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError('Use a new output file')
    result = replay(args.model, args.scores, args.initial, args.fraction, args.tail_fraction,
                    args.exploration, args.top, args.seed, args.max_conformers)
    with args.output.open('x', encoding='utf-8') as stream:
        json.dump(result, stream, indent=2, allow_nan=False)


if __name__ == '__main__':
    main()
