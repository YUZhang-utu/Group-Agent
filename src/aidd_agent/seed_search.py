"""Explicit experimental seed budgets with bounded physical-mask batches."""
from itertools import chain, islice
import time

import numpy as np

from .gaussian_overlay import iter_pair_alignment_seeds, principal_axis_seeds


def policy(value=None):
    value = {} if value is None else dict(value)
    if set(value)-{'max_pair_seeds','survivor_target','seed_batch','backend'}:
        raise ValueError('Unknown seed search fields')
    result = dict(max_pair_seeds=512, survivor_target=0, seed_batch=64, backend='reference')
    result.update(value)
    for key in ('max_pair_seeds','survivor_target','seed_batch'):
        number = result[key]
        if type(number) is not int or number < (1 if key == 'seed_batch' else 0):
            raise ValueError(f'Invalid seed search {key}')
    if result['backend'] not in ('reference','batched'):
        raise ValueError('Unknown seed backend')
    return result


def prepare(candidate, query, guided, settings, pose_parameters):
    """Return the generated prefix and its exact physical mask, including PCA.

    The survivor target counts seeds, not distinct poses. No target changes the
    generated cap. Targets cannot cause unbounded search of a clashing conformer.
    """
    settings = policy(settings)
    if guided is None:
        raise ValueError('Experimental seed search requires a physical pocket mask')
    started = time.perf_counter()
    pca = principal_axis_seeds(candidate.shape_points, query['shape_points'])
    indices = query['anchor_feature_indices']
    pairs = iter_pair_alignment_seeds(candidate.feature_points, candidate.feature_types,
        query['feature_points'][indices], query['feature_types'][indices],
        tolerance=pose_parameters['pair_tolerance'], axial_samples=pose_parameters['axial_samples'],
        max_seeds=settings['max_pair_seeds'], backend=settings['backend'])
    stream = chain(pca,pairs)
    generation = time.perf_counter()-started
    pocket_seconds = 0.
    seeds, masks = [], []
    survived = 0
    target = settings['survivor_target']
    while not target or survived < target:
        width = min(settings['seed_batch'], target-survived) if target else settings['seed_batch']
        started = time.perf_counter()
        batch = list(islice(stream,width))
        generation += time.perf_counter()-started
        if not batch:
            break
        started = time.perf_counter()
        matrices = np.asarray([s.transform_matrix for s in batch]).reshape(-1,4,4)
        mask = np.asarray(guided.pocket_mask(candidate.shape_points,matrices),dtype=bool)
        if mask.shape != (len(batch),):
            raise ValueError('Physical mask shape mismatch')
        pocket_seconds += time.perf_counter()-started
        seeds.extend(batch); masks.extend(mask.tolist()); survived += int(mask.sum())
    pair_count = sum(s.axial_sample >= 0 for s in seeds)
    reason = ('survivor_target' if target and survived >= target else
              'generated_cap' if pair_count >= settings['max_pair_seeds'] else 'exhausted_pairs')
    stats = dict(generated=len(seeds),pair_seeds=pair_count,surviving=survived,
                 stop_reason=reason,target_reached=bool(target and survived >= target))
    return seeds,np.asarray(masks,dtype=bool),stats,dict(seed_generation=generation,pocket_exclusion=pocket_seconds)
