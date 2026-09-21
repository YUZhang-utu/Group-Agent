"""Conservative necessary bounds for positive typed rigid-pose feature matches."""
from __future__ import annotations

import numpy as np


def distinct_assignment(domains):
    """Whether the domain graph admits an injective assignment (not a pose)."""
    owners = {}
    def augment(i, seen):
        for j in np.flatnonzero(domains[i]):
            j = int(j)
            if j in seen: continue
            seen.add(j)
            if j not in owners or augment(owners[j], seen):
                owners[j] = i
                return True
        return False
    return all(augment(i, set()) for i in sorted(range(len(domains)), key=lambda k: domains[k].sum()))


class NecessaryConditions:
    """Reject only when a requested rigid feature condition is impossible."""
    def __init__(self, query, indices, mode, threshold, *, sigma=1., cutoff=4.5):
        if mode not in ('all','any') or not 0 < threshold <= 1 or sigma <= 0:
            raise ValueError('Invalid necessary-condition rule')
        indices = sorted(set(indices))
        if not indices: raise ValueError('No requested features')
        self.points = np.asarray(query['feature_points'][indices], dtype=float)
        self.types = np.asarray(query['feature_types'][indices])
        self.kinds = np.asarray(query['feature_direction_kinds'][indices])
        if not np.isfinite(self.points).all(): raise ValueError('Nonfinite query features')
        self.distances = np.linalg.norm(self.points[:,None]-self.points[None,:],axis=2)
        radius = sigma * np.sqrt(-2*np.log(threshold))
        if cutoff is not None: radius = min(radius,cutoff)
        self.pair_slack = 2*radius
        self.mode = mode

    def check(self, candidate):
        points = np.asarray(candidate.feature_points,dtype=float)
        if points.ndim != 2 or points.shape[1] != 3 or not np.isfinite(points).all():
            raise ValueError('Invalid candidate features')
        types, kinds = np.asarray(candidate.feature_types), np.asarray(candidate.feature_kinds)
        if types.shape != (len(points),) or kinds.shape != types.shape:
            raise ValueError('Invalid feature metadata')
        domains = (self.types[:,None] == types[None,:]) & (
            (self.kinds[:,None] == 0) | (self.kinds[:,None] == kinds[None,:]))
        if self.mode == 'any': return (bool(domains.any()), 'possible' if domains.any() else 'feature_compatibility')
        if not domains.any(axis=1).all(): return False, 'feature_compatibility'
        if not distinct_assignment(domains): return False, 'distinct_assignment'
        if len(domains) == 1: return True, 'possible'
        distances = np.linalg.norm(points[:,None]-points[None,:],axis=2)
        changed = True
        while changed:
            changed = False
            for i in range(len(domains)):
                for j in range(len(domains)):
                    if i == j: continue
                    # An expanded margin makes floating-point uncertainty retain a row.
                    margin = 1e-5*(1+self.distances[i,j]+distances)
                    compatible = np.abs(distances-self.distances[i,j]) <= self.pair_slack+margin
                    np.fill_diagonal(compatible,False)
                    supported = (compatible & domains[j][None,:]).any(axis=1)
                    reduced = domains[i] & supported
                    if not reduced.any(): return False, 'pair_distance'
                    if not np.array_equal(reduced,domains[i]):
                        domains[i] = reduced
                        changed = True
        if not distinct_assignment(domains): return False, 'distinct_assignment'
        return True, 'possible'


def pose_mask(arrays, columns, policy, objective):
    hits = (arrays[objective+'__anchor_scores'][:,columns] >= policy['minimum_score']) & (
        arrays[objective+'__anchor_assignments'][:,columns] >= 0)
    return hits.all(axis=1) if policy['match_mode']=='all' else hits.any(axis=1)
