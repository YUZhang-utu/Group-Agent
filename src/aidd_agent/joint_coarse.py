"""Explicit pose-independent eligibility; not bounds on Gaussian scores."""
from __future__ import annotations

import numpy as np


def validate_constraints(rule):
    fields = {'heavy_atom_ratio', 'maximum_extent_distance', 'minimum_feature_coverage'}
    if not isinstance(rule, dict) or set(rule) != fields:
        raise ValueError('Joint coarse constraints require heavy_atom_ratio, maximum_extent_distance and minimum_feature_coverage')
    ratio = rule['heavy_atom_ratio']
    def number(x):
        return type(x) in (int, float) and np.isfinite(x)
    if not isinstance(ratio, list) or len(ratio) != 2 or not all(number(x) for x in ratio) or not 0 < ratio[0] <= ratio[1]:
        raise ValueError('heavy_atom_ratio needs positive inclusive [lower, upper] bounds')
    if not number(rule['maximum_extent_distance']) or rule['maximum_extent_distance'] < 0:
        raise ValueError('maximum_extent_distance must be finite and nonnegative')
    if not number(rule['minimum_feature_coverage']) or not 0 < rule['minimum_feature_coverage'] <= 1:
        raise ValueError('minimum_feature_coverage must be in (0,1]')
    return rule


def extents(points):
    points = np.asarray(points, dtype=float)
    if points.ndim != 2 or points.shape[1] != 3 or not len(points) or not np.isfinite(points).all():
        raise ValueError('Finite nonempty heavy-atom coordinates required')
    centered = points - points.mean(axis=0)
    return np.sqrt(np.maximum(np.linalg.eigvalsh(centered.T @ centered / len(points)), 0))


class JointCoarse:
    """Size AND principal extents AND query typed-feature count coverage.

    Extent distance is ||candidate_extents-query_extents|| / ||query_extents||.
    Coverage is sum_type min(candidate_count, query_count) / query_feature_count.
    Neither quantity is a Gaussian score or a candidate-protein contact score.
    """
    def __init__(self, query, rule):
        self.rule = validate_constraints(rule)
        self.count = len(query['shape_points'])
        self.extent = extents(query['shape_points'])
        self.norm = np.linalg.norm(self.extent)
        if self.norm <= 1e-12:
            raise ValueError('Query extent is degenerate')
        self.types, self.counts = np.unique(query['feature_types'], return_counts=True)
        if not self.counts.sum():
            raise ValueError('Query chemical features required')

    def check(self, candidate):
        ratio = len(candidate.shape_points) / self.count
        low, high = self.rule['heavy_atom_ratio']
        if not low <= ratio <= high:
            return False, 'joint_heavy_atom_ratio'
        types, counts = np.unique(candidate.feature_types, return_counts=True)
        table = dict(zip(types.tolist(), counts.tolist()))
        coverage = sum(min(int(n), table.get(t, 0)) for t, n in zip(self.types, self.counts)) / self.counts.sum()
        if coverage < self.rule['minimum_feature_coverage']:
            return False, 'joint_feature_coverage'
        distance = np.linalg.norm(extents(candidate.shape_points) - self.extent) / self.norm
        # A fixed numerical guard retains geometrically ambiguous boundary rows.
        if distance > self.rule['maximum_extent_distance'] + 1e-7:
            return False, 'joint_shape_extent'
        return True, 'possible'


def eligibility(query, reader, ids, rule):
    bound = JointCoarse(query, rule)
    return np.asarray([bound.check(reader.get(int(gid)))[0] for gid in ids], dtype=bool)
