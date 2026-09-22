"""Explicit anchor geometry and receptor-frame exclusion before Gaussian work."""
from pathlib import Path

import numpy as np

from .necessary_conditions import NecessaryConditions
from .pose_feasibility import possible_seed_mask

DEFAULTS = dict(minimum_score=.5,coarse_constraints=dict(heavy_atom_ratio=[.7,1.3],
    maximum_extent_distance=.35,minimum_feature_coverage=.7),
    pocket=dict(minimum_heavy_atom_distance=1.2,maximum_clashing_fraction=0.0))


def rule_passes(mask, anchor_order, design):
    present={a for i,a in enumerate(anchor_order) if mask & (1 << i)}
    return (set(design['mandatory_anchors']) <= present and
            all(present.intersection(group) for group in design.get('alternative_groups',[])))


class GuidedFilter:
    def __init__(self, query, expanded, design):
        from scipy.spatial import cKDTree
        self.design=design
        lookup={a['anchor_id']:a['feature_index'] for a in query['anchors']}
        self.bounds=[]
        if design['mandatory_anchors']:
            self.bounds.append(NecessaryConditions(expanded,[lookup[a] for a in design['mandatory_anchors']],
                                                    'all',design['minimum_score']))
        for group in design.get('alternative_groups',[]):
            self.bounds.append(NecessaryConditions(expanded,[lookup[a] for a in group],
                                                    'any',design['minimum_score']))
        with np.load(Path(design['receptor_npz']),allow_pickle=False) as z:
            points=z['points']
        if points.ndim!=2 or points.shape[1]!=3 or not len(points) or not np.isfinite(points).all():
            raise ValueError('Invalid receptor coordinates')
        self.tree=cKDTree(points)
        self.cutoff=design['pocket']['minimum_heavy_atom_distance']
        self.fraction=design['pocket']['maximum_clashing_fraction']
        self.allowed_reference_clash_fraction=None

    def geometry_possible(self, features):
        return all(bound.check(features)[0] for bound in self.bounds)

    def geometry_seeds(self, features, seeds, possible):
        result=possible.copy()
        for bound in self.bounds:
            result &= possible_seed_mask(bound,features,seeds)
        return result

    def pocket_mask(self, points, transforms):
        matrices=np.asarray(transforms).reshape(-1,4,4)
        answer=[]
        for matrix in matrices:
            moved=points @ matrix[:3,:3].T + matrix[:3,3]
            distances,_=self.tree.query(moved,k=1)
            answer.append(float(np.mean(distances < self.cutoff-1e-8)) <= self.fraction)
        return np.asarray(answer,dtype=bool)
