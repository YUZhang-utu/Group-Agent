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
            physical=float(np.mean(distances < self.cutoff-1e-8)) <= self.fraction
            excluded=any(region['mode']=='hard' and np.any(np.linalg.norm(moved-np.array(region['center']),axis=1)<region['radius'])
                         for region in self.design.get('exclusions',[]))
            answer.append(physical and not excluded)
        return np.asarray(answer,dtype=bool)


def adaptive_stage(candidate, rule):
    from .joint_coarse import extents
    if not rule['heavy_atoms'][0]<=len(candidate.shape_points)<=rule['heavy_atoms'][1]:return 0
    if any(np.count_nonzero(candidate.feature_types==int(t))<n for t,n in rule['feature_minimum'].items()):return 1
    extent=extents(candidate.shape_points)
    if np.any(extent<np.array(rule['extent_lower'])-1e-7) or np.any(extent>np.array(rule['extent_upper'])+1e-7):return 2
    return 3


def same_pose_gaussian(query, candidate, seeds):
    """Normalized shape/color Tanimoto for each actual seed, with cached self terms."""
    from .pose_acceleration import cross_overlaps
    from .gaussian_batch import _query_self_overlaps
    from .gaussian_overlay import gaussian_self_overlap
    sigma=1.;cutoff=None
    cross=cross_overlaps(query,candidate,seeds,sigma=sigma,cutoff=cutoff,backend='numpy')
    qself=_query_self_overlaps(query,sigma,cutoff)
    cs=gaussian_self_overlap(candidate.shape_points,sigma=sigma,cutoff=cutoff)
    cf=gaussian_self_overlap(candidate.feature_points,types=candidate.feature_types,sigma=sigma,cutoff=cutoff)
    shape=cross[:,0]/np.maximum(qself['shape']+cs-cross[:,0],1e-12)
    color=cross[:,1]/np.maximum(qself['color_unweighted']+cf-cross[:,1],1e-12)
    return np.clip(.5*(shape+color),0,1)


def pose_rank(values, assignments, order, gaussian, points, transform, design):
    weights=np.array([design.get('optional_weights',{}).get(a,0.) for a in order])
    optional=float(np.dot(weights,np.where(assignments>=0,values,0))/weights.sum()) if weights.sum() else 0.
    moved=points@transform[:3,:3].T+transform[:3,3]
    extensions={}
    if design.get('optional_groups') or design.get('spatial_groups'):
        from .contact_groups import grouped_terms
        extensions=grouped_terms(values,assignments,order,moved,design)
        families=design.get('optional_groups',[])
        total=weights.sum()+sum(g['weight'] for g in families)
        numerator=float(np.dot(weights,np.where(assignments>=0,values,0)))+sum(
            g['weight']*extensions['optional_group_scores'][g['id']] for g in families)
        optional=numerator/total if total else 0.
    if design.get('optional_normalization')=='fixed_budget':
        numerator=float(np.dot(weights,np.where(assignments>=0,values,0)))+sum(
            g['weight']*extensions['optional_group_scores'][g['id']] for g in design.get('optional_groups',[]))
        optional=numerator/design['optional_budget']
        extensions.update(optional_numerator=numerator,optional_denominator=design['optional_budget'],
                          optional_normalization='fixed_budget')
    penalty=sum(r['weight']*float(np.mean(np.linalg.norm(moved-np.array(r['center']),axis=1)<r['radius']))
                for r in design.get('exclusions',[]) if r['mode']=='soft')
    gw=design['gaussian_weight'];ow=design['optional_weight']
    sw=design.get('occupancy_weight',0.)
    combined=(gw*float(gaussian)+ow*optional+sw*extensions.get('occupancy_score',0.))/(gw+ow+sw)-penalty
    return dict(gaussian_same_pose=float(gaussian),optional_score=optional,exclusion_penalty=penalty,composite_score=combined,**extensions)
