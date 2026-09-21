"""Optimistic same-seed rule checks before whole-ligand Gaussian evaluation."""
from __future__ import annotations

import time
import numpy as np

from .gaussian_batch import prepare_seeds


def possible_seed_mask(bound, candidate, seeds, *, batch_size=32):
    """A false row cannot satisfy the rule under any E031 assignment at that seed.

    Ignore assignment competition to obtain an upper bound per required anchor.
    Conditions are combined inside each pose, never across poses. This implements
    the current E031 sigma=1, cutoff=4.5, angular_power=2 protocol only.
    """
    if batch_size < 1: raise ValueError('Positive seed batch size required')
    qp,qt,qk = bound.points,bound.types,bound.kinds
    qd = bound.directions
    if qd is None:
        if np.any(qk): raise ValueError('Directional query needs directions')
        qd=np.zeros_like(qp)
    points=np.asarray(candidate.feature_points,dtype=np.float64)
    kinds=np.asarray(candidate.feature_kinds)
    directions=np.asarray(candidate.feature_directions,dtype=np.float64)
    if points.shape != directions.shape or not np.isfinite(directions).all():
        raise ValueError('Invalid candidate feature directions')
    if not np.isfinite(qd).all(): raise ValueError('Invalid query directions')
    typed=(qt[:,None]==candidate.feature_types[None,:]) & ((qk[:,None]==0)|(qk[:,None]==kinds[None,:]))
    result=np.zeros(len(seeds),dtype=bool)
    if not len(points): return result
    for start in range(0,len(seeds),batch_size):
        matrices=np.array([s.transform_matrix for s in seeds[start:start+batch_size]],dtype=np.float64).reshape(-1,4,4)
        rotations=matrices[:,:3,:3].transpose(0,2,1)
        cp=points@rotations+matrices[:,None,:3,3]
        cd=directions@rotations
        lengths=np.linalg.norm(cd,axis=2);valid=lengths>1e-12
        cd[valid]/=lengths[valid,None];cd[~valid]=0
        if not np.allclose(np.linalg.norm(cd[:,kinds!=0],axis=2),1.,atol=5e-3):
            raise ValueError('Invalid directional vectors')
        delta=qp[None,:,None,:]-cp[:,None,:,:]
        squared=np.einsum('bijk,bijk->bij',delta,delta)
        cosine=np.clip(qd[None]@cd.transpose(0,2,1),-1.,1.)
        agreement=np.ones_like(cosine)
        agreement=np.where(qk[None,:,None]==1,np.maximum(cosine,0.),agreement)
        agreement=np.where(qk[None,:,None]==2,np.abs(cosine),agreement)
        # Retain uncertainty near a cutoff/threshold instead of rejecting it.
        values=np.where(typed[None] & (squared <= 4.5**2+1e-8),
                        np.exp(-squared/2.)*agreement**2,0.)
        hits=values.max(axis=2) >= bound.threshold-1e-10
        result[start:start+len(matrices)]=hits.all(axis=1) if bound.mode=='all' else hits.any(axis=1)
    return result


def prepare_survivors(reader, query, bound, records, ids, invariant_keep, parameters):
    """Return a subset of conformers, with all original seeds cached for survivors."""
    keep=np.array(invariant_keep,dtype=bool,copy=True)
    prepared={}
    stats=dict(seed_generation_seconds=0.,pose_feasibility_seconds=0.,
               tested_seeds=0,possible_seeds=0,pose_feasibility_rejected=0)
    for position in np.flatnonzero(keep):
        gid=int(ids[position]);candidate=reader.get(gid);features=records[position]
        if candidate.molecule_id != features.molecule_id or candidate.conformer_id != features.conformer_id:
            raise ValueError('Pose feasibility artifact/chemical identity mismatch')
        if not np.array_equal(candidate.feature_types,features.feature_types) or not np.array_equal(candidate.feature_points,features.feature_points):
            raise ValueError('Pose feasibility artifact/chemical feature mismatch')
        started=time.perf_counter()
        seeds,pair_count=prepare_seeds(candidate,query,**parameters)
        stats['seed_generation_seconds']+=time.perf_counter()-started
        started=time.perf_counter()
        mask=possible_seed_mask(bound,features,seeds)
        stats['pose_feasibility_seconds']+=time.perf_counter()-started
        stats['tested_seeds']+=len(seeds);stats['possible_seeds']+=int(mask.sum())
        if not mask.any():
            keep[position]=False;stats['pose_feasibility_rejected']+=1
        else:
            # Never drop a seed from a surviving conformer's Gaussian competition.
            prepared[gid]=(candidate,seeds,pair_count)
    return keep,prepared,stats
