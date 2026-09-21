"""Bounded seed batches for unchanged Gaussian overlaps; CUDA is opt-in."""
from __future__ import annotations

import os
import numpy as np


def hardware_options(backend=None, workers=None):
    backend = backend or os.environ.get('AIDD_POSE_BACKEND', 'numpy')
    if backend not in ('reference', 'numpy', 'cupy'):
        raise ValueError('AIDD_POSE_BACKEND must be reference, numpy or cupy')
    available = len(os.sched_getaffinity(0)) if hasattr(os, 'sched_getaffinity') else (os.cpu_count() or 1)
    default = 1 if backend == 'cupy' else max(1, min(24, available - 2))
    workers = int(workers if workers is not None else os.environ.get('AIDD_SCREEN_WORKERS', default))
    if workers < 1 or backend == 'cupy' and workers != 1:
        raise ValueError('Workers must be positive; CuPy uses exactly one GPU owner')
    return backend, workers


def array_module(backend):
    if backend == 'numpy':
        return np
    if backend != 'cupy':
        raise ValueError('Unknown accelerated pose backend')
    try:
        import cupy as cp
        if cp.cuda.runtime.getDeviceCount() < 1:
            raise RuntimeError('No CUDA device')
        return cp
    except (ImportError, RuntimeError) as exc:
        raise RuntimeError('CuPy/CUDA unavailable; install a compatible CuPy build or select numpy') from exc


def cross_overlaps(query, candidate, seeds, *, sigma, cutoff, backend='numpy', batch_size=32):
    """Return shape, color and anchored color for every original ordered seed.

    Float64 throughout; no TF32, reduced seed budget or changed cutoff. Output
    transfers synchronize CUDA. Temporary pair arrays are bounded by seed batches
    and 2048 query points, matching the reference reduction block size.
    """
    if batch_size < 1 or sigma <= 0 or not np.isfinite(sigma):
        raise ValueError('Invalid overlap batch size or sigma')
    if cutoff is not None and (cutoff <= 0 or not np.isfinite(cutoff)):
        raise ValueError('Invalid overlap cutoff')
    xp = array_module(backend)
    transforms = np.asarray([s.transform_matrix for s in seeds], dtype=np.float64).reshape(-1,4,4)
    result = np.empty((len(seeds),3), dtype=np.float64)
    qs, qf = [xp.asarray(query[k], dtype=xp.float64) for k in ('shape_points','feature_points')]
    cs, cf = [xp.asarray(p, dtype=xp.float64) for p in (candidate.shape_points,candidate.feature_points)]
    weights = xp.asarray(query['anchored_weights'], dtype=xp.float64)
    compatible = xp.asarray(query['feature_types'])[:,None] == xp.asarray(candidate.feature_types)[None,:]

    def kernel(a, moved, mask=None):
        delta = a[None,:,None,:] - moved[:,None,:,:]
        squared = xp.einsum('bijk,bijk->bij',delta,delta)
        valid = xp.ones(squared.shape, dtype=bool)
        if cutoff is not None:
            valid &= squared <= cutoff*cutoff
        if mask is not None:
            valid &= mask[None]
        return xp.where(valid, xp.exp(-squared/(2.0*sigma*sigma)), 0.0)

    # Bound allocation also for unusually large conformers; no geometry truncation.
    pairs = max(min(len(qs),2048)*len(cs), min(len(qf),2048)*len(cf), 1)
    batch_size = min(batch_size, max(1, 2_000_000//pairs))
    for start in range(0,len(seeds),batch_size):
        matrices = xp.asarray(transforms[start:start+batch_size])
        moved_shape = cs @ matrices[:,:3,:3].transpose(0,2,1) + matrices[:,None,:3,3]
        moved_features = cf @ matrices[:,:3,:3].transpose(0,2,1) + matrices[:,None,:3,3]
        values = xp.zeros((len(matrices),3), dtype=xp.float64)
        for offset in range(0,len(qs),2048):
            v = kernel(qs[offset:offset+2048],moved_shape)
            values[:,0] += v.reshape(len(matrices),-1).sum(axis=1)
        for offset in range(0,len(qf),2048):
            v = kernel(qf[offset:offset+2048],moved_features,compatible[offset:offset+2048])
            values[:,1] += v.reshape(len(matrices),-1).sum(axis=1)
            values[:,2] += (v*weights[None,offset:offset+2048,None]).reshape(len(matrices),-1).sum(axis=1)
        result[start:start+len(matrices)] = xp.asnumpy(values) if backend == 'cupy' else values
    if not np.isfinite(result).all():
        raise ValueError('Nonfinite accelerated overlaps')
    return result


def winner_indices(cross, query_self, candidate_self, *, contenders=False):
    """Stable first maximum, as in the reference strict-greater update."""
    denominators = np.array([.95*query_self['shape']+.05*candidate_self['shape'],
        .95*query_self['color_unweighted']+.05*candidate_self['color'],
        .95*query_self['color_anchored']+.05*candidate_self['color']])
    ratios = np.divide(cross,denominators,out=np.zeros_like(cross),where=denominators>0)
    objectives = np.column_stack((ratios[:,0],.5*(ratios[:,0]+ratios[:,1]),.5*(ratios[:,0]+ratios[:,2])))
    if contenders:
        # Re-evaluate near ties with the reference CPU kernel, including first ties.
        # This is a numerical guard, not a proof of CUDA equivalence on all inputs.
        tolerance = 1e-10*(1+np.abs(objectives).max(axis=0))
        return np.flatnonzero((objectives >= objectives.max(axis=0)-tolerance).any(axis=1))
    return objectives.argmax(axis=0)
