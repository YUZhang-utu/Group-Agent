from __future__ import annotations

from dataclasses import asdict, dataclass
import math
from typing import Iterable, Sequence

import numpy as np


@dataclass(frozen=True)
class RigidSeed:
    seed_id: str
    candidate_pair: tuple[int, int]
    query_pair: tuple[int, int]
    axial_sample: int
    transform_matrix: tuple[float, ...]

    def __post_init__(self) -> None:
        matrix = np.asarray(self.transform_matrix, dtype=np.float64)
        if matrix.shape != (16,) or not np.isfinite(matrix).all():
            raise ValueError("seed requires a finite row-major 4x4 transform")
        rotation = matrix.reshape(4, 4)[:3, :3]
        if not np.isclose(np.linalg.det(rotation), 1.0, atol=1e-6):
            raise ValueError("seed rotation must be proper")


def _points(value: Sequence[Sequence[float]], name: str) -> np.ndarray:
    result = np.asarray(value, dtype=np.float64)
    if result.ndim != 2 or result.shape[1:] != (3,) or not np.isfinite(result).all():
        raise ValueError(f"{name} must have finite shape (N,3)")
    return result


def _weights(value: Sequence[float] | None, count: int, name: str) -> np.ndarray:
    result = np.ones(count, dtype=np.float64) if value is None else np.asarray(
        value, dtype=np.float64)
    if result.shape != (count,) or not np.isfinite(result).all() or np.any(result < 0):
        raise ValueError(f"{name} must contain one finite non-negative value per point")
    return result


def apply_transform(points: Sequence[Sequence[float]],
                    transform_matrix: Sequence[float]) -> np.ndarray:
    coordinates = _points(points, "points")
    matrix = np.asarray(transform_matrix, dtype=np.float64)
    if matrix.size != 16:
        raise ValueError("transform must contain 16 values")
    matrix = matrix.reshape(4, 4)
    if not np.isfinite(matrix).all() or not np.allclose(matrix[3], [0, 0, 0, 1]):
        raise ValueError("transform must be a finite homogeneous 4x4 matrix")
    return coordinates @ matrix[:3, :3].T + matrix[:3, 3]


def gaussian_overlap(points_a: Sequence[Sequence[float]],
                     points_b: Sequence[Sequence[float]], *, sigma: float = 1.0,
                     weights_a: Sequence[float] | None = None,
                     weights_b: Sequence[float] | None = None,
                     types_a: Sequence[int] | None = None,
                     types_b: Sequence[int] | None = None,
                     cutoff: float | None = 4.5,
                     block_size: int = 2048) -> float:
    """Weighted exact-order Gaussian cross overlap; optional types require equality."""
    if sigma <= 0 or not math.isfinite(sigma):
        raise ValueError("sigma must be finite and positive")
    if cutoff is not None and (cutoff <= 0 or not math.isfinite(cutoff)):
        raise ValueError("cutoff must be finite and positive or None")
    if block_size <= 0:
        raise ValueError("block_size must be positive")
    a, b = _points(points_a, "points_a"), _points(points_b, "points_b")
    wa, wb = _weights(weights_a, len(a), "weights_a"), _weights(
        weights_b, len(b), "weights_b")
    ta = tb = None
    if (types_a is None) != (types_b is None):
        raise ValueError("both type arrays must be supplied together")
    if types_a is not None:
        ta, tb = np.asarray(types_a), np.asarray(types_b)
        if ta.shape != (len(a),) or tb.shape != (len(b),):
            raise ValueError("type arrays must contain one value per point")
    total = 0.0
    cutoff_squared = None if cutoff is None else cutoff * cutoff
    denominator = 2.0 * sigma * sigma
    for start in range(0, len(a), block_size):
        stop = min(len(a), start + block_size)
        delta = a[start:stop, None, :] - b[None, :, :]
        squared = np.einsum("ijk,ijk->ij", delta, delta)
        mask = np.ones(squared.shape, dtype=bool)
        if cutoff_squared is not None:
            mask &= squared <= cutoff_squared
        if ta is not None:
            mask &= ta[start:stop, None] == tb[None, :]
        values = np.zeros(squared.shape, dtype=np.float64)
        np.exp(-squared / denominator, out=values, where=mask)
        values[~mask] = 0.0
        total += float(np.sum(values * wa[start:stop, None] * wb[None, :]))
    return total


def gaussian_self_overlap(points: Sequence[Sequence[float]], **kwargs) -> float:
    weights = kwargs.pop("weights", None)
    types = kwargs.pop("types", None)
    return gaussian_overlap(points, points, weights_a=weights, weights_b=weights,
                            types_a=types, types_b=types, **kwargs)


def tanimoto(cross: float, self_a: float, self_b: float) -> float:
    values = np.asarray([cross, self_a, self_b], dtype=np.float64)
    if not np.isfinite(values).all() or np.any(values < 0):
        raise ValueError("overlaps must be finite and non-negative")
    denominator = self_a + self_b - cross
    return 0.0 if denominator <= 0 else float(cross / denominator)


def query_biased_tversky(cross: float, query_self: float, candidate_self: float,
                         *, alpha: float = 0.95, beta: float = 0.05) -> float:
    values = np.asarray([cross, query_self, candidate_self, alpha, beta])
    if not np.isfinite(values).all() or np.any(values[:3] < 0):
        raise ValueError("overlaps and coefficients must be finite; overlaps non-negative")
    if alpha < 0 or beta < 0 or not np.isclose(alpha + beta, 1.0):
        raise ValueError("alpha and beta must be non-negative and sum to one")
    denominator = alpha * query_self + beta * candidate_self
    return 0.0 if denominator <= 0 else float(cross / denominator)


def split_parent_feature_weights(parent_ids: Sequence[int],
                                 parent_weights: Sequence[float]) -> np.ndarray:
    """Divide each functional parent weight across its alternative points."""
    parents = np.asarray(parent_ids)
    weights = np.asarray(parent_weights, dtype=np.float64)
    if parents.ndim != 1 or weights.ndim != 1 or len(parents) != len(weights):
        raise ValueError("parent IDs and weights require equal one-dimensional arrays")
    if not np.isfinite(weights).all() or np.any(weights < 0):
        raise ValueError("parent weights must be finite and non-negative")
    result = np.empty(len(weights), dtype=np.float64)
    for parent in np.unique(parents):
        positions = np.flatnonzero(parents == parent)
        declared = weights[positions]
        if not np.allclose(declared, declared[0]):
            raise ValueError("all alternative points for one parent need one declared weight")
        result[positions] = declared[0] / len(positions)
    return result


def kabsch_transform(source: Sequence[Sequence[float]],
                     target: Sequence[Sequence[float]],
                     weights: Sequence[float] | None = None) -> tuple[float, ...]:
    """Return a proper homogeneous transform mapping source into target."""
    source_points, target_points = _points(source, "source"), _points(target, "target")
    if source_points.shape != target_points.shape or len(source_points) < 2:
        raise ValueError("source and target must have the same shape and at least two points")
    point_weights = _weights(weights, len(source_points), "weights")
    if point_weights.sum() <= 0:
        raise ValueError("Kabsch weights must have positive total")
    point_weights = point_weights / point_weights.sum()
    source_center = np.sum(source_points * point_weights[:, None], axis=0)
    target_center = np.sum(target_points * point_weights[:, None], axis=0)
    covariance = (source_points - source_center).T @ (
        (target_points - target_center) * point_weights[:, None])
    left, _, right_t = np.linalg.svd(covariance)
    rotation = right_t.T @ left.T
    if np.linalg.det(rotation) < 0:
        right_t[-1] *= -1
        rotation = right_t.T @ left.T
    translation = target_center - rotation @ source_center
    matrix = np.eye(4); matrix[:3, :3] = rotation; matrix[:3, 3] = translation
    return tuple(map(float, matrix.ravel()))


def _rotation_axis_angle(axis: np.ndarray, angle: float) -> np.ndarray:
    axis = axis / np.linalg.norm(axis)
    cross = np.asarray([[0, -axis[2], axis[1]],
                        [axis[2], 0, -axis[0]],
                        [-axis[1], axis[0], 0]])
    return np.eye(3) + math.sin(angle) * cross + (1.0 - math.cos(angle)) * (cross @ cross)


def _rotation_between(source: np.ndarray, target: np.ndarray) -> np.ndarray:
    source, target = source / np.linalg.norm(source), target / np.linalg.norm(target)
    cross, dot = np.cross(source, target), float(np.clip(np.dot(source, target), -1, 1))
    norm = np.linalg.norm(cross)
    if norm > 1e-12:
        axis = cross / norm
        return _rotation_axis_angle(axis, math.atan2(norm, dot))
    if dot > 0:
        return np.eye(3)
    helper = np.asarray([1.0, 0.0, 0.0]) if abs(source[0]) < 0.9 else np.asarray([0.0, 1.0, 0.0])
    return _rotation_axis_angle(np.cross(source, helper), math.pi)


def pair_alignment_seeds(candidate_points: Sequence[Sequence[float]],
                         candidate_types: Sequence[int],
                         query_points: Sequence[Sequence[float]],
                         query_types: Sequence[int], *, tolerance: float = 2.0,
                         axial_samples: int = 6,
                         max_seeds: int | None = None, backend: str = 'reference') -> tuple[RigidSeed, ...]:
    """Generate a bounded ordered prefix; reference remains the default kernel."""
    return tuple(iter_pair_alignment_seeds(candidate_points, candidate_types, query_points,
        query_types, tolerance=tolerance, axial_samples=axial_samples,
        max_seeds=max_seeds, backend=backend))


def iter_pair_alignment_seeds(candidate_points, candidate_types, query_points, query_types,
                              *, tolerance=2.0, axial_samples=6, max_seeds=None,
                              backend='reference'):
    """Stream unique seeds without expanding the full Cartesian product."""
    if backend not in ('reference', 'batched'):
        raise ValueError('Unknown seed backend')
    candidate = _points(candidate_points, "candidate_points")
    query = _points(query_points, "query_points")
    candidate_types, query_types = np.asarray(candidate_types), np.asarray(query_types)
    if candidate_types.shape != (len(candidate),) or query_types.shape != (len(query),):
        raise ValueError("feature types require one value per point")
    if tolerance < 0 or axial_samples <= 0:
        raise ValueError("tolerance must be non-negative and axial_samples positive")
    if max_seeds is not None and (isinstance(max_seeds, bool) or not isinstance(max_seeds, (int, np.integer)) or max_seeds < 0):
        raise ValueError("max_seeds must be a non-negative integer or None")
    if max_seeds == 0:
        return ()
    count, seen = 0, set()
    candidate_distances={}
    for qi in range(len(query)):
        for qj in range(qi + 1, len(query)):
            query_vector = query[qj] - query[qi]
            query_distance = np.linalg.norm(query_vector)
            if query_distance <= 1e-12:
                continue
            axials=None
            target_mid = (query[qi] + query[qj]) / 2.0
            for ci in range(len(candidate)):
                for cj in range(ci + 1, len(candidate)):
                    orientations = []
                    if candidate_types[ci] == query_types[qi] and candidate_types[cj] == query_types[qj]:
                        orientations.append((ci, cj))
                    if candidate_types[cj] == query_types[qi] and candidate_types[ci] == query_types[qj]:
                        orientations.append((cj, ci))
                    if not orientations:continue
                    if (ci,cj) not in candidate_distances:
                        candidate_distances[ci,cj]=np.linalg.norm(candidate[cj]-candidate[ci])
                    candidate_distance=candidate_distances[ci,cj]
                    if candidate_distance <= 1e-12 or abs(candidate_distance-query_distance)>tolerance:continue
                    if axials is None:
                        axis=query_vector/query_distance
                        axials=[_rotation_axis_angle(axis,2.0*math.pi*sample/axial_samples) for sample in range(axial_samples)]
                    for first, second in orientations:
                        vector = candidate[second] - candidate[first]
                        base = _rotation_between(vector, query_vector)
                        source_mid = (candidate[first] + candidate[second]) / 2.0
                        if backend == 'batched':
                            rotations = np.asarray(axials) @ base
                            matrices = np.broadcast_to(np.eye(4), (axial_samples,4,4)).copy()
                            matrices[:,:3,:3] = rotations
                            matrices[:,:3,3] = target_mid - (rotations @ source_mid)
                            rounded = np.round(matrices.reshape(-1,16),10)
                            # Validate the entire numeric batch before constructing records.
                            if not np.isfinite(matrices).all() or not np.allclose(
                                    np.linalg.det(rotations),1.,atol=1e-6):
                                raise ValueError('Invalid batched rigid transform')
                        for sample,axial in enumerate(axials):
                            if backend == 'reference':
                                rotation = axial @ base
                                translation = target_mid - rotation @ source_mid
                                matrix = np.eye(4); matrix[:3, :3] = rotation; matrix[:3, 3] = translation
                                key = tuple(np.round(matrix.ravel(), 10))
                            else:
                                matrix = matrices[sample]
                                key = tuple(rounded[sample])
                            if key in seen:
                                continue
                            seen.add(key)
                            values = (f"pair-q{qi}-{qj}-c{first}-{second}-r{sample}",
                                      (first,second),(qi,qj),sample,tuple(map(float,matrix.ravel())))
                            if backend == 'reference':
                                seed = RigidSeed(*values)
                            else:
                                # The same finite/proper-rotation checks were applied in bulk above.
                                seed = object.__new__(RigidSeed)
                                for name,value in zip(('seed_id','candidate_pair','query_pair',
                                                       'axial_sample','transform_matrix'),values):
                                    object.__setattr__(seed,name,value)
                            count += 1
                            yield seed
                            if max_seeds is not None and count >= max_seeds:
                                return



def principal_axis_seeds(candidate_points: Sequence[Sequence[float]],
                         query_points: Sequence[Sequence[float]]) -> tuple[RigidSeed, ...]:
    """Return deterministic proper-rotation PCA seeds plus centroid translation."""
    candidate = _points(candidate_points, "candidate_points")
    query = _points(query_points, "query_points")
    if not len(candidate) or not len(query):
        raise ValueError("principal-axis seeds require non-empty point clouds")
    candidate_center, query_center = candidate.mean(axis=0), query.mean(axis=0)

    def axes(points: np.ndarray, center: np.ndarray) -> np.ndarray:
        if len(points) < 2:
            return np.eye(3)
        covariance = (points - center).T @ (points - center)
        values, vectors = np.linalg.eigh(covariance)
        basis = vectors[:, np.argsort(values)[::-1]]
        if np.linalg.det(basis) < 0:
            basis[:, -1] *= -1
        return basis

    candidate_axes = axes(candidate, candidate_center)
    query_axes = axes(query, query_center)
    rotations = [np.eye(3)]
    for signs in ((1, 1, 1), (1, -1, -1), (-1, 1, -1), (-1, -1, 1)):
        rotations.append(query_axes @ np.diag(signs) @ candidate_axes.T)
    seeds, seen = [], set()
    for index, rotation in enumerate(rotations):
        translation = query_center - rotation @ candidate_center
        matrix = np.eye(4); matrix[:3, :3] = rotation; matrix[:3, 3] = translation
        key = tuple(np.round(matrix.ravel(), 10))
        if key in seen:
            continue
        seen.add(key)
        seeds.append(RigidSeed(
            "centroid" if index == 0 else f"pca-{index}", (-1, -1), (-1, -1),
            -1, tuple(map(float, matrix.ravel()))))
    return tuple(seeds)


def score_overlay(query_shape: Sequence[Sequence[float]],
                  candidate_shape: Sequence[Sequence[float]],
                  query_features: Sequence[Sequence[float]],
                  query_feature_types: Sequence[int],
                  candidate_features: Sequence[Sequence[float]],
                  candidate_feature_types: Sequence[int],
                  transform_matrix: Sequence[float], *, sigma: float = 1.0,
                  cutoff: float | None = 4.5,
                  query_feature_weights: Sequence[float] | None = None,
                  candidate_feature_weights: Sequence[float] | None = None,
                  alpha: float = 0.95, beta: float = 0.05) -> dict:
    """Score one explicit candidate-to-query pose and retain all primitives."""
    query_shape = _points(query_shape, "query_shape")
    candidate_shape = apply_transform(candidate_shape, transform_matrix)
    query_features = _points(query_features, "query_features")
    candidate_features = apply_transform(candidate_features, transform_matrix)
    qfw = _weights(query_feature_weights, len(query_features), "query_feature_weights")
    cfw = _weights(candidate_feature_weights, len(candidate_features), "candidate_feature_weights")

    shape_cross = gaussian_overlap(query_shape, candidate_shape, sigma=sigma, cutoff=cutoff)
    shape_query_self = gaussian_self_overlap(query_shape, sigma=sigma, cutoff=cutoff)
    shape_candidate_self = gaussian_self_overlap(candidate_shape, sigma=sigma, cutoff=cutoff)
    color_cross = gaussian_overlap(
        query_features, candidate_features, sigma=sigma, cutoff=cutoff,
        weights_a=qfw, weights_b=cfw, types_a=query_feature_types,
        types_b=candidate_feature_types)
    color_query_self = gaussian_self_overlap(
        query_features, sigma=sigma, cutoff=cutoff, weights=qfw,
        types=query_feature_types)
    color_candidate_self = gaussian_self_overlap(
        candidate_features, sigma=sigma, cutoff=cutoff, weights=cfw,
        types=candidate_feature_types)
    return {
        "operand_direction": {"A": "query", "B": "candidate"},
        "parameters": {"sigma_angstrom": sigma, "cutoff_angstrom": cutoff,
                       "tversky_alpha_query": alpha, "tversky_beta_candidate": beta},
        "shape": {"cross_overlap_raw": shape_cross,
                  "query_self_overlap_raw": shape_query_self,
                  "candidate_self_overlap_raw": shape_candidate_self,
                  "tanimoto": tanimoto(shape_cross, shape_query_self, shape_candidate_self),
                  "query_biased_tversky": query_biased_tversky(
                      shape_cross, shape_query_self, shape_candidate_self,
                      alpha=alpha, beta=beta)},
        "color": {"cross_overlap_raw": color_cross,
                  "query_self_overlap_raw": color_query_self,
                  "candidate_self_overlap_raw": color_candidate_self,
                  "tanimoto": tanimoto(color_cross, color_query_self, color_candidate_self),
                  "query_biased_tversky": query_biased_tversky(
                      color_cross, color_query_self, color_candidate_self,
                      alpha=alpha, beta=beta)},
        "transform_matrix": list(map(float, transform_matrix)),
    }
