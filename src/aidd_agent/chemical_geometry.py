from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Callable, Sequence

import numpy as np


DIRECTION_NONE = 0
DIRECTION_SIGNED = 1
DIRECTION_AXIAL = 2

VDW_RADII = {
    6: 1.70, 7: 1.55, 8: 1.52, 9: 1.47, 15: 1.80, 16: 1.80,
    17: 1.75, 35: 1.85, 53: 1.98,
}


def _points(values: Sequence[Sequence[float]], name: str) -> np.ndarray:
    result = np.asarray(values, dtype=np.float64)
    if result.ndim != 2 or result.shape[1:] != (3,) or not np.isfinite(result).all():
        raise ValueError(f"{name} must have finite shape (N,3)")
    return result


def unit_vector(vector: Sequence[float]) -> np.ndarray | None:
    value = np.asarray(vector, dtype=np.float64)
    if value.shape != (3,) or not np.isfinite(value).all():
        raise ValueError("direction must be a finite 3-vector")
    length = float(np.linalg.norm(value))
    return None if length <= 1e-12 else value / length


def transform_directions(directions: Sequence[Sequence[float]],
                         transform_matrix: Sequence[float]) -> np.ndarray:
    """Rotate unit directions by a homogeneous transform; never translate."""
    values = _points(directions, "directions")
    matrix = np.asarray(transform_matrix, dtype=np.float64)
    if matrix.size != 16:
        raise ValueError("transform must contain 16 values")
    matrix = matrix.reshape(4, 4)
    if (not np.isfinite(matrix).all()
            or not np.allclose(matrix[3], [0, 0, 0, 1])):
        raise ValueError("transform must be a finite homogeneous 4x4 matrix")
    moved = values @ matrix[:3, :3].T
    lengths = np.linalg.norm(moved, axis=1)
    valid = lengths > 1e-12
    moved[valid] /= lengths[valid, None]
    moved[~valid] = 0.0
    return moved


def directional_gaussian_overlap(
        points_a: Sequence[Sequence[float]], points_b: Sequence[Sequence[float]],
        types_a: Sequence[int], types_b: Sequence[int],
        directions_a: Sequence[Sequence[float]],
        directions_b: Sequence[Sequence[float]],
        kinds_a: Sequence[int], kinds_b: Sequence[int], *,
        sigma: float = 1.0, cutoff: float | None = 4.5,
        angular_power: float = 2.0,
        weights_a: Sequence[float] | None = None,
        weights_b: Sequence[float] | None = None) -> float:
    """Typed Gaussian overlap weighted by signed/axial direction agreement."""
    if sigma <= 0 or not math.isfinite(sigma) or angular_power <= 0:
        raise ValueError("sigma and angular_power must be finite and positive")
    if cutoff is not None and (cutoff <= 0 or not math.isfinite(cutoff)):
        raise ValueError("cutoff must be finite and positive or None")
    a, b = _points(points_a, "points_a"), _points(points_b, "points_b")
    da, db = _points(directions_a, "directions_a"), _points(
        directions_b, "directions_b")
    ta, tb = np.asarray(types_a), np.asarray(types_b)
    ka, kb = np.asarray(kinds_a, dtype=np.uint8), np.asarray(kinds_b, dtype=np.uint8)
    if not (da.shape == a.shape and db.shape == b.shape
            and ta.shape == ka.shape == (len(a),)
            and tb.shape == kb.shape == (len(b),)):
        raise ValueError("directional feature arrays have inconsistent shapes")
    if (not set(map(int, ka)).issubset({DIRECTION_NONE, DIRECTION_SIGNED, DIRECTION_AXIAL})
            or not set(map(int, kb)).issubset(
                {DIRECTION_NONE, DIRECTION_SIGNED, DIRECTION_AXIAL})):
        raise ValueError("unknown direction kind")
    for directions, kinds in ((da, ka), (db, kb)):
        valid = kinds != DIRECTION_NONE
        if np.any(valid) and not np.allclose(
                np.linalg.norm(directions[valid], axis=1), 1.0, atol=5e-3):
            raise ValueError("valid feature directions must be unit vectors")
    wa = np.ones(len(a)) if weights_a is None else np.asarray(weights_a, dtype=float)
    wb = np.ones(len(b)) if weights_b is None else np.asarray(weights_b, dtype=float)
    if (wa.shape != (len(a),) or wb.shape != (len(b),)
            or np.any(wa < 0) or np.any(wb < 0)
            or not np.isfinite(wa).all() or not np.isfinite(wb).all()):
        raise ValueError("weights must be finite and non-negative")
    delta = a[:, None, :] - b[None, :, :]
    squared = np.einsum("ijk,ijk->ij", delta, delta)
    mask = (ta[:, None] == tb[None, :]) & (ka[:, None] != DIRECTION_NONE)
    mask &= ka[:, None] == kb[None, :]
    if cutoff is not None:
        mask &= squared <= cutoff * cutoff
    cosine = np.clip(da @ db.T, -1.0, 1.0)
    axial = ka[:, None] == DIRECTION_AXIAL
    agreement = np.where(axial, np.abs(cosine), np.maximum(cosine, 0.0))
    agreement = agreement ** angular_power
    spatial = np.exp(-squared / (2.0 * sigma * sigma))
    return float(np.sum(np.where(mask, spatial * agreement, 0.0)
                        * wa[:, None] * wb[None, :]))


def vdw_exclusion_metrics(
        candidate_points: Sequence[Sequence[float]],
        candidate_atomic_numbers: Sequence[int],
        protein_points: Sequence[Sequence[float]],
        protein_atomic_numbers: Sequence[int], *,
        tolerance: float = 0.4, severe_penetration: float = 0.75) -> dict:
    """Element-aware soft penetration metrics; returns annotations, not a filter."""
    if tolerance < 0 or severe_penetration <= tolerance:
        raise ValueError("require 0 <= tolerance < severe_penetration")
    candidate = _points(candidate_points, "candidate_points")
    protein = _points(protein_points, "protein_points")
    ca = np.asarray(candidate_atomic_numbers, dtype=np.int16)
    pa = np.asarray(protein_atomic_numbers, dtype=np.int16)
    if ca.shape != (len(candidate),) or pa.shape != (len(protein),):
        raise ValueError("atomic number arrays have inconsistent lengths")
    cr = np.asarray([VDW_RADII.get(int(value), 1.70) for value in ca])
    pr = np.asarray([VDW_RADII.get(int(value), 1.70) for value in pa])
    delta = candidate[:, None, :] - protein[None, :, :]
    distances = np.sqrt(np.einsum("ijk,ijk->ij", delta, delta))
    penetration = cr[:, None] + pr[None, :] - distances
    maximum_by_atom = np.maximum(0.0, penetration.max(axis=1))
    excess = np.maximum(0.0, maximum_by_atom - tolerance)
    return {
        "minimum_center_distance_angstrom": float(distances.min()),
        "maximum_vdw_penetration_angstrom": float(maximum_by_atom.max()),
        "clashing_candidate_atoms": int(np.count_nonzero(maximum_by_atom > tolerance)),
        "severe_candidate_atoms": int(np.count_nonzero(
            maximum_by_atom > severe_penetration)),
        "clashing_candidate_fraction": float(np.mean(maximum_by_atom > tolerance)),
        "soft_exclusion_penalty": float(np.dot(excess, excess)),
    }


@dataclass(frozen=True)
class TerminalTorsion:
    atom_a: int
    atom_b: int
    atom_c: int
    atom_d: int
    moving_atoms: tuple[int, ...]


@dataclass(frozen=True)
class RigidMicroSeed:
    seed_id: str
    transform_matrix: np.ndarray


def rigid_micro_seeds(
        points: Sequence[Sequence[float]], *,
        translation_angstrom: float = 0.25,
        rotation_degrees: float = 5.0) -> tuple[RigidMicroSeed, ...]:
    """Return identity plus locked +/- axis translations and centroid rotations."""
    coordinates = _points(points, "points")
    if not len(coordinates):
        raise ValueError("points must not be empty")
    if (translation_angstrom <= 0 or not math.isfinite(translation_angstrom)
            or rotation_degrees <= 0 or not math.isfinite(rotation_degrees)):
        raise ValueError("micro-seed translation and rotation must be positive")
    seeds = [RigidMicroSeed("identity", np.eye(4, dtype=np.float64))]
    axes = (("x", np.asarray([1.0, 0.0, 0.0])),
            ("y", np.asarray([0.0, 1.0, 0.0])),
            ("z", np.asarray([0.0, 0.0, 1.0])))
    for name, axis in axes:
        for label, sign in (("minus", -1.0), ("plus", 1.0)):
            matrix = np.eye(4, dtype=np.float64)
            matrix[:3, 3] = sign * translation_angstrom * axis
            seeds.append(RigidMicroSeed(f"translate_{name}_{label}", matrix))
    center = coordinates.mean(axis=0)
    for name, axis in axes:
        cross = np.asarray([
            [0.0, -axis[2], axis[1]],
            [axis[2], 0.0, -axis[0]],
            [-axis[1], axis[0], 0.0],
        ])
        for label, sign in (("minus", -1.0), ("plus", 1.0)):
            angle = math.radians(sign * rotation_degrees)
            rotation = (np.eye(3) * math.cos(angle)
                        + (1.0 - math.cos(angle)) * np.outer(axis, axis)
                        + math.sin(angle) * cross)
            matrix = np.eye(4, dtype=np.float64)
            matrix[:3, :3] = rotation
            matrix[:3, 3] = center - rotation @ center
            seeds.append(RigidMicroSeed(f"rotate_{name}_{label}", matrix))
    return tuple(seeds)


def rotate_terminal_torsion(points: Sequence[Sequence[float]],
                            torsion: TerminalTorsion,
                            angle_degrees: float) -> np.ndarray:
    coordinates = _points(points, "points").copy()
    count = len(coordinates)
    indices = np.asarray(torsion.moving_atoms, dtype=np.int64)
    referenced = (torsion.atom_a, torsion.atom_b, torsion.atom_c, torsion.atom_d)
    if (any(index < 0 or index >= count for index in referenced)
            or indices.ndim != 1 or not len(indices)
            or np.any(indices < 0) or np.any(indices >= count)
            or torsion.atom_b in indices or torsion.atom_c not in indices):
        raise ValueError("invalid terminal torsion atom indices")
    axis = coordinates[torsion.atom_c] - coordinates[torsion.atom_b]
    unit = unit_vector(axis)
    if unit is None:
        raise ValueError("torsion axis has zero length")
    angle = math.radians(float(angle_degrees))
    cross = np.asarray([[0, -unit[2], unit[1]], [unit[2], 0, -unit[0]],
                        [-unit[1], unit[0], 0]], dtype=np.float64)
    rotation = (np.eye(3) * math.cos(angle)
                + (1.0 - math.cos(angle)) * np.outer(unit, unit)
                + math.sin(angle) * cross)
    origin = coordinates[torsion.atom_b]
    coordinates[indices] = (coordinates[indices] - origin) @ rotation.T + origin
    return coordinates


def refine_terminal_torsions(
        points: Sequence[Sequence[float]], torsions: Sequence[TerminalTorsion],
        score_function: Callable[[np.ndarray], float], *,
        angle_offsets: Sequence[float] = (-60, -30, 0, 30, 60),
        max_torsions: int = 2, beam_width: int = 8) -> dict:
    """Deterministic bounded beam refinement; zero-angle keeps baseline eligible."""
    if max_torsions < 0 or beam_width <= 0 or not angle_offsets:
        raise ValueError("invalid torsion refinement limits")
    if 0.0 not in {float(value) for value in angle_offsets}:
        raise ValueError("angle offsets must include zero to preserve baseline")
    initial = _points(points, "points")
    baseline = float(score_function(initial.copy()))
    if not math.isfinite(baseline):
        raise ValueError("score function returned a non-finite baseline")
    chosen = sorted(torsions, key=lambda value: (
        len(value.moving_atoms), value.atom_b, value.atom_c,
        value.moving_atoms))[:max_torsions]
    beam = [(baseline, tuple(), initial.copy())]
    for torsion in chosen:
        expanded = []
        for _, angles, coordinates in beam:
            for angle in sorted({float(value) for value in angle_offsets}):
                moved = rotate_terminal_torsion(coordinates, torsion, angle)
                score = float(score_function(moved))
                if not math.isfinite(score):
                    raise ValueError("score function returned a non-finite value")
                expanded.append((score, angles + (angle,), moved))
        expanded.sort(key=lambda item: (-item[0], item[1]))
        beam = expanded[:beam_width]
    best = beam[0]
    return {
        "baseline_score": baseline,
        "refined_score": float(best[0]),
        "score_delta": float(best[0] - baseline),
        "angles_degrees": list(best[1]),
        "coordinates": best[2],
        "torsions_considered": len(chosen),
        "evaluations": 1 + sum(
            min(beam_width, len(angle_offsets) ** step) * len(angle_offsets)
            for step in range(len(chosen))),
    }
