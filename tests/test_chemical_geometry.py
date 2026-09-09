import numpy as np
import pytest

from aidd_agent.chemical_geometry import (
    DIRECTION_AXIAL, DIRECTION_SIGNED, TerminalTorsion,
    directional_gaussian_overlap, refine_terminal_torsions,
    rigid_micro_seeds, rotate_terminal_torsion, transform_directions,
    vdw_exclusion_metrics,
)
from aidd_agent.gaussian_overlay import apply_transform


def test_directions_rotate_without_translation_and_score_signed_vs_axial():
    transform = np.asarray([[0, -1, 0, 100], [1, 0, 0, -50],
                            [0, 0, 1, 25], [0, 0, 0, 1]], dtype=float)
    moved = transform_directions([[1, 0, 0]], transform)
    assert np.allclose(moved, [[0, 1, 0]])
    common = dict(points_a=[[0, 0, 0]], points_b=[[0, 0, 0]],
                  types_a=[1], types_b=[1], directions_a=[[1, 0, 0]])
    aligned = directional_gaussian_overlap(
        **common, directions_b=[[1, 0, 0]], kinds_a=[DIRECTION_SIGNED],
        kinds_b=[DIRECTION_SIGNED])
    reversed_signed = directional_gaussian_overlap(
        **common, directions_b=[[-1, 0, 0]], kinds_a=[DIRECTION_SIGNED],
        kinds_b=[DIRECTION_SIGNED])
    reversed_axial = directional_gaussian_overlap(
        **common, directions_b=[[-1, 0, 0]], kinds_a=[DIRECTION_AXIAL],
        kinds_b=[DIRECTION_AXIAL])
    assert aligned == pytest.approx(1.0)
    assert reversed_signed == pytest.approx(0.0)
    assert reversed_axial == pytest.approx(1.0)


def test_vdw_exclusion_orders_separated_close_and_overlapping_controls():
    far = vdw_exclusion_metrics([[0, 0, 0]], [6], [[10, 0, 0]], [6])
    close = vdw_exclusion_metrics([[0, 0, 0]], [6], [[3, 0, 0]], [6])
    overlap = vdw_exclusion_metrics([[0, 0, 0]], [6], [[1, 0, 0]], [6])
    assert far["soft_exclusion_penalty"] == 0.0
    assert close["maximum_vdw_penetration_angstrom"] == pytest.approx(0.4)
    assert overlap["soft_exclusion_penalty"] > close["soft_exclusion_penalty"]
    assert overlap["severe_candidate_atoms"] == 1


def test_terminal_rotation_preserves_axis_bonds_and_nonmoving_atoms():
    points = np.asarray([[-1, 0, 0], [0, 0, 0], [1, 0, 0],
                         [2, 1, 0], [2, 2, 0]], dtype=float)
    torsion = TerminalTorsion(0, 1, 2, 3, (2, 3, 4))
    moved = rotate_terminal_torsion(points, torsion, 90)
    assert moved[:2] == pytest.approx(points[:2])
    assert np.linalg.norm(moved[2] - moved[1]) == pytest.approx(1.0)
    assert np.linalg.norm(moved[3] - moved[2]) == pytest.approx(
        np.linalg.norm(points[3] - points[2]))
    assert moved[3] == pytest.approx([2, 0, 1])


def test_bounded_torsion_refinement_never_loses_zero_angle_baseline():
    points = np.asarray([[-1, 0, 0], [0, 0, 0], [1, 0, 0], [2, 1, 0]], dtype=float)
    torsion = TerminalTorsion(0, 1, 2, 3, (2, 3))
    target = np.asarray([2, 0, 1], dtype=float)
    result = refine_terminal_torsions(
        points, [torsion], lambda xyz: -float(np.linalg.norm(xyz[3] - target)),
        angle_offsets=(-90, 0, 90), max_torsions=1, beam_width=3)
    assert result["refined_score"] >= result["baseline_score"]
    assert result["angles_degrees"] == [90.0]
    assert result["coordinates"][3] == pytest.approx(target)


def test_rigid_micro_seeds_are_fixed_deterministic_and_preserve_geometry():
    points = np.asarray([[10, 0, 0], [11, 1, 0], [12, 0, 1]], dtype=float)
    first = rigid_micro_seeds(points)
    second = rigid_micro_seeds(points)

    assert len(first) == 13
    assert [seed.seed_id for seed in first] == [seed.seed_id for seed in second]
    assert np.array_equal(apply_transform(points, first[0].transform_matrix), points)
    assert apply_transform(points, first[1].transform_matrix) == pytest.approx(
        points + [-0.25, 0, 0])

    baseline_distances = np.linalg.norm(points[:, None] - points[None, :], axis=2)
    for left, right in zip(first, second):
        assert np.array_equal(left.transform_matrix, right.transform_matrix)
        moved = apply_transform(points, left.transform_matrix)
        distances = np.linalg.norm(moved[:, None] - moved[None, :], axis=2)
        assert distances == pytest.approx(baseline_distances)
    rotated = apply_transform(points, first[7].transform_matrix)
    assert rotated.mean(axis=0) == pytest.approx(points.mean(axis=0))
