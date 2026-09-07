import numpy as np
import pytest

from aidd_agent.gaussian_overlay import (
    apply_transform, gaussian_overlap, gaussian_self_overlap, kabsch_transform,
    pair_alignment_seeds, query_biased_tversky, score_overlay,
    split_parent_feature_weights, tanimoto,
)


def test_identical_cloud_scores_one_and_self_is_rigid_invariant():
    query = np.asarray([[0., 0., 0.], [1., 0., 0.], [0., 2., 0.]])
    self_query = gaussian_self_overlap(query)
    moved = query @ np.asarray([[0., -1., 0.], [1., 0., 0.], [0., 0., 1.]]) + [5, 6, 7]
    assert gaussian_self_overlap(moved) == pytest.approx(self_query)
    cross = gaussian_overlap(query, query)
    assert tanimoto(cross, self_query, self_query) == pytest.approx(1.)
    assert query_biased_tversky(cross, self_query, self_query) == pytest.approx(1.)


def test_color_requires_compatible_type_and_parent_weight_is_conserved():
    point = [[0., 0., 0.]]
    assert gaussian_overlap(point, point, types_a=[1], types_b=[2]) == 0
    weights = split_parent_feature_weights([4, 4, 7], [1.5, 1.5, 0.5])
    assert weights.tolist() == [.75, .75, .5]
    assert weights[:2].sum() == pytest.approx(1.5)


def test_kabsch_recovers_known_proper_rigid_transform():
    source = np.asarray([[0., 0., 0.], [1., 0., 0.], [0., 2., 0.]])
    rotation = np.asarray([[0., -1., 0.], [1., 0., 0.], [0., 0., 1.]])
    target = source @ rotation.T + [4., 5., 6.]
    transform = kabsch_transform(source, target)
    assert np.allclose(apply_transform(source, transform), target)
    assert np.linalg.det(np.asarray(transform).reshape(4, 4)[:3, :3]) == pytest.approx(1.)


def test_pair_seed_recovers_typed_pair_and_all_rotations_are_proper():
    candidate = np.asarray([[0., 0., 0.], [4., 0., 0.], [2., 1., 0.]])
    query = np.asarray([[10., 3., 1.], [10., 7., 1.]])
    seeds = pair_alignment_seeds(candidate, [2, 1, 5], query, [2, 1], tolerance=.01)
    assert len(seeds) == 6
    mapped = apply_transform(candidate[[0, 1]], seeds[0].transform_matrix)
    assert np.allclose(mapped, query)
    for seed in seeds:
        assert np.linalg.det(np.asarray(seed.transform_matrix).reshape(4, 4)[:3, :3]) == pytest.approx(1.)


def test_overlay_schema_keeps_primitives_direction_and_pose():
    points = np.asarray([[0., 0., 0.], [1., 0., 0.]])
    identity = tuple(np.eye(4).ravel())
    result = score_overlay(points, points, points, [1, 2], points, [1, 2], identity)
    assert result["operand_direction"] == {"A": "query", "B": "candidate"}
    assert result["shape"]["tanimoto"] == pytest.approx(1.)
    assert result["color"]["query_biased_tversky"] == pytest.approx(1.)
    assert len(result["transform_matrix"]) == 16


@pytest.mark.parametrize("call", [
    lambda: gaussian_overlap([[0, 0, 0]], [[0, 0, 0]], sigma=0),
    lambda: gaussian_overlap([[0, 0, 0]], [[0, 0, 0]], weights_a=[-1]),
    lambda: query_biased_tversky(1, 1, 1, alpha=.8, beta=.3),
    lambda: apply_transform([[0, 0]], np.eye(4).ravel()),
])
def test_invalid_inputs_fail_explicitly(call):
    with pytest.raises(ValueError):
        call()
