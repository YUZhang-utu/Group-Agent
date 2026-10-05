"""Scientific join and missingness invariants for the reserved scoring interface."""
from dataclasses import replace

import pytest

from aidd_agent.rescoring_contract import ModelSpec, PoseInput, PoseScore, validate_batch


SPEC = ModelSpec("fixture", "test", "a" * 64, "b" * 64,
                 "interaction_score", "model_native", "lower")
POSE = PoseInput("request-1", "mol-1", "cid-1", "receptor-1", "rec.pdb",
                 "c" * 64, "poses.sdf", 0, "d" * 64, "state-1", "e" * 64)
SCORE = PoseScore("request-1", "c" * 64, "d" * 64, 0, "ok", -7.0)


def test_preserves_valid_scores_and_explicit_failures():
    other = replace(POSE, request_id="request-2", conformer_id="cid-2")
    failure = replace(SCORE, request_id="request-2", status="unsupported_input",
                      value=None, reason="Unsupported atom type")
    result = validate_batch(SPEC, [POSE, other], [failure, SCORE])
    assert result["request-1"].value == -7.0
    assert result["request-2"].value is None


@pytest.mark.parametrize("change", [
    {"receptor_sha256": "f" * 64}, {"pose_file_sha256": "f" * 64},
    {"pose_index": 1}, {"pose_index": False}, {"request_id": "unknown"},
    {"value": float("nan")}, {"value": float("inf")}, {"value": True},
    {"status": "model_failed", "value": 0, "reason": "Failure"},
    {"status": "not_run", "value": None},
])
def test_rejects_misattributed_or_misleading_scores(change):
    with pytest.raises(ValueError):
        validate_batch(SPEC, [POSE], [replace(SCORE, **change)])


def test_no_dropped_or_duplicated_candidates():
    for poses, results in [([POSE], []), ([POSE], [SCORE, SCORE]),
                           ([POSE, POSE], [SCORE])]:
        with pytest.raises(ValueError):
            validate_batch(SPEC, poses, results)
