import hashlib
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

import aidd_agent.interaction_matching as matching
from aidd_agent.chemical_geometry import DIRECTION_SIGNED
from aidd_agent.interaction_matching import (
    _maximum_weight_assignment, interaction_match, score_interaction_matches,
)


def _hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def test_assignment_finds_global_optimum_instead_of_greedy_pairing():
    values = np.asarray([[0.90, 0.80], [0.85, 0.10]])
    result = _maximum_weight_assignment(values, np.ones(2))
    assert result.tolist() == [1, 0]


def test_interaction_match_is_typed_directional_and_one_to_one():
    common = {
        "query_points": [[0, 0, 0], [1, 0, 0]],
        "query_types": [1, 1],
        "query_directions": [[1, 0, 0], [1, 0, 0]],
        "query_kinds": [DIRECTION_SIGNED, DIRECTION_SIGNED],
        "query_weights": [2, 1],
        "candidate_types": [1, 1],
        "candidate_kinds": [DIRECTION_SIGNED, DIRECTION_SIGNED],
    }
    exact = interaction_match(
        **common, candidate_points=[[0, 0, 0], [1, 0, 0]],
        candidate_directions=[[1, 0, 0], [1, 0, 0]])
    assert exact["interaction_match_score"] == pytest.approx(1.0)
    assert exact["assignments"].tolist() == [0, 1]

    reversed_direction = interaction_match(
        **common, candidate_points=[[0, 0, 0], [1, 0, 0]],
        candidate_directions=[[-1, 0, 0], [1, 0, 0]])
    assert reversed_direction["interaction_match_score"] == pytest.approx(
        2 * np.exp(-0.5) / 3)
    assert reversed_direction["assignments"].tolist() == [1, -1]


def test_sidecar_reuses_rigid_pose_and_preserves_source(monkeypatch, tmp_path: Path):
    artifact_catalog = tmp_path / "artifact-catalog.json"
    companion_catalog = tmp_path / "companion-catalog.json"
    artifact_catalog.write_text("{}", encoding="utf-8")
    companion_catalog.write_text(json.dumps({
        "artifact_v1_catalog_sha256": _hash(artifact_catalog)}), encoding="utf-8")

    query = tmp_path / "query.npz"
    np.savez(
        query,
        feature_points=np.asarray([[0, 0, 0], [1, 0, 0]], dtype=float),
        feature_types=np.asarray([1, 1], dtype=np.uint8),
        feature_directions=np.asarray([[1, 0, 0], [1, 0, 0]], dtype=float),
        feature_direction_kinds=np.asarray(
            [DIRECTION_SIGNED, DIRECTION_SIGNED], dtype=np.uint8),
        anchored_weights=np.asarray([2, 1], dtype=float),
        anchor_feature_indices=np.asarray([0, 1], dtype=np.int64),
    )
    rigid = tmp_path / "rigid.npz"
    identity = np.eye(4).reshape(1, 16)
    np.savez(
        rigid, global_ids=np.asarray([7], dtype=np.int64),
        molecule_ids=np.asarray(["MOL-7"]),
        conformer_ids=np.asarray(["CNF-7"]),
        objective_names=np.asarray(["shape_only"]),
        shape_only__transform=identity,
        shape_only__objective=np.asarray([0.75]),
    )
    source_hash = _hash(rigid)
    rigid.with_suffix(".manifest.json").write_text(json.dumps({
        "result_sha256": source_hash}), encoding="utf-8")

    artifact = SimpleNamespace(
        molecule_id="MOL-7", conformer_id="CNF-7",
        feature_points=np.asarray([[0, 0, 0], [1, 0, 0]], dtype=float),
        feature_types=np.asarray([1, 1], dtype=np.uint8))
    chemistry = SimpleNamespace(
        molecule_id="MOL-7", conformer_id="CNF-7",
        feature_directions=np.asarray([[1, 0, 0], [1, 0, 0]], dtype=float),
        feature_kinds=np.asarray(
            [DIRECTION_SIGNED, DIRECTION_SIGNED], dtype=np.uint8))

    class FakeArtifactReader:
        def __init__(self, _path):
            pass

        def get(self, global_id):
            assert global_id == 7
            return artifact

    class FakeChemistryReader:
        def __init__(self, _path):
            pass

        def get(self, global_id):
            assert global_id == 7
            return chemistry

    monkeypatch.setattr(matching, "ArtifactCatalogReader", FakeArtifactReader)
    monkeypatch.setattr(matching, "ChemicalCompanionReader", FakeChemistryReader)

    output = tmp_path / "interaction-matches.npz"
    manifest = score_interaction_matches(
        artifact_catalog, companion_catalog, query, rigid, output)
    assert manifest["candidates"] == 1
    assert manifest["anchors"] == 2
    assert manifest["analysis"]["shape_only"]["interaction_min_median_max"] == [
        1.0, 1.0, 1.0]
    assert _hash(rigid) == source_hash
    with np.load(output, allow_pickle=False) as result:
        assert result["global_ids"].tolist() == [7]
        assert result["shape_only__interaction_match_score"].tolist() == [1.0]
        assert result["shape_only__anchor_assignments"].tolist() == [[0, 1]]
        assert result["shape_only__anchor_scores"].tolist() == [[1.0, 1.0]]
    document = json.loads(output.with_suffix(".manifest.json").read_text())
    assert document["inputs"]["rigid_result"]["sha256"] == source_hash
    assert document["output_sha256"] == _hash(output)
