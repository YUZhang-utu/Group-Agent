import hashlib
import json
from pathlib import Path

import numpy as np
import pytest

from aidd_agent.conformer_artifacts import FEATURE_DTYPE, META_DTYPE
from aidd_agent.gaussian_batch import (
    ArtifactCatalogReader, load_candidate_ids, score_gaussian_candidates,
    write_gaussian_query,
)
from aidd_agent.gaussian_overlay import score_overlay


def _hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _artifact(tmp_path: Path) -> Path:
    shard = tmp_path / "shard-a"
    shard.mkdir()
    origin = np.asarray([10.0, -2.0, 4.0], dtype=np.float32)
    points = np.asarray([[0, 0, 0], [100, 0, 0], [0, 200, 0], [0, 0, 300]], dtype="<i2")
    points.tofile(shard / "coords.bin")
    features = np.empty(2, dtype=FEATURE_DTYPE)
    features["xyz"] = np.asarray([[0, 0, 0], [0, 0, 300]], dtype="<i2")
    features["type"] = [1, 2]
    features["parent"] = [0, 1]
    features.tofile(shard / "feats.bin")
    meta = np.asarray([(0, 0, 0, 4, 2, origin, [1, 2, 3], [1, 2, 3],
                        [1, 1, 0, 0, 0, 0], [0, 0])], dtype=META_DTYPE)
    meta.tofile(shard / "meta.bin")
    np.asarray([b"CONF-1"], dtype="S16").tofile(shard / "conformer_ids.bin")
    np.asarray([b"MOL-1"], dtype="S16").tofile(shard / "molecule_ids.bin")
    manifest = {"format": "aidd-conformer-artifact-shard", "version": 1,
                "library_id": "LIB-X", "global_id_start": 0, "conformers": 1,
                "coordinate_scale": 100.0}
    (shard / "manifest.json").write_text(json.dumps(manifest), encoding="utf-8")
    catalog = {"format": "aidd-conformer-artifact-catalog", "version": 1,
               "library_id": "LIB-X", "conformers": 1,
               "shards": [{"name": "shard-a", "path": "Z:/stale/windows/path",
                           "global_id_start": 0, "conformers": 1}]}
    path = tmp_path / "catalog.json"
    path.write_text(json.dumps(catalog), encoding="utf-8")
    return path


def test_artifact_reader_dequantizes_and_uses_relocated_shard(tmp_path: Path):
    record = ArtifactCatalogReader(_artifact(tmp_path)).get(0)
    assert record.molecule_id == "MOL-1"
    assert record.conformer_id == "CONF-1"
    assert np.allclose(record.shape_points,
                       [[10, -2, 4], [11, -2, 4], [10, 0, 4], [10, -2, 7]])
    assert np.allclose(record.feature_points, [[10, -2, 4], [10, -2, 7]])
    with pytest.raises(IndexError, match="outside"):
        ArtifactCatalogReader(tmp_path / "catalog.json").get(1)


def test_candidate_ids_preserve_order_and_reject_duplicates(tmp_path: Path):
    path = tmp_path / "ids.npy"
    np.save(path, np.asarray([7, 2, 9], dtype=np.int64))
    assert load_candidate_ids(path).tolist() == [7, 2, 9]
    np.save(path, np.asarray([7, 7], dtype=np.int64))
    with pytest.raises(ValueError, match="duplicates"):
        load_candidate_ids(path)
    np.save(path, np.asarray([1.5]))
    with pytest.raises(ValueError, match="integer"):
        load_candidate_ids(path)
    text = tmp_path / "ids.txt"; text.write_text("9\n2\n", encoding="utf-8")
    assert load_candidate_ids(text).tolist() == [9, 2]


def test_batch_retains_ids_and_recovers_locked_rigid_pose(tmp_path: Path):
    catalog = _artifact(tmp_path)
    candidate = ArtifactCatalogReader(catalog).get(0)
    angle = np.pi / 3
    rotation = np.asarray([[np.cos(angle), -np.sin(angle), 0],
                           [np.sin(angle), np.cos(angle), 0], [0, 0, 1]])
    translation = np.asarray([4.0, 6.0, -3.0])
    transform = lambda x: x @ rotation.T + translation
    query = tmp_path / "query.npz"
    write_gaussian_query(
        query, shape_points=transform(candidate.shape_points),
        feature_points=transform(candidate.feature_points), feature_types=[1, 2],
        anchored_weights=[1.5, 1.5], anchor_feature_indices=[0, 1],
        source={"query_id": "synthetic"})
    ids = tmp_path / "ids.npy"; np.save(ids, np.asarray([0], dtype=np.int64))
    output = tmp_path / "scores.npz"
    manifest = score_gaussian_candidates(catalog, query, ids, output)
    assert manifest["candidates_input"] == manifest["candidates_output"] == 1
    with np.load(output, allow_pickle=False) as result:
        assert result["global_ids"].tolist() == [0]
        assert result["best_seed_ids"].shape == (1, 3)
        for objective in result["objective_names"]:
            prefix = str(objective)
            assert result[f"{prefix}__transform"].shape == (1, 16)
            assert result[f"{prefix}__shape_tanimoto"][0] == pytest.approx(1.0, abs=1e-10)
        assert result["shape_only__color_tanimoto"][0] == pytest.approx(1.0, abs=1e-10)
        assert result["atomcentered_unweighted_joint__color_tanimoto"][0] == pytest.approx(
            1.0, abs=1e-10)
        assert result["atomcentered_anchored_joint__color_tanimoto"][0] == pytest.approx(
            6.0 / 7.0, abs=1e-10)
        for objective in result["objective_names"]:
            prefix = str(objective)
            rescored = score_overlay(
                transform(candidate.shape_points), candidate.shape_points,
                transform(candidate.feature_points), [1, 2], candidate.feature_points,
                [1, 2], result[f"{prefix}__transform"][0],
                query_feature_weights=([1.5, 1.5] if "anchored" in prefix else None))
            assert rescored["shape"]["tanimoto"] == pytest.approx(
                result[f"{prefix}__shape_tanimoto"][0])
            assert rescored["color"]["tanimoto"] == pytest.approx(
                result[f"{prefix}__color_tanimoto"][0])
    first = {k: v.copy() for k, v in np.load(output, allow_pickle=False).items()}
    score_gaussian_candidates(catalog, query, ids, output)
    with np.load(output, allow_pickle=False) as second:
        assert all(np.array_equal(first[key], second[key]) for key in first)
