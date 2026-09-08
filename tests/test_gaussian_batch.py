import hashlib
import json
from pathlib import Path

import numpy as np
import pytest

from aidd_agent.conformer_artifacts import FEATURE_DTYPE, META_DTYPE
from aidd_agent.gaussian_batch import (
    ArtifactCatalogReader, load_candidate_ids, run_scaled_gaussian_reranking,
    run_staged_gaussian_reranking,
    score_gaussian_candidates, write_gaussian_query,
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


def _multi_artifact(tmp_path: Path, count: int = 6) -> Path:
    shard = tmp_path / "shard-multi"
    shard.mkdir()
    coords, feature_records, meta_rows = [], [], []
    for index in range(count):
        origin = np.asarray([10.0 + index, -2.0, 4.0], dtype=np.float32)
        points = np.asarray([
            [0, 0, 0], [100 + 5 * index, 0, 0],
            [0, 200 - 3 * index, 0], [0, 0, 300 + 2 * index]], dtype="<i2")
        features = np.empty(2, dtype=FEATURE_DTYPE)
        features["xyz"] = points[[0, 3]]
        features["type"] = [1, 2]
        features["parent"] = [0, 1]
        coords.append(points); feature_records.append(features)
        meta_rows.append((index, index * 4, index * 2, 4, 2, origin,
                          [1, 2, 3], [1, 2, 3], [1, 1, 0, 0, 0, 0], [0, 0]))
    np.concatenate(coords).tofile(shard / "coords.bin")
    np.concatenate(feature_records).tofile(shard / "feats.bin")
    np.asarray(meta_rows, dtype=META_DTYPE).tofile(shard / "meta.bin")
    np.asarray([f"CONF-{i}".encode() for i in range(count)], dtype="S16").tofile(
        shard / "conformer_ids.bin")
    np.asarray([f"MOL-{i}".encode() for i in range(count)], dtype="S16").tofile(
        shard / "molecule_ids.bin")
    manifest = {"format": "aidd-conformer-artifact-shard", "version": 1,
                "library_id": "LIB-M", "global_id_start": 0,
                "conformers": count, "coordinate_scale": 100.0}
    (shard / "manifest.json").write_text(json.dumps(manifest), encoding="utf-8")
    catalog = {"format": "aidd-conformer-artifact-catalog", "version": 1,
               "library_id": "LIB-M", "conformers": count,
               "shards": [{"name": shard.name, "path": str(shard),
                           "global_id_start": 0, "conformers": count}]}
    path = tmp_path / "catalog-multi.json"
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


def test_gaussian_query_retains_optional_directional_features(tmp_path: Path):
    output = tmp_path / "directional-query.npz"
    manifest = write_gaussian_query(
        output, shape_points=[[0, 0, 0]],
        feature_points=[[0, 0, 0], [1, 0, 0]], feature_types=[1, 6],
        anchored_weights=[1.5, 1.0], anchor_feature_indices=[0],
        feature_directions=[[1, 0, 0], [0, 0, 1]],
        feature_direction_kinds=[1, 2])
    assert manifest["projected_color_available"] is True
    assert manifest["directional_features"] == 2
    with np.load(output, allow_pickle=False) as archive:
        assert archive["feature_direction_kinds"].tolist() == [1, 2]
        assert np.allclose(archive["feature_directions"],
                           [[1, 0, 0], [0, 0, 1]])


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


def _staged_inputs(tmp_path: Path):
    catalog = _multi_artifact(tmp_path)
    reference = ArtifactCatalogReader(catalog).get(3)
    query = tmp_path / "staged-query.npz"
    write_gaussian_query(
        query, shape_points=reference.shape_points,
        feature_points=reference.feature_points, feature_types=[1, 2],
        anchored_weights=[1.5, 1.5], anchor_feature_indices=[0, 1],
        source={"query_id": "staged-synthetic"})
    ids = tmp_path / "staged-ids.npy"
    np.save(ids, np.arange(6, dtype=np.int64))
    return catalog, query, ids


def _load_arrays(path: Path) -> dict[str, np.ndarray]:
    with np.load(path, allow_pickle=False) as archive:
        return {name: archive[name].copy() for name in archive.files}


def test_staged_workers_match_and_selection_is_exact_union(tmp_path: Path):
    catalog, query, ids = _staged_inputs(tmp_path)
    common = dict(stage="all", chunk_size=2, top_n_per_objective=2,
                  max_pair_seeds=6, progress_every=0)
    one = tmp_path / "run-one"
    two = tmp_path / "run-two"
    run_staged_gaussian_reranking(catalog, query, ids, one, workers=1, **common)
    run_staged_gaussian_reranking(catalog, query, ids, two, workers=2, **common)
    for relative in ("coarse/merged-scores.npz", "refine/selected-candidates.npz",
                     "refine/merged-scores.npz"):
        left, right = _load_arrays(one / relative), _load_arrays(two / relative)
        assert left.keys() == right.keys()
        assert all(np.array_equal(left[key], right[key]) for key in left)
    coarse = _load_arrays(one / "coarse/merged-scores.npz")
    selection = _load_arrays(one / "refine/selected-candidates.npz")
    refined = _load_arrays(one / "refine/merged-scores.npz")
    assert coarse["global_ids"].tolist() == list(range(6))
    expected = set()
    for name in coarse["objective_names"].astype(str):
        expected.update(np.argsort(-coarse[f"{name}__objective"], kind="stable")[:2])
    assert set(selection["coarse_indices"]) == expected
    assert len(selection["global_ids"]) <= 6
    assert np.array_equal(refined["global_ids"], selection["global_ids"])
    assert not list(one.rglob("*.partial"))


def test_staged_resume_reuses_valid_chunk_and_repairs_corrupt_chunk(tmp_path: Path):
    catalog, query, ids = _staged_inputs(tmp_path)
    output = tmp_path / "resume-run"
    kwargs = dict(stage="coarse", workers=1, chunk_size=2,
                  top_n_per_objective=2, max_pair_seeds=6, progress_every=0)
    run_staged_gaussian_reranking(catalog, query, ids, output, **kwargs)
    chunk0 = output / "coarse/chunks/chunk-000000.npz"
    chunk1 = output / "coarse/chunks/chunk-000001.npz"
    chunk0_hash, chunk0_mtime = _hash(chunk0), chunk0.stat().st_mtime_ns
    chunk1.write_bytes(b"interrupted-or-corrupt")
    run_staged_gaussian_reranking(catalog, query, ids, output, **kwargs)
    assert _hash(chunk0) == chunk0_hash
    assert chunk0.stat().st_mtime_ns == chunk0_mtime
    with np.load(chunk1, allow_pickle=False) as repaired:
        assert repaired["global_ids"].tolist() == [2, 3]
    stage_manifest = json.loads(
        (output / "coarse/merged-scores.manifest.json").read_text(encoding="utf-8"))
    assert stage_manifest["chunks_reused"] == 2
    assert stage_manifest["chunks_computed"] == 1
    with pytest.raises(ValueError, match="different inputs or parameters"):
        run_staged_gaussian_reranking(
            catalog, query, ids, output, **dict(kwargs, top_n_per_objective=3))


def test_scaled_slim_streaming_selection_matches_staged_and_resumes(tmp_path: Path):
    catalog, query, ids = _staged_inputs(tmp_path)
    staged = tmp_path / "staged-reference"
    run_staged_gaussian_reranking(
        catalog, query, ids, staged, stage="all", workers=1, chunk_size=2,
        top_n_per_objective=2, max_pair_seeds=6, progress_every=0)
    schedule = tmp_path / "schedule.npz"
    np.savez(schedule, global_ids=np.arange(6, dtype=np.int64))
    scaled = tmp_path / "scaled"
    first = run_scaled_gaussian_reranking(
        catalog, query, schedule, scaled, workers=2, coarse_chunk_size=2,
        refine_chunk_size=1, top_n_per_objective=2,
        max_pair_seeds=6, progress_every=0)
    assert first["status"] == "complete"
    assert not (scaled / "coarse/merged-scores.npz").exists()
    for chunk in sorted((scaled / "coarse/chunks").glob("*.npz")):
        with np.load(chunk, allow_pickle=False) as result:
            assert set(result.files) == {
                "global_ids", "objective_names", "objective_scores"}
            assert result["objective_scores"].dtype == np.float32
    reference = _load_arrays(staged / "refine/selected-candidates.npz")
    streamed = _load_arrays(scaled / "refine/selected-candidates.npz")
    assert np.array_equal(reference["global_ids"], streamed["global_ids"])
    assert np.array_equal(reference["selected_by_objective"],
                          streamed["selected_by_objective"])
    chunk0 = scaled / "coarse/chunks/chunk-000000.npz"
    before = (_hash(chunk0), chunk0.stat().st_mtime_ns)
    second = run_scaled_gaussian_reranking(
        catalog, query, schedule, scaled, workers=2, coarse_chunk_size=2,
        refine_chunk_size=1, top_n_per_objective=2,
        max_pair_seeds=6, progress_every=0)
    assert second["stages"]["coarse"]["chunks_reused"] == 3
    assert second["stages"]["coarse"]["chunks_computed"] == 0
    assert (_hash(chunk0), chunk0.stat().st_mtime_ns) == before
