import hashlib
import json

import numpy as np

from aidd_agent.conformer_artifacts import FEATURE_DTYPE, META_DTYPE
from aidd_agent.pharmacophore_index import (
    build_pharmacophore_index, compile_pharmacophore_query,
    search_pharmacophore_index,
)


def _sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _artifact_catalog(tmp_path):
    shard = tmp_path / "artifacts" / "s0"; shard.mkdir(parents=True)
    # gid 0: HBA--HBD at 4.0 A and HBA--HBA at 8.0 A (strict three-anchor match)
    # gid 1: only an HBA--HBD pair at 5.5 A (loose partial match)
    rows = np.zeros(2, dtype=META_DTYPE)
    rows["global_id"] = [0, 1]; rows["feature_offset"] = [0, 3]
    rows["features"] = [3, 2]; rows["origin"] = [[0, 0, 0], [0, 0, 0]]
    rows.tofile(shard / "meta.bin")
    features = np.zeros(5, dtype=FEATURE_DTYPE)
    features["xyz"] = np.asarray([[0, 0, 0], [400, 0, 0], [800, 0, 0],
                                  [0, 0, 0], [550, 0, 0]], dtype=np.int16)
    features["type"] = [2, 1, 2, 2, 1]
    features.tofile(shard / "feats.bin")
    source = {
        "format": "aidd-conformer-artifact-shard", "version": 1,
        "library_id": "LIB-X", "global_id_start": 0, "conformers": 2,
        "coordinate_scale": 100.0,
        "feature_types": {"Donor": 1, "Acceptor": 2}, "files": {},
    }
    (shard / "manifest.json").write_text(json.dumps(source))
    catalog = {
        "format": "aidd-conformer-artifact-catalog", "version": 1,
        "library_id": "LIB-X", "conformers": 2,
        "shards": [{"name": "s0", "path": str(shard.resolve()),
                    "global_id_start": 0, "conformers": 2}],
    }
    path = tmp_path / "artifacts" / "catalog.json"
    path.write_text(json.dumps(catalog)); return path


def _query(tmp_path):
    query = {
        "query_id": "Q", "anchors": [
            {"anchor_id": "a", "feature_type": "HBA", "atom_center": [10, 0, 0]},
            {"anchor_id": "b", "feature_type": "HBD", "atom_center": [14, 0, 0]},
            {"anchor_id": "c", "feature_type": "HBA", "atom_center": [18, 0, 0]},
        ]}
    path = tmp_path / "query.json"; path.write_text(json.dumps(query)); return path


def test_one_time_shard_index_and_nested_query_tiers(tmp_path):
    artifact_catalog = _artifact_catalog(tmp_path)
    index_root = tmp_path / "index"
    first = build_pharmacophore_index(artifact_catalog, index_root)
    manifest_before = (index_root / "shards" / "s0" / "manifest.json").read_bytes()
    second = build_pharmacophore_index(artifact_catalog, index_root)
    assert first["conformers"] == second["conformers"] == 2
    assert (index_root / "shards" / "s0" / "manifest.json").read_bytes() == manifest_before

    query_plan = tmp_path / "query-plan.json"
    compiled = compile_pharmacophore_query(_query(tmp_path), query_plan)
    assert len(compiled["pairs"]) == 3
    output = tmp_path / "hits.npz"
    result = search_pharmacophore_index(index_root / "catalog.json", query_plan,
                                         output, external_l1_ids=[99])
    assert result["profile_counts"] == {"loose": 2, "balanced": 1, "strict": 1}
    arrays = np.load(output)
    assert arrays["global_ids"].tolist() == [0, 1, 99]
    assert arrays["highest_tier"].tolist() == [3, 1, 0]
    assert arrays["source_flags"].tolist() == [1, 1, 2]


def test_query_compilation_is_library_independent_and_hashed(tmp_path):
    plan = tmp_path / "plan.json"
    result = compile_pharmacophore_query(_query(tmp_path), plan)
    assert result["supported_anchors"] == 3
    assert len(result["query_hash"]) == 64
    assert "library" not in json.dumps(result).lower()
    copied = tmp_path / "elsewhere" / "same-query.json"
    copied.parent.mkdir(); copied.write_bytes(_query(tmp_path).read_bytes())
    second = compile_pharmacophore_query(copied, tmp_path / "second-plan.json")
    assert second["query_hash"] == result["query_hash"]


def test_one_anchor_query_does_not_fabricate_spatial_hits(tmp_path):
    query = {"query_id": "one", "anchors": [
        {"anchor_id": "a", "feature_type": "HBA", "atom_center": [0, 0, 0]}]}
    source = tmp_path / "one.json"; source.write_text(json.dumps(query))
    plan = compile_pharmacophore_query(source, tmp_path / "one-plan.json")
    assert plan["supported_anchors"] == 1
    assert plan["pairs"] == []
