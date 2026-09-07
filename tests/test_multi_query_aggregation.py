import json
from pathlib import Path

import numpy as np
import pytest

from aidd_agent.multi_query_aggregation import aggregate_multi_cocrystal_results


OBJECTIVES = ("shape_only", "atomcentered_unweighted_joint",
              "atomcentered_anchored_joint")


def _result(path: Path, rows, primary_scores, shape_scores=None):
    count = len(rows)
    arrays = {
        "global_ids": np.asarray([row[0] for row in rows], dtype=np.int64),
        "molecule_ids": np.asarray([row[1] for row in rows], dtype="U16"),
        "conformer_ids": np.asarray([row[2] for row in rows], dtype="U16"),
        "objective_names": np.asarray(OBJECTIVES, dtype="U40"),
        "best_seed_ids": np.full((count, 3), "pca:0", dtype="U96"),
    }
    shape_scores = primary_scores if shape_scores is None else shape_scores
    for column, objective in enumerate(OBJECTIVES):
        values = shape_scores if column == 0 else primary_scores
        arrays[f"{objective}__objective"] = np.asarray(values, dtype=np.float64)
        transforms = np.tile(np.eye(4).reshape(1, 16), (count, 1))
        transforms[:, 3] = np.arange(count) + column
        arrays[f"{objective}__transform"] = transforms
    np.savez(path, **arrays)
    return path


def _plan(path: Path, q1: Path, q2: Path):
    document = {
        "format": "aidd-multi-cocrystal-query-plan", "version": 1,
        "library_id": "LIB-X",
        "queries": [
            {"query_id": "Q1", "receptor_id": "R1", "site_id": "SITE-A",
             "result": str(q1)},
            {"query_id": "Q2", "receptor_id": "R2", "site_id": "SITE-A",
             "result": str(q2)},
        ],
    }
    path.write_text(json.dumps(document), encoding="utf-8")
    return path


def _jsonl(path: Path):
    return [json.loads(line) for line in path.read_text(encoding="utf-8").splitlines()]


def test_multi_query_union_conformer_collapse_and_docking_tasks(tmp_path: Path):
    q1 = _result(
        tmp_path / "q1.npz",
        [(1, "M1", "C1"), (2, "M1", "C2"), (3, "M2", "C3")],
        [.80, .90, .85], shape_scores=[.95, .70, .80])
    q2 = _result(
        tmp_path / "q2.npz", [(4, "M2", "C4"), (5, "M3", "C5")],
        [.95, .80])
    plan = _plan(tmp_path / "plan.json", q1, q2)
    output = tmp_path / "aggregated"
    manifest = aggregate_multi_cocrystal_results(
        plan, output, top_conformers_per_query=2, per_query_quota=1,
        consensus_quota=3, global_limit_per_site=4)
    assert manifest["status"] == "complete"
    evidence = _jsonl(output / "query-molecule-evidence.jsonl")
    m1 = next(row for row in evidence if row["query_id"] == "Q1"
              and row["molecule_id"] == "M1")
    assert m1["objective_best"]["atomcentered_anchored_joint"]["global_id"] == 2
    assert m1["objective_best"]["shape_only"]["global_id"] == 1
    assert [row["global_id"] for row in m1["primary_conformers"]] == [2, 1]
    summaries = _jsonl(output / "molecule-summary.jsonl")
    assert {row["molecule_id"] for row in summaries} == {"M1", "M2", "M3"}
    assert next(row for row in summaries if row["molecule_id"] == "M2")[
        "query_support_count"] == 2
    admissions = _jsonl(output / "docking-admission.jsonl")
    assert {row["molecule_id"] for row in admissions} == {"M1", "M2", "M3"}
    tasks = _jsonl(output / "docking-tasks.jsonl")
    assert len(tasks) == 4
    assert {(row["molecule_id"], row["receptor_id"]) for row in tasks} == {
        ("M1", "R1"), ("M2", "R1"), ("M2", "R2"), ("M3", "R2")}
    first_hashes = {name: row["sha256"] for name, row in manifest["outputs"].items()}
    repeated = aggregate_multi_cocrystal_results(
        plan, output, top_conformers_per_query=2, per_query_quota=1,
        consensus_quota=3, global_limit_per_site=4)
    assert first_hashes == {name: row["sha256"]
                            for name, row in repeated["outputs"].items()}


def test_rank_fusion_is_invariant_to_monotonic_score_rescaling(tmp_path: Path):
    rows1 = [(1, "M1", "C1"), (2, "M2", "C2")]
    rows2 = [(3, "M2", "C3"), (4, "M3", "C4")]
    q1 = _result(tmp_path / "q1.npz", rows1, [.9, .8])
    q2 = _result(tmp_path / "q2.npz", rows2, [.7, .6])
    first_plan = _plan(tmp_path / "plan1.json", q1, q2)
    aggregate_multi_cocrystal_results(
        first_plan, tmp_path / "out1", per_query_quota=1,
        consensus_quota=3, global_limit_per_site=3)
    q1_scaled = _result(tmp_path / "q1-scaled.npz", rows1, [90.0, 80.0])
    second_plan = _plan(tmp_path / "plan2.json", q1_scaled, q2)
    aggregate_multi_cocrystal_results(
        second_plan, tmp_path / "out2", per_query_quota=1,
        consensus_quota=3, global_limit_per_site=3)
    order1 = [row["molecule_id"] for row in _jsonl(tmp_path / "out1/molecule-summary.jsonl")]
    order2 = [row["molecule_id"] for row in _jsonl(tmp_path / "out2/molecule-summary.jsonl")]
    assert order1 == order2


def test_query_protection_rejects_an_impossible_site_limit(tmp_path: Path):
    q1 = _result(tmp_path / "q1.npz", [(1, "M1", "C1")], [.9])
    q2 = _result(tmp_path / "q2.npz", [(2, "M2", "C2")], [.8])
    plan = _plan(tmp_path / "plan.json", q1, q2)
    with pytest.raises(ValueError, match="cannot guarantee"):
        aggregate_multi_cocrystal_results(
            plan, tmp_path / "out", per_query_quota=2,
            consensus_quota=1, global_limit_per_site=3)


def test_different_sites_are_aggregated_independently(tmp_path: Path):
    q1 = _result(tmp_path / "q1.npz", [(1, "M1", "C1")], [.9])
    q2 = _result(tmp_path / "q2.npz", [(2, "M1", "C2")], [.8])
    plan = {
        "format": "aidd-multi-cocrystal-query-plan", "version": 1,
        "library_id": "LIB-X", "queries": [
            {"query_id": "Q1", "receptor_id": "R1", "site_id": "SITE-A",
             "result": str(q1)},
            {"query_id": "Q2", "receptor_id": "R2", "site_id": "SITE-B",
             "result": str(q2)},
        ]}
    plan_path = tmp_path / "plan.json"
    plan_path.write_text(json.dumps(plan), encoding="utf-8")
    output = tmp_path / "out"
    aggregate_multi_cocrystal_results(
        plan_path, output, per_query_quota=1,
        consensus_quota=1, global_limit_per_site=1)
    summaries = _jsonl(output / "molecule-summary.jsonl")
    assert [(row["site_id"], row["molecule_id"]) for row in summaries] == [
        ("SITE-A", "M1"), ("SITE-B", "M1")]


def test_declared_result_checksum_is_enforced(tmp_path: Path):
    q1 = _result(tmp_path / "q1.npz", [(1, "M1", "C1")], [.9])
    bad_manifest = tmp_path / "q1.manifest.json"
    bad_manifest.write_text(json.dumps({"result_sha256": "0" * 64}), encoding="utf-8")
    plan = {
        "format": "aidd-multi-cocrystal-query-plan", "version": 1,
        "library_id": "LIB-X", "queries": [{
            "query_id": "Q1", "receptor_id": "R1", "site_id": "SITE-A",
            "result": str(q1), "result_manifest": str(bad_manifest)}]}
    plan_path = tmp_path / "plan.json"
    plan_path.write_text(json.dumps(plan), encoding="utf-8")
    with pytest.raises(ValueError, match="checksum mismatch"):
        aggregate_multi_cocrystal_results(
            plan_path, tmp_path / "out", per_query_quota=1,
            consensus_quota=1, global_limit_per_site=1)
