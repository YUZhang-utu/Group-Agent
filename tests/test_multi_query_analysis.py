import json
from pathlib import Path

import numpy as np
import pytest

from aidd_agent.multi_query_aggregation import aggregate_multi_cocrystal_results
from aidd_agent.multi_query_analysis import analyze_multi_cocrystal_result


OBJECTIVES = ("shape_only", "atomcentered_unweighted_joint",
              "atomcentered_anchored_joint")


def _result(path: Path, rows, scores, offsets=(0, 0, 1)) -> Path:
    count = len(rows)
    arrays = {
        "global_ids": np.asarray([row[0] for row in rows], dtype=np.int64),
        "molecule_ids": np.asarray([row[1] for row in rows], dtype="U16"),
        "conformer_ids": np.asarray([row[2] for row in rows], dtype="U16"),
        "objective_names": np.asarray(OBJECTIVES, dtype="U40"),
        "best_seed_ids": np.full((count, 3), "pca:0", dtype="U96"),
    }
    for column, objective in enumerate(OBJECTIVES):
        arrays[f"{objective}__objective"] = np.asarray(scores, dtype=np.float64)
        transforms = np.tile(np.eye(4).reshape(1, 16), (count, 1))
        transforms[:, 3] = offsets[column]
        arrays[f"{objective}__transform"] = transforms
    np.savez(path, **arrays)
    return path


def _aggregation(tmp_path: Path) -> Path:
    q1 = _result(tmp_path / "q1.npz",
                 [(1, "M1", "C1"), (2, "M2", "C2")], [.9, .8])
    q2 = _result(tmp_path / "q2.npz",
                 [(3, "M2", "C3"), (4, "M3", "C4")], [.95, .7])
    plan = {"format": "aidd-multi-cocrystal-query-plan", "version": 1,
            "library_id": "LIB-X", "queries": [
                {"query_id": "Q1", "receptor_id": "R1", "site_id": "SITE",
                 "result": str(q1)},
                {"query_id": "Q2", "receptor_id": "R2", "site_id": "SITE",
                 "result": str(q2)}]}
    plan_path = tmp_path / "plan.json"
    plan_path.write_text(json.dumps(plan), encoding="utf-8")
    result = tmp_path / "aggregation"
    aggregate_multi_cocrystal_results(
        plan_path, result, per_query_quota=1,
        consensus_quota=3, global_limit_per_site=3)
    return result


def test_overlap_analysis_recomputes_union_and_queue(tmp_path: Path):
    result = _aggregation(tmp_path)
    report = analyze_multi_cocrystal_result(
        result, tmp_path / "analysis", top_k=(1, 2))

    assert report["accepted"] is True
    assert report["molecules"] == {
        "union": 3, "shared_all": 1, "shared_all_fraction": 1 / 3}
    assert report["queries"]["Q1"]["exclusive_molecules"] == 1
    assert report["queries"]["Q2"]["exclusive_molecules"] == 1
    assert [row["intersection"] for row in report["pairwise"][0]["top_k"]] == [0, 1]
    assert report["docking_queue"]["admitted_molecules"] == 3
    assert report["docking_queue"]["docking_tasks"] == 4
    assert report["docking_queue"]["multi_query_admitted"] == 1
    assert report["docking_queue"]["task_counts_by_query"] == {"Q1": 2, "Q2": 2}
    assert (tmp_path / "analysis/analysis.json").is_file()
    assert "Pairwise Top-K overlap" in (
        tmp_path / "analysis/analysis.md").read_text(encoding="utf-8")


def test_overlap_analysis_rejects_modified_aggregation_output(tmp_path: Path):
    result = _aggregation(tmp_path)
    with (result / "docking-tasks.jsonl").open("a", encoding="utf-8") as stream:
        stream.write("{}\n")
    with pytest.raises(ValueError, match="checksum mismatch"):
        analyze_multi_cocrystal_result(result, tmp_path / "analysis")
