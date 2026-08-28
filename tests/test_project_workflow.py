import json
from pathlib import Path

import pytest

from aidd_agent.campaign import create_campaign
from aidd_agent.context import AccessDeniedError, activate_task, create_task, create_user
from aidd_agent.project_context import (
    activate_project, create_scientific_project, list_owned_projects,
    project_summary,
)
from aidd_agent.provenance import (
    create_workflow_run, import_result_manifest, record_project_decision,
    run_summary,
)
from aidd_agent.registry import connect, initialize


def make_project(connection, tmp_path: Path, username: str, name: str):
    user_id = create_user(connection, username)
    project_id = create_scientific_project(
        connection, user_id, name, f"Objective for {name}", tmp_path / "aidd-storage"
    )
    return user_id, project_id


def test_project_is_the_active_isolation_boundary(tmp_path: Path) -> None:
    database = tmp_path / "aidd.sqlite3"
    initialize(database)
    with connect(database) as connection:
        alice, wee1 = make_project(connection, tmp_path, "alice", "WEE1 Macrocycles")
        bob, kras = make_project(connection, tmp_path, "bob", "KRAS G12D")
        assert [row["project_id"] for row in list_owned_projects(connection, alice)] == [wee1]
        with pytest.raises(AccessDeniedError):
            activate_project(connection, alice, kras)
        activate_project(connection, alice, wee1)
        task_id = create_task(connection, alice, wee1, "Structure Selection",
                              "Select structures for the ATP site")
        activate_task(connection, alice, task_id)
        campaign_id = create_campaign(
            connection, alice, wee1, task_id, "ATP Site", "Screen the ATP site")
        assert campaign_id.startswith("CAM-")
        with pytest.raises(AccessDeniedError):
            project_summary(connection, bob, wee1)


def test_only_executed_methods_appear_in_project_summary(tmp_path: Path) -> None:
    database = tmp_path / "aidd.sqlite3"
    initialize(database)
    with connect(database) as connection:
        user_id, project_id = make_project(connection, tmp_path, "alice", "Selective Screen")
        activate_project(connection, user_id, project_id)
        create_workflow_run(
            connection, user_id, project_id, run_type="docking", tool_name="PLANTS",
            tool_version="1.2", parameters={"poses": 5}, execution_backend="slurm",
        )
        decision_id = record_project_decision(
            connection, user_id, project_id, decision_type="method_skipped",
            subject="md_simulation", rationale="No stable pose passed the prespecified gate.",
        )
        summary = project_summary(connection, user_id, project_id)
    assert summary["executed_methods"] == ["docking"]
    assert decision_id > 0
    assert "md_simulation" not in summary["executed_methods"]


def test_result_manifest_registers_verified_lineage_and_metrics(tmp_path: Path) -> None:
    database = tmp_path / "aidd.sqlite3"
    initialize(database)
    with connect(database) as connection:
        user_id, project_id = make_project(connection, tmp_path, "alice", "Property Prediction")
        active = activate_project(connection, user_id, project_id)
        run_id = create_workflow_run(
            connection, user_id, project_id, run_type="model_prediction",
            tool_name="ResAddNet", tool_version="v9",
            parameters={"property": "predicted_pKi"}, execution_backend="local",
        )
        run_path = Path(active["project_root"]) / "runs" / run_id
        output = run_path / "outputs" / "predictions.csv"
        output.write_text("molecule_id,predicted_pKi\nMOL-1,8.2\n", encoding="utf-8")
        manifest = run_path / "result_manifest.json"
        manifest.write_text(json.dumps({
            "schema_version": 1, "run_id": run_id, "project_id": project_id,
            "outputs": [{"role": "predictions", "artifact_type": "prediction_table",
                         "path": "outputs/predictions.csv", "format": "csv"}],
            "metrics": {"successful_predictions": 1, "failed_predictions": 0},
        }), encoding="utf-8")
        imported = import_result_manifest(connection, user_id, manifest)
        summary = run_summary(connection, user_id, run_id)
    assert imported["status"] == "completed"
    assert summary["run"]["status"] == "completed"
    assert summary["outputs"][0]["role"] == "predictions"
    assert {row["name"] for row in summary["metrics"]} == {
        "successful_predictions", "failed_predictions"
    }


def test_manifest_cannot_register_output_outside_project(tmp_path: Path) -> None:
    database = tmp_path / "aidd.sqlite3"
    initialize(database)
    with connect(database) as connection:
        user_id, project_id = make_project(connection, tmp_path, "alice", "Path Safety")
        active = activate_project(connection, user_id, project_id)
        run_id = create_workflow_run(
            connection, user_id, project_id, run_type="rescoring", tool_name="test",
            tool_version="1", parameters={}, execution_backend="local",
        )
        outside = tmp_path / "secret.csv"
        outside.write_text("private", encoding="utf-8")
        manifest = Path(active["project_root"]) / "runs" / run_id / "result_manifest.json"
        manifest.write_text(json.dumps({
            "schema_version": 1, "run_id": run_id, "project_id": project_id,
            "outputs": [{"role": "scores", "path": str(outside), "format": "csv"}],
        }), encoding="utf-8")
        with pytest.raises(ValueError, match="escapes allowed root"):
            import_result_manifest(connection, user_id, manifest)
