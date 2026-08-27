from __future__ import annotations

import hashlib
import json
from pathlib import Path
import sqlite3
from typing import Any

from .project_context import ensure_within, require_active_project, require_project_owner
from .registry import stable_id, utc_now


def sha256_file(path: Path, chunk_size: int = 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(chunk_size):
            digest.update(chunk)
    return digest.hexdigest()


def canonical_hash(value: Any) -> str:
    return hashlib.sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode("utf-8")
    ).hexdigest()


def _project_root(connection: sqlite3.Connection, user_id: str, project_id: str) -> Path:
    project = require_project_owner(connection, user_id, project_id)
    require_active_project(connection, user_id, project_id)
    return Path(project["workspace_path"]).resolve()


def create_workflow_run(connection: sqlite3.Connection, user_id: str, project_id: str,
                        *, run_type: str, tool_name: str, tool_version: str,
                        parameters: dict[str, Any], execution_backend: str,
                        campaign_id: str | None = None,
                        input_artifact_ids: list[str] | None = None) -> str:
    root = _project_root(connection, user_id, project_id)
    if campaign_id:
        campaign = connection.execute(
            "SELECT id FROM campaign WHERE id=? AND project_id=?", (campaign_id, project_id)
        ).fetchone()
        if not campaign:
            raise ValueError(f"Campaign {campaign_id} does not belong to Project {project_id}")
    inputs = input_artifact_ids or []
    input_rows = []
    for artifact_id in inputs:
        row = connection.execute("SELECT * FROM artifact WHERE id=?", (artifact_id,)).fetchone()
        if not row or (row["scope"] != "shared" and row["project_id"] != project_id):
            raise ValueError(f"Artifact {artifact_id} is not accessible to Project {project_id}")
        input_rows.append(row)
    compute_payload = {
        "run_type": run_type, "tool_name": tool_name, "tool_version": tool_version,
        "parameters": parameters,
        "inputs": sorted((row["id"], row["sha256"]) for row in input_rows),
    }
    run_id = stable_id("RUN")
    run_path = ensure_within(root / "runs" / run_id, root)
    run_path.mkdir(parents=True, exist_ok=False)
    (run_path / "outputs").mkdir()
    (run_path / "parameters.json").write_text(
        json.dumps(parameters, indent=2, sort_keys=True), encoding="utf-8"
    )
    compute_key = canonical_hash(compute_payload)
    connection.execute(
        """INSERT INTO workflow_run(id, project_id, campaign_id, owner_user_id,
           run_type, status, tool_name, tool_version, execution_backend,
           parameters_json, compute_key, run_path, created_at)
           VALUES (?, ?, ?, ?, ?, 'created', ?, ?, ?, ?, ?, ?, ?)""",
        (run_id, project_id, campaign_id, user_id, run_type, tool_name, tool_version,
         execution_backend, json.dumps(parameters, sort_keys=True), compute_key,
         str(run_path), utc_now()),
    )
    for row in input_rows:
        connection.execute(
            "INSERT INTO run_input(run_id, artifact_id, role) VALUES (?, ?, ?)",
            (run_id, row["id"], row["role"]),
        )
    run_metadata = compute_payload | {"run_id": run_id, "project_id": project_id,
                                      "campaign_id": campaign_id, "status": "created"}
    (run_path / "run.json").write_text(
        json.dumps(run_metadata, indent=2, sort_keys=True), encoding="utf-8"
    )
    return run_id


def set_run_status(connection: sqlite3.Connection, user_id: str, run_id: str,
                   status: str, *, external_job_id: str | None = None,
                   error_message: str | None = None) -> None:
    row = connection.execute("SELECT * FROM workflow_run WHERE id=?", (run_id,)).fetchone()
    if not row:
        raise KeyError(f"Unknown Run: {run_id}")
    _project_root(connection, user_id, row["project_id"])
    allowed = {
        "created": {"submitted", "running", "cancelled", "failed"},
        "submitted": {"running", "cancelled", "failed"},
        "running": {"completed", "cancelled", "failed"},
    }
    if status not in allowed.get(row["status"], set()):
        raise ValueError(f"Invalid Run transition: {row['status']} -> {status}")
    now = utc_now()
    started_at = now if status == "running" and not row["started_at"] else row["started_at"]
    completed_at = now if status in {"completed", "failed", "cancelled"} else None
    connection.execute(
        """UPDATE workflow_run SET status=?, external_job_id=COALESCE(?, external_job_id),
           error_message=?, started_at=?, completed_at=? WHERE id=?""",
        (status, external_job_id, error_message, started_at, completed_at, run_id),
    )


def import_result_manifest(connection: sqlite3.Connection, user_id: str,
                           manifest_path: Path) -> dict[str, Any]:
    manifest_path = manifest_path.resolve()
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if manifest.get("schema_version") != 1:
        raise ValueError("Unsupported result manifest schema_version")
    run_id = manifest["run_id"]
    run = connection.execute("SELECT * FROM workflow_run WHERE id=?", (run_id,)).fetchone()
    if not run:
        raise KeyError(f"Unknown Run: {run_id}")
    root = _project_root(connection, user_id, run["project_id"])
    ensure_within(manifest_path, root)
    if manifest.get("project_id") != run["project_id"]:
        raise ValueError("Manifest project_id does not match the Run")
    inserted = []
    for output in manifest.get("outputs", []):
        output_path = Path(output["path"])
        if not output_path.is_absolute():
            output_path = manifest_path.parent / output_path
        output_path = ensure_within(output_path, root)
        if not output_path.is_file():
            raise FileNotFoundError(output_path)
        actual_hash = sha256_file(output_path)
        declared_hash = output.get("sha256")
        if declared_hash and declared_hash.lower() != actual_hash:
            raise ValueError(f"SHA-256 mismatch for {output_path}")
        artifact_id = stable_id("ART")
        role = output["role"]
        now = utc_now()
        existing = connection.execute(
            """SELECT id FROM artifact WHERE scope='project' AND project_id=?
               AND sha256=? AND role=?""", (run["project_id"], actual_hash, role),
        ).fetchone()
        if existing:
            artifact_id = existing["id"]
        else:
            connection.execute(
                """INSERT INTO artifact(id, project_id, run_id, scope, artifact_type,
                   role, format, storage_path, sha256, size_bytes, metadata_json,
                   created_at, verified_at) VALUES (?, ?, ?, 'project', ?, ?, ?, ?, ?, ?, ?, ?, ?)""",
                (artifact_id, run["project_id"], run_id,
                 output.get("artifact_type", role), role, output.get("format", "unknown"),
                 str(output_path), actual_hash, output_path.stat().st_size,
                 json.dumps(output.get("metadata", {}), sort_keys=True), now, now),
            )
        connection.execute(
            "INSERT OR IGNORE INTO run_output(run_id, artifact_id, role) VALUES (?, ?, ?)",
            (run_id, artifact_id, role),
        )
        inserted.append({"artifact_id": artifact_id, "role": role, "sha256": actual_hash})
    for name, value in manifest.get("metrics", {}).items():
        numeric = float(value) if isinstance(value, (int, float)) else None
        text = None if numeric is not None else str(value)
        connection.execute(
            """INSERT OR REPLACE INTO run_metric(run_id, name, numeric_value, text_value, unit)
               VALUES (?, ?, ?, ?, ?)""", (run_id, name, numeric, text, None),
        )
    current = run["status"]
    if current == "created":
        set_run_status(connection, user_id, run_id, "running")
    refreshed = connection.execute("SELECT status FROM workflow_run WHERE id=?", (run_id,)).fetchone()
    if refreshed["status"] == "running":
        set_run_status(connection, user_id, run_id, "completed")
    return {"run_id": run_id, "registered_outputs": inserted,
            "metric_count": len(manifest.get("metrics", {})), "status": "completed"}


def record_project_decision(connection: sqlite3.Connection, user_id: str,
                            project_id: str, *, decision_type: str, subject: str,
                            rationale: str, campaign_id: str | None = None,
                            payload: dict[str, Any] | None = None) -> int:
    _project_root(connection, user_id, project_id)
    if not rationale.strip():
        raise ValueError("A scientific rationale is required")
    cursor = connection.execute(
        """INSERT INTO project_decision(project_id, campaign_id, user_id,
           decision_type, subject, rationale, payload_json, created_at)
           VALUES (?, ?, ?, ?, ?, ?, ?, ?)""",
        (project_id, campaign_id, user_id, decision_type, subject, rationale.strip(),
         json.dumps(payload or {}, sort_keys=True), utc_now()),
    )
    return int(cursor.lastrowid)


def run_summary(connection: sqlite3.Connection, user_id: str, run_id: str) -> dict[str, Any]:
    run = connection.execute("SELECT * FROM workflow_run WHERE id=?", (run_id,)).fetchone()
    if not run:
        raise KeyError(f"Unknown Run: {run_id}")
    _project_root(connection, user_id, run["project_id"])
    outputs = connection.execute(
        """SELECT a.id AS artifact_id, o.role, a.format, a.storage_path, a.sha256,
           a.size_bytes FROM run_output o JOIN artifact a ON a.id=o.artifact_id
           WHERE o.run_id=? ORDER BY o.role""", (run_id,),
    ).fetchall()
    metrics = connection.execute(
        "SELECT name, numeric_value, text_value, unit FROM run_metric WHERE run_id=? ORDER BY name",
        (run_id,),
    ).fetchall()
    return {"run": dict(run), "outputs": [dict(row) for row in outputs],
            "metrics": [dict(row) for row in metrics]}
