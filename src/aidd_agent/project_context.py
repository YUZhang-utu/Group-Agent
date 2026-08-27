from __future__ import annotations

import json
from pathlib import Path
import sqlite3
from typing import Any

from .context import AccessDeniedError, _require_user, slugify
from .registry import stable_id, utc_now


def ensure_within(path: Path, root: Path) -> Path:
    resolved = path.resolve()
    resolved_root = root.resolve()
    if resolved != resolved_root and resolved_root not in resolved.parents:
        raise ValueError(f"Path escapes allowed root: {resolved}")
    return resolved


def initialize_storage_root(storage_root: Path) -> dict[str, str]:
    root = storage_root.resolve()
    paths = {
        "root": root,
        "registry": root / "registry",
        "shared": root / "shared",
        "artifacts": root / "artifacts",
        "users": root / "users",
    }
    for path in paths.values():
        path.mkdir(parents=True, exist_ok=True)
    for relative in ("libraries", "models", "software", "environments", "receptors"):
        (paths["shared"] / relative).mkdir(exist_ok=True)
    return {name: str(path) for name, path in paths.items()}


def _write_project_workspace(path: Path, metadata: dict[str, str]) -> None:
    path.mkdir(parents=True, exist_ok=False)
    for directory in (
        "target", "target/pdb_candidates", "target/prepared", "campaigns",
        "runs", "results", "reports", "data", "tmp",
    ):
        (path / directory).mkdir(parents=True)
    (path / "project.yaml").write_text(
        "\n".join([
            "schema_version: 1", f"project_id: {metadata['project_id']}",
            f"project_name: {metadata['name']}", f"owner_user_id: {metadata['owner_user_id']}",
            f"objective: {metadata['objective']}", "status: created",
            "active_campaign_id: null", "current_stage: project_initialization", "",
        ]), encoding="utf-8",
    )
    templates = {
        "research-state.yaml": "stage: project_initialization\ncurrent_experiment: null\nnext_actions:\n  - define the biological target and desired state\n",
        "research_log.md": "# Research Log\n\nAppend-only record for this Project.\n",
        "findings.md": "# Findings\n\n## Current Understanding\n\nNo Project-specific findings yet.\n",
        "experiments.md": "# Experiments\n\nOnly executed or explicitly preregistered experiments are recorded here.\n",
        "benchmarks.md": "# Benchmarks\n\nNo Project-specific benchmarks yet.\n",
        "decisions.md": "# Decisions\n\nOnly scientifically meaningful decisions are recorded here.\n",
        "target/structure_selection.md": "# Target Structure Selection\n\nNo PDB structure has been selected.\n",
    }
    for relative, content in templates.items():
        (path / relative).write_text(content, encoding="utf-8")


def create_scientific_project(connection: sqlite3.Connection, user_id: str, name: str,
                              objective: str, storage_root: Path,
                              description: str = "") -> str:
    user = _require_user(connection, user_id)
    if not name.strip() or not objective.strip():
        raise ValueError("Project name and objective are required")
    roots = initialize_storage_root(storage_root)
    project_id = stable_id("PRJ")
    project_slug = slugify(name)
    workspace = ensure_within(
        Path(roots["users"]) / user["username"] / "projects" /
        f"{project_id.lower()}-{project_slug}", Path(roots["users"]),
    )
    _write_project_workspace(workspace, {
        "project_id": project_id, "name": name.strip(), "owner_user_id": user_id,
        "objective": objective.strip(),
    })
    now = utc_now()
    connection.execute(
        """INSERT INTO project(id, owner_user_id, name, slug, description, objective, status,
           workspace_path, updated_at, created_at) VALUES (?, ?, ?, ?, ?, ?, 'created', ?, ?, ?)""",
        (project_id, user_id, name.strip(), project_slug, description.strip(), objective.strip(),
         str(workspace), now, now),
    )
    return project_id


def require_project_owner(connection: sqlite3.Connection, user_id: str,
                          project_id: str) -> sqlite3.Row:
    row = connection.execute(
        "SELECT * FROM project WHERE id=? AND owner_user_id=?",
        (project_id, user_id),
    ).fetchone()
    if not row:
        raise AccessDeniedError(f"User {user_id} does not own Project {project_id}")
    if not row["workspace_path"]:
        raise ValueError(f"Project {project_id} has no project-centric workspace")
    return row


def activate_project(connection: sqlite3.Connection, user_id: str,
                     project_id: str) -> dict[str, Any]:
    project = require_project_owner(connection, user_id, project_id)
    now = utc_now()
    connection.execute(
        """INSERT INTO active_project_context(user_id, project_id, activated_at)
           VALUES (?, ?, ?) ON CONFLICT(user_id) DO UPDATE SET
           project_id=excluded.project_id, activated_at=excluded.activated_at""",
        (user_id, project_id, now),
    )
    connection.execute(
        "UPDATE project SET status='active', updated_at=? WHERE id=?",
        (now, project_id),
    )
    return {"user_id": user_id, "project_id": project_id,
            "project_root": project["workspace_path"]}


def active_project(connection: sqlite3.Connection, user_id: str) -> dict[str, Any] | None:
    row = connection.execute(
        """SELECT p.id AS project_id, p.name, p.status, p.workspace_path AS project_root
           FROM active_project_context c JOIN project p ON p.id=c.project_id
           WHERE c.user_id=? AND p.owner_user_id=?""", (user_id, user_id),
    ).fetchone()
    return dict(row) if row else None


def require_active_project(connection: sqlite3.Connection, user_id: str,
                           expected_project_id: str | None = None) -> dict[str, Any]:
    current = active_project(connection, user_id)
    if not current:
        raise AccessDeniedError(
            f"User {user_id} has no active Project; activate a Project before scientific work"
        )
    if expected_project_id and current["project_id"] != expected_project_id:
        raise AccessDeniedError(
            f"Active Project {current['project_id']} does not match required Project {expected_project_id}"
        )
    return current


def deactivate_project(connection: sqlite3.Connection, user_id: str) -> None:
    connection.execute("DELETE FROM active_project_context WHERE user_id=?", (user_id,))


def list_owned_projects(connection: sqlite3.Connection, user_id: str) -> list[dict[str, Any]]:
    _require_user(connection, user_id)
    rows = connection.execute(
        """SELECT id AS project_id, name, objective, status, workspace_path AS project_root,
           created_at, updated_at FROM project WHERE owner_user_id=? AND workspace_path IS NOT NULL
           ORDER BY COALESCE(updated_at, created_at) DESC, id""", (user_id,),
    ).fetchall()
    return [dict(row) for row in rows]


def project_summary(connection: sqlite3.Connection, user_id: str,
                    project_id: str) -> dict[str, Any]:
    project = require_project_owner(connection, user_id, project_id)
    runs = connection.execute(
        """SELECT id AS run_id, run_type, status, tool_name, tool_version,
           execution_backend, created_at, completed_at FROM workflow_run
           WHERE project_id=? ORDER BY created_at, id""", (project_id,),
    ).fetchall()
    campaigns = connection.execute(
        "SELECT id AS campaign_id, name, state FROM campaign WHERE project_id=? ORDER BY created_at",
        (project_id,),
    ).fetchall()
    return {
        "project_id": project_id, "name": project["name"], "objective": project["objective"],
        "status": project["status"], "project_root": project["workspace_path"],
        "executed_methods": sorted({row["run_type"] for row in runs}),
        "campaigns": [dict(row) for row in campaigns], "runs": [dict(row) for row in runs],
    }
