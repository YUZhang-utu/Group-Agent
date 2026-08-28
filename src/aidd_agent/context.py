from __future__ import annotations

import re
import sqlite3
from pathlib import Path
from typing import Any

from .registry import stable_id, utc_now


class AccessDeniedError(PermissionError):
    pass


def slugify(value: str) -> str:
    slug = re.sub(r"[^a-z0-9]+", "-", value.strip().lower()).strip("-")
    if not slug:
        raise ValueError("A name must contain at least one ASCII letter or digit")
    return slug


def create_user(connection: sqlite3.Connection, username: str,
                display_name: str | None = None) -> str:
    username = slugify(username)
    user_id = stable_id("USR")
    connection.execute(
        "INSERT INTO app_user(id, username, display_name, created_at) VALUES (?, ?, ?, ?)",
        (user_id, username, (display_name or username).strip(), utc_now()),
    )
    return user_id


def create_project(connection: sqlite3.Connection, owner_user_id: str, name: str,
                   description: str = "") -> str:
    _require_user(connection, owner_user_id)
    project_id = stable_id("PRJ")
    connection.execute(
        """INSERT INTO project(id, owner_user_id, name, slug, description, created_at)
           VALUES (?, ?, ?, ?, ?, ?)""",
        (project_id, owner_user_id, name.strip(), slugify(name), description.strip(), utc_now()),
    )
    return project_id


def _require_user(connection: sqlite3.Connection, user_id: str) -> sqlite3.Row:
    row = connection.execute("SELECT * FROM app_user WHERE id=?", (user_id,)).fetchone()
    if not row:
        raise KeyError(f"Unknown user: {user_id}")
    return row


def require_task_access(connection: sqlite3.Connection, user_id: str, task_id: str,
                        *, write: bool = False) -> sqlite3.Row:
    row = connection.execute(
        """SELECT t.*, a.role FROM task t JOIN task_access a ON a.task_id=t.id
           WHERE t.id=? AND a.user_id=?""", (task_id, user_id),
    ).fetchone()
    if not row:
        raise AccessDeniedError(f"User {user_id} cannot access task {task_id}")
    if write and row["role"] == "viewer":
        raise AccessDeniedError(f"User {user_id} has read-only access to task {task_id}")
    return row


def _write_task_workspace(path: Path, task: dict[str, str]) -> None:
    path.mkdir(parents=True, exist_ok=False)
    for directory in ("structures/candidates", "structures/prepared", "docking",
                      "models", "md", "data", "reports"):
        (path / directory).mkdir(parents=True)
    (path / "task.yaml").write_text(
        "\n".join([f"task_id: {task['task_id']}", f"task_name: {task['name']}",
                    f"project_id: {task['project_id']}", f"owner_id: {task['owner_id']}",
                    f"objective: {task['objective']}", "status: created",
                    "active_campaign_id: null", "current_stage: task_initialization", ""]),
        encoding="utf-8",
    )
    templates = {
        "research-state.yaml": "stage: task_initialization\ncurrent_experiment: null\nnext_actions:\n  - define the first campaign and biological target\n",
        "research_log.md": "# Research Log\n\nAppend-only record for this task.\n",
        "findings.md": "# Findings\n\n## Current Understanding\n\nNo task-specific findings yet.\n",
        "experiments.md": "# Experiments\n\nNo task-specific experiments yet.\n",
        "benchmarks.md": "# Benchmarks\n\nNo task-specific benchmarks yet.\n",
        "structures/selection.md": "# Target Structure Selection\n\nNo PDB structure has been selected.\n",
    }
    for relative, content in templates.items():
        (path / relative).write_text(content, encoding="utf-8")


def create_task(connection: sqlite3.Connection, user_id: str, project_id: str,
                name: str, objective: str, workspace_root: Path | None = None) -> str:
    from .project_context import require_active_project, require_project_owner

    project = require_project_owner(connection, user_id, project_id)
    require_active_project(connection, user_id, project_id)
    if not name.strip() or not objective.strip():
        raise ValueError("Task name and objective are required")
    task_id = stable_id("TSK")
    task_slug = slugify(name)
    if workspace_root is None:
        workspace = (
            Path(project["workspace_path"]) / "tasks" /
            f"{task_id.lower()}-{task_slug}"
        ).resolve()
    else:
        user = _require_user(connection, user_id)
        workspace = (workspace_root / user["username"] / project["slug"] / task_slug).resolve()
    _write_task_workspace(workspace, {"task_id": task_id, "name": name.strip(),
                                     "project_id": project_id, "owner_id": user_id,
                                     "objective": objective.strip()})
    now = utc_now()
    try:
        connection.execute(
            """INSERT INTO task(id, project_id, name, slug, objective, status,
               workspace_path, created_at, updated_at) VALUES (?, ?, ?, ?, ?, 'created', ?, ?, ?)""",
            (task_id, project_id, name.strip(), task_slug, objective.strip(), str(workspace), now, now),
        )
        connection.execute(
            "INSERT INTO task_access(task_id, user_id, role, created_at) VALUES (?, ?, 'owner', ?)",
            (task_id, user_id, now),
        )
    except Exception:
        # The just-created empty workspace is intentionally left visible for recovery.
        raise
    return task_id


def list_tasks(connection: sqlite3.Connection, user_id: str) -> list[dict[str, Any]]:
    _require_user(connection, user_id)
    rows = connection.execute(
        """SELECT t.id AS task_id, t.name, t.status, t.objective, t.workspace_path,
           p.id AS project_id, p.name AS project_name, a.role
           FROM task_access a JOIN task t ON t.id=a.task_id
           JOIN project p ON p.id=t.project_id WHERE a.user_id=?
           ORDER BY t.updated_at DESC, t.id""", (user_id,),
    ).fetchall()
    return [dict(row) for row in rows]


def activate_task(connection: sqlite3.Connection, user_id: str, task_id: str) -> dict[str, Any]:
    from .project_context import require_active_project

    task = require_task_access(connection, user_id, task_id)
    require_active_project(connection, user_id, task["project_id"])
    connection.execute(
        """INSERT INTO active_task_context(user_id, task_id, activated_at) VALUES (?, ?, ?)
           ON CONFLICT(user_id) DO UPDATE SET task_id=excluded.task_id,
           activated_at=excluded.activated_at""", (user_id, task_id, utc_now()),
    )
    if task["status"] == "created":
        connection.execute("UPDATE task SET status='active', updated_at=? WHERE id=?",
                           (utc_now(), task_id))
    return {"user_id": user_id, "task_id": task_id,
            "workspace_path": task["workspace_path"], "role": task["role"]}


def active_task(connection: sqlite3.Connection, user_id: str) -> dict[str, Any] | None:
    row = connection.execute(
        """SELECT t.id AS task_id, t.name, t.status, t.workspace_path, a.role
           FROM active_task_context c JOIN task t ON t.id=c.task_id
           JOIN task_access a ON a.task_id=t.id AND a.user_id=c.user_id
           WHERE c.user_id=?""", (user_id,),
    ).fetchone()
    return dict(row) if row else None


def require_active_task(connection: sqlite3.Connection, user_id: str,
                        expected_task_id: str | None = None) -> dict[str, Any]:
    current = active_task(connection, user_id)
    if not current:
        raise AccessDeniedError(
            f"User {user_id} has no active task; activate a task before scientific work"
        )
    if expected_task_id and current["task_id"] != expected_task_id:
        raise AccessDeniedError(
            f"Active task {current['task_id']} does not match required task {expected_task_id}"
        )
    return current


def deactivate_task(connection: sqlite3.Connection, user_id: str) -> None:
    connection.execute("DELETE FROM active_task_context WHERE user_id=?", (user_id,))
