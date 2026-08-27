from pathlib import Path

import pytest

from aidd_agent.context import (
    AccessDeniedError, activate_task, active_task, create_project, create_task,
    create_user, list_tasks, require_task_access,
)
from aidd_agent.registry import connect, initialize


def test_user_tasks_are_isolated_and_workspaces_are_separate(tmp_path: Path) -> None:
    database = tmp_path / "aidd.sqlite3"
    initialize(database)
    with connect(database) as connection:
        alice = create_user(connection, "alice", "Alice")
        bob = create_user(connection, "bob", "Bob")
        alice_project = create_project(connection, alice, "Kinase Program")
        bob_project = create_project(connection, bob, "Oncology Program")
        task_a = create_task(connection, alice, alice_project, "WEE1 ATP Site",
                             "Screen macrocycles against WEE1", tmp_path / "workspaces")
        task_b = create_task(connection, bob, bob_project, "KRAS Switch II",
                             "Screen macrocycles against KRAS G12D", tmp_path / "workspaces")
        assert [row["task_id"] for row in list_tasks(connection, alice)] == [task_a]
        assert [row["task_id"] for row in list_tasks(connection, bob)] == [task_b]
        with pytest.raises(AccessDeniedError):
            activate_task(connection, alice, task_b)
        activated = activate_task(connection, alice, task_a)
        assert active_task(connection, alice)["task_id"] == task_a
    workspace_a = Path(activated["workspace_path"])
    assert workspace_a != tmp_path / "workspaces" / "bob" / "oncology-program" / "kras-switch-ii"
    assert (workspace_a / "task.yaml").is_file()
    assert (workspace_a / "structures" / "selection.md").is_file()
    assert "No task-specific findings yet" in (workspace_a / "findings.md").read_text()


def test_schema_migration_preserves_library_data(tmp_path: Path) -> None:
    database = tmp_path / "aidd.sqlite3"
    initialize(database)
    with connect(database) as connection:
        connection.execute("INSERT INTO library VALUES ('LIB-X', 'shared-library', 'now')")
    initialize(database)
    with connect(database) as connection:
        assert connection.execute("SELECT COUNT(*) FROM library").fetchone()[0] == 1
        assert connection.execute("SELECT COUNT(*) FROM app_user").fetchone()[0] == 0
