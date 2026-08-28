import json
from pathlib import Path

import pytest

from aidd_agent.context import (
    AccessDeniedError, activate_task, active_task, create_task,
    create_user, list_tasks, require_task_access,
)
from aidd_agent.project_context import activate_project, create_scientific_project
from aidd_agent.registry import connect, initialize
from aidd_agent.cli import main


def test_user_tasks_are_isolated_and_workspaces_are_separate(tmp_path: Path) -> None:
    database = tmp_path / "aidd.sqlite3"
    initialize(database)
    with connect(database) as connection:
        alice = create_user(connection, "alice", "Alice")
        bob = create_user(connection, "bob", "Bob")
        alice_project = create_scientific_project(
            connection, alice, "Kinase Program", "Discover kinase inhibitors", tmp_path / "storage")
        bob_project = create_scientific_project(
            connection, bob, "Oncology Program", "Discover oncology compounds", tmp_path / "storage")
        activate_project(connection, alice, alice_project)
        activate_project(connection, bob, bob_project)
        task_a = create_task(connection, alice, alice_project, "WEE1 ATP Site",
                             "Screen macrocycles against WEE1")
        task_b = create_task(connection, bob, bob_project, "KRAS Switch II",
                             "Screen macrocycles against KRAS G12D")
        assert [row["task_id"] for row in list_tasks(connection, alice)] == [task_a]
        assert [row["task_id"] for row in list_tasks(connection, bob)] == [task_b]
        with pytest.raises(AccessDeniedError):
            activate_task(connection, alice, task_b)
        activated = activate_task(connection, alice, task_a)
        assert active_task(connection, alice)["task_id"] == task_a
    workspace_a = Path(activated["workspace_path"])
    assert workspace_a.name.startswith(f"{task_a.lower()}-")
    assert workspace_a.parent.name == "tasks"
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


def test_cli_creates_active_task_scoped_campaign(tmp_path: Path, capsys: pytest.CaptureFixture) -> None:
    database = tmp_path / "aidd.sqlite3"
    storage = tmp_path / "storage"

    assert main(["create-user", "--db", str(database), "--username", "alice"]) == 0
    user_id = json.loads(capsys.readouterr().out)["user_id"]
    assert main([
        "create-project", "--db", str(database), "--user", user_id,
        "--name", "Kinase Program", "--objective", "Discover inhibitors",
        "--storage-root", str(storage),
    ]) == 0
    project_id = json.loads(capsys.readouterr().out)["project_id"]
    assert main([
        "activate-project", "--db", str(database), "--user", user_id,
        "--project", project_id,
    ]) == 0
    capsys.readouterr()
    assert main([
        "create-task", "--db", str(database), "--user", user_id,
        "--project", project_id, "--name", "Structure Selection",
        "--objective", "Select experimental structures",
    ]) == 0
    task_id = json.loads(capsys.readouterr().out)["task_id"]
    assert main([
        "activate-task", "--db", str(database), "--user", user_id, "--task", task_id,
    ]) == 0
    capsys.readouterr()
    assert main([
        "create-campaign", "--db", str(database), "--user", user_id,
        "--project", project_id, "--task", task_id, "--name", "pdb-review",
        "--objective", "Review PDB candidates",
    ]) == 0
    campaign_id = json.loads(capsys.readouterr().out)["campaign_id"]

    with connect(database) as connection:
        campaign = connection.execute(
            "SELECT project_id, task_id FROM campaign WHERE id=?", (campaign_id,)
        ).fetchone()
    assert dict(campaign) == {"project_id": project_id, "task_id": task_id}
