import json
from pathlib import Path

import pytest

from aidd_agent.cluster import (
    export_slurm_bundle, import_submission_receipt, register_cluster,
    remote_within, set_user_cluster_identity, transport_commands, validate_bundle,
)
from aidd_agent.context import create_user
from aidd_agent.project_context import activate_project, create_scientific_project
from aidd_agent.provenance import create_workflow_run
from aidd_agent.registry import connect, initialize


def setup_run(connection, tmp_path: Path):
    user_id = create_user(connection, "alice")
    project_id = create_scientific_project(
        connection, user_id, "Cluster Project", "Run docking", tmp_path / "storage"
    )
    activate_project(connection, user_id, project_id)
    run_id = create_workflow_run(
        connection, user_id, project_id, run_type="docking", tool_name="PLANTS",
        tool_version="1.2", parameters={"library_id": "LIB-AFA68EE6888C"},
        execution_backend="slurm",
    )
    cluster_id = register_cluster(
        connection, name="institutional-hpc", host="login.hpc.example.edu", port=22,
        transport="auto", ssh_profile="institutional-hpc",
        expected_host_key="SHA256:test", project_root="/project/group/aidd",
        inbox_root="/project/group/aidd/inbox", run_root="/project/group/aidd/runs",
        library_root="/project/group/aidd/libraries",
        software_root="/project/group/aidd/software",
        container_root="/project/group/aidd/containers",
        scratch_template="/scratch/{cluster_username}/aidd",
    )
    set_user_cluster_identity(connection, user_id, cluster_id, "alice123")
    return user_id, project_id, run_id, cluster_id


def test_openssh_and_putty_commands_are_argument_arrays(tmp_path: Path) -> None:
    database = tmp_path / "db.sqlite3"
    initialize(database)
    with connect(database) as connection:
        user_id, _, _, cluster_id = setup_run(connection, tmp_path)
        profile = connection.execute("SELECT * FROM cluster_profile WHERE id=?", (cluster_id,)).fetchone()
        identity = connection.execute("SELECT * FROM user_cluster_identity").fetchone()
        openssh = transport_commands(
            profile, identity, tmp_path / "bundle", "/project/group/aidd/inbox/RUN.bundle",
            transport="openssh",
        )
        connection.execute(
            "UPDATE user_cluster_identity SET putty_saved_session='MyHPC' WHERE user_id=?",
            (user_id,),
        )
        identity = connection.execute("SELECT * FROM user_cluster_identity").fetchone()
        putty = transport_commands(
            profile, identity, tmp_path / "bundle", "/project/group/aidd/inbox/RUN.bundle",
            transport="putty",
        )
    assert openssh["upload"][0] == "scp"
    assert openssh["submit"][0] == "ssh"
    assert putty["upload"][:4] == ["pscp", "-load", "MyHPC", "-share"]
    assert putty["submit"][:4] == ["plink", "-load", "MyHPC", "-share"]


def test_remote_path_escape_is_rejected() -> None:
    assert remote_within("/project/aidd/inbox/run", "/project/aidd/inbox").endswith("/run")
    with pytest.raises(ValueError):
        remote_within("/project/aidd/inbox/../secret", "/project/aidd/inbox")
    with pytest.raises(ValueError):
        remote_within("/other/run", "/project/aidd/inbox")


def test_bundle_is_deterministic_and_detects_tampering(tmp_path: Path) -> None:
    database = tmp_path / "db.sqlite3"
    initialize(database)
    with connect(database) as connection:
        user_id, project_id, run_id, cluster_id = setup_run(connection, tmp_path)
        result = export_slurm_bundle(
            connection, user_id, run_id, cluster_id, tmp_path / "exports",
            {"partition": "cpu", "cpus_per_task": 16, "memory_gb": 32,
             "time": "08:00:00", "array": {"start": 1, "end": 100, "max_concurrent": 20}},
        )
    bundle = Path(result["bundle_path"])
    manifest = validate_bundle(bundle)
    script = (bundle / "submit.slurm").read_text(encoding="utf-8")
    assert manifest["run_id"] == run_id
    assert "#SBATCH --array=1-100%20" in script
    assert "aidd-cluster-runner execute run.json" in script
    assert result["authentication_required"] is True
    (bundle / "run.json").write_text("tampered", encoding="utf-8")
    with pytest.raises(ValueError, match="verification failed"):
        validate_bundle(bundle)


def test_receipt_updates_only_matching_active_project(tmp_path: Path) -> None:
    database = tmp_path / "db.sqlite3"
    initialize(database)
    with connect(database) as connection:
        user_id, project_id, run_id, cluster_id = setup_run(connection, tmp_path)
        exported = export_slurm_bundle(
            connection, user_id, run_id, cluster_id, tmp_path / "exports",
            {"partition": "cpu", "cpus_per_task": 4, "memory_gb": 8,
             "time": "01:00:00"},
        )
        receipt = tmp_path / "receipt.json"
        receipt.write_text(json.dumps({
            "schema_version": 1, "run_id": run_id, "project_id": project_id,
            "cluster_id": cluster_id, "slurm_job_id": "1837462",
            "submitted_by": "alice123", "submitted_at": "2026-08-27T14:30:00Z",
            "bundle_sha256": exported["bundle_sha256"],
        }), encoding="utf-8")
        bad = json.loads(receipt.read_text(encoding="utf-8"))
        bad["bundle_sha256"] = "wrong"
        receipt.write_text(json.dumps(bad), encoding="utf-8")
        with pytest.raises(ValueError, match="exported Cluster Bundle"):
            import_submission_receipt(connection, user_id, receipt)
        bad["bundle_sha256"] = exported["bundle_sha256"]
        receipt.write_text(json.dumps(bad), encoding="utf-8")
        result = import_submission_receipt(connection, user_id, receipt)
        run = connection.execute("SELECT status, external_job_id FROM workflow_run WHERE id=?", (run_id,)).fetchone()
    assert result["status"] == "submitted"
    assert (run["status"], run["external_job_id"]) == ("submitted", "1837462")
