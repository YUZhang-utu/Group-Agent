from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import platform
import re
import shutil
import sqlite3
from typing import Any

from .project_context import require_active_project, require_project_owner
from .registry import stable_id, utc_now

SAFE_ID = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]{0,127}$")
SAFE_SLURM = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._:-]{0,127}$")


def _safe(value: str, label: str, pattern: re.Pattern[str] = SAFE_ID) -> str:
    if not pattern.fullmatch(value):
        raise ValueError(f"Invalid {label}: {value!r}")
    return value


def _remote_root(value: str) -> str:
    path = PurePosixPath(value)
    if not path.is_absolute() or ".." in path.parts:
        raise ValueError(f"Remote root must be an absolute normalized POSIX path: {value}")
    return str(path)


def remote_within(path: str, root: str) -> str:
    candidate, allowed = PurePosixPath(path), PurePosixPath(root)
    if not candidate.is_absolute() or ".." in candidate.parts:
        raise ValueError(f"Unsafe remote path: {path}")
    if candidate != allowed and allowed not in candidate.parents:
        raise ValueError(f"Remote path escapes registered root: {path}")
    return str(candidate)


def register_cluster(connection: sqlite3.Connection, *, name: str, host: str,
                     port: int, transport: str, project_root: str,
                     inbox_root: str, run_root: str, library_root: str,
                     software_root: str, container_root: str,
                     scratch_template: str, ssh_profile: str | None = None,
                     expected_host_key: str | None = None) -> str:
    if transport not in {"auto", "openssh", "putty"}:
        raise ValueError("Unsupported transport")
    if not 1 <= port <= 65535:
        raise ValueError("Invalid SSH port")
    cluster_id = stable_id("CLU")
    roots = [_remote_root(x) for x in (project_root, inbox_root, run_root,
                                       library_root, software_root, container_root)]
    if not scratch_template.startswith("/") or ".." in PurePosixPath(scratch_template).parts:
        raise ValueError("Invalid scratch_template")
    connection.execute(
        """INSERT INTO cluster_profile VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)""",
        (cluster_id, name.strip(), host.strip(), port, transport, ssh_profile,
         expected_host_key, *roots, scratch_template, utc_now()),
    )
    return cluster_id


def set_user_cluster_identity(connection: sqlite3.Connection, user_id: str,
                              cluster_id: str, cluster_username: str,
                              transport_preference: str = "auto",
                              putty_saved_session: str | None = None) -> None:
    _safe(cluster_username, "cluster username")
    if transport_preference not in {"auto", "openssh", "putty"}:
        raise ValueError("Unsupported transport preference")
    if not connection.execute("SELECT 1 FROM app_user WHERE id=?", (user_id,)).fetchone():
        raise KeyError(f"Unknown user: {user_id}")
    if not connection.execute("SELECT 1 FROM cluster_profile WHERE id=?", (cluster_id,)).fetchone():
        raise KeyError(f"Unknown cluster: {cluster_id}")
    connection.execute(
        """INSERT OR REPLACE INTO user_cluster_identity VALUES (?, ?, ?, ?, ?, ?)""",
        (user_id, cluster_id, cluster_username, transport_preference,
         putty_saved_session, utc_now()),
    )


def resolve_transport(profile: sqlite3.Row, identity: sqlite3.Row,
                      system: str | None = None) -> str:
    requested = identity["transport_preference"]
    if requested == "auto":
        requested = profile["transport"]
    if requested != "auto":
        return requested
    return "openssh" if (system or platform.system()).lower() in {
        "darwin", "linux", "windows"
    } else "openssh"


def transport_commands(profile: sqlite3.Row, identity: sqlite3.Row,
                       local_bundle: Path, remote_bundle: str,
                       *, transport: str | None = None) -> dict[str, list[str]]:
    selected = transport or resolve_transport(profile, identity)
    remote_bundle = remote_within(remote_bundle, profile["inbox_root"])
    username = _safe(identity["cluster_username"], "cluster username")
    if selected == "openssh":
        target = profile["ssh_profile"] or f"{username}@{profile['host']}"
        common = [] if profile["ssh_profile"] else ["-p", str(profile["port"])]
        upload = ["scp", *common, str(local_bundle), f"{target}:{remote_bundle}"]
        execute = ["ssh", *common, target,
                   "aidd-cluster-runner", "submit", remote_bundle]
    elif selected == "putty":
        session = identity["putty_saved_session"] or profile["ssh_profile"]
        if not session:
            raise ValueError("PuTTY transport requires a saved session")
        upload = ["pscp", "-load", session, "-share", str(local_bundle),
                  f"{username}@{profile['host']}:{remote_bundle}"]
        execute = ["plink", "-load", session, "-share",
                   "aidd-cluster-runner", "submit", remote_bundle]
    else:
        raise ValueError(f"Unsupported transport: {selected}")
    return {"transport": [selected], "upload": upload, "submit": execute}


def _file_hash(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def _slurm_script(run_id: str, spec: dict[str, Any]) -> str:
    _safe(run_id, "Run ID")
    partition = _safe(str(spec["partition"]), "partition", SAFE_SLURM)
    time_limit = str(spec["time"])
    if not re.fullmatch(r"(?:\d+-)?\d{1,2}:\d{2}:\d{2}", time_limit):
        raise ValueError("Invalid Slurm time")
    cpus, memory = int(spec["cpus_per_task"]), int(spec["memory_gb"])
    if cpus < 1 or memory < 1:
        raise ValueError("Slurm resources must be positive")
    lines = ["#!/bin/bash", f"#SBATCH --job-name={run_id}",
             f"#SBATCH --partition={partition}", f"#SBATCH --cpus-per-task={cpus}",
             f"#SBATCH --mem={memory}G", f"#SBATCH --time={time_limit}",
             "#SBATCH --output=logs/%A_%a.out", "#SBATCH --error=logs/%A_%a.err"]
    array = spec.get("array")
    if array:
        start, end, maximum = int(array["start"]), int(array["end"]), int(array["max_concurrent"])
        if start < 0 or end < start or maximum < 1:
            raise ValueError("Invalid Slurm array")
        lines.append(f"#SBATCH --array={start}-{end}%{maximum}")
    lines += ["", "set -euo pipefail", "mkdir -p logs outputs",
              "aidd-cluster-runner execute run.json", ""]
    return "\n".join(lines)


def export_slurm_bundle(connection: sqlite3.Connection, user_id: str, run_id: str,
                        cluster_id: str, output_dir: Path,
                        slurm_spec: dict[str, Any]) -> dict[str, Any]:
    run = connection.execute("SELECT * FROM workflow_run WHERE id=?", (run_id,)).fetchone()
    if not run:
        raise KeyError(f"Unknown Run: {run_id}")
    project = require_project_owner(connection, user_id, run["project_id"])
    require_active_project(connection, user_id, run["project_id"])
    profile = connection.execute("SELECT * FROM cluster_profile WHERE id=?", (cluster_id,)).fetchone()
    identity = connection.execute(
        "SELECT * FROM user_cluster_identity WHERE user_id=? AND cluster_id=?",
        (user_id, cluster_id),
    ).fetchone()
    if not profile or not identity:
        raise ValueError("Cluster Profile and user identity are required")
    bundle = output_dir.resolve() / f"{run_id}.bundle"
    bundle.mkdir(parents=True, exist_ok=False)
    run_payload = {"schema_version": 1, "run_id": run_id,
                   "project_id": run["project_id"], "campaign_id": run["campaign_id"],
                   "run_type": run["run_type"], "tool": {"name": run["tool_name"],
                   "version": run["tool_version"]}, "parameters": json.loads(run["parameters_json"]),
                   "cluster_id": cluster_id}
    (bundle / "run.json").write_text(json.dumps(run_payload, indent=2, sort_keys=True), encoding="utf-8")
    (bundle / "submit.slurm").write_text(_slurm_script(run_id, slurm_spec), encoding="utf-8", newline="\n")
    files = []
    for path in sorted(bundle.iterdir(), key=lambda p: p.name):
        files.append({"path": path.name, "sha256": _file_hash(path), "size_bytes": path.stat().st_size})
    manifest = {"schema_version": 1, "run_id": run_id, "project_id": run["project_id"],
                "cluster_id": cluster_id, "files": files}
    (bundle / "bundle_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True), encoding="utf-8"
    )
    bundle_hash = _file_hash(bundle / "bundle_manifest.json")
    remote = str(PurePosixPath(profile["inbox_root"]) / f"{run_id}.bundle")
    commands = transport_commands(profile, identity, bundle, remote)
    connection.execute(
        """INSERT INTO cluster_bundle(id, run_id, cluster_id, bundle_path,
           bundle_sha256, remote_bundle_path, created_at) VALUES (?, ?, ?, ?, ?, ?, ?)""",
        (stable_id("BND"), run_id, cluster_id, str(bundle), bundle_hash, remote, utc_now()),
    )
    return {"bundle_path": str(bundle), "bundle_sha256": bundle_hash,
            "remote_bundle_path": remote, "commands": commands,
            "authentication_required": True, "project_root": project["workspace_path"]}


def validate_bundle(bundle: Path) -> dict[str, Any]:
    manifest_path = bundle / "bundle_manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    for item in manifest["files"]:
        path = bundle / item["path"]
        if not path.is_file() or _file_hash(path) != item["sha256"]:
            raise ValueError(f"Bundle file verification failed: {item['path']}")
    return manifest


def import_submission_receipt(connection: sqlite3.Connection, user_id: str,
                              receipt_path: Path) -> dict[str, Any]:
    receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
    required = {"schema_version", "run_id", "project_id", "cluster_id", "slurm_job_id",
                "submitted_by", "submitted_at", "bundle_sha256"}
    if receipt.get("schema_version") != 1 or not required.issubset(receipt):
        raise ValueError("Invalid submission receipt")
    run = connection.execute("SELECT * FROM workflow_run WHERE id=?", (receipt["run_id"],)).fetchone()
    if not run or run["project_id"] != receipt["project_id"]:
        raise ValueError("Receipt does not match its Run and Project")
    require_project_owner(connection, user_id, run["project_id"])
    require_active_project(connection, user_id, run["project_id"])
    _safe(str(receipt["slurm_job_id"]), "Slurm job ID", SAFE_SLURM)
    bundle = connection.execute(
        """SELECT 1 FROM cluster_bundle WHERE run_id=? AND cluster_id=?
           AND bundle_sha256=?""",
        (run["id"], receipt["cluster_id"], receipt["bundle_sha256"]),
    ).fetchone()
    if not bundle:
        raise ValueError("Receipt does not match an exported Cluster Bundle")
    receipt_id = stable_id("SUB")
    connection.execute(
        """INSERT INTO submission_receipt VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)""",
        (receipt_id, run["id"], receipt["cluster_id"], str(receipt["slurm_job_id"]),
         receipt["submitted_by"], receipt["submitted_at"], receipt["bundle_sha256"],
         json.dumps(receipt, sort_keys=True), utc_now()),
    )
    connection.execute(
        """UPDATE workflow_run SET status='submitted', external_job_id=? WHERE id=?""",
        (str(receipt["slurm_job_id"]), run["id"]),
    )
    return {"receipt_id": receipt_id, "run_id": run["id"], "status": "submitted",
            "slurm_job_id": str(receipt["slurm_job_id"])}
