from __future__ import annotations

from contextlib import contextmanager
from datetime import datetime, timezone
import json
from pathlib import Path
import sqlite3
from typing import Iterator
from uuid import uuid4

from .mol2 import Mol2Record


SCHEMA = """
PRAGMA foreign_keys = ON;
CREATE TABLE IF NOT EXISTS library (
    id TEXT PRIMARY KEY,
    name TEXT NOT NULL UNIQUE,
    created_at TEXT NOT NULL
);
CREATE TABLE IF NOT EXISTS molecule (
    id TEXT PRIMARY KEY,
    library_id TEXT NOT NULL REFERENCES library(id),
    source_name TEXT NOT NULL,
    created_at TEXT NOT NULL,
    UNIQUE(library_id, source_name)
);
CREATE TABLE IF NOT EXISTS conformer (
    id TEXT PRIMARY KEY,
    molecule_id TEXT NOT NULL REFERENCES molecule(id),
    conformer_index INTEGER NOT NULL,
    source_record_name TEXT NOT NULL,
    atom_count INTEGER NOT NULL,
    bond_count INTEGER NOT NULL,
    content_sha256 TEXT NOT NULL,
    topology_sha256 TEXT NOT NULL,
    source_path TEXT NOT NULL,
    source_record_index INTEGER NOT NULL,
    warnings_json TEXT NOT NULL,
    created_at TEXT NOT NULL,
    UNIQUE(molecule_id, conformer_index),
    UNIQUE(molecule_id, content_sha256)
);
CREATE INDEX IF NOT EXISTS idx_conformer_content_hash ON conformer(content_sha256);
CREATE INDEX IF NOT EXISTS idx_conformer_topology_hash ON conformer(topology_sha256);
CREATE TABLE IF NOT EXISTS app_user (
    id TEXT PRIMARY KEY,
    username TEXT NOT NULL UNIQUE,
    display_name TEXT NOT NULL,
    created_at TEXT NOT NULL
);
CREATE TABLE IF NOT EXISTS project (
    id TEXT PRIMARY KEY,
    owner_user_id TEXT NOT NULL REFERENCES app_user(id),
    name TEXT NOT NULL,
    slug TEXT NOT NULL,
    description TEXT NOT NULL,
    objective TEXT,
    status TEXT NOT NULL DEFAULT 'created',
    workspace_path TEXT,
    updated_at TEXT,
    created_at TEXT NOT NULL,
    UNIQUE(owner_user_id, slug)
);
CREATE TABLE IF NOT EXISTS task (
    id TEXT PRIMARY KEY,
    project_id TEXT NOT NULL REFERENCES project(id),
    name TEXT NOT NULL,
    slug TEXT NOT NULL,
    objective TEXT NOT NULL,
    status TEXT NOT NULL CHECK(status IN ('created', 'active', 'paused', 'completed', 'archived')),
    workspace_path TEXT NOT NULL UNIQUE,
    created_at TEXT NOT NULL,
    updated_at TEXT NOT NULL,
    UNIQUE(project_id, slug)
);
CREATE TABLE IF NOT EXISTS task_access (
    task_id TEXT NOT NULL REFERENCES task(id),
    user_id TEXT NOT NULL REFERENCES app_user(id),
    role TEXT NOT NULL CHECK(role IN ('owner', 'collaborator', 'viewer')),
    created_at TEXT NOT NULL,
    PRIMARY KEY(task_id, user_id)
);
CREATE TABLE IF NOT EXISTS active_task_context (
    user_id TEXT PRIMARY KEY REFERENCES app_user(id),
    task_id TEXT NOT NULL REFERENCES task(id),
    activated_at TEXT NOT NULL
);
CREATE TABLE IF NOT EXISTS active_project_context (
    user_id TEXT PRIMARY KEY REFERENCES app_user(id),
    project_id TEXT NOT NULL REFERENCES project(id),
    activated_at TEXT NOT NULL
);
CREATE TABLE IF NOT EXISTS campaign (
    id TEXT PRIMARY KEY,
    task_id TEXT REFERENCES task(id),
    project_id TEXT REFERENCES project(id),
    name TEXT NOT NULL UNIQUE,
    objective TEXT NOT NULL,
    state TEXT NOT NULL CHECK(state IN ('draft', 'target_set', 'structures_review', 'structure_selected')),
    selected_structure_id TEXT,
    created_at TEXT NOT NULL,
    updated_at TEXT NOT NULL
);
CREATE TABLE IF NOT EXISTS target (
    id TEXT PRIMARY KEY,
    campaign_id TEXT NOT NULL UNIQUE REFERENCES campaign(id),
    name TEXT NOT NULL,
    organism TEXT NOT NULL,
    uniprot_id TEXT,
    gene_name TEXT,
    notes TEXT NOT NULL,
    created_at TEXT NOT NULL
);
CREATE TABLE IF NOT EXISTS structure_candidate (
    id TEXT PRIMARY KEY,
    campaign_id TEXT NOT NULL REFERENCES campaign(id),
    pdb_id TEXT NOT NULL,
    title TEXT NOT NULL,
    experimental_method TEXT,
    resolution_angstrom REAL,
    deposition_date TEXT,
    ligand_ids_json TEXT NOT NULL,
    metadata_json TEXT NOT NULL,
    query_json TEXT NOT NULL,
    created_at TEXT NOT NULL,
    UNIQUE(campaign_id, pdb_id)
);
CREATE TABLE IF NOT EXISTS decision_event (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    campaign_id TEXT NOT NULL REFERENCES campaign(id),
    event_type TEXT NOT NULL,
    from_state TEXT,
    to_state TEXT,
    rationale TEXT NOT NULL,
    payload_json TEXT NOT NULL,
    created_at TEXT NOT NULL
);
CREATE INDEX IF NOT EXISTS idx_structure_candidate_campaign ON structure_candidate(campaign_id);
CREATE INDEX IF NOT EXISTS idx_decision_event_campaign ON decision_event(campaign_id, id);
CREATE TABLE IF NOT EXISTS predicted_structure_candidate (
    id TEXT PRIMARY KEY,
    campaign_id TEXT NOT NULL REFERENCES campaign(id),
    backend TEXT NOT NULL,
    model_name TEXT NOT NULL,
    model_version TEXT NOT NULL,
    construct_name TEXT NOT NULL,
    chain_ids_json TEXT NOT NULL,
    structure_path TEXT NOT NULL,
    structure_sha256 TEXT NOT NULL,
    ranking_score REAL,
    confidence_json TEXT NOT NULL,
    provenance_json TEXT NOT NULL,
    created_at TEXT NOT NULL,
    UNIQUE(campaign_id, structure_sha256)
);
CREATE INDEX IF NOT EXISTS idx_predicted_structure_campaign
    ON predicted_structure_candidate(campaign_id, created_at);
CREATE TABLE IF NOT EXISTS receptor_selection_lock (
    campaign_id TEXT PRIMARY KEY REFERENCES campaign(id),
    candidate_kind TEXT NOT NULL CHECK(candidate_kind IN ('experimental', 'predicted')),
    candidate_id TEXT NOT NULL,
    snapshot_json TEXT NOT NULL,
    rationale TEXT NOT NULL,
    selected_by TEXT NOT NULL REFERENCES app_user(id),
    created_at TEXT NOT NULL
);
CREATE TABLE IF NOT EXISTS receptor_ensemble (
    id TEXT PRIMARY KEY,
    campaign_id TEXT NOT NULL REFERENCES campaign(id),
    name TEXT NOT NULL,
    reference_candidate_id TEXT NOT NULL,
    pocket_residues_json TEXT NOT NULL,
    rationale TEXT NOT NULL,
    created_at TEXT NOT NULL,
    UNIQUE(campaign_id, name)
);
CREATE TABLE IF NOT EXISTS receptor_ensemble_member (
    ensemble_id TEXT NOT NULL REFERENCES receptor_ensemble(id),
    candidate_kind TEXT NOT NULL CHECK(candidate_kind IN ('experimental', 'predicted')),
    candidate_id TEXT NOT NULL,
    chain_id TEXT NOT NULL,
    residue_offset INTEGER NOT NULL DEFAULT 0,
    PRIMARY KEY(ensemble_id, candidate_id)
);
CREATE TABLE IF NOT EXISTS receptor_ensemble_exclusion (
    ensemble_id TEXT NOT NULL REFERENCES receptor_ensemble(id),
    candidate_id TEXT NOT NULL,
    candidate_kind TEXT NOT NULL CHECK(candidate_kind IN ('experimental', 'predicted')),
    reason TEXT NOT NULL,
    PRIMARY KEY(ensemble_id, candidate_id)
);
CREATE INDEX IF NOT EXISTS idx_task_project ON task(project_id);
CREATE INDEX IF NOT EXISTS idx_task_access_user ON task_access(user_id, task_id);
CREATE INDEX IF NOT EXISTS idx_project_owner ON project(owner_user_id, id);
CREATE TABLE IF NOT EXISTS workflow_run (
    id TEXT PRIMARY KEY,
    project_id TEXT NOT NULL REFERENCES project(id),
    campaign_id TEXT REFERENCES campaign(id),
    owner_user_id TEXT NOT NULL REFERENCES app_user(id),
    run_type TEXT NOT NULL,
    status TEXT NOT NULL CHECK(status IN ('created', 'submitted', 'running', 'completed', 'failed', 'cancelled')),
    tool_name TEXT NOT NULL,
    tool_version TEXT NOT NULL,
    execution_backend TEXT NOT NULL,
    parameters_json TEXT NOT NULL,
    compute_key TEXT NOT NULL,
    run_path TEXT NOT NULL UNIQUE,
    external_job_id TEXT,
    error_message TEXT,
    created_at TEXT NOT NULL,
    started_at TEXT,
    completed_at TEXT,
    UNIQUE(project_id, compute_key, id)
);
CREATE TABLE IF NOT EXISTS artifact (
    id TEXT PRIMARY KEY,
    project_id TEXT REFERENCES project(id),
    run_id TEXT REFERENCES workflow_run(id),
    scope TEXT NOT NULL CHECK(scope IN ('project', 'shared')),
    artifact_type TEXT NOT NULL,
    role TEXT NOT NULL,
    format TEXT NOT NULL,
    storage_path TEXT NOT NULL,
    sha256 TEXT NOT NULL,
    size_bytes INTEGER NOT NULL,
    metadata_json TEXT NOT NULL,
    created_at TEXT NOT NULL,
    verified_at TEXT NOT NULL,
    UNIQUE(scope, project_id, sha256, role)
);
CREATE TABLE IF NOT EXISTS run_input (
    run_id TEXT NOT NULL REFERENCES workflow_run(id),
    artifact_id TEXT NOT NULL REFERENCES artifact(id),
    role TEXT NOT NULL,
    PRIMARY KEY(run_id, artifact_id, role)
);
CREATE TABLE IF NOT EXISTS run_output (
    run_id TEXT NOT NULL REFERENCES workflow_run(id),
    artifact_id TEXT NOT NULL REFERENCES artifact(id),
    role TEXT NOT NULL,
    PRIMARY KEY(run_id, artifact_id, role)
);
CREATE TABLE IF NOT EXISTS run_metric (
    run_id TEXT NOT NULL REFERENCES workflow_run(id),
    name TEXT NOT NULL,
    numeric_value REAL,
    text_value TEXT,
    unit TEXT,
    PRIMARY KEY(run_id, name)
);
CREATE TABLE IF NOT EXISTS project_decision (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    project_id TEXT NOT NULL REFERENCES project(id),
    campaign_id TEXT REFERENCES campaign(id),
    user_id TEXT NOT NULL REFERENCES app_user(id),
    decision_type TEXT NOT NULL,
    subject TEXT NOT NULL,
    rationale TEXT NOT NULL,
    payload_json TEXT NOT NULL,
    created_at TEXT NOT NULL
);
CREATE INDEX IF NOT EXISTS idx_run_project ON workflow_run(project_id, created_at);
CREATE INDEX IF NOT EXISTS idx_artifact_project ON artifact(project_id, created_at);
CREATE INDEX IF NOT EXISTS idx_artifact_hash ON artifact(sha256);
CREATE INDEX IF NOT EXISTS idx_project_decision ON project_decision(project_id, id);
CREATE TABLE IF NOT EXISTS cluster_profile (
    id TEXT PRIMARY KEY,
    name TEXT NOT NULL UNIQUE,
    host TEXT NOT NULL,
    port INTEGER NOT NULL,
    transport TEXT NOT NULL CHECK(transport IN ('auto', 'openssh', 'putty')),
    ssh_profile TEXT,
    expected_host_key TEXT,
    project_root TEXT NOT NULL,
    inbox_root TEXT NOT NULL,
    run_root TEXT NOT NULL,
    library_root TEXT NOT NULL,
    software_root TEXT NOT NULL,
    container_root TEXT NOT NULL,
    scratch_template TEXT NOT NULL,
    created_at TEXT NOT NULL
);
CREATE TABLE IF NOT EXISTS user_cluster_identity (
    user_id TEXT NOT NULL REFERENCES app_user(id),
    cluster_id TEXT NOT NULL REFERENCES cluster_profile(id),
    cluster_username TEXT NOT NULL,
    transport_preference TEXT NOT NULL CHECK(transport_preference IN ('auto', 'openssh', 'putty')),
    putty_saved_session TEXT,
    created_at TEXT NOT NULL,
    PRIMARY KEY(user_id, cluster_id)
);
CREATE TABLE IF NOT EXISTS submission_receipt (
    id TEXT PRIMARY KEY,
    run_id TEXT NOT NULL UNIQUE REFERENCES workflow_run(id),
    cluster_id TEXT NOT NULL REFERENCES cluster_profile(id),
    slurm_job_id TEXT NOT NULL,
    submitted_by TEXT NOT NULL,
    submitted_at TEXT NOT NULL,
    bundle_sha256 TEXT NOT NULL,
    receipt_json TEXT NOT NULL,
    imported_at TEXT NOT NULL
);
CREATE TABLE IF NOT EXISTS cluster_bundle (
    id TEXT PRIMARY KEY,
    run_id TEXT NOT NULL REFERENCES workflow_run(id),
    cluster_id TEXT NOT NULL REFERENCES cluster_profile(id),
    bundle_path TEXT NOT NULL,
    bundle_sha256 TEXT NOT NULL,
    remote_bundle_path TEXT NOT NULL,
    created_at TEXT NOT NULL,
    UNIQUE(run_id, cluster_id, bundle_sha256)
);
CREATE INDEX IF NOT EXISTS idx_cluster_identity_user ON user_cluster_identity(user_id, cluster_id);
CREATE INDEX IF NOT EXISTS idx_cluster_bundle_run ON cluster_bundle(run_id, cluster_id);
CREATE TABLE IF NOT EXISTS structure_comparison_set (
    id TEXT PRIMARY KEY,
    campaign_id TEXT NOT NULL REFERENCES campaign(id),
    name TEXT NOT NULL,
    reference_candidate_id TEXT NOT NULL REFERENCES structure_candidate(id),
    pocket_residues_json TEXT NOT NULL,
    rationale TEXT NOT NULL,
    created_at TEXT NOT NULL,
    UNIQUE(campaign_id, name)
);
CREATE TABLE IF NOT EXISTS structure_comparison_member (
    comparison_set_id TEXT NOT NULL REFERENCES structure_comparison_set(id),
    candidate_id TEXT NOT NULL REFERENCES structure_candidate(id),
    chain_id TEXT NOT NULL,
    ligand_ids_json TEXT NOT NULL,
    PRIMARY KEY(comparison_set_id, candidate_id)
);
CREATE TABLE IF NOT EXISTS similarity_index (
    id TEXT PRIMARY KEY,
    library_id TEXT NOT NULL REFERENCES library(id),
    index_type TEXT NOT NULL CHECK(index_type IN ('morgan2d', 'usrcat3d')),
    parameters_json TEXT NOT NULL,
    artifact_id TEXT REFERENCES artifact(id),
    molecule_count INTEGER NOT NULL,
    conformer_count INTEGER NOT NULL,
    status TEXT NOT NULL CHECK(status IN ('planned', 'building', 'ready', 'failed')),
    created_at TEXT NOT NULL,
    UNIQUE(library_id, index_type, parameters_json)
);
CREATE TABLE IF NOT EXISTS ai_review_request (
    id TEXT PRIMARY KEY,
    project_id TEXT NOT NULL REFERENCES project(id),
    campaign_id TEXT REFERENCES campaign(id),
    user_id TEXT NOT NULL REFERENCES app_user(id),
    task_type TEXT NOT NULL,
    subject_type TEXT NOT NULL,
    subject_id TEXT NOT NULL,
    prompt_version TEXT NOT NULL,
    prompt_text TEXT NOT NULL,
    evidence_json TEXT NOT NULL,
    data_classes_json TEXT NOT NULL,
    privacy_policy_json TEXT NOT NULL,
    status TEXT NOT NULL CHECK(status IN ('pending', 'completed')),
    created_at TEXT NOT NULL
);
CREATE TABLE IF NOT EXISTS ai_recommendation (
    id TEXT PRIMARY KEY,
    request_id TEXT NOT NULL UNIQUE REFERENCES ai_review_request(id),
    provider TEXT NOT NULL,
    model_name TEXT NOT NULL,
    model_version TEXT,
    action TEXT NOT NULL,
    rationale TEXT NOT NULL,
    confidence REAL NOT NULL CHECK(confidence >= 0 AND confidence <= 1),
    evidence_refs_json TEXT NOT NULL,
    uncertainties_json TEXT NOT NULL,
    requires_human_review INTEGER NOT NULL CHECK(requires_human_review IN (0, 1)),
    response_json TEXT NOT NULL,
    created_at TEXT NOT NULL
);
CREATE INDEX IF NOT EXISTS idx_ai_review_project ON ai_review_request(project_id, created_at);
"""


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def stable_id(prefix: str) -> str:
    return f"{prefix}-{uuid4().hex[:12].upper()}"


@contextmanager
def connect(db_path: Path) -> Iterator[sqlite3.Connection]:
    db_path.parent.mkdir(parents=True, exist_ok=True)
    connection = sqlite3.connect(db_path)
    connection.row_factory = sqlite3.Row
    connection.execute("PRAGMA foreign_keys = ON")
    try:
        yield connection
        connection.commit()
    except Exception:
        connection.rollback()
        raise
    finally:
        connection.close()


def initialize(db_path: Path) -> None:
    with connect(db_path) as connection:
        connection.executescript(SCHEMA)
        columns = {row["name"] for row in connection.execute("PRAGMA table_info(campaign)")}
        if "task_id" not in columns:
            connection.execute("ALTER TABLE campaign ADD COLUMN task_id TEXT REFERENCES task(id)")
        if "project_id" not in columns:
            connection.execute("ALTER TABLE campaign ADD COLUMN project_id TEXT REFERENCES project(id)")
        project_columns = {row["name"] for row in connection.execute("PRAGMA table_info(project)")}
        for name, definition in (
            ("status", "TEXT NOT NULL DEFAULT 'created'"),
            ("workspace_path", "TEXT"),
            ("updated_at", "TEXT"),
            ("objective", "TEXT"),
        ):
            if name not in project_columns:
                connection.execute(f"ALTER TABLE project ADD COLUMN {name} {definition}")


def ensure_library(connection: sqlite3.Connection, name: str) -> str:
    row = connection.execute("SELECT id FROM library WHERE name = ?", (name,)).fetchone()
    if row:
        return str(row["id"])
    library_id = stable_id("LIB")
    connection.execute(
        "INSERT INTO library(id, name, created_at) VALUES (?, ?, ?)",
        (library_id, name, utc_now()),
    )
    return library_id


def register_record(connection: sqlite3.Connection, library_id: str, record: Mol2Record) -> str:
    row = connection.execute(
        "SELECT id FROM molecule WHERE library_id = ? AND source_name = ?",
        (library_id, record.molecule_name),
    ).fetchone()
    if row:
        molecule_id = str(row["id"])
    else:
        molecule_id = stable_id("MOL")
        connection.execute(
            "INSERT INTO molecule(id, library_id, source_name, created_at) VALUES (?, ?, ?, ?)",
            (molecule_id, library_id, record.molecule_name, utc_now()),
        )

    duplicate = connection.execute(
        "SELECT id FROM conformer WHERE molecule_id = ? AND content_sha256 = ?",
        (molecule_id, record.content_sha256),
    ).fetchone()
    if duplicate:
        return "duplicate"

    occupied = connection.execute(
        "SELECT id FROM conformer WHERE molecule_id = ? AND conformer_index = ?",
        (molecule_id, record.conformer_index),
    ).fetchone()
    if occupied:
        raise ValueError(
            f"Conformer index conflict for {record.molecule_name} conf{record.conformer_index}; "
            "rename the source records or provide a corrected library"
        )

    connection.execute(
        """
        INSERT INTO conformer(
            id, molecule_id, conformer_index, source_record_name, atom_count,
            bond_count, content_sha256, topology_sha256, source_path,
            source_record_index, warnings_json, created_at
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        """,
        (
            stable_id("CNF"), molecule_id, record.conformer_index, record.name,
            record.atom_count, record.bond_count, record.content_sha256,
            record.topology_sha256, str(record.source_path), record.record_index,
            json.dumps(record.warnings), utc_now(),
        ),
    )
    return "inserted"


def library_summary(connection: sqlite3.Connection, library_name: str) -> dict[str, object]:
    library = connection.execute(
        "SELECT id, name, created_at FROM library WHERE name = ?", (library_name,)
    ).fetchone()
    if not library:
        raise KeyError(f"Unknown library: {library_name}")
    counts = connection.execute(
        """
        SELECT COUNT(DISTINCT m.id) AS molecules, COUNT(c.id) AS conformers,
               COALESCE(SUM(c.atom_count), 0) AS atoms
        FROM molecule m LEFT JOIN conformer c ON c.molecule_id = m.id
        WHERE m.library_id = ?
        """,
        (library["id"],),
    ).fetchone()
    distribution = connection.execute(
        """
        SELECT conformer_count, COUNT(*) AS molecule_count FROM (
            SELECT m.id, COUNT(c.id) AS conformer_count
            FROM molecule m LEFT JOIN conformer c ON c.molecule_id = m.id
            WHERE m.library_id = ? GROUP BY m.id
        ) GROUP BY conformer_count ORDER BY conformer_count
        """,
        (library["id"],),
    ).fetchall()
    return {
        "library_id": library["id"],
        "library_name": library["name"],
        "created_at": library["created_at"],
        "molecule_count": counts["molecules"],
        "conformer_count": counts["conformers"],
        "atom_count": counts["atoms"],
        "conformers_per_molecule": {
            str(row["conformer_count"]): row["molecule_count"] for row in distribution
        },
    }
