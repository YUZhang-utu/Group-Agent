from __future__ import annotations

import json
import sqlite3
from typing import Any, Iterable

from .registry import stable_id, utc_now
from .project_context import require_active_project, require_project_owner
from .context import require_active_task, require_task_access


class CampaignStateError(ValueError):
    pass


def _campaign(connection: sqlite3.Connection, campaign_id: str) -> sqlite3.Row:
    row = connection.execute("SELECT * FROM campaign WHERE id = ?", (campaign_id,)).fetchone()
    if not row:
        raise KeyError(f"Unknown campaign: {campaign_id}")
    return row


def _campaign_for_user(connection: sqlite3.Connection, campaign_id: str,
                       user_id: str, *, write: bool = True) -> sqlite3.Row:
    row = _campaign(connection, campaign_id)
    if not row["project_id"]:
        raise CampaignStateError("Legacy Campaign has no Project context and must be migrated")
    if not row["task_id"]:
        raise CampaignStateError("Legacy Campaign has no Task context and must be migrated")
    require_project_owner(connection, user_id, row["project_id"])
    require_active_project(connection, user_id, row["project_id"])
    require_task_access(connection, user_id, row["task_id"], write=write)
    require_active_task(connection, user_id, row["task_id"])
    return row


def _event(connection: sqlite3.Connection, campaign_id: str, event_type: str,
           from_state: str | None, to_state: str | None, rationale: str,
           payload: dict[str, Any]) -> None:
    connection.execute(
        """INSERT INTO decision_event(campaign_id, event_type, from_state, to_state,
           rationale, payload_json, created_at) VALUES (?, ?, ?, ?, ?, ?, ?)""",
        (campaign_id, event_type, from_state, to_state, rationale,
         json.dumps(payload, sort_keys=True), utc_now()),
    )


def create_campaign(connection: sqlite3.Connection, user_id: str, project_id: str,
                    task_id: str, name: str, objective: str) -> str:
    require_project_owner(connection, user_id, project_id)
    require_active_project(connection, user_id, project_id)
    task = require_task_access(connection, user_id, task_id, write=True)
    require_active_task(connection, user_id, task_id)
    if task["project_id"] != project_id:
        raise ValueError(f"Task {task_id} does not belong to Project {project_id}")
    if not name.strip() or not objective.strip():
        raise ValueError("Campaign name and objective are required")
    campaign_id = stable_id("CAM")
    now = utc_now()
    connection.execute(
        """INSERT INTO campaign(id, task_id, project_id, name, objective, state,
           selected_structure_id, created_at, updated_at)
           VALUES (?, ?, ?, ?, ?, 'draft', NULL, ?, ?)""",
        (campaign_id, task_id, project_id, name.strip(), objective.strip(), now, now),
    )
    _event(connection, campaign_id, "campaign_created", None, "draft", objective.strip(), {})
    return campaign_id


def set_target(connection: sqlite3.Connection, user_id: str, campaign_id: str, *, name: str,
               organism: str, uniprot_id: str | None = None,
               gene_name: str | None = None, notes: str = "") -> str:
    campaign = _campaign_for_user(connection, campaign_id, user_id)
    if campaign["state"] != "draft":
        raise CampaignStateError("Target can only be set while campaign is draft")
    if not name.strip() or not organism.strip():
        raise ValueError("Target name and organism are required")
    target_id = stable_id("TGT")
    connection.execute(
        """INSERT INTO target(id, campaign_id, name, organism, uniprot_id,
           gene_name, notes, created_at) VALUES (?, ?, ?, ?, ?, ?, ?, ?)""",
        (target_id, campaign_id, name.strip(), organism.strip(), uniprot_id,
         gene_name, notes, utc_now()),
    )
    connection.execute(
        "UPDATE campaign SET state='target_set', updated_at=? WHERE id=?",
        (utc_now(), campaign_id),
    )
    _event(connection, campaign_id, "target_set", "draft", "target_set", notes,
           {"target_id": target_id, "uniprot_id": uniprot_id, "gene_name": gene_name})
    return target_id


def register_structure_candidates(connection: sqlite3.Connection, user_id: str,
                                  campaign_id: str,
                                  candidates: Iterable[dict[str, Any]],
                                  query: dict[str, Any]) -> int:
    campaign = _campaign_for_user(connection, campaign_id, user_id)
    if campaign["state"] not in {"target_set", "structures_review"}:
        raise CampaignStateError("Structure candidates require a target_set campaign")
    inserted = 0
    for item in candidates:
        pdb_id = str(item["pdb_id"]).upper()
        cursor = connection.execute(
            """INSERT OR IGNORE INTO structure_candidate(
               id, campaign_id, pdb_id, title, experimental_method,
               resolution_angstrom, deposition_date, ligand_ids_json,
               metadata_json, query_json, created_at
               ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)""",
            (stable_id("PDB"), campaign_id, pdb_id, item.get("title", ""),
             item.get("experimental_method"), item.get("resolution_angstrom"),
             item.get("deposition_date"), json.dumps(item.get("ligand_ids", [])),
             json.dumps(item, sort_keys=True), json.dumps(query, sort_keys=True), utc_now()),
        )
        inserted += cursor.rowcount
    connection.execute(
        "UPDATE campaign SET state='structures_review', updated_at=? WHERE id=?",
        (utc_now(), campaign_id),
    )
    _event(connection, campaign_id, "structures_registered", campaign["state"],
           "structures_review", "RCSB candidates registered for review",
           {"inserted": inserted, "query": query})
    return inserted


def select_structure(connection: sqlite3.Connection, user_id: str, campaign_id: str, pdb_id: str,
                     rationale: str) -> None:
    campaign = _campaign_for_user(connection, campaign_id, user_id)
    if campaign["state"] != "structures_review":
        raise CampaignStateError("Structure selection requires structures_review state")
    if not rationale.strip():
        raise ValueError("A scientific selection rationale is required")
    candidate = connection.execute(
        "SELECT id FROM structure_candidate WHERE campaign_id=? AND pdb_id=?",
        (campaign_id, pdb_id.upper()),
    ).fetchone()
    if not candidate:
        raise KeyError(f"PDB {pdb_id.upper()} is not a candidate for {campaign_id}")
    connection.execute(
        """UPDATE campaign SET state='structure_selected', selected_structure_id=?,
           updated_at=? WHERE id=?""",
        (candidate["id"], utc_now(), campaign_id),
    )
    _event(connection, campaign_id, "structure_selected", "structures_review",
           "structure_selected", rationale.strip(), {"pdb_id": pdb_id.upper()})


def campaign_status(connection: sqlite3.Connection, user_id: str,
                    campaign_id: str) -> dict[str, Any]:
    campaign = _campaign_for_user(connection, campaign_id, user_id, write=False)
    target = connection.execute(
        "SELECT * FROM target WHERE campaign_id=?", (campaign_id,)
    ).fetchone()
    candidates = connection.execute(
        """SELECT pdb_id, title, experimental_method, resolution_angstrom,
           deposition_date, ligand_ids_json FROM structure_candidate
           WHERE campaign_id=? ORDER BY resolution_angstrom IS NULL, resolution_angstrom, pdb_id""",
        (campaign_id,),
    ).fetchall()
    selected = None
    if campaign["selected_structure_id"]:
        row = connection.execute(
            "SELECT pdb_id FROM structure_candidate WHERE id=?",
            (campaign["selected_structure_id"],),
        ).fetchone()
        selected = row["pdb_id"] if row else None
    return {
        "campaign_id": campaign_id, "project_id": campaign["project_id"],
        "task_id": campaign["task_id"], "name": campaign["name"],
        "objective": campaign["objective"], "state": campaign["state"],
        "target": dict(target) if target else None, "selected_pdb_id": selected,
        "candidates": [dict(row) | {"ligand_ids": json.loads(row["ligand_ids_json"])}
                       for row in candidates],
    }
