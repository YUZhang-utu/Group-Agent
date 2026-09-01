from __future__ import annotations

import json
import hashlib
from pathlib import Path
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
        existing = connection.execute(
            "SELECT id FROM structure_candidate WHERE campaign_id=? AND pdb_id=?",
            (campaign_id, pdb_id),
        ).fetchone()
        if existing:
            connection.execute(
                """UPDATE structure_candidate SET title=?, experimental_method=?,
                   resolution_angstrom=?, deposition_date=?, ligand_ids_json=?,
                   metadata_json=?, query_json=? WHERE id=?""",
                (item.get("title", ""), item.get("experimental_method"),
                 item.get("resolution_angstrom"), item.get("deposition_date"),
                 json.dumps(item.get("ligand_ids", [])), json.dumps(item, sort_keys=True),
                 json.dumps(query, sort_keys=True), existing["id"]),
            )
            continue
        connection.execute(
            """INSERT INTO structure_candidate(
               id, campaign_id, pdb_id, title, experimental_method,
               resolution_angstrom, deposition_date, ligand_ids_json,
               metadata_json, query_json, created_at
               ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)""",
            (stable_id("PDB"), campaign_id, pdb_id, item.get("title", ""),
             item.get("experimental_method"), item.get("resolution_angstrom"),
             item.get("deposition_date"), json.dumps(item.get("ligand_ids", [])),
             json.dumps(item, sort_keys=True), json.dumps(query, sort_keys=True), utc_now()),
        )
        inserted += 1
    connection.execute(
        "UPDATE campaign SET state='structures_review', updated_at=? WHERE id=?",
        (utc_now(), campaign_id),
    )
    _event(connection, campaign_id, "structures_registered", campaign["state"],
           "structures_review", "RCSB candidates registered for review",
           {"inserted": inserted, "query": query})
    return inserted


def register_predicted_structure(
    connection: sqlite3.Connection, user_id: str, campaign_id: str,
    manifest: dict[str, Any],
) -> str:
    campaign = _campaign_for_user(connection, campaign_id, user_id)
    if campaign["state"] not in {"target_set", "structures_review"}:
        raise CampaignStateError("Predicted structures require a target_set campaign")
    required = {
        "backend", "model_name", "model_version", "construct_name",
        "chain_ids", "structure_path", "structure_sha256", "provenance",
    }
    missing = sorted(required - manifest.keys())
    if missing:
        raise ValueError(f"Prediction manifest is missing: {', '.join(missing)}")
    structure_path = Path(str(manifest["structure_path"])).resolve()
    if not structure_path.is_file():
        raise FileNotFoundError(f"Predicted structure does not exist: {structure_path}")
    digest = hashlib.sha256(structure_path.read_bytes()).hexdigest()
    if digest != str(manifest["structure_sha256"]).lower():
        raise ValueError("Predicted structure SHA-256 does not match the manifest")
    chain_ids = manifest["chain_ids"]
    if not isinstance(chain_ids, list) or not chain_ids or any(
            not isinstance(chain, str) or not chain.strip() for chain in chain_ids):
        raise ValueError("Prediction manifest chain_ids must be a non-empty string list")
    existing = connection.execute(
        "SELECT id FROM predicted_structure_candidate WHERE campaign_id=? AND structure_sha256=?",
        (campaign_id, digest),
    ).fetchone()
    if existing:
        return str(existing["id"])
    candidate_id = stable_id("PRD")
    connection.execute(
        """INSERT INTO predicted_structure_candidate(
           id, campaign_id, backend, model_name, model_version, construct_name,
           chain_ids_json, structure_path, structure_sha256, ranking_score,
           confidence_json, provenance_json, created_at
           ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)""",
        (candidate_id, campaign_id, str(manifest["backend"]),
         str(manifest["model_name"]), str(manifest["model_version"]),
         str(manifest["construct_name"]), json.dumps(chain_ids), str(structure_path),
         digest, manifest.get("ranking_score"),
         json.dumps(manifest.get("confidence", {}), sort_keys=True),
         json.dumps(manifest["provenance"], sort_keys=True), utc_now()),
    )
    connection.execute(
        "UPDATE campaign SET state='structures_review', updated_at=? WHERE id=?",
        (utc_now(), campaign_id),
    )
    _event(connection, campaign_id, "predicted_structure_registered", campaign["state"],
           "structures_review", "Versioned prediction output registered for review",
           {"candidate_id": candidate_id, "backend": manifest["backend"],
            "structure_sha256": digest})
    return candidate_id


def select_receptor_candidate(
    connection: sqlite3.Connection, user_id: str, campaign_id: str,
    candidate_id: str, rationale: str,
) -> None:
    campaign = _campaign_for_user(connection, campaign_id, user_id)
    if campaign["state"] != "structures_review":
        raise CampaignStateError("Receptor selection requires structures_review state")
    if not rationale.strip():
        raise ValueError("A scientific selection rationale is required")
    experimental = connection.execute(
        "SELECT * FROM structure_candidate WHERE campaign_id=? AND id=?",
        (campaign_id, candidate_id),
    ).fetchone()
    predicted = connection.execute(
        "SELECT * FROM predicted_structure_candidate WHERE campaign_id=? AND id=?",
        (campaign_id, candidate_id),
    ).fetchone()
    if not experimental and not predicted:
        raise KeyError(f"Unknown receptor candidate for {campaign_id}: {candidate_id}")
    kind = "experimental" if experimental else "predicted"
    row = experimental or predicted
    snapshot = dict(row)
    connection.execute(
        """INSERT INTO receptor_selection_lock(
           campaign_id, candidate_kind, candidate_id, snapshot_json, rationale,
           selected_by, created_at) VALUES (?, ?, ?, ?, ?, ?, ?)""",
        (campaign_id, kind, candidate_id, json.dumps(snapshot, sort_keys=True),
         rationale.strip(), user_id, utc_now()),
    )
    connection.execute(
        """UPDATE campaign SET state='structure_selected', selected_structure_id=?,
           updated_at=? WHERE id=?""",
        (candidate_id, utc_now(), campaign_id),
    )
    event_type = "structure_selected" if kind == "experimental" else "predicted_structure_selected"
    _event(connection, campaign_id, event_type, "structures_review",
           "structure_selected", rationale.strip(),
           {"candidate_id": candidate_id, "candidate_kind": kind})


def select_structure(connection: sqlite3.Connection, user_id: str, campaign_id: str, pdb_id: str,
                     rationale: str) -> None:
    candidate = connection.execute(
        "SELECT id FROM structure_candidate WHERE campaign_id=? AND pdb_id=?",
        (campaign_id, pdb_id.upper()),
    ).fetchone()
    if not candidate:
        raise KeyError(f"PDB {pdb_id.upper()} is not a candidate for {campaign_id}")
    select_receptor_candidate(
        connection, user_id, campaign_id, str(candidate["id"]), rationale)


def campaign_status(connection: sqlite3.Connection, user_id: str,
                    campaign_id: str) -> dict[str, Any]:
    campaign = _campaign_for_user(connection, campaign_id, user_id, write=False)
    target = connection.execute(
        "SELECT * FROM target WHERE campaign_id=?", (campaign_id,)
    ).fetchone()
    candidates = connection.execute(
        """SELECT id AS candidate_id, pdb_id, title, experimental_method, resolution_angstrom,
           deposition_date, ligand_ids_json, metadata_json FROM structure_candidate
           WHERE campaign_id=? ORDER BY resolution_angstrom IS NULL, resolution_angstrom, pdb_id""",
        (campaign_id,),
    ).fetchall()
    predicted = connection.execute(
        """SELECT id AS candidate_id, backend, model_name, model_version,
           construct_name, chain_ids_json, structure_path, structure_sha256,
           ranking_score, confidence_json, provenance_json, created_at
           FROM predicted_structure_candidate WHERE campaign_id=? ORDER BY created_at""",
        (campaign_id,),
    ).fetchall()
    selection_lock = connection.execute(
        "SELECT * FROM receptor_selection_lock WHERE campaign_id=?", (campaign_id,)
    ).fetchone()
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
        "selected_receptor": ({
            "candidate_kind": selection_lock["candidate_kind"],
            "candidate_id": selection_lock["candidate_id"],
            "snapshot": json.loads(selection_lock["snapshot_json"]),
            "rationale": selection_lock["rationale"],
            "selected_by": selection_lock["selected_by"],
            "created_at": selection_lock["created_at"],
        } if selection_lock else None),
        "candidates": [
            {key: value for key, value in dict(row).items()
             if key not in {"ligand_ids_json", "metadata_json"}} | {
                "ligand_ids": json.loads(row["ligand_ids_json"]),
                "nonpolymer_ids": json.loads(row["metadata_json"]).get(
                    "nonpolymer_ids", json.loads(row["ligand_ids_json"])),
                "excluded_nonpolymer_components": json.loads(row["metadata_json"]).get(
                    "excluded_nonpolymer_components", []),
            }
            for row in candidates
        ],
        "predicted_candidates": [
            {key: value for key, value in dict(row).items()
             if key not in {"chain_ids_json", "confidence_json", "provenance_json"}} | {
                "candidate_kind": "predicted",
                "chain_ids": json.loads(row["chain_ids_json"]),
                "confidence": json.loads(row["confidence_json"]),
                "provenance": json.loads(row["provenance_json"]),
            }
            for row in predicted
        ],
    }
