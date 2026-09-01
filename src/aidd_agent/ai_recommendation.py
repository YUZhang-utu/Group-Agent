from __future__ import annotations

import json
import sqlite3
from typing import Any, Sequence

from .campaign import _campaign_for_user
from .project_context import require_active_project, require_project_owner
from .registry import stable_id, utc_now


FORBIDDEN_MODEL_DATA_CLASSES = {"experimental_raw", "experimental_measurements"}
ALLOWED_MODEL_DATA_CLASSES = {
    "public_target_metadata", "public_structure_metadata", "literature",
    "computed_summary", "project_configuration", "user_prompt",
}


def _authorize(connection: sqlite3.Connection, user_id: str, project_id: str,
               campaign_id: str | None = None) -> None:
    require_project_owner(connection, user_id, project_id)
    require_active_project(connection, user_id, project_id)
    if campaign_id:
        campaign = _campaign_for_user(connection, campaign_id, user_id, write=False)
        if campaign["project_id"] != project_id:
            raise ValueError("Campaign does not belong to the Project")


def create_ai_review_request(
    connection: sqlite3.Connection, user_id: str, project_id: str, *,
    task_type: str, subject_type: str, subject_id: str, prompt_version: str,
    prompt_text: str, evidence: dict[str, Any], data_classes: Sequence[str],
    campaign_id: str | None = None,
) -> str:
    _authorize(connection, user_id, project_id, campaign_id)
    required = (task_type, subject_type, subject_id, prompt_version, prompt_text)
    if any(not str(value).strip() for value in required):
        raise ValueError("AI review task, subject, prompt version, and prompt are required")
    classes = sorted(set(data_classes))
    forbidden = sorted(set(classes) & FORBIDDEN_MODEL_DATA_CLASSES)
    if forbidden:
        raise ValueError("Experimental data cannot enter AI context: " + ", ".join(forbidden))
    unknown = sorted(set(classes) - ALLOWED_MODEL_DATA_CLASSES)
    if unknown:
        raise ValueError("Unknown AI data classes: " + ", ".join(unknown))
    request_id = stable_id("AIR")
    privacy = {
        "experimental_data_visible_to_model": False,
        "forbidden_data_classes": sorted(FORBIDDEN_MODEL_DATA_CLASSES),
    }
    connection.execute(
        """INSERT INTO ai_review_request(
           id, project_id, campaign_id, user_id, task_type, subject_type,
           subject_id, prompt_version, prompt_text, evidence_json,
           data_classes_json, privacy_policy_json, status, created_at)
           VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, 'pending', ?)""",
        (request_id, project_id, campaign_id, user_id, task_type.strip(),
         subject_type.strip(), subject_id.strip(), prompt_version.strip(),
         prompt_text.strip(), json.dumps(evidence, sort_keys=True),
         json.dumps(classes), json.dumps(privacy, sort_keys=True), utc_now()),
    )
    return request_id


def import_ai_recommendation(
    connection: sqlite3.Connection, user_id: str, request_id: str, *,
    provider: str, model_name: str, model_version: str | None,
    response: dict[str, Any],
) -> str:
    request = connection.execute(
        "SELECT * FROM ai_review_request WHERE id=?", (request_id,),
    ).fetchone()
    if not request:
        raise KeyError(f"Unknown AI review request: {request_id}")
    _authorize(connection, user_id, request["project_id"], request["campaign_id"])
    if request["status"] != "pending":
        raise ValueError("AI review request already has a recommendation")
    action = str(response.get("action", "")).strip()
    rationale = str(response.get("rationale", "")).strip()
    if not action or not rationale:
        raise ValueError("AI recommendation action and rationale are required")
    confidence = response.get("confidence")
    if not isinstance(confidence, (int, float)) or isinstance(confidence, bool) \
            or not 0 <= float(confidence) <= 1:
        raise ValueError("AI recommendation confidence must be between 0 and 1")
    evidence_refs = response.get("evidence_refs", [])
    if not isinstance(evidence_refs, list) or not all(isinstance(x, str) for x in evidence_refs):
        raise ValueError("evidence_refs must be a list of strings")
    evidence = json.loads(request["evidence_json"])
    unknown_refs = sorted(set(evidence_refs) - set(evidence))
    if unknown_refs:
        raise ValueError("Recommendation cites unknown evidence: " + ", ".join(unknown_refs))
    uncertainties = response.get("uncertainties", [])
    if not isinstance(uncertainties, list) or not all(isinstance(x, str) for x in uncertainties):
        raise ValueError("uncertainties must be a list of strings")
    human_review = response.get("requires_human_review", True)
    if not isinstance(human_review, bool):
        raise ValueError("requires_human_review must be boolean")
    if not provider.strip() or not model_name.strip():
        raise ValueError("AI provider and model name are required")
    recommendation_id = stable_id("REC")
    connection.execute(
        """INSERT INTO ai_recommendation(
           id, request_id, provider, model_name, model_version, action,
           rationale, confidence, evidence_refs_json, uncertainties_json,
           requires_human_review, response_json, created_at)
           VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)""",
        (recommendation_id, request_id, provider.strip(), model_name.strip(),
         model_version, action, rationale, float(confidence),
         json.dumps(evidence_refs), json.dumps(uncertainties), int(human_review),
         json.dumps(response, sort_keys=True), utc_now()),
    )
    connection.execute(
        "UPDATE ai_review_request SET status='completed' WHERE id=?", (request_id,),
    )
    return recommendation_id


def ai_review_status(connection: sqlite3.Connection, user_id: str,
                     request_id: str) -> dict[str, Any]:
    request = connection.execute(
        "SELECT * FROM ai_review_request WHERE id=?", (request_id,),
    ).fetchone()
    if not request:
        raise KeyError(f"Unknown AI review request: {request_id}")
    _authorize(connection, user_id, request["project_id"], request["campaign_id"])
    recommendation = connection.execute(
        "SELECT * FROM ai_recommendation WHERE request_id=?", (request_id,),
    ).fetchone()
    request_result = dict(request)
    for key in ("evidence_json", "data_classes_json", "privacy_policy_json"):
        request_result[key.removesuffix("_json")] = json.loads(request_result.pop(key))
    recommendation_result = None
    if recommendation:
        recommendation_result = dict(recommendation)
        for key in ("evidence_refs_json", "uncertainties_json", "response_json"):
            recommendation_result[key.removesuffix("_json")] = json.loads(
                recommendation_result.pop(key))
        recommendation_result["requires_human_review"] = bool(
            recommendation_result["requires_human_review"])
    return {"request": request_result, "recommendation": recommendation_result}
