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


def create_receptor_eligibility_review_request(
    connection: sqlite3.Connection, user_id: str, campaign_id: str, *,
    requirements: dict[str, Any], computed_evidence: dict[str, Any] | None = None,
    prompt_version: str = "receptor-eligibility-v1",
) -> str:
    campaign = _campaign_for_user(connection, campaign_id, user_id, write=False)
    target = connection.execute(
        "SELECT * FROM target WHERE campaign_id=?", (campaign_id,)
    ).fetchone()
    if not target:
        raise ValueError("Campaign target metadata is required")
    evidence: dict[str, Any] = {
        "target_requirement": requirements,
        "target_metadata": dict(target),
    }
    experimental = connection.execute(
        """SELECT id, pdb_id, title, experimental_method, resolution_angstrom,
           deposition_date, ligand_ids_json, metadata_json
           FROM structure_candidate WHERE campaign_id=? ORDER BY pdb_id""",
        (campaign_id,),
    ).fetchall()
    for row in experimental:
        evidence[f"candidate:{row['id']}"] = {
            "candidate_id": row["id"], "candidate_kind": "experimental",
            "pdb_id": row["pdb_id"], "title": row["title"],
            "experimental_method": row["experimental_method"],
            "resolution_angstrom": row["resolution_angstrom"],
            "deposition_date": row["deposition_date"],
            "ligand_ids": json.loads(row["ligand_ids_json"]),
            "public_metadata": json.loads(row["metadata_json"]),
            "computed": (computed_evidence or {}).get(row["id"], {}),
        }
    predicted = connection.execute(
        """SELECT id, backend, model_name, model_version, construct_name,
           chain_ids_json, ranking_score, confidence_json, provenance_json
           FROM predicted_structure_candidate WHERE campaign_id=? ORDER BY created_at""",
        (campaign_id,),
    ).fetchall()
    for row in predicted:
        evidence[f"candidate:{row['id']}"] = {
            "candidate_id": row["id"], "candidate_kind": "predicted",
            "backend": row["backend"], "model_name": row["model_name"],
            "model_version": row["model_version"],
            "construct_name": row["construct_name"],
            "chain_ids": json.loads(row["chain_ids_json"]),
            "ranking_score": row["ranking_score"],
            "confidence": json.loads(row["confidence_json"]),
            "provenance": json.loads(row["provenance_json"]),
            "computed": (computed_evidence or {}).get(row["id"], {}),
        }
    if not any(key.startswith("candidate:") for key in evidence):
        raise ValueError("At least one receptor candidate is required")
    prompt = (
        "Review every supplied receptor candidate for the stated structural purpose. "
        "Use only supplied evidence. A UniProt match alone does not prove that the "
        "required domain or chain is present. Return JSON with top-level action, "
        "rationale, confidence, evidence_refs, uncertainties, "
        "requires_human_review=true, and decisions. decisions must contain exactly "
        "one object per candidate with candidate_id, verdict "
        "(include|exclude|manual_review), rationale, evidence_refs, and "
        "uncertainties. Do not select a final receptor or infer missing data."
    )
    return create_ai_review_request(
        connection, user_id, campaign["project_id"], campaign_id=campaign_id,
        task_type="receptor_eligibility", subject_type="campaign",
        subject_id=campaign_id, prompt_version=prompt_version, prompt_text=prompt,
        evidence=evidence,
        data_classes=["public_target_metadata", "public_structure_metadata",
                      "computed_summary", "project_configuration"],
    )


def _validate_receptor_eligibility_response(
    request: sqlite3.Row, response: dict[str, Any],
) -> None:
    evidence = json.loads(request["evidence_json"])
    expected = {key.removeprefix("candidate:") for key in evidence
                if key.startswith("candidate:")}
    decisions = response.get("decisions")
    if not isinstance(decisions, list):
        raise ValueError("Receptor eligibility response requires a decisions list")
    received: set[str] = set()
    for decision in decisions:
        if not isinstance(decision, dict):
            raise ValueError("Each receptor eligibility decision must be an object")
        candidate_id = str(decision.get("candidate_id", ""))
        if candidate_id in received:
            raise ValueError(f"Duplicate receptor eligibility decision: {candidate_id}")
        received.add(candidate_id)
        if decision.get("verdict") not in {"include", "exclude", "manual_review"}:
            raise ValueError(f"Invalid receptor eligibility verdict for {candidate_id}")
        if not str(decision.get("rationale", "")).strip():
            raise ValueError(f"Missing receptor eligibility rationale for {candidate_id}")
        refs = decision.get("evidence_refs")
        if not isinstance(refs, list) or not refs:
            raise ValueError(f"Missing receptor eligibility evidence for {candidate_id}")
        unknown = sorted(set(refs) - set(evidence))
        if unknown:
            raise ValueError("Eligibility decision cites unknown evidence: " + ", ".join(unknown))
        uncertainties = decision.get("uncertainties", [])
        if not isinstance(uncertainties, list) or not all(
                isinstance(value, str) for value in uncertainties):
            raise ValueError("Eligibility decision uncertainties must be strings")
    if received != expected:
        raise ValueError("Eligibility decisions must cover exactly all Campaign candidates")
    if response.get("requires_human_review") is not True:
        raise ValueError("Receptor eligibility recommendations require human review")


def create_ligand_query_review_request(
    connection: sqlite3.Connection, user_id: str, campaign_id: str, *,
    selection_requirements: dict[str, Any],
    comparison_summary: dict[str, Any] | None = None,
    prompt_version: str = "ligand-query-selection-v1",
) -> str:
    campaign = _campaign_for_user(connection, campaign_id, user_id, write=False)
    if connection.execute(
            "SELECT 1 FROM query_ligand_lock WHERE campaign_id=?", (campaign_id,)).fetchone():
        raise ValueError("Campaign query ligand is already locked")
    rows = connection.execute(
        """SELECT cl.*, sc.pdb_id FROM campaign_ligand cl
           JOIN structure_candidate sc ON sc.id=cl.structure_candidate_id
           WHERE cl.campaign_id=? ORDER BY sc.pdb_id, cl.ccd_id, cl.id""",
        (campaign_id,),
    ).fetchall()
    if not rows:
        raise ValueError("At least one registered Campaign ligand is required")
    evidence: dict[str, Any] = {
        "selection_requirements": selection_requirements,
        "comparison_summary": comparison_summary or {},
    }
    for row in rows:
        evidence[f"ligand:{row['id']}"] = {
            "ligand_id": row["id"], "pdb_id": row["pdb_id"],
            "ccd_id": row["ccd_id"], "chain_id": row["chain_id"],
            "residue_number": row["residue_number"], "altloc": row["altloc"],
            "standardized_smiles": row["standardized_smiles"],
            "has_3d_descriptor": row["usrcat_json"] is not None,
            "metadata": json.loads(row["metadata_json"]),
        }
    prompt = (
        "Assess every registered co-crystal ligand as a possible similarity-search "
        "query using only supplied evidence. Return JSON with action, rationale, "
        "confidence, evidence_refs, uncertainties, requires_human_review=true, "
        "recommended_ligand_id (a supplied ID or null), and assessments containing "
        "exactly one object per ligand with ligand_id, verdict "
        "(recommend|alternative|exclude|manual_review), rationale, evidence_refs, "
        "and uncertainties. Do not lock a query or start a search."
    )
    return create_ai_review_request(
        connection, user_id, campaign["project_id"], campaign_id=campaign_id,
        task_type="ligand_query_selection", subject_type="campaign",
        subject_id=campaign_id, prompt_version=prompt_version, prompt_text=prompt,
        evidence=evidence,
        data_classes=["public_structure_metadata", "computed_summary",
                      "project_configuration"],
    )


def _validate_ligand_query_response(request: sqlite3.Row,
                                    response: dict[str, Any]) -> None:
    evidence = json.loads(request["evidence_json"])
    expected = {key.removeprefix("ligand:") for key in evidence
                if key.startswith("ligand:")}
    recommended = response.get("recommended_ligand_id")
    if recommended is not None and recommended not in expected:
        raise ValueError("Recommended query ligand is not in the Campaign evidence")
    assessments = response.get("assessments")
    if not isinstance(assessments, list):
        raise ValueError("Ligand query response requires an assessments list")
    received: set[str] = set()
    recommend_verdicts: set[str] = set()
    for assessment in assessments:
        if not isinstance(assessment, dict):
            raise ValueError("Each ligand assessment must be an object")
        ligand_id = str(assessment.get("ligand_id", ""))
        if ligand_id in received:
            raise ValueError(f"Duplicate ligand assessment: {ligand_id}")
        received.add(ligand_id)
        verdict = assessment.get("verdict")
        if verdict not in {"recommend", "alternative", "exclude", "manual_review"}:
            raise ValueError(f"Invalid ligand assessment verdict for {ligand_id}")
        if verdict == "recommend":
            recommend_verdicts.add(ligand_id)
        if not str(assessment.get("rationale", "")).strip():
            raise ValueError(f"Missing ligand assessment rationale for {ligand_id}")
        refs = assessment.get("evidence_refs")
        if not isinstance(refs, list) or not refs or set(refs) - set(evidence):
            raise ValueError(f"Invalid ligand assessment evidence for {ligand_id}")
        uncertainties = assessment.get("uncertainties", [])
        if not isinstance(uncertainties, list) or not all(
                isinstance(value, str) for value in uncertainties):
            raise ValueError("Ligand assessment uncertainties must be strings")
    if received != expected:
        raise ValueError("Ligand assessments must cover exactly all Campaign ligands")
    if recommended is not None and recommend_verdicts != {recommended}:
        raise ValueError("Recommended ligand must be the single recommend verdict")
    if recommended is None and recommend_verdicts:
        raise ValueError("A recommend verdict requires recommended_ligand_id")
    if response.get("requires_human_review") is not True:
        raise ValueError("Ligand query recommendations require human review")


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
    if request["task_type"] == "receptor_eligibility":
        _validate_receptor_eligibility_response(request, response)
    if request["task_type"] == "ligand_query_selection":
        _validate_ligand_query_response(request, response)
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
