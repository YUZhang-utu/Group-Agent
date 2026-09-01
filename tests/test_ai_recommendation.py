from pathlib import Path

import pytest

from aidd_agent.ai_recommendation import (
    ai_review_status, create_ai_review_request, import_ai_recommendation,
)
from aidd_agent.campaign import campaign_status
from aidd_agent.registry import connect
from test_campaign import setup_campaign


def test_ai_recommendation_round_trip_does_not_mutate_campaign(tmp_path: Path) -> None:
    database = tmp_path / "aidd.sqlite3"
    user_id, campaign_id = setup_campaign(database, tmp_path / "workspaces")
    with connect(database) as connection:
        campaign = campaign_status(connection, user_id, campaign_id)
        request_id = create_ai_review_request(
            connection, user_id, campaign["project_id"], campaign_id=campaign_id,
            task_type="classify_nonpolymer", subject_type="pdb_component",
            subject_id="1ABC:GOL", prompt_version="component-review-v1",
            prompt_text="Classify this component using only the supplied evidence.",
            evidence={"ccd": {"name": "GLYCEROL"}, "pocket": {"distance": 14.2}},
            data_classes=["public_structure_metadata", "computed_summary"],
        )
        recommendation_id = import_ai_recommendation(
            connection, user_id, request_id, provider="test-provider",
            model_name="test-model", model_version="1",
            response={
                "action": "exclude_from_candidate_ligands",
                "rationale": "Glycerol is distant from the reviewed pocket.",
                "confidence": 0.96, "evidence_refs": ["ccd", "pocket"],
                "uncertainties": ["No occupancy evidence supplied"],
                "requires_human_review": True,
            },
        )
        status = ai_review_status(connection, user_id, request_id)
        campaign_after = campaign_status(connection, user_id, campaign_id)
    assert recommendation_id.startswith("REC-")
    assert status["request"]["status"] == "completed"
    assert status["request"]["privacy_policy"]["experimental_data_visible_to_model"] is False
    assert status["recommendation"]["confidence"] == 0.96
    assert status["recommendation"]["requires_human_review"] is True
    assert campaign_after["state"] == "target_set"


def test_ai_review_rejects_experimental_data_and_invalid_evidence(tmp_path: Path) -> None:
    database = tmp_path / "aidd.sqlite3"
    user_id, campaign_id = setup_campaign(database, tmp_path / "workspaces")
    with connect(database) as connection:
        campaign = campaign_status(connection, user_id, campaign_id)
        with pytest.raises(ValueError, match="Experimental data cannot enter AI context"):
            create_ai_review_request(
                connection, user_id, campaign["project_id"], campaign_id=campaign_id,
                task_type="rank", subject_type="campaign", subject_id=campaign_id,
                prompt_version="v1", prompt_text="Rank candidates", evidence={},
                data_classes=["experimental_raw"],
            )
        request_id = create_ai_review_request(
            connection, user_id, campaign["project_id"], campaign_id=campaign_id,
            task_type="rank", subject_type="campaign", subject_id=campaign_id,
            prompt_version="v1", prompt_text="Rank candidates",
            evidence={"structures": []}, data_classes=["public_structure_metadata"],
        )
        with pytest.raises(ValueError, match="unknown evidence"):
            import_ai_recommendation(
                connection, user_id, request_id, provider="p", model_name="m",
                model_version=None,
                response={"action": "review", "rationale": "Reason",
                          "confidence": 0.5, "evidence_refs": ["missing"],
                          "uncertainties": [], "requires_human_review": True},
            )


def test_ai_review_rejects_invalid_confidence(tmp_path: Path) -> None:
    database = tmp_path / "aidd.sqlite3"
    user_id, campaign_id = setup_campaign(database, tmp_path / "workspaces")
    with connect(database) as connection:
        campaign = campaign_status(connection, user_id, campaign_id)
        request_id = create_ai_review_request(
            connection, user_id, campaign["project_id"], campaign_id=campaign_id,
            task_type="review", subject_type="campaign", subject_id=campaign_id,
            prompt_version="v1", prompt_text="Review", evidence={},
            data_classes=["user_prompt"],
        )
        with pytest.raises(ValueError, match="between 0 and 1"):
            import_ai_recommendation(
                connection, user_id, request_id, provider="p", model_name="m",
                model_version=None,
                response={"action": "review", "rationale": "Reason",
                          "confidence": 1.2, "evidence_refs": [],
                          "uncertainties": [], "requires_human_review": True},
            )
