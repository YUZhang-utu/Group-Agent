from pathlib import Path

import pytest

from aidd_agent.ai_recommendation import (
    ai_review_status, create_ai_review_request,
    create_receptor_eligibility_review_request, import_ai_recommendation,
)
from aidd_agent.campaign import campaign_status, register_structure_candidates
from aidd_agent.registry import connect
from test_campaign import candidate, setup_campaign


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


def test_receptor_eligibility_review_is_complete_and_advisory(tmp_path: Path) -> None:
    database = tmp_path / "aidd.sqlite3"
    user_id, campaign_id = setup_campaign(database, tmp_path / "workspaces")
    kinase = candidate("8BJU")
    kinase["title"] = "WEE1 kinase domain with inhibitor"
    degron = candidate("9TG7")
    degron["title"] = "beta-TrCP bound by a WEE1 degron peptide"
    with connect(database) as connection:
        register_structure_candidates(
            connection, user_id, campaign_id, [kinase, degron], {"source": "RCSB"})
        campaign = campaign_status(connection, user_id, campaign_id)
        ids = {item["pdb_id"]: item["candidate_id"] for item in campaign["candidates"]}
        request_id = create_receptor_eligibility_review_request(
            connection, user_id, campaign_id,
            requirements={"purpose": "kinase ATP-pocket receptor",
                          "required_domain": {"start": 299, "end": 569}},
            computed_evidence={ids["9TG7"]: {"target_chain_length": 12}},
        )
        request_status = ai_review_status(connection, user_id, request_id)
        response = {
            "action": "review_receptor_eligibility",
            "rationale": "Each candidate was assessed against the required domain.",
            "confidence": 0.95,
            "evidence_refs": ["target_requirement", f"candidate:{ids['8BJU']}",
                              f"candidate:{ids['9TG7']}"],
            "uncertainties": [], "requires_human_review": True,
            "decisions": [
                {"candidate_id": ids["8BJU"], "verdict": "include",
                 "rationale": "The kinase construct matches the requested purpose.",
                 "evidence_refs": ["target_requirement", f"candidate:{ids['8BJU']}"],
                 "uncertainties": []},
                {"candidate_id": ids["9TG7"], "verdict": "exclude",
                 "rationale": "Only a short degron is present, not the kinase domain.",
                 "evidence_refs": ["target_requirement", f"candidate:{ids['9TG7']}"],
                 "uncertainties": []},
            ],
        }
        import_ai_recommendation(
            connection, user_id, request_id, provider="test-provider",
            model_name="test-model", model_version="1", response=response)
        final_status = ai_review_status(connection, user_id, request_id)
        campaign_after = campaign_status(connection, user_id, campaign_id)

    assert request_status["request"]["task_type"] == "receptor_eligibility"
    assert request_status["request"]["evidence"][f"candidate:{ids['9TG7']}"][
        "computed"]["target_chain_length"] == 12
    assert final_status["recommendation"]["response"]["decisions"][1]["verdict"] == "exclude"
    assert campaign_after["state"] == "structures_review"


def test_receptor_eligibility_rejects_partial_or_automatic_decisions(tmp_path: Path) -> None:
    database = tmp_path / "aidd.sqlite3"
    user_id, campaign_id = setup_campaign(database, tmp_path / "workspaces")
    with connect(database) as connection:
        register_structure_candidates(
            connection, user_id, campaign_id, [candidate("8BJU"), candidate("9TG7")], {})
        campaign = campaign_status(connection, user_id, campaign_id)
        ids = [item["candidate_id"] for item in campaign["candidates"]]

        def make_request() -> str:
            return create_receptor_eligibility_review_request(
                connection, user_id, campaign_id,
                requirements={"purpose": "kinase receptor"})

        base = {
            "action": "review", "rationale": "Reviewed supplied evidence.",
            "confidence": 0.7, "evidence_refs": ["target_requirement"],
            "uncertainties": [], "requires_human_review": True,
            "decisions": [{"candidate_id": ids[0], "verdict": "include",
                           "rationale": "Matches purpose.",
                           "evidence_refs": [f"candidate:{ids[0]}"],
                           "uncertainties": []}],
        }
        with pytest.raises(ValueError, match="cover exactly all"):
            import_ai_recommendation(
                connection, user_id, make_request(), provider="p", model_name="m",
                model_version=None, response=base)
        complete = dict(base)
        complete["requires_human_review"] = False
        complete["decisions"] = [
            {"candidate_id": item_id, "verdict": "manual_review",
             "rationale": "Insufficient supplied evidence.",
             "evidence_refs": [f"candidate:{item_id}"], "uncertainties": []}
            for item_id in ids
        ]
        with pytest.raises(ValueError, match="require human review"):
            import_ai_recommendation(
                connection, user_id, make_request(), provider="p", model_name="m",
                model_version=None, response=complete)
