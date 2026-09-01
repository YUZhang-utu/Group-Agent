from pathlib import Path
import hashlib

import pytest

from aidd_agent.campaign import (
    CampaignStateError, campaign_status, create_campaign,
    register_predicted_structure, register_structure_candidates,
    select_receptor_candidate, select_structure, set_target,
)
from aidd_agent.rcsb import search_structures
from aidd_agent.registry import connect, initialize
from aidd_agent.context import activate_task, create_task, create_user
from aidd_agent.project_context import activate_project, create_scientific_project


def candidate(pdb_id: str = "1ABC") -> dict:
    return {"pdb_id": pdb_id, "title": "Target with inhibitor",
            "experimental_method": "X-RAY DIFFRACTION", "resolution_angstrom": 1.8,
            "deposition_date": "2025-01-01", "ligand_ids": ["LIG"]}


def setup_campaign(path: Path, workspace: Path) -> tuple[str, str]:
    initialize(path)
    with connect(path) as connection:
        user_id = create_user(connection, "alice")
        project_id = create_scientific_project(
            connection, user_id, "Kinase Discovery", "Find selective binders", workspace)
        activate_project(connection, user_id, project_id)
        task_id = create_task(connection, user_id, project_id, "Structure Selection",
                              "Select a target structure")
        activate_task(connection, user_id, task_id)
        campaign_id = create_campaign(connection, user_id, project_id, task_id,
                                      "macrocycle-campaign", "Find selective binders")
        set_target(connection, user_id, campaign_id, name="Example kinase", organism="Homo sapiens",
                   uniprot_id="P00001", gene_name="EXK")
    return user_id, campaign_id


def test_campaign_transition_and_selection_audit(tmp_path: Path) -> None:
    database = tmp_path / "aidd.sqlite3"
    user_id, campaign_id = setup_campaign(database, tmp_path / "workspaces")
    with connect(database) as connection:
        assert register_structure_candidates(connection, user_id, campaign_id, [candidate()], {"q": 1}) == 1
        assert register_structure_candidates(connection, user_id, campaign_id, [candidate()], {"q": 1}) == 0
        select_structure(connection, user_id, campaign_id, "1abc", "Best resolution and relevant bound ligand")
        status = campaign_status(connection, user_id, campaign_id)
        events = connection.execute(
            "SELECT event_type, rationale FROM decision_event WHERE campaign_id=? ORDER BY id",
            (campaign_id,),
        ).fetchall()
    assert status["state"] == "structure_selected"
    assert status["selected_pdb_id"] == "1ABC"
    assert [row["event_type"] for row in events] == [
        "campaign_created", "target_set", "structures_registered",
        "structures_registered", "structure_selected"]
    assert events[-1]["rationale"] == "Best resolution and relevant bound ligand"


def test_invalid_transitions_and_rationale_are_rejected(tmp_path: Path) -> None:
    database = tmp_path / "aidd.sqlite3"
    initialize(database)
    with connect(database) as connection:
        user_id = create_user(connection, "alice")
        project_id = create_scientific_project(
            connection, user_id, "Project", "Objective", tmp_path / "storage")
        activate_project(connection, user_id, project_id)
        task_id = create_task(connection, user_id, project_id, "Task", "Objective")
        activate_task(connection, user_id, task_id)
        campaign_id = create_campaign(connection, user_id, project_id, task_id, "c", "o")
        with pytest.raises(CampaignStateError):
            register_structure_candidates(connection, user_id, campaign_id, [candidate()], {})
        set_target(connection, user_id, campaign_id, name="T", organism="human")
        register_structure_candidates(connection, user_id, campaign_id, [candidate()], {})
        with pytest.raises(ValueError, match="rationale"):
            select_structure(connection, user_id, campaign_id, "1ABC", "")
        with pytest.raises(KeyError):
            select_structure(connection, user_id, campaign_id, "9ZZZ", "looks good")


def test_mocked_rcsb_search_preserves_selection_metadata() -> None:
    def post(url: str, payload: dict) -> dict:
        if "search.rcsb" in url:
            return {"result_set": [{"identifier": "2xyz"}]}
        assert payload["variables"]["ids"] == ["2XYZ"]
        return {"data": {"entries": [{
            "rcsb_id": "2XYZ", "struct": {"title": "Protein-ligand complex"},
            "exptl": [{"method": "X-RAY DIFFRACTION"}],
            "rcsb_entry_info": {"resolution_combined": [2.1]},
            "rcsb_accession_info": {"deposit_date": "2024-02-03"},
            "nonpolymer_entities": [
                {"pdbx_entity_nonpoly": {"comp_id": "ATP"}},
                {"pdbx_entity_nonpoly": {"comp_id": "GOL"}},
                {"pdbx_entity_nonpoly": {"comp_id": "CL"}},
                {"pdbx_entity_nonpoly": {"comp_id": "NA"}},
                {"pdbx_entity_nonpoly": {"comp_id": "EDO"}},
                {"pdbx_entity_nonpoly": {"comp_id": "PO4"}},
                {"pdbx_entity_nonpoly": {"comp_id": "MG"}},
            ],
        }]}}

    query, candidates = search_structures("P00001", max_resolution=2.5, post_json=post)
    assert query["return_type"] == "entry"
    assert candidates == [{"pdb_id": "2XYZ", "title": "Protein-ligand complex",
                           "experimental_method": "X-RAY DIFFRACTION",
                           "resolution_angstrom": 2.1, "deposition_date": "2024-02-03",
                           "ligand_ids": ["ATP"],
                           "nonpolymer_ids": ["ATP", "CL", "EDO", "GOL", "MG", "NA", "PO4"],
                           "excluded_nonpolymer_components": [
                               {"component_id": "CL", "reason": "inorganic_ion"},
                               {"component_id": "EDO", "reason": "solvent_or_cryoprotectant"},
                               {"component_id": "GOL", "reason": "solvent_or_cryoprotectant"},
                               {"component_id": "MG", "reason": "inorganic_ion"},
                               {"component_id": "NA", "reason": "inorganic_ion"},
                               {"component_id": "PO4", "reason": "buffer_or_salt"},
                           ]}]


def test_rcsb_search_accepts_null_optional_lists() -> None:
    def post(url: str, payload: dict) -> dict:
        if "search.rcsb" in url:
            return {"result_set": [{"identifier": "3abc"}]}
        return {"data": {"entries": [None, {
            "rcsb_id": "3ABC", "struct": None, "exptl": None,
            "rcsb_entry_info": {"resolution_combined": [2.4]},
            "rcsb_accession_info": None, "nonpolymer_entities": None,
        }]}}

    _, candidates = search_structures("P00001", post_json=post)
    assert candidates == [{
        "pdb_id": "3ABC", "title": "", "experimental_method": None,
        "resolution_angstrom": 2.4, "deposition_date": None, "ligand_ids": [],
        "nonpolymer_ids": [], "excluded_nonpolymer_components": [],
    }]


def test_candidate_reregistration_refreshes_ligand_classification(tmp_path: Path) -> None:
    database = tmp_path / "aidd.sqlite3"
    user_id, campaign_id = setup_campaign(database, tmp_path / "workspaces")
    old = candidate()
    old["ligand_ids"] = ["GOL", "LIG"]
    refreshed = candidate()
    refreshed.update({
        "ligand_ids": ["LIG"], "nonpolymer_ids": ["GOL", "LIG"],
        "excluded_nonpolymer_components": [
            {"component_id": "GOL", "reason": "solvent_or_cryoprotectant"}
        ],
    })
    with connect(database) as connection:
        assert register_structure_candidates(
            connection, user_id, campaign_id, [old], {"version": 1}) == 1
        assert register_structure_candidates(
            connection, user_id, campaign_id, [refreshed], {"version": 2}) == 0
        status = campaign_status(connection, user_id, campaign_id)
    item = status["candidates"][0]
    assert item["ligand_ids"] == ["LIG"]
    assert item["nonpolymer_ids"] == ["GOL", "LIG"]
    assert item["excluded_nonpolymer_components"][0]["component_id"] == "GOL"


def test_predicted_structure_registration_and_version_lock(tmp_path: Path) -> None:
    database = tmp_path / "aidd.sqlite3"
    user_id, campaign_id = setup_campaign(database, tmp_path / "workspaces")
    structure = tmp_path / "wee1_af3.cif"
    structure.write_text("data_wee1\n#\n", encoding="utf-8")
    digest = hashlib.sha256(structure.read_bytes()).hexdigest()
    manifest = {
        "backend": "alphafold3", "model_name": "AlphaFold 3",
        "model_version": "source-abc123", "construct_name": "WEE1_299_569",
        "chain_ids": ["A"], "structure_path": str(structure),
        "structure_sha256": digest, "ranking_score": 0.85,
        "confidence": {"ranking_score": 0.85},
        "provenance": {"input_sha256": "input-hash", "sif_sha256": "sif-hash"},
    }
    with connect(database) as connection:
        candidate_id = register_predicted_structure(
            connection, user_id, campaign_id, manifest)
        assert register_predicted_structure(
            connection, user_id, campaign_id, manifest) == candidate_id
        status = campaign_status(connection, user_id, campaign_id)
        assert status["state"] == "structures_review"
        assert status["predicted_candidates"][0]["structure_sha256"] == digest
        select_receptor_candidate(
            connection, user_id, campaign_id, candidate_id,
            "AF3 kinase core agrees with the reviewed experimental ensemble")
        selected = campaign_status(connection, user_id, campaign_id)
    assert selected["selected_pdb_id"] is None
    assert selected["selected_receptor"]["candidate_kind"] == "predicted"
    assert selected["selected_receptor"]["snapshot"]["structure_sha256"] == digest


def test_predicted_structure_rejects_changed_file(tmp_path: Path) -> None:
    database = tmp_path / "aidd.sqlite3"
    user_id, campaign_id = setup_campaign(database, tmp_path / "workspaces")
    structure = tmp_path / "changed.cif"
    structure.write_text("data_changed\n", encoding="utf-8")
    manifest = {
        "backend": "alphafold3", "model_name": "AlphaFold 3",
        "model_version": "v1", "construct_name": "target_domain",
        "chain_ids": ["A"], "structure_path": str(structure),
        "structure_sha256": "0" * 64, "provenance": {},
    }
    with connect(database) as connection:
        with pytest.raises(ValueError, match="SHA-256"):
            register_predicted_structure(connection, user_id, campaign_id, manifest)
