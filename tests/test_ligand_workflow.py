from pathlib import Path

import pytest

from aidd_agent.campaign import register_structure_candidates
from aidd_agent.ai_recommendation import (
    ai_review_status, create_ligand_query_review_request, import_ai_recommendation,
)
from aidd_agent.ligand_workflow import (
    campaign_ligand_status, register_campaign_ligands,
    run_hierarchical_library_search, select_query_ligand,
)
from aidd_agent.registry import connect, ensure_library
from aidd_agent.similarity import (
    HierarchicalSimilarityHit, SimilarityHit, hierarchical_similarity_search,
)
from test_campaign import candidate, setup_campaign


def test_campaign_ligand_registration_filter_and_query_lock(tmp_path: Path) -> None:
    database = tmp_path / "aidd.sqlite3"
    user_id, campaign_id = setup_campaign(database, tmp_path / "workspaces")
    pdb = candidate("8BJU")
    pdb["ligand_ids"] = ["ATP", "LIG"]
    with connect(database) as connection:
        register_structure_candidates(connection, user_id, campaign_id, [pdb], {})
        with pytest.raises(ValueError, match="not a retained ligand"):
            register_campaign_ligands(connection, user_id, campaign_id, [{
                "pdb_id": "8BJU", "ccd_id": "GOL", "chain_id": "A",
                "residue_number": "901", "standardized_smiles": "OCC(O)CO",
            }])
        ids = register_campaign_ligands(connection, user_id, campaign_id, [{
            "pdb_id": "8BJU", "ccd_id": "LIG", "chain_id": "A",
            "residue_number": "901", "standardized_smiles": "CC1=CC=CC=C1",
            "metadata": {"standardization": "neutral-parent-v1"},
        }])
        assert len(ids) == 1
        assert register_campaign_ligands(connection, user_id, campaign_id, [{
            "pdb_id": "8BJU", "ccd_id": "LIG", "chain_id": "A",
            "residue_number": "901", "standardized_smiles": "CC1=CC=CC=C1",
        }]) == []
        select_query_ligand(
            connection, user_id, campaign_id, ids[0], "Representative ATP-site chemotype")
        with pytest.raises(ValueError, match="already locked"):
            select_query_ligand(connection, user_id, campaign_id, ids[0], "Again")
        status = campaign_ligand_status(connection, user_id, campaign_id)
    assert status["ligands"][0]["ccd_id"] == "LIG"
    assert status["query_lock"]["ligand_id"] == ids[0]
    assert status["query_lock"]["snapshot"]["standardized_smiles"] == "CC1=CC=CC=C1"


def test_hierarchical_search_keeps_best_conformer_and_2d_fallback(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path,
) -> None:
    monkeypatch.setattr(
        "aidd_agent.similarity.query_morgan_index",
        lambda *_args, **_kwargs: [SimilarityHit("M1", 0.8), SimilarityHit("M2", 0.7)])
    monkeypatch.setattr(
        "aidd_agent.similarity.query_usrcat_index",
        lambda *_args, **_kwargs: [SimilarityHit("C1", 0.9), SimilarityHit("C2", 0.6)])
    hits = hierarchical_similarity_search(
        "CC", [0.0] * 60, tmp_path / "morgan.json.gz", tmp_path / "shape.npz",
        conformer_to_molecule={"C1": "M1", "C2": "M1"}, two_d_pool=2, limit=2)
    assert [hit.molecule_id for hit in hits] == ["M1", "M2"]
    assert hits[0].conformer_id == "C1"
    assert hits[0].stages_completed == ("morgan2d", "usrcat3d")
    assert hits[1].usrcat_similarity is None
    assert hits[1].stages_completed == ("morgan2d",)


def test_hierarchical_search_validates_limits(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="limit"):
        hierarchical_similarity_search(
            "CC", None, tmp_path / "index", two_d_pool=5, limit=6)


def test_ligand_query_ai_review_is_advisory_and_complete(tmp_path: Path) -> None:
    database = tmp_path / "aidd.sqlite3"
    user_id, campaign_id = setup_campaign(database, tmp_path / "workspaces")
    pdb = candidate("8BJU")
    pdb["ligand_ids"] = ["L1", "L2"]
    with connect(database) as connection:
        register_structure_candidates(connection, user_id, campaign_id, [pdb], {})
        ids = register_campaign_ligands(connection, user_id, campaign_id, [
            {"pdb_id": "8BJU", "ccd_id": "L1", "chain_id": "A",
             "residue_number": "901", "standardized_smiles": "CC"},
            {"pdb_id": "8BJU", "ccd_id": "L2", "chain_id": "A",
             "residue_number": "902", "standardized_smiles": "CCC"},
        ])
        request_id = create_ligand_query_review_request(
            connection, user_id, campaign_id,
            selection_requirements={"purpose": "macrocycle similarity search"})
        response = {
            "action": "recommend_query", "rationale": "Reviewed both ligands.",
            "confidence": 0.8, "evidence_refs": [f"ligand:{ids[0]}"],
            "uncertainties": [], "requires_human_review": True,
            "recommended_ligand_id": ids[0],
            "assessments": [
                {"ligand_id": ids[0], "verdict": "recommend",
                 "rationale": "Best-supported query.",
                 "evidence_refs": [f"ligand:{ids[0]}"], "uncertainties": []},
                {"ligand_id": ids[1], "verdict": "alternative",
                 "rationale": "Useful alternative chemotype.",
                 "evidence_refs": [f"ligand:{ids[1]}"], "uncertainties": []},
            ],
        }
        import_ai_recommendation(
            connection, user_id, request_id, provider="test", model_name="model",
            model_version="1", response=response)
        status = ai_review_status(connection, user_id, request_id)
        ligand_status = campaign_ligand_status(connection, user_id, campaign_id)
    assert status["recommendation"]["response"]["recommended_ligand_id"] == ids[0]
    assert ligand_status["query_lock"] is None


def test_campaign_hierarchical_search_requires_lock_and_records_provenance(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path,
) -> None:
    database = tmp_path / "aidd.sqlite3"
    user_id, campaign_id = setup_campaign(database, tmp_path / "workspaces")
    index = tmp_path / "morgan.json.gz"
    index.write_bytes(b"index")
    with connect(database) as connection:
        library_id = ensure_library(connection, "search-library")
        with pytest.raises(ValueError, match="locked"):
            run_hierarchical_library_search(
                connection, user_id, campaign_id, library_id, index)
        pdb = candidate("8BJU")
        pdb["ligand_ids"] = ["LIG"]
        register_structure_candidates(connection, user_id, campaign_id, [pdb], {})
        ligand_id = register_campaign_ligands(connection, user_id, campaign_id, [{
            "pdb_id": "8BJU", "ccd_id": "LIG", "chain_id": "A",
            "residue_number": "901", "standardized_smiles": "CC",
        }])[0]
        select_query_ligand(connection, user_id, campaign_id, ligand_id, "Reviewed query")
        monkeypatch.setattr(
            "aidd_agent.ligand_workflow.hierarchical_similarity_search",
            lambda *_args, **_kwargs: [HierarchicalSimilarityHit(
                "M1", "C1", 1, 0.8, 0.9, 0.86, ("morgan2d", "usrcat3d"))])
        result = run_hierarchical_library_search(
            connection, user_id, campaign_id, library_id, index, two_d_pool=10, limit=1)
        stored = connection.execute(
            "SELECT * FROM ligand_similarity_search WHERE id=?", (result["search_id"],)
        ).fetchone()
    assert result["hits"][0]["molecule_id"] == "M1"
    assert stored["query_ligand_id"] == ligand_id
    assert result["parameters"]["morgan_sha256"]
