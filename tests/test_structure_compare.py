from pathlib import Path
import hashlib
import sys

import numpy as np
import pytest

from aidd_agent.similarity import ChemistryDependencyError, search_morgan_2d
from aidd_agent.structure_compare import (
    compare_pocket_coordinates, comparison_set_status,
    create_campaign_comparison_set, create_receptor_ensemble,
    generate_comparison_pymol_review, generate_receptor_ensemble_pymol_review,
    receptor_ensemble_status,
    generate_pymol_review, kabsch_align,
)
from test_campaign import candidate, setup_campaign
from aidd_agent.campaign import (
    campaign_status, register_predicted_structure, register_structure_candidates,
)
from aidd_agent.registry import connect


def test_kabsch_recovers_rigid_transform() -> None:
    reference = np.array([[0, 0, 0], [2, 0, 0], [0, 3, 0], [0, 0, 4]], dtype=float)
    rotation = np.array([[0, -1, 0], [1, 0, 0], [0, 0, 1]], dtype=float)
    mobile = reference @ rotation + np.array([10, -4, 7])
    _, _, rmsd = kabsch_align(mobile, reference)
    assert rmsd < 1e-10


def test_pocket_comparison_reports_deviation_and_missing_residues() -> None:
    reference = {1: [0, 0, 0], 2: [2, 0, 0], 3: [0, 3, 0], 4: [0, 0, 4]}
    shifted = {key: np.asarray(value) + [5, 6, 7] for key, value in reference.items()}
    shifted[4] = shifted[4] + [0.5, 0, 0]
    missing = {1: [0, 0, 0], 2: [2, 0, 0]}
    result = compare_pocket_coordinates(
        {"REF": reference, "MOBILE": shifted, "MISSING": missing}, "REF", [1, 2, 3, 4]
    )
    assert result["structures"]["MOBILE"]["status"] == "ok"
    assert result["structures"]["MOBILE"]["per_residue_deviation_angstrom"]["4"] > 0
    assert result["structures"]["MISSING"]["status"] == "insufficient_overlap"
    assert result["structures"]["MISSING"]["missing_residues"] == [3, 4]


def test_pymol_review_script_is_project_scoped(tmp_path: Path) -> None:
    project = tmp_path / "project"
    structures = project / "target" / "review"
    structures.mkdir(parents=True)
    first, second = structures / "1ABC.cif", structures / "2XYZ.cif"
    first.write_text("data_1ABC", encoding="utf-8")
    second.write_text("data_2XYZ", encoding="utf-8")
    output = generate_pymol_review(
        [{"structure_id": "PDB_1ABC", "path": str(first), "chain_id": "A"},
         {"structure_id": "PDB_2XYZ", "path": str(second), "chain_id": "B"}],
        "PDB_1ABC", [12, 62, 95], structures / "review.pml", project,
    )
    script = output.read_text(encoding="utf-8")
    assert "align PDB_2XYZ and chain B, PDB_1ABC" in script
    assert "resi 12+62+95" in script
    assert "organic" in script
    with pytest.raises(ValueError, match="escapes allowed root"):
        generate_pymol_review(
            [{"structure_id": "REF", "path": str(tmp_path / "outside.cif")}],
            "REF", [1], project / "review.pml", project,
        )


def test_similarity_dependency_error_is_actionable(monkeypatch: pytest.MonkeyPatch) -> None:
    # Simulate a machine without RDKit even when the chemistry test environment
    # has it installed. This keeps the optional-dependency test deterministic.
    monkeypatch.setitem(sys.modules, "rdkit", None)
    with pytest.raises(ChemistryDependencyError, match="RDKit is required"):
        search_morgan_2d("CCO", [("MOL-1", "CCO")])


def test_campaign_comparison_and_pymol_review_preserve_review_state(tmp_path: Path) -> None:
    database = tmp_path / "aidd.sqlite3"
    project_root = tmp_path / "project"
    user_id, campaign_id = setup_campaign(database, tmp_path / "workspaces")
    with connect(database) as connection:
        register_structure_candidates(
            connection, user_id, campaign_id,
            [candidate("1ABC"), candidate("2XYZ")], {"source": "test"},
        )
        comparison_id = create_campaign_comparison_set(
            connection, user_id, campaign_id, "ATP-pocket-review",
            ["1abc", "2xyz"], "1abc", [41, 54, 104],
            {"1ABC": "A", "2XYZ": "B"}, "Compare ATP-pocket conformations",
        )
        status = comparison_set_status(connection, user_id, comparison_id)
        assert campaign_status(connection, user_id, campaign_id)["state"] == "structures_review"
    structure_dir = project_root / "inputs" / "structures"
    structure_dir.mkdir(parents=True)
    (structure_dir / "1ABC.cif").write_text("data_1ABC", encoding="utf-8")
    (structure_dir / "2XYZ.cif").write_text("data_2XYZ", encoding="utf-8")
    with connect(database) as connection:
        result = generate_comparison_pymol_review(
            connection, user_id, comparison_id, project_root)
    assert [item["pdb_id"] for item in status["members"]] == ["1ABC", "2XYZ"]
    assert status["members"][0]["is_reference"] is True
    script = Path(result["pymol_script"]).read_text(encoding="utf-8")
    assert "align PDB_2XYZ and chain B, PDB_1ABC" in script
    assert "resi 41+54+104" in script


def test_campaign_comparison_rejects_foreign_pdb_and_missing_file(tmp_path: Path) -> None:
    database = tmp_path / "aidd.sqlite3"
    user_id, campaign_id = setup_campaign(database, tmp_path / "workspaces")
    with connect(database) as connection:
        register_structure_candidates(
            connection, user_id, campaign_id,
            [candidate("1ABC"), candidate("2XYZ")], {},
        )
        with pytest.raises(ValueError, match="do not belong"):
            create_campaign_comparison_set(
                connection, user_id, campaign_id, "bad", ["1ABC", "9ZZZ"],
                "1ABC", [1, 2, 3], {}, "Invalid foreign candidate",
            )
        comparison_id = create_campaign_comparison_set(
            connection, user_id, campaign_id, "valid", ["1ABC", "2XYZ"],
            "1ABC", [1, 2, 3], {}, "Valid comparison",
        )
        with pytest.raises(FileNotFoundError, match="Missing downloaded structure files"):
            generate_comparison_pymol_review(
                connection, user_id, comparison_id, tmp_path / "project")


def test_all_candidate_receptor_ensemble_includes_prediction_and_ligands(tmp_path: Path) -> None:
    database = tmp_path / "aidd.sqlite3"
    project_root = tmp_path / "project"
    user_id, campaign_id = setup_campaign(database, tmp_path / "workspaces")
    first, second, fragment = candidate("8BJU"), candidate("8ZZZ"), candidate("9TG7")
    first["ligand_ids"] = ["ATP"]
    second["ligand_ids"] = ["LIG"]
    predicted_path = tmp_path / "af3_wee1.cif"
    predicted_path.write_text("data_af3\n", encoding="utf-8")
    manifest = {
        "backend": "alphafold3", "model_name": "AlphaFold 3",
        "model_version": "commit", "construct_name": "WEE1_P30291_299_569",
        "chain_ids": ["A"], "structure_path": str(predicted_path),
        "structure_sha256": hashlib.sha256(predicted_path.read_bytes()).hexdigest(),
        "provenance": {},
    }
    with connect(database) as connection:
        register_structure_candidates(
            connection, user_id, campaign_id, [first, second, fragment], {})
        predicted_id = register_predicted_structure(
            connection, user_id, campaign_id, manifest)
        ensemble_id = create_receptor_ensemble(
            connection, user_id, campaign_id, "all-target-structures", "8BJU",
            [320, 337, 463], {"8BJU": "A", predicted_id: "A"},
            "Observe the complete experimental and predicted ensemble",
            exclusions={"9TG7": "WEE1 is only a 12-residue degron peptide"})
        status = receptor_ensemble_status(connection, user_id, ensemble_id)
    assert len(status["members"]) == 3
    assert status["members"][0]["display_id"] == "8BJU"
    assert any(item["candidate_kind"] == "predicted" for item in status["members"])
    assert status["excluded_candidates"] == [{
        "candidate_id": status["excluded_candidates"][0]["candidate_id"],
        "candidate_kind": "experimental",
        "reason": "WEE1 is only a 12-residue degron peptide",
        "pdb_id": "9TG7",
    }]
    structure_dir = project_root / "inputs" / "structures"
    structure_dir.mkdir(parents=True)
    (structure_dir / "8BJU.cif").write_text("data_8BJU\n", encoding="utf-8")
    (structure_dir / "8ZZZ.cif").write_text("data_8ZZZ\n", encoding="utf-8")
    with connect(database) as connection:
        result = generate_receptor_ensemble_pymol_review(
            connection, user_id, ensemble_id, project_root,
            project_root / "target" / "reviews" / "all.pml")
    script = Path(result["pymol_script"]).read_text(encoding="utf-8")
    assert "cealign PDB_8BJU" in script
    assert "resn ATP" in script
    assert "resn LIG" in script
    assert predicted_id.replace("-", "_") in script
    assert "within 6 of PDB_8BJU_ligands" in script
    assert "resi 22+39+165" in script
