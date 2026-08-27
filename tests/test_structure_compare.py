from pathlib import Path

import numpy as np
import pytest

from aidd_agent.similarity import ChemistryDependencyError, search_morgan_2d
from aidd_agent.structure_compare import (
    compare_pocket_coordinates, generate_pymol_review, kabsch_align,
)


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


def test_similarity_dependency_error_is_actionable() -> None:
    with pytest.raises(ChemistryDependencyError, match="RDKit is required"):
        search_morgan_2d("CCO", [("MOL-1", "CCO")])
