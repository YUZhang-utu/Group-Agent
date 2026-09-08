import hashlib
import json
from pathlib import Path

import numpy as np
import pytest

from aidd_agent.conformer_artifacts import FEATURE_DTYPE, META_DTYPE
import aidd_agent.predocking_qc as predocking_qc
from aidd_agent.predocking_qc import run_predocking_pocket_qc


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _artifact_catalog(root: Path) -> Path:
    shard = root / "shard-0"
    shard.mkdir(parents=True)
    np.asarray([[0, 0, 0], [100, 0, 0]], dtype="<i2").tofile(shard / "coords.bin")
    features = np.zeros(1, dtype=FEATURE_DTYPE)
    features.tofile(shard / "feats.bin")
    meta = np.asarray([
        (0, 0, 0, 2, 1, (0, 0, 0), (1, 0, 0), (0, 0, 0),
         (0, 0, 0, 0, 0, 0), (0, 0))
    ], dtype=META_DTYPE)
    meta.tofile(shard / "meta.bin")
    np.asarray([b"C1"], dtype="S16").tofile(shard / "conformer_ids.bin")
    np.asarray([b"M1"], dtype="S16").tofile(shard / "molecule_ids.bin")
    (shard / "manifest.json").write_text(json.dumps({
        "format": "aidd-conformer-artifact-shard", "version": 1,
        "global_id_start": 0, "conformers": 1, "coordinate_scale": 100.0,
    }), encoding="utf-8")
    catalog = root / "catalog.json"
    catalog.write_text(json.dumps({
        "format": "aidd-conformer-artifact-catalog", "version": 1,
        "shards": [{"name": "shard-0", "path": str(shard),
                    "global_id_start": 0, "conformers": 1}],
    }), encoding="utf-8")
    return catalog


def _receptor(path: Path) -> Path:
    path.write_text("""data_TEST
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.pdbx_formal_charge
_atom_site.auth_seq_id
_atom_site.auth_comp_id
_atom_site.auth_asym_id
_atom_site.auth_atom_id
_atom_site.pdbx_PDB_model_num
HETATM 1 C C1 . LIG A 2 . ? 10.000 0.000 0.000 1.00 10.0 ? 901 LIG A C1 1
HETATM 2 O O1 . LIG A 2 . ? 11.000 0.000 0.000 1.00 10.0 ? 901 LIG A O1 1
ATOM 3 C CA . ALA A 1 1 ? 11.000 0.500 0.000 1.00 10.0 ? 1 ALA A CA 1
ATOM 4 N N . ALA A 1 1 ? 11.000 2.500 0.000 1.00 10.0 ? 1 ALA A N 1
#
""", encoding="utf-8")
    return path


def _aggregation(root: Path, transform=None) -> Path:
    root.mkdir()
    transform = (np.asarray(transform, dtype=float) if transform is not None
                 else np.asarray([[1, 0, 0, 10], [0, 1, 0, 0],
                                  [0, 0, 1, 0], [0, 0, 0, 1]], dtype=float))
    task = {
        "docking_task_id": "DT-1", "query_id": "TEST:LIG:A:901",
        "receptor_id": "R1", "molecule_id": "M1", "conformer_id": "C1",
        "global_id": 0, "query_molecule_rank": 1,
        "objective": "atomcentered_anchored_joint", "search_score": 0.9,
        "candidate_to_query_transform": transform.reshape(-1).tolist(),
    }
    files = {
        "query_molecule_evidence": root / "query-molecule-evidence.jsonl",
        "molecule_summary": root / "molecule-summary.jsonl",
        "docking_admission": root / "docking-admission.jsonl",
        "docking_tasks": root / "docking-tasks.jsonl",
    }
    for name, path in files.items():
        row = task if name == "docking_tasks" else {"placeholder": name}
        path.write_text(json.dumps(row) + "\n", encoding="utf-8")
    (root / "manifest.json").write_text(json.dumps({
        "format": "aidd-multi-cocrystal-molecule-aggregation",
        "version": 1, "status": "complete", "config_hash": "abc",
        "outputs": {name: {"path": str(path), "sha256": _sha256(path)}
                    for name, path in files.items()},
    }), encoding="utf-8")
    return root


def _mock_structure(monkeypatch):
    monkeypatch.setattr(
        predocking_qc, "_query_instance",
        lambda path, query_id: (
            np.asarray([[10, 0, 0], [11, 0, 0]], dtype=float),
            "1", ("LIG", "A", "901")))
    monkeypatch.setattr(
        predocking_qc, "_protein_atoms",
        lambda path, model, excluded: [
            {"chain": "A", "xyz": np.asarray([11, .5, 0], dtype=float)},
            {"chain": "A", "xyz": np.asarray([11, 2.5, 0], dtype=float)},
            {"chain": "B", "xyz": np.asarray([10, 0, 0], dtype=float)},
        ])


def test_predocking_qc_exports_transformed_point_cloud_and_metrics(
        tmp_path: Path, monkeypatch):
    _mock_structure(monkeypatch)
    catalog = _artifact_catalog(tmp_path / "artifacts")
    aggregation = _aggregation(tmp_path / "aggregation")
    receptor = _receptor(tmp_path / "receptor.cif")
    output = tmp_path / "qc"
    first = run_predocking_pocket_qc(
        catalog, aggregation, output, receptors={"R1": receptor}, top_n_per_query=5)
    row = json.loads((output / "pose-qc.jsonl").read_text(encoding="utf-8"))
    assert first["status"] == "complete"
    assert first["counts"]["selected_by_query"] == {"TEST:LIG:A:901": 1}
    assert row["centroid_distance_angstrom"] == pytest.approx(0.0)
    assert row["candidate_points_near_query_fraction"] == 1.0
    assert row["query_points_covered_fraction"] == 1.0
    assert row["protein_chain"] == "A"
    assert row["minimum_protein_distance_angstrom"] == pytest.approx(0.5)
    assert row["severe_protein_points"] == 2
    pose = Path(row["point_cloud_pdb"])
    assert row["point_cloud_pdb_sha256"] == _sha256(pose)
    assert "GENERIC X ATOMS" in pose.read_text(encoding="utf-8")
    assert "10.000   0.000   0.000" in pose.read_text(encoding="utf-8")
    hashes = {name: value["sha256"] for name, value in first["outputs"].items()}
    manifest_hash = _sha256(output / "manifest.json")
    second = run_predocking_pocket_qc(
        catalog, aggregation, output, receptors={"R1": receptor}, top_n_per_query=5)
    assert hashes == {name: value["sha256"] for name, value in second["outputs"].items()}
    assert manifest_hash == _sha256(output / "manifest.json")


def test_predocking_qc_rejects_missing_receptor_and_non_affine_transform(
        tmp_path: Path, monkeypatch):
    _mock_structure(monkeypatch)
    catalog = _artifact_catalog(tmp_path / "artifacts")
    receptor = _receptor(tmp_path / "receptor.cif")
    aggregation = _aggregation(tmp_path / "aggregation")
    with pytest.raises(ValueError, match="missing receptor"):
        run_predocking_pocket_qc(
            catalog, aggregation, tmp_path / "missing", receptors={"OTHER": receptor})
    bad = np.eye(4); bad[3, 3] = 2
    aggregation = _aggregation(tmp_path / "bad-aggregation", bad)
    with pytest.raises(ValueError, match="affine"):
        run_predocking_pocket_qc(
            catalog, aggregation, tmp_path / "bad", receptors={"R1": receptor})
