from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import re
from typing import Mapping, Sequence

import numpy as np

from .anchor_extraction import _protein_atoms
from .chemistry_prep import enumerate_ligand_instances
from .gaussian_batch import ArtifactCatalogReader
from .gaussian_overlay import apply_transform
from .chemical_companion import ChemicalCompanionReader
from .chemical_geometry import vdw_exclusion_metrics


RESULT_FORMAT = "aidd-predocking-pocket-qc"
RESULT_VERSION = 1
ELEMENTS = {5: "B", 6: "C", 7: "N", 8: "O", 9: "F", 14: "Si",
            15: "P", 16: "S", 17: "Cl", 34: "Se", 35: "Br", 53: "I"}
ATOMIC_NUMBERS = {value.upper(): key for key, value in ELEMENTS.items()}


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _atomic_text(path: Path, value: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    partial = path.with_name(path.name + ".partial")
    with partial.open("w", encoding="utf-8", newline="\n") as stream:
        stream.write(value)
        stream.flush()
        os.fsync(stream.fileno())
    os.replace(partial, path)


def _atomic_json(path: Path, document: Mapping) -> None:
    _atomic_text(path, json.dumps(
        document, indent=2, sort_keys=True, ensure_ascii=False) + "\n")


def _atomic_jsonl(path: Path, rows: Sequence[Mapping]) -> None:
    _atomic_text(path, "".join(
        json.dumps(row, sort_keys=True, ensure_ascii=False,
                   separators=(",", ":")) + "\n" for row in rows))


def _jsonl(path: Path) -> list[dict]:
    return [json.loads(line) for line in path.read_text(
        encoding="utf-8").splitlines() if line.strip()]


def _safe_name(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", value).strip("._") or "item"


def _validate_aggregation(directory: Path) -> tuple[dict, list[dict]]:
    manifest_path = directory / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if (manifest.get("format") != "aidd-multi-cocrystal-molecule-aggregation"
            or manifest.get("status") != "complete"):
        raise ValueError("aggregation input is not a complete supported result")
    for output in manifest.get("outputs", {}).values():
        path = Path(output["path"])
        if not path.is_file():
            path = directory / path.name
        if _sha256(path) != output["sha256"]:
            raise ValueError(f"aggregation output checksum mismatch: {path}")
    tasks_path = directory / "docking-tasks.jsonl"
    if not tasks_path.is_file():
        raise FileNotFoundError(tasks_path)
    return manifest, _jsonl(tasks_path)


def _query_instance(path: Path, query_id: str) -> tuple[np.ndarray, str,
                                                         tuple[str, str, str]]:
    fields = query_id.split(":")
    if len(fields) != 4:
        raise ValueError(f"query ID must use PDB:CCD:chain:residue: {query_id}")
    _, ccd_id, chain_id, residue = fields
    matches = [row for row in enumerate_ligand_instances(path, [ccd_id])
               if row["chain_id"] == chain_id
               and row["residue_number"] == residue]
    if len(matches) != 1:
        raise ValueError(
            f"expected one {query_id} ligand instance in {path}, found {len(matches)}")
    instance = matches[0]
    xyz = np.asarray([[atom["x"], atom["y"], atom["z"]]
                      for atom in instance["atoms"]
                      if str(atom["element"]).upper() != "H"], dtype=np.float64)
    if xyz.ndim != 2 or xyz.shape[1:] != (3,) or not len(xyz):
        raise ValueError(f"query ligand has no heavy coordinates: {query_id}")
    return xyz, str(instance["model"]), (ccd_id, chain_id, residue)


def _distances(left: np.ndarray, right: np.ndarray) -> np.ndarray:
    if not len(left) or not len(right):
        raise ValueError("distance operands must not be empty")
    delta = left[:, None, :] - right[None, :, :]
    return np.sqrt(np.einsum("ijk,ijk->ij", delta, delta))


def _point_cloud_pdb(points: np.ndarray, molecule_id: str,
                     conformer_id: str, query_id: str) -> str:
    lines = [
        "REMARK 950 AIDD ARTIFACT V1 HEAVY-ATOM POINT CLOUD",
        "REMARK 950 GENERIC X ATOMS; NOT A CHEMICAL OR DOCKING STRUCTURE",
        f"REMARK 950 MOLECULE {molecule_id}",
        f"REMARK 950 CONFORMER {conformer_id}",
        f"REMARK 950 QUERY {query_id}",
    ]
    for serial, (x, y, z) in enumerate(points, start=1):
        atom_name = f"X{serial % 1000:03d}"
        lines.append(
            f"HETATM{serial:5d} {atom_name:<4s} PCD Z   1    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           X")
    lines.extend(("TER", "END"))
    return "\n".join(lines) + "\n"


def _chemical_sdf(points: np.ndarray, chemistry, query_id: str) -> str:
    if len(points) != len(chemistry.atomic_numbers):
        raise ValueError("chemical companion and coordinate atom counts differ")
    bonds = chemistry.bonds
    lines = [chemistry.molecule_id, "  AIDD E029", "",
             f"{len(points):3d}{len(bonds):3d}  0  0  0  0            3D V2000"]
    for xyz, atomic_number in zip(points, chemistry.atomic_numbers):
        symbol = ELEMENTS.get(int(atomic_number))
        if symbol is None:
            raise ValueError(f"unsupported SDF atomic number: {atomic_number}")
        lines.append(f"{xyz[0]:10.4f}{xyz[1]:10.4f}{xyz[2]:10.4f} {symbol:<3s}"
                     " 0  0  0  0  0  0  0  0  0  0  0  0")
    for bond in bonds:
        lines.append(f"{int(bond['begin'])+1:3d}{int(bond['end'])+1:3d}"
                     f"{int(bond['order']):3d}  0  0  0  0")
    charged = [(index + 1, int(charge)) for index, charge in
               enumerate(chemistry.formal_charges) if int(charge)]
    for start in range(0, len(charged), 8):
        group = charged[start:start + 8]
        lines.append(f"M  CHG{len(group):3d}" + "".join(
            f"{index:4d}{charge:4d}" for index, charge in group))
    lines.extend(("M  END", ">  <AIDD_GLOBAL_ID>", str(chemistry.global_id), "",
                  ">  <AIDD_CONFORMER_ID>", chemistry.conformer_id, "",
                  ">  <AIDD_QUERY_ID>", query_id, "", "$$$$"))
    return "\n".join(lines) + "\n"


def _parse_receptor_assignments(values: Sequence[str]) -> dict[str, Path]:
    result: dict[str, Path] = {}
    for value in values:
        receptor_id, separator, path = value.partition("=")
        receptor_id, path = receptor_id.strip(), path.strip()
        if not separator or not receptor_id or not path:
            raise ValueError("receptors must use RECEPTOR_ID=/path/structure.cif")
        if receptor_id in result:
            raise ValueError(f"duplicate receptor assignment: {receptor_id}")
        resolved = Path(path).resolve()
        if not resolved.is_file():
            raise FileNotFoundError(resolved)
        result[receptor_id] = resolved
    if not result:
        raise ValueError("at least one receptor assignment is required")
    return result


def run_predocking_pocket_qc(
        artifact_catalog: Path, aggregation_dir: Path, output_dir: Path, *,
        receptors: Mapping[str, Path] | Sequence[str],
        chemical_companion: Path | None = None,
        top_n_per_query: int = 100,
        query_neighborhood: float = 4.0,
        close_distance: float = 2.0,
        severe_distance: float = 1.5) -> dict:
    """Export transformed point clouds and annotate coordinate-only pocket QC."""
    if top_n_per_query <= 0:
        raise ValueError("top_n_per_query must be positive")
    if not (0 < severe_distance < close_distance < query_neighborhood):
        raise ValueError(
            "distances must satisfy 0 < severe < close < query neighborhood")
    artifact_catalog = artifact_catalog.resolve()
    aggregation_dir = aggregation_dir.resolve()
    output_dir = output_dir.resolve()
    receptor_paths = (_parse_receptor_assignments(receptors)
                      if not isinstance(receptors, Mapping)
                      else {str(key): Path(value).resolve()
                            for key, value in receptors.items()})
    for receptor_id, path in receptor_paths.items():
        if not receptor_id or not path.is_file():
            raise FileNotFoundError(path)

    aggregation, tasks = _validate_aggregation(aggregation_dir)
    expected_receptors = {str(row["receptor_id"]) for row in tasks}
    missing = sorted(expected_receptors - set(receptor_paths))
    if missing:
        raise ValueError("missing receptor assignments: " + ", ".join(missing))

    by_query: dict[str, list[dict]] = {}
    for task in tasks:
        by_query.setdefault(str(task["query_id"]), []).append(task)
    selected: list[dict] = []
    for query_id in sorted(by_query):
        ordered = sorted(by_query[query_id], key=lambda row: (
            int(row["query_molecule_rank"]), str(row["molecule_id"]),
            int(row["global_id"]), str(row["docking_task_id"])))
        selected.extend(ordered[:top_n_per_query])

    reader = ArtifactCatalogReader(artifact_catalog)
    chemistry_reader = None
    if chemical_companion is not None:
        chemistry_reader = ChemicalCompanionReader(Path(chemical_companion))
    receptor_cache: dict[tuple[str, str], tuple[np.ndarray, np.ndarray, np.ndarray]] = {}
    rows: list[dict] = []
    pml_lines = [
        "reinitialize",
        "set auto_zoom, off",
        "bg_color white",
    ]
    loaded_receptors: set[str] = set()
    for task in selected:
        query_id, receptor_id = str(task["query_id"]), str(task["receptor_id"])
        cache_key = (receptor_id, query_id)
        if cache_key not in receptor_cache:
            receptor_path = receptor_paths[receptor_id]
            query_xyz, model, excluded = _query_instance(receptor_path, query_id)
            protein = [atom for atom in _protein_atoms(receptor_path, model, excluded)
                       if atom["chain"] == excluded[1]]
            protein_xyz = np.asarray([atom["xyz"] for atom in protein], dtype=np.float64)
            protein_atomic_numbers = np.asarray([
                ATOMIC_NUMBERS.get(str(atom["element"]).upper(), 6) for atom in protein],
                dtype=np.uint8)
            if not len(protein_xyz):
                raise ValueError(f"receptor has no protein heavy atoms: {receptor_path}")
            receptor_cache[cache_key] = query_xyz, protein_xyz, protein_atomic_numbers
        query_xyz, protein_xyz, protein_atomic_numbers = receptor_cache[cache_key]

        candidate = reader.get(int(task["global_id"]))
        if (candidate.molecule_id != str(task["molecule_id"])
                or candidate.conformer_id != str(task["conformer_id"])):
            raise ValueError(
                f"artifact identity mismatch for global ID {task['global_id']}")
        transform = np.asarray(task["candidate_to_query_transform"], dtype=np.float64)
        if transform.size != 16:
            raise ValueError("candidate transform must contain 16 values")
        transform = transform.reshape(4, 4)
        if (not np.isfinite(transform).all()
                or not np.allclose(transform[3], [0, 0, 0, 1], atol=1e-7)):
            raise ValueError("candidate transform is not a finite affine 4x4 matrix")
        moved = apply_transform(candidate.shape_points, transform)
        chemistry = chemistry_reader.get(candidate.global_id) if chemistry_reader else None
        if chemistry is not None and (chemistry.molecule_id != candidate.molecule_id
                                      or chemistry.conformer_id != candidate.conformer_id):
            raise ValueError(f"chemical companion identity mismatch for {candidate.global_id}")
        query_distances = _distances(moved, query_xyz)
        protein_distances = _distances(moved, protein_xyz)
        candidate_query_nearest = query_distances.min(axis=1)
        query_candidate_nearest = query_distances.min(axis=0)
        candidate_protein_nearest = protein_distances.min(axis=1)
        close_count = int(np.count_nonzero(candidate_protein_nearest < close_distance))
        severe_count = int(np.count_nonzero(candidate_protein_nearest < severe_distance))

        query_rank = int(task["query_molecule_rank"])
        basename = _safe_name(
            f"{query_id}.rank-{query_rank:05d}.{candidate.molecule_id}.gid-{candidate.global_id}")
        pose_path = output_dir / "poses" / f"{basename}.pdb"
        _atomic_text(pose_path, _point_cloud_pdb(
            moved, candidate.molecule_id, candidate.conformer_id, query_id))
        sdf_path = None
        vdw = None
        if chemistry is not None:
            sdf_path = output_dir / "chemical-poses" / f"{basename}.sdf"
            _atomic_text(sdf_path, _chemical_sdf(moved, chemistry, query_id))
            vdw = vdw_exclusion_metrics(
                moved, chemistry.atomic_numbers, protein_xyz, protein_atomic_numbers)
        record = {
            "docking_task_id": str(task["docking_task_id"]),
            "query_id": query_id,
            "receptor_id": receptor_id,
            "molecule_id": candidate.molecule_id,
            "conformer_id": candidate.conformer_id,
            "global_id": candidate.global_id,
            "query_molecule_rank": query_rank,
            "objective": str(task["objective"]),
            "search_score": float(task["search_score"]),
            "candidate_heavy_points": len(moved),
            "query_heavy_atoms": len(query_xyz),
            "protein_chain": query_id.split(":")[2],
            "centroid_distance_angstrom": float(np.linalg.norm(
                moved.mean(axis=0) - query_xyz.mean(axis=0))),
            "candidate_points_near_query_fraction": float(np.mean(
                candidate_query_nearest <= query_neighborhood)),
            "query_points_covered_fraction": float(np.mean(
                query_candidate_nearest <= query_neighborhood)),
            "minimum_protein_distance_angstrom": float(candidate_protein_nearest.min()),
            "close_protein_points": close_count,
            "close_protein_fraction": close_count / len(moved),
            "severe_protein_points": severe_count,
            "severe_protein_fraction": severe_count / len(moved),
            "point_cloud_pdb": str(pose_path),
            "point_cloud_pdb_sha256": _sha256(pose_path),
            "chemical_sdf": str(sdf_path) if sdf_path else None,
            "chemical_sdf_sha256": _sha256(sdf_path) if sdf_path else None,
            "vdw_exclusion": vdw,
            "representation": "artifact-v1-heavy-atom-centers-generic-element",
        }
        rows.append(record)

        receptor_object = _safe_name(receptor_id)
        pose_object = _safe_name(f"pose_{query_id}_{query_rank:05d}")
        if receptor_id not in loaded_receptors:
            pml_lines.extend((
                f"load {receptor_paths[receptor_id].as_posix()}, {receptor_object}",
                f"hide everything, {receptor_object}",
                f"show cartoon, {receptor_object} and polymer.protein",
                f"color gray80, {receptor_object} and polymer.protein",
            ))
            loaded_receptors.add(receptor_id)
        pml_lines.extend((
            f"load {pose_path.as_posix()}, {pose_object}",
            f"show spheres, {pose_object}",
            f"set sphere_scale, 0.22, {pose_object}",
            f"color magenta, {pose_object}",
            f"disable {pose_object}",
        ))

    rows.sort(key=lambda row: (
        row["query_id"], row["query_molecule_rank"], row["molecule_id"],
        row["global_id"]))
    qc_path = output_dir / "pose-qc.jsonl"
    _atomic_jsonl(qc_path, rows)
    if rows:
        first_object = _safe_name(
            f"pose_{rows[0]['query_id']}_{rows[0]['query_molecule_rank']:05d}")
        pml_lines.extend((f"enable {first_object}", f"zoom {first_object}, 8"))
    pml_path = output_dir / "review.pml"
    _atomic_text(pml_path, "\n".join(pml_lines) + "\n")

    receptor_inputs = {key: {"path": str(path), "sha256": _sha256(path)}
                       for key, path in sorted(receptor_paths.items())}
    manifest = {
        "format": RESULT_FORMAT,
        "version": RESULT_VERSION,
        "status": "complete",
        "artifact_catalog": {"path": str(artifact_catalog),
                             "sha256": _sha256(artifact_catalog)},
        "chemical_companion": ({"path": str(Path(chemical_companion).resolve()),
                                "sha256": _sha256(Path(chemical_companion).resolve())}
                               if chemical_companion is not None else None),
        "aggregation": {"path": str(aggregation_dir / "manifest.json"),
                        "sha256": _sha256(aggregation_dir / "manifest.json"),
                        "config_hash": aggregation.get("config_hash")},
        "receptors": receptor_inputs,
        "configuration": {
            "top_n_per_query": top_n_per_query,
            "query_neighborhood_angstrom": query_neighborhood,
            "close_distance_angstrom": close_distance,
            "severe_distance_angstrom": severe_distance,
        },
        "counts": {
            "queries": len(by_query),
            "selected_tasks": len(rows),
            "selected_by_query": {
                query_id: sum(row["query_id"] == query_id for row in rows)
                for query_id in sorted(by_query)},
            "poses_with_close_points": sum(row["close_protein_points"] > 0 for row in rows),
            "poses_with_severe_points": sum(row["severe_protein_points"] > 0 for row in rows),
        },
        "outputs": {
            "pose_qc": {"path": str(qc_path), "sha256": _sha256(qc_path)},
            "pymol_review": {"path": str(pml_path), "sha256": _sha256(pml_path)},
        },
        "invariants": {
            "retrieval_membership_changed": False,
            "admission_order_changed": False,
            "coordinate_source": "immutable-artifact-v1",
            "stored_transform_applied_without_optimization": True,
            "chemical_topology_available": chemistry_reader is not None,
            "chemical_sdf_available": chemistry_reader is not None,
            "point_clouds_are_docking_inputs": False,
        },
        "interpretation_boundary": (
            "Coordinate-only geometry diagnostics; not a chemical structure, "
            "retrieval filter, docking score, binding prediction, or activity claim."),
    }
    _atomic_json(output_dir / "manifest.json", manifest)
    return manifest
