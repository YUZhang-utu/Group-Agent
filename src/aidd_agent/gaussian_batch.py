from __future__ import annotations

from dataclasses import dataclass
from concurrent.futures import ProcessPoolExecutor, as_completed
import hashlib
import heapq
import json
import os
from pathlib import Path
import sys
import time
from typing import Mapping, Sequence

import numpy as np

from .conformer_artifacts import FEATURE_DTYPE, FEATURE_TYPES, META_DTYPE
from .gaussian_overlay import (
    apply_transform, gaussian_overlap, gaussian_self_overlap, pair_alignment_seeds,
    principal_axis_seeds, query_biased_tversky, tanimoto,
)


QUERY_FORMAT = "aidd-gaussian-query"
RESULT_FORMAT = "aidd-gaussian-batch-scores"
SCHEMA_VERSION = 1
STAGED_FORMAT = "aidd-staged-gaussian-reranking"
STAGED_VERSION = 1
SCALED_FORMAT = "aidd-scaled-gaussian-reranking"
SCALED_VERSION = 1
OBJECTIVES = (
    "shape_only", "atomcentered_unweighted_joint", "atomcentered_anchored_joint")
ANCHOR_TO_ARTIFACT_TYPE = {"HBD": FEATURE_TYPES["Donor"],
                           "HBA": FEATURE_TYPES["Acceptor"]}


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


@dataclass(frozen=True)
class ArtifactConformer:
    global_id: int
    conformer_id: str
    molecule_id: str
    shape_points: np.ndarray
    feature_points: np.ndarray
    feature_types: np.ndarray


class ArtifactCatalogReader:
    """Resolve stable IDs and expose conformer arrays from read-only shard maps."""

    def __init__(self, catalog_path: Path):
        self.catalog_path = catalog_path.resolve()
        self.catalog = json.loads(self.catalog_path.read_text(encoding="utf-8"))
        if self.catalog.get("format") != "aidd-conformer-artifact-catalog":
            raise ValueError("not an AIDD conformer artifact catalog")
        self._shards = []
        for record in self.catalog["shards"]:
            directory = Path(record["path"])
            if not directory.is_dir():
                directory = self.catalog_path.parent / record["name"]
            manifest_path = directory / "manifest.json"
            manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
            if manifest.get("format") != "aidd-conformer-artifact-shard":
                raise ValueError(f"not an artifact shard: {directory}")
            start, count = int(record["global_id_start"]), int(record["conformers"])
            if (int(manifest["global_id_start"]), int(manifest["conformers"])) != (start, count):
                raise ValueError(f"catalog/shard range mismatch: {directory}")
            scale = float(manifest["coordinate_scale"])
            if not np.isfinite(scale) or scale <= 0:
                raise ValueError(f"invalid coordinate scale: {directory}")
            meta = np.memmap(directory / "meta.bin", dtype=META_DTYPE, mode="r")
            coords = np.memmap(directory / "coords.bin", dtype="<i2", mode="r").reshape(-1, 3)
            features = np.memmap(directory / "feats.bin", dtype=FEATURE_DTYPE, mode="r")
            conformer_ids = np.memmap(directory / "conformer_ids.bin", dtype="S16", mode="r")
            molecule_ids = np.memmap(directory / "molecule_ids.bin", dtype="S16", mode="r")
            if not (len(meta) == len(conformer_ids) == len(molecule_ids) == count):
                raise ValueError(f"artifact row-count mismatch: {directory}")
            self._shards.append((start, start + count, scale, meta, coords, features,
                                 conformer_ids, molecule_ids))
        self._shards.sort(key=lambda item: item[0])

    def get(self, global_id: int) -> ArtifactConformer:
        if isinstance(global_id, (bool, np.bool_)) or int(global_id) != global_id:
            raise ValueError("global ID must be an integer")
        gid = int(global_id)
        for start, stop, scale, meta, coords, features, conformer_ids, molecule_ids in self._shards:
            if start <= gid < stop:
                local = gid - start
                row = meta[local]
                if int(row["global_id"]) != gid:
                    raise ValueError(f"artifact global ID mismatch at {gid}")
                co, cn = int(row["coord_offset"]), int(row["heavy_atoms"])
                fo, fn = int(row["feature_offset"]), int(row["features"])
                if co + cn > len(coords) or fo + fn > len(features):
                    raise ValueError(f"artifact offset out of bounds at {gid}")
                origin = np.asarray(row["origin"], dtype=np.float64)
                shape = origin + np.asarray(coords[co:co + cn], dtype=np.float64) / scale
                feature_rows = features[fo:fo + fn]
                feature_points = origin + np.asarray(
                    feature_rows["xyz"], dtype=np.float64) / scale
                decode = lambda value: bytes(value).split(b"\0", 1)[0].decode("utf-8")
                return ArtifactConformer(
                    gid, decode(conformer_ids[local]), decode(molecule_ids[local]), shape,
                    feature_points, np.asarray(feature_rows["type"], dtype=np.uint8))
        raise IndexError(f"global ID outside artifact catalog: {gid}")


def write_gaussian_query(output_path: Path, *, shape_points: Sequence[Sequence[float]],
                         feature_points: Sequence[Sequence[float]],
                         feature_types: Sequence[int],
                         anchored_weights: Sequence[float],
                         anchor_feature_indices: Sequence[int],
                         source: Mapping | None = None) -> dict:
    output_path = output_path.resolve()
    if output_path.suffix.lower() != ".npz":
        raise ValueError("Gaussian query output must use .npz")
    shape = np.asarray(shape_points, dtype=np.float64)
    points = np.asarray(feature_points, dtype=np.float64)
    types = np.asarray(feature_types, dtype=np.uint8)
    weights = np.asarray(anchored_weights, dtype=np.float64)
    anchors = np.asarray(anchor_feature_indices, dtype=np.int64)
    if shape.ndim != 2 or shape.shape[1:] != (3,) or not len(shape):
        raise ValueError("query shape points require non-empty shape (N,3)")
    if points.ndim != 2 or points.shape[1:] != (3,) or not len(points):
        raise ValueError("query features require non-empty shape (N,3)")
    if types.shape != (len(points),) or weights.shape != (len(points),):
        raise ValueError("query feature arrays have inconsistent lengths")
    if not np.isfinite(shape).all() or not np.isfinite(points).all():
        raise ValueError("query coordinates must be finite")
    if not np.isfinite(weights).all() or np.any(weights < 0):
        raise ValueError("query weights must be finite and non-negative")
    if anchors.ndim != 1 or np.any(anchors < 0) or np.any(anchors >= len(points)):
        raise ValueError("anchor feature index is out of range")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(output_path, shape_points=shape, feature_points=points,
                        feature_types=types, anchored_weights=weights,
                        anchor_feature_indices=np.unique(anchors))
    manifest = {
        "format": QUERY_FORMAT, "version": SCHEMA_VERSION,
        "query_npz": output_path.name, "query_npz_sha256": _sha256(output_path),
        "shape_points": len(shape), "feature_points": len(points),
        "anchor_feature_indices": np.unique(anchors).tolist(),
        "color_model": "atom-centered-feature-v1",
        "projected_color_available": False,
        "source": dict(source or {}),
    }
    manifest_path = output_path.with_suffix(".manifest.json")
    manifest_path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    return manifest


def prepare_gaussian_query(mmcif: Path, ccd: Path, query_manifest: Path,
                           output_path: Path) -> dict:
    """Build a query package from one locked co-crystal ligand instance."""
    from rdkit import RDConfig
    from rdkit.Chem import ChemicalFeatures
    from .chemistry_prep import _ccd_molecule, enumerate_ligand_instances

    document = json.loads(query_manifest.read_text(encoding="utf-8"))
    fields = str(document["query_id"]).split(":")
    if len(fields) != 4:
        raise ValueError("query_id must use PDB:CCD:chain:residue syntax")
    _, ccd_id, chain_id, residue = fields
    instances = [item for item in enumerate_ligand_instances(mmcif, [ccd_id])
                 if item["chain_id"] == chain_id and item["residue_number"] == residue]
    if len(instances) != 1:
        raise ValueError("query manifest must resolve to exactly one ligand instance")
    molecule = _ccd_molecule(ccd, instances[0]["atoms"])
    conformer = molecule.GetConformer()
    heavy = [atom.GetIdx() for atom in molecule.GetAtoms() if atom.GetAtomicNum() > 1]
    shape = np.asarray([[conformer.GetAtomPosition(i).x, conformer.GetAtomPosition(i).y,
                         conformer.GetAtomPosition(i).z] for i in heavy])
    factory = ChemicalFeatures.BuildFeatureFactory(str(Path(RDConfig.RDDataDir) / "BaseFeatures.fdef"))
    feature_rows = []
    for feature in factory.GetFeaturesForMol(molecule):
        type_id = FEATURE_TYPES.get(feature.GetFamily())
        if type_id is None:
            continue
        atom_ids = tuple(feature.GetAtomIds())
        xyz = np.asarray([[conformer.GetAtomPosition(i).x, conformer.GetAtomPosition(i).y,
                           conformer.GetAtomPosition(i).z] for i in atom_ids]).mean(axis=0)
        feature_rows.append((xyz, type_id, atom_ids))
    if not feature_rows:
        raise ValueError("query ligand has no supported atom-centered features")
    points = np.asarray([row[0] for row in feature_rows])
    types = np.asarray([row[1] for row in feature_rows], dtype=np.uint8)
    ordinary = float(document.get("weights", {}).get("ordinary", 1.0))
    weights = np.full(len(points), ordinary, dtype=np.float64)
    anchor_indices = []
    for anchor in document.get("anchors", []):
        expected = ANCHOR_TO_ARTIFACT_TYPE.get(anchor.get("feature_type"))
        atom_ids = set(map(int, anchor.get("ligand_atom_indices", [])))
        compatible = [i for i, row in enumerate(feature_rows)
                      if row[1] == expected and atom_ids.intersection(row[2])]
        if not compatible:
            raise ValueError(f"anchor does not map to a query feature: {anchor.get('anchor_id')}")
        best = min(compatible, key=lambda i: (np.linalg.norm(
            points[i] - np.asarray(anchor["atom_center"], dtype=np.float64)), i))
        weights[best] = max(weights[best], float(anchor.get("weight", ordinary)))
        anchor_indices.append(best)
    return write_gaussian_query(
        output_path, shape_points=shape, feature_points=points, feature_types=types,
        anchored_weights=weights, anchor_feature_indices=anchor_indices,
        source={"mmcif": str(mmcif.resolve()), "mmcif_sha256": _sha256(mmcif),
                "ccd": str(ccd.resolve()), "ccd_sha256": _sha256(ccd),
                "query_manifest": str(query_manifest.resolve()),
                "query_manifest_sha256": _sha256(query_manifest),
                "query_id": document["query_id"]})


def _load_query(path: Path) -> tuple[dict, dict[str, np.ndarray]]:
    manifest_path = path.with_suffix(".manifest.json")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if manifest.get("format") != QUERY_FORMAT or manifest.get("version") != SCHEMA_VERSION:
        raise ValueError("not a supported Gaussian query package")
    if manifest.get("query_npz_sha256") != _sha256(path):
        raise ValueError("Gaussian query checksum mismatch")
    with np.load(path, allow_pickle=False) as archive:
        arrays = {name: archive[name] for name in archive.files}
    required = {"shape_points", "feature_points", "feature_types",
                "anchored_weights", "anchor_feature_indices"}
    if not required.issubset(arrays):
        raise ValueError("Gaussian query package is missing required arrays")
    return manifest, arrays


def load_candidate_ids(path: Path) -> np.ndarray:
    suffix = path.suffix.lower()
    if suffix == ".npy":
        payload = np.load(path, allow_pickle=False)
    elif suffix == ".npz":
        with np.load(path, allow_pickle=False) as archive:
            key = next((name for name in ("global_ids", "ids") if name in archive), None)
            if key is None:
                raise ValueError("candidate NPZ needs global_ids or ids")
            payload = archive[key]
    elif suffix == ".json":
        payload = json.loads(path.read_text(encoding="utf-8"))
        if isinstance(payload, dict):
            payload = payload.get("global_ids", payload.get("ids"))
    else:
        try:
            payload = [int(line.strip()) for line in path.read_text(encoding="utf-8").splitlines()
                       if line.strip()]
        except ValueError as exc:
            raise ValueError("text candidate IDs must contain one integer per line") from exc
    values = np.asarray(payload)
    if values.ndim != 1 or values.dtype.kind not in "iu" or np.any(values < 0):
        raise ValueError("candidate IDs must be a one-dimensional non-negative integer array")
    result = values.astype(np.int64)
    if len(np.unique(result)) != len(result):
        raise ValueError("candidate IDs must not contain duplicates")
    if not len(result):
        raise ValueError("candidate ID input must not be empty")
    return result


def _objective(score: Mapping, name: str) -> float:
    shape = float(score["shape"]["query_biased_tversky"])
    if name == "shape_only":
        return shape
    return 0.5 * (shape + float(score["color"]["query_biased_tversky"]))


def _overlap_record(cross: float, query_self: float, candidate_self: float) -> dict:
    return {
        "cross_overlap_raw": cross,
        "query_self_overlap_raw": query_self,
        "candidate_self_overlap_raw": candidate_self,
        "tanimoto": tanimoto(cross, query_self, candidate_self),
        "query_biased_tversky": query_biased_tversky(
            cross, query_self, candidate_self),
    }


def _query_self_overlaps(query: Mapping[str, np.ndarray], sigma: float,
                         cutoff: float | None) -> dict[str, float]:
    return {
        "shape": gaussian_self_overlap(
            query["shape_points"], sigma=sigma, cutoff=cutoff),
        "color_unweighted": gaussian_self_overlap(
            query["feature_points"], sigma=sigma, cutoff=cutoff,
            types=query["feature_types"]),
        "color_anchored": gaussian_self_overlap(
            query["feature_points"], sigma=sigma, cutoff=cutoff,
            weights=query["anchored_weights"], types=query["feature_types"]),
    }


def _score_ids(reader: ArtifactCatalogReader, query: Mapping[str, np.ndarray],
               ids: np.ndarray, *, sigma: float, cutoff: float | None,
               pair_tolerance: float, axial_samples: int,
               max_pair_seeds: int) -> dict[str, np.ndarray]:
    """Score a locked ID slice while caching every rigid-invariant self overlap."""
    query_shape = query["shape_points"]
    query_features = query["feature_points"]
    query_types = query["feature_types"]
    query_weights = query["anchored_weights"]
    anchor_indices = query["anchor_feature_indices"]
    query_self = _query_self_overlaps(query, sigma, cutoff)
    seed_ids, molecule_ids, conformer_ids = [], [], []
    seed_counts, pair_seed_counts = [], []
    values = {name: {"objective": [], "shape_tanimoto": [], "shape_tversky": [],
                     "color_tanimoto": [], "color_tversky": [], "shape_cross": [],
                     "shape_query_self": [], "shape_candidate_self": [], "color_cross": [],
                     "color_query_self": [], "color_candidate_self": [], "transform": []}
              for name in OBJECTIVES}
    for gid in ids:
        candidate = reader.get(int(gid))
        molecule_ids.append(candidate.molecule_id)
        conformer_ids.append(candidate.conformer_id)
        seeds = list(principal_axis_seeds(candidate.shape_points, query_shape))
        pair_count = 0
        if (max_pair_seeds and len(anchor_indices) >= 2
                and len(candidate.feature_points) >= 2):
            pair_seeds = pair_alignment_seeds(
                candidate.feature_points, candidate.feature_types,
                query_features[anchor_indices], query_types[anchor_indices],
                tolerance=pair_tolerance, axial_samples=axial_samples)
            pair_seeds = pair_seeds[:max_pair_seeds]
            pair_count = len(pair_seeds)
            seeds.extend(pair_seeds)
        seed_counts.append(len(seeds)); pair_seed_counts.append(pair_count)
        candidate_self = {
            "shape": gaussian_self_overlap(
                candidate.shape_points, sigma=sigma, cutoff=cutoff),
            "color": gaussian_self_overlap(
                candidate.feature_points, sigma=sigma, cutoff=cutoff,
                types=candidate.feature_types),
        }
        best = {name: None for name in OBJECTIVES}
        for seed_index, seed in enumerate(seeds):
            moved_shape = apply_transform(candidate.shape_points, seed.transform_matrix)
            moved_features = apply_transform(candidate.feature_points, seed.transform_matrix)
            shape_cross = gaussian_overlap(
                query_shape, moved_shape, sigma=sigma, cutoff=cutoff)
            color_cross = gaussian_overlap(
                query_features, moved_features, sigma=sigma, cutoff=cutoff,
                types_a=query_types, types_b=candidate.feature_types)
            anchored_cross = gaussian_overlap(
                query_features, moved_features, sigma=sigma, cutoff=cutoff,
                weights_a=query_weights, types_a=query_types,
                types_b=candidate.feature_types)
            shape_record = _overlap_record(
                shape_cross, query_self["shape"], candidate_self["shape"])
            unweighted = {
                "shape": shape_record,
                "color": _overlap_record(
                    color_cross, query_self["color_unweighted"],
                    candidate_self["color"]),
                "transform_matrix": list(seed.transform_matrix),
            }
            anchored = {
                "shape": shape_record,
                "color": _overlap_record(
                    anchored_cross, query_self["color_anchored"],
                    candidate_self["color"]),
                "transform_matrix": list(seed.transform_matrix),
            }
            scored = {"shape_only": unweighted,
                      "atomcentered_unweighted_joint": unweighted,
                      "atomcentered_anchored_joint": anchored}
            for name, score in scored.items():
                objective = _objective(score, name)
                if best[name] is None or objective > best[name][0]:
                    best[name] = (objective, seed_index, seed.seed_id, score)
        seed_ids.append([best[name][2] for name in OBJECTIVES])
        for name in OBJECTIVES:
            objective, _, _, score = best[name]
            shape, color = score["shape"], score["color"]
            target = values[name]
            target["objective"].append(objective)
            for prefix, record in (("shape", shape), ("color", color)):
                target[f"{prefix}_tanimoto"].append(record["tanimoto"])
                target[f"{prefix}_tversky"].append(record["query_biased_tversky"])
                target[f"{prefix}_cross"].append(record["cross_overlap_raw"])
                target[f"{prefix}_query_self"].append(record["query_self_overlap_raw"])
                target[f"{prefix}_candidate_self"].append(record["candidate_self_overlap_raw"])
            target["transform"].append(score["transform_matrix"])
    arrays = {
        "global_ids": np.asarray(ids, dtype=np.int64),
        "molecule_ids": np.asarray(molecule_ids, dtype="U16"),
        "conformer_ids": np.asarray(conformer_ids, dtype="U16"),
        "best_seed_ids": np.asarray(seed_ids, dtype="U96"),
        "seed_counts": np.asarray(seed_counts, dtype=np.int32),
        "pair_seed_counts": np.asarray(pair_seed_counts, dtype=np.int32),
        "objective_names": np.asarray(OBJECTIVES, dtype="U40"),
    }
    for name, records in values.items():
        for metric, payload in records.items():
            arrays[f"{name}__{metric}"] = np.asarray(payload, dtype=np.float64)
    return arrays


def _atomic_savez(path: Path, arrays: Mapping[str, np.ndarray]) -> None:
    temporary = path.with_name(path.name + ".partial")
    with temporary.open("wb") as stream:
        np.savez_compressed(stream, **arrays)
        stream.flush()
        os.fsync(stream.fileno())
    os.replace(temporary, path)


def _atomic_json(path: Path, document: Mapping) -> None:
    temporary = path.with_name(path.name + ".partial")
    with temporary.open("w", encoding="utf-8") as stream:
        json.dump(document, stream, indent=2)
        stream.flush()
        os.fsync(stream.fileno())
    os.replace(temporary, path)


def _result_manifest(artifact_catalog: Path, query_path: Path,
                     candidate_ids_path: Path, output_path: Path, ids: np.ndarray,
                     *, sigma: float, cutoff: float | None, pair_tolerance: float,
                     axial_samples: int, max_pair_seeds: int) -> dict:
    return {
        "format": RESULT_FORMAT, "version": SCHEMA_VERSION,
        "result_npz": output_path.name, "result_npz_sha256": _sha256(output_path),
        "artifact_catalog": str(artifact_catalog.resolve()),
        "artifact_catalog_sha256": _sha256(artifact_catalog.resolve()),
        "query": str(query_path.resolve()), "query_sha256": _sha256(query_path.resolve()),
        "candidate_ids": str(candidate_ids_path.resolve()),
        "candidate_ids_sha256": _sha256(candidate_ids_path.resolve()),
        "candidates_input": len(ids), "candidates_output": len(ids),
        "operand_direction": {"A": "query", "B": "candidate"},
        "objectives": {
            "shape_only": "shape query-biased Tversky",
            "atomcentered_unweighted_joint":
                "mean(shape query-biased Tversky, unit-weight atom-centered color query-biased Tversky)",
            "atomcentered_anchored_joint":
                "mean(shape query-biased Tversky, anchor-weighted atom-centered color query-biased Tversky)",
        },
        "parameters": {"sigma_angstrom": sigma, "cutoff_angstrom": cutoff,
                       "pair_tolerance_angstrom": pair_tolerance,
                       "axial_samples": axial_samples, "max_pair_seeds": max_pair_seeds,
                       "tversky_alpha_query": 0.95, "tversky_beta_candidate": 0.05},
        "limitations": ["artifact-v1 has atom-centered features only",
                        "projected pharmacophore color was not scored",
                        "scores rerank and never alter L1 admission"],
    }


def score_gaussian_candidates(artifact_catalog: Path, query_path: Path,
                              candidate_ids_path: Path, output_path: Path, *,
                              sigma: float = 1.0, cutoff: float | None = 4.5,
                              pair_tolerance: float = 2.0, axial_samples: int = 6,
                              max_pair_seeds: int = 512) -> dict:
    if max_pair_seeds < 0:
        raise ValueError("max_pair_seeds must be non-negative")
    reader = ArtifactCatalogReader(artifact_catalog)
    _, query = _load_query(query_path.resolve())
    ids = load_candidate_ids(candidate_ids_path.resolve())
    arrays = _score_ids(
        reader, query, ids, sigma=sigma, cutoff=cutoff,
        pair_tolerance=pair_tolerance, axial_samples=axial_samples,
        max_pair_seeds=max_pair_seeds)
    output_path = output_path.resolve()
    if output_path.suffix.lower() != ".npz":
        raise ValueError("Gaussian score output must use .npz")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    _atomic_savez(output_path, arrays)
    manifest = _result_manifest(
        artifact_catalog, query_path, candidate_ids_path, output_path, ids,
        sigma=sigma, cutoff=cutoff, pair_tolerance=pair_tolerance,
        axial_samples=axial_samples, max_pair_seeds=max_pair_seeds)
    _atomic_json(output_path.with_suffix(".manifest.json"), manifest)
    return manifest


_WORKER_READER: ArtifactCatalogReader | None = None
_WORKER_QUERY: dict[str, np.ndarray] | None = None
_WORKER_PARAMETERS: dict | None = None


def _worker_initialize(artifact_catalog: str, query_path: str, parameters: dict) -> None:
    global _WORKER_READER, _WORKER_QUERY, _WORKER_PARAMETERS
    _WORKER_READER = ArtifactCatalogReader(Path(artifact_catalog))
    _, _WORKER_QUERY = _load_query(Path(query_path))
    _WORKER_PARAMETERS = parameters


def _ids_sha256(ids: np.ndarray) -> str:
    values = np.ascontiguousarray(ids, dtype="<i8")
    return hashlib.sha256(values.tobytes()).hexdigest()


def _chunk_paths(stage_directory: Path, chunk_index: int) -> tuple[Path, Path]:
    stem = f"chunk-{chunk_index:06d}"
    return stage_directory / "chunks" / f"{stem}.npz", \
        stage_directory / "chunks" / f"{stem}.manifest.json"


def _write_chunk(stage_directory: Path, stage: str, chunk_index: int,
                 start: int, stop: int, ids: np.ndarray, arrays: Mapping[str, np.ndarray],
                 config_hash: str, seconds: float, storage: str = "detailed") -> dict:
    result_path, manifest_path = _chunk_paths(stage_directory, chunk_index)
    _atomic_savez(result_path, arrays)
    manifest = {
        "format": "aidd-gaussian-score-chunk", "version": 1,
        "stage": stage, "chunk_index": chunk_index, "start": start, "stop": stop,
        "candidates": len(ids), "ids_sha256": _ids_sha256(ids),
        "config_hash": config_hash, "result": result_path.name,
        "result_sha256": _sha256(result_path), "seconds": seconds,
        "storage": storage,
    }
    _atomic_json(manifest_path, manifest)
    return manifest


def _validate_chunk(stage_directory: Path, stage: str, chunk_index: int,
                    start: int, stop: int, ids: np.ndarray,
                    config_hash: str, storage: str = "detailed") -> dict | None:
    result_path, manifest_path = _chunk_paths(stage_directory, chunk_index)
    if not result_path.is_file() or not manifest_path.is_file():
        return None
    try:
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        expected = {
            "format": "aidd-gaussian-score-chunk", "version": 1,
            "stage": stage, "chunk_index": chunk_index, "start": start,
            "stop": stop, "candidates": len(ids), "ids_sha256": _ids_sha256(ids),
            "config_hash": config_hash,
        }
        if any(manifest.get(key) != value for key, value in expected.items()):
            return None
        if manifest.get("storage", "detailed") != storage:
            return None
        if manifest.get("result_sha256") != _sha256(result_path):
            return None
        with np.load(result_path, allow_pickle=False) as archive:
            if not np.array_equal(archive["global_ids"], ids):
                return None
            if tuple(map(str, archive["objective_names"])) != OBJECTIVES:
                return None
        return manifest
    except (OSError, ValueError, KeyError, json.JSONDecodeError):
        return None


def _worker_run_chunk(task: tuple) -> dict:
    if _WORKER_READER is None or _WORKER_QUERY is None or _WORKER_PARAMETERS is None:
        raise RuntimeError("Gaussian worker was not initialized")
    stage_directory, stage, chunk_index, start, stop, ids, config_hash, storage = task
    started = time.perf_counter()
    arrays = _score_ids(
        _WORKER_READER, _WORKER_QUERY, ids, **_WORKER_PARAMETERS)
    if storage == "slim":
        arrays = {
            "global_ids": arrays["global_ids"],
            "objective_names": arrays["objective_names"],
            "objective_scores": np.column_stack([
                arrays[f"{name}__objective"] for name in OBJECTIVES
            ]).astype(np.float32),
        }
    return _write_chunk(
        Path(stage_directory), stage, chunk_index, start, stop, ids, arrays,
        config_hash, time.perf_counter() - started, storage)


def _merge_chunks(stage_directory: Path, stage: str, ids: np.ndarray,
                  chunk_size: int, config_hash: str) -> tuple[Path, dict]:
    pieces: dict[str, list[np.ndarray]] = {}
    objective_names = None
    chunk_count = (len(ids) + chunk_size - 1) // chunk_size
    for chunk_index in range(chunk_count):
        result_path, _ = _chunk_paths(stage_directory, chunk_index)
        with np.load(result_path, allow_pickle=False) as archive:
            current_names = archive["objective_names"].copy()
            if objective_names is None:
                objective_names = current_names
            elif not np.array_equal(objective_names, current_names):
                raise ValueError("objective names differ between Gaussian chunks")
            for name in archive.files:
                if name != "objective_names":
                    pieces.setdefault(name, []).append(archive[name].copy())
    arrays = {name: np.concatenate(records, axis=0) for name, records in pieces.items()}
    arrays["objective_names"] = objective_names
    if not np.array_equal(arrays["global_ids"], ids):
        raise ValueError(f"{stage} merged IDs do not reproduce locked input order")
    output_path = stage_directory / "merged-scores.npz"
    _atomic_savez(output_path, arrays)
    manifest = {
        "format": "aidd-gaussian-merged-stage", "version": 1,
        "stage": stage, "config_hash": config_hash, "chunks": chunk_count,
        "candidates": len(ids), "ids_sha256": _ids_sha256(ids),
        "result": output_path.name, "result_sha256": _sha256(output_path),
    }
    _atomic_json(output_path.with_suffix(".manifest.json"), manifest)
    return output_path, manifest


def _execute_stage(stage_directory: Path, stage: str, ids: np.ndarray, *,
                   artifact_catalog: Path, query_path: Path, parameters: dict,
                   config_hash: str, workers: int, chunk_size: int,
                   resume: bool, progress_every: int, storage: str = "detailed",
                   merge: bool = True) -> tuple[Path, dict]:
    if storage not in {"detailed", "slim"}:
        raise ValueError("storage must be detailed or slim")
    stage_directory.mkdir(parents=True, exist_ok=True)
    (stage_directory / "chunks").mkdir(exist_ok=True)
    chunks = []
    reused = 0
    for chunk_index, start in enumerate(range(0, len(ids), chunk_size)):
        stop = min(len(ids), start + chunk_size)
        chunk_ids = ids[start:stop]
        valid = _validate_chunk(
            stage_directory, stage, chunk_index, start, stop, chunk_ids,
            config_hash, storage)
        if resume and valid is not None:
            reused += 1
        else:
            chunks.append((str(stage_directory), stage, chunk_index, start, stop,
                           chunk_ids, config_hash, storage))
    total_chunks = reused + len(chunks)
    finished = reused
    started = time.perf_counter()

    def report(manifest: Mapping) -> None:
        nonlocal finished
        finished += 1
        if progress_every and (finished % progress_every == 0 or finished == total_chunks):
            elapsed = time.perf_counter() - started
            newly_done = max(1, finished - reused)
            remaining = max(0, total_chunks - finished)
            eta = elapsed / newly_done * remaining
            print(f"[{stage}] chunks {finished}/{total_chunks}; "
                  f"last={manifest['candidates']} in {manifest['seconds']:.2f}s; "
                  f"resume_reused={reused}; ETA={eta/60:.1f} min",
                  file=sys.stderr, flush=True)

    initializer_args = (str(artifact_catalog), str(query_path), parameters)
    if workers == 1:
        _worker_initialize(*initializer_args)
        for task in chunks:
            report(_worker_run_chunk(task))
    elif chunks:
        for variable in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
            os.environ.setdefault(variable, "1")
        with ProcessPoolExecutor(
                max_workers=workers, initializer=_worker_initialize,
                initargs=initializer_args) as executor:
            futures = [executor.submit(_worker_run_chunk, task) for task in chunks]
            for future in as_completed(futures):
                report(future.result())
    for chunk_index, start in enumerate(range(0, len(ids), chunk_size)):
        stop = min(len(ids), start + chunk_size)
        if _validate_chunk(stage_directory, stage, chunk_index, start, stop,
                           ids[start:stop], config_hash, storage) is None:
            raise ValueError(f"{stage} chunk {chunk_index} failed final validation")
    if merge:
        output_path, merged = _merge_chunks(
            stage_directory, stage, ids, chunk_size, config_hash)
    else:
        output_path = stage_directory / "stage-manifest.json"
        chunk_hashes = []
        for chunk_index in range(total_chunks):
            _, manifest_path = _chunk_paths(stage_directory, chunk_index)
            chunk_hashes.append(json.loads(
                manifest_path.read_text(encoding="utf-8"))["result_sha256"])
        merged = {
            "format": "aidd-gaussian-sharded-stage", "version": 1,
            "stage": stage, "config_hash": config_hash, "chunks": total_chunks,
            "candidates": len(ids), "ids_sha256": _ids_sha256(ids),
            "storage": storage, "logical_bytes_per_row": 20 if storage == "slim" else None,
            "chunk_result_hashes_sha256": hashlib.sha256(
                "".join(chunk_hashes).encode()).hexdigest(),
        }
    merged["chunks_reused"] = reused
    merged["chunks_computed"] = len(chunks)
    merged["workers"] = workers
    merged["wall_seconds_this_invocation"] = time.perf_counter() - started
    if merge:
        _atomic_json(output_path.with_suffix(".manifest.json"), merged)
    else:
        _atomic_json(output_path, merged)
    return output_path, merged


def _select_refinement_candidates_streaming(
        coarse_directory: Path, ids: np.ndarray, chunk_size: int,
        output_path: Path, top_n_per_objective: int,
        config_hash: str) -> tuple[np.ndarray, dict]:
    """Select exact stable per-objective Top-N without a full coarse merge."""
    limit = min(top_n_per_objective, len(ids))
    heaps: list[list[tuple[float, int, int]]] = [[] for _ in OBJECTIVES]
    chunk_count = (len(ids) + chunk_size - 1) // chunk_size
    chunk_hashes = []
    for chunk_index in range(chunk_count):
        result_path, manifest_path = _chunk_paths(coarse_directory, chunk_index)
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        chunk_hashes.append(manifest["result_sha256"])
        with np.load(result_path, allow_pickle=False) as archive:
            chunk_ids = archive["global_ids"]
            values = archive["objective_scores"]
        start = chunk_index * chunk_size
        if values.shape != (len(chunk_ids), len(OBJECTIVES)):
            raise ValueError(f"invalid slim coarse score shape in chunk {chunk_index}")
        for local, gid in enumerate(chunk_ids):
            position = start + local
            for column in range(len(OBJECTIVES)):
                item = (float(values[local, column]), -position, int(gid))
                if len(heaps[column]) < limit:
                    heapq.heappush(heaps[column], item)
                elif item > heaps[column][0]:
                    heapq.heapreplace(heaps[column], item)
    selected_records: dict[int, dict] = {}
    for column, heap in enumerate(heaps):
        ordered = sorted(heap, key=lambda item: (-item[0], -item[1]))
        for rank, (_, negative_position, gid) in enumerate(ordered):
            position = -negative_position
            record = selected_records.setdefault(position, {
                "gid": gid, "selected": [False] * len(OBJECTIVES),
                "ranks": [-1] * len(OBJECTIVES), "best": len(ids) + 1})
            record["selected"][column] = True
            record["ranks"][column] = rank
            record["best"] = min(record["best"], rank)
    ordered_records = sorted(selected_records.items(), key=lambda item: (item[1]["best"], item[0]))
    positions = np.asarray([position for position, _ in ordered_records], dtype=np.int64)
    selected_ids = ids[positions]
    arrays = {
        "global_ids": selected_ids, "coarse_indices": positions,
        "selected_by_objective": np.asarray(
            [record["selected"] for _, record in ordered_records], dtype=bool),
        "coarse_ranks": np.asarray(
            [record["ranks"] for _, record in ordered_records], dtype=np.int64),
        "best_coarse_rank": np.asarray(
            [record["best"] for _, record in ordered_records], dtype=np.int64),
        "objective_names": np.asarray(OBJECTIVES, dtype="U40"),
    }
    _atomic_savez(output_path, arrays)
    manifest = {
        "format": "aidd-streaming-gaussian-refinement-selection", "version": 1,
        "config_hash": config_hash, "coarse_chunks": chunk_count,
        "coarse_chunk_hashes_sha256": hashlib.sha256(
            "".join(chunk_hashes).encode()).hexdigest(),
        "top_n_per_objective": top_n_per_objective,
        "selected_candidates": len(selected_ids),
        "maximum_candidates": min(len(ids), top_n_per_objective * len(OBJECTIVES)),
        "ids_sha256": _ids_sha256(selected_ids), "selection": output_path.name,
        "selection_sha256": _sha256(output_path),
        "policy": "streaming exact union of stable per-objective Top-N; compute allocation only",
    }
    _atomic_json(output_path.with_suffix(".manifest.json"), manifest)
    return selected_ids, manifest


def _select_refinement_candidates(coarse_path: Path, output_path: Path,
                                  top_n_per_objective: int,
                                  config_hash: str) -> tuple[np.ndarray, dict]:
    with np.load(coarse_path, allow_pickle=False) as coarse:
        ids = coarse["global_ids"].copy()
        objective_names = coarse["objective_names"].astype(str)
        count = len(ids)
        limit = min(top_n_per_objective, count)
        selected_by = np.zeros((count, len(objective_names)), dtype=bool)
        best_rank = np.full(count, count + 1, dtype=np.int64)
        ranks = np.full((count, len(objective_names)), -1, dtype=np.int64)
        for column, name in enumerate(objective_names):
            order = np.argsort(-coarse[f"{name}__objective"], kind="stable")
            chosen = order[:limit]
            selected_by[chosen, column] = True
            ranks[order, column] = np.arange(count, dtype=np.int64)
            best_rank[chosen] = np.minimum(best_rank[chosen], ranks[chosen, column])
    indices = np.flatnonzero(selected_by.any(axis=1))
    indices = indices[np.lexsort((indices, best_rank[indices]))]
    selected_ids = ids[indices]
    arrays = {
        "global_ids": selected_ids,
        "coarse_indices": indices.astype(np.int64),
        "selected_by_objective": selected_by[indices],
        "coarse_ranks": ranks[indices],
        "best_coarse_rank": best_rank[indices],
        "objective_names": objective_names.astype("U40"),
    }
    _atomic_savez(output_path, arrays)
    manifest = {
        "format": "aidd-gaussian-refinement-selection", "version": 1,
        "config_hash": config_hash, "coarse_result": str(coarse_path.resolve()),
        "coarse_result_sha256": _sha256(coarse_path),
        "top_n_per_objective": top_n_per_objective,
        "selected_candidates": len(selected_ids),
        "maximum_candidates": min(count, top_n_per_objective * len(OBJECTIVES)),
        "ids_sha256": _ids_sha256(selected_ids), "selection": output_path.name,
        "selection_sha256": _sha256(output_path),
        "policy": "union of stable per-objective Top-N; compute allocation only",
    }
    _atomic_json(output_path.with_suffix(".manifest.json"), manifest)
    return selected_ids, manifest


def run_staged_gaussian_reranking(artifact_catalog: Path, query_path: Path,
                                  candidate_ids_path: Path, output_dir: Path, *,
                                  stage: str = "all", workers: int = 16,
                                  chunk_size: int = 1000,
                                  top_n_per_objective: int = 5000,
                                  sigma: float = 1.0, cutoff: float | None = 4.5,
                                  pair_tolerance: float = 2.0,
                                  axial_samples: int = 6,
                                  max_pair_seeds: int = 512,
                                  resume: bool = True,
                                  progress_every: int = 1) -> dict:
    """Run all-candidate PCA coarse scoring then Top-N-union pair refinement."""
    if stage not in {"coarse", "refine", "all"}:
        raise ValueError("stage must be coarse, refine, or all")
    if workers <= 0 or chunk_size <= 0 or top_n_per_objective <= 0:
        raise ValueError("workers, chunk_size, and top_n_per_objective must be positive")
    if max_pair_seeds <= 0:
        raise ValueError("staged refinement requires max_pair_seeds to be positive")
    if progress_every < 0:
        raise ValueError("progress_every must be non-negative")
    artifact_catalog = artifact_catalog.resolve()
    query_path = query_path.resolve()
    candidate_ids_path = candidate_ids_path.resolve()
    output_dir = output_dir.resolve()
    ids = load_candidate_ids(candidate_ids_path)
    config = {
        "artifact_catalog": str(artifact_catalog),
        "artifact_catalog_sha256": _sha256(artifact_catalog),
        "query": str(query_path), "query_sha256": _sha256(query_path),
        "candidate_ids": str(candidate_ids_path),
        "candidate_ids_sha256": _sha256(candidate_ids_path),
        "candidate_count": len(ids), "candidate_ids_sha256_canonical": _ids_sha256(ids),
        "chunk_size": chunk_size, "top_n_per_objective": top_n_per_objective,
        "sigma_angstrom": sigma, "cutoff_angstrom": cutoff,
        "pair_tolerance_angstrom": pair_tolerance,
        "axial_samples": axial_samples, "max_pair_seeds": max_pair_seeds,
        "objectives": list(OBJECTIVES),
    }
    canonical = json.dumps(config, sort_keys=True, separators=(",", ":")).encode()
    config_hash = hashlib.sha256(canonical).hexdigest()
    output_dir.mkdir(parents=True, exist_ok=True)
    run_manifest_path = output_dir / "run-manifest.json"
    if run_manifest_path.is_file():
        run_manifest = json.loads(run_manifest_path.read_text(encoding="utf-8"))
        if (run_manifest.get("format") != STAGED_FORMAT
                or run_manifest.get("config_hash") != config_hash):
            raise ValueError(
                "existing staged Gaussian run has different inputs or parameters; "
                "use a new output directory")
    else:
        run_manifest = {
            "format": STAGED_FORMAT, "version": STAGED_VERSION,
            "config_hash": config_hash, "config": config,
            "status": "initialized", "stages": {},
        }
        _atomic_json(run_manifest_path, run_manifest)

    coarse_parameters = {
        "sigma": sigma, "cutoff": cutoff, "pair_tolerance": pair_tolerance,
        "axial_samples": axial_samples, "max_pair_seeds": 0,
    }
    refine_parameters = dict(coarse_parameters, max_pair_seeds=max_pair_seeds)
    coarse_dir = output_dir / "coarse"
    coarse_path = coarse_dir / "merged-scores.npz"
    if stage in {"coarse", "all"}:
        coarse_path, coarse_manifest = _execute_stage(
            coarse_dir, "coarse", ids, artifact_catalog=artifact_catalog,
            query_path=query_path, parameters=coarse_parameters,
            config_hash=config_hash, workers=workers, chunk_size=chunk_size,
            resume=resume, progress_every=progress_every)
        run_manifest["stages"]["coarse"] = coarse_manifest
        run_manifest["status"] = "coarse_complete"
        _atomic_json(run_manifest_path, run_manifest)
    elif not coarse_path.is_file():
        raise ValueError("refine stage requires a completed coarse/merged-scores.npz")

    if stage == "coarse":
        return run_manifest

    selection_path = output_dir / "refine" / "selected-candidates.npz"
    selection_path.parent.mkdir(parents=True, exist_ok=True)
    selected_ids, selection_manifest = _select_refinement_candidates(
        coarse_path, selection_path, top_n_per_objective, config_hash)
    run_manifest["stages"]["selection"] = selection_manifest
    _atomic_json(run_manifest_path, run_manifest)
    refine_path, refine_manifest = _execute_stage(
        output_dir / "refine", "refine", selected_ids,
        artifact_catalog=artifact_catalog, query_path=query_path,
        parameters=refine_parameters, config_hash=config_hash, workers=workers,
        chunk_size=chunk_size, resume=resume, progress_every=progress_every)
    run_manifest["stages"]["refine"] = refine_manifest
    run_manifest["status"] = "complete"
    run_manifest["final_result"] = str(refine_path)
    run_manifest["final_result_sha256"] = _sha256(refine_path)
    _atomic_json(run_manifest_path, run_manifest)
    return run_manifest


def run_scaled_gaussian_reranking(artifact_catalog: Path, query_path: Path,
                                  candidate_schedule: Path, output_dir: Path, *,
                                  stage: str = "all", workers: int = 16,
                                  coarse_chunk_size: int = 2000,
                                  refine_chunk_size: int = 250,
                                  top_n_per_objective: int = 5000,
                                  sigma: float = 1.0,
                                  cutoff: float | None = 4.5,
                                  pair_tolerance: float = 2.0,
                                  axial_samples: int = 6,
                                  max_pair_seeds: int = 512,
                                  resume: bool = True,
                                  progress_every: int = 1) -> dict:
    """Run fixed-budget Gaussian scoring with retained slim coarse shards."""
    if stage not in {"coarse", "refine", "all"}:
        raise ValueError("stage must be coarse, refine, or all")
    if min(workers, coarse_chunk_size, refine_chunk_size,
           top_n_per_objective, max_pair_seeds) <= 0:
        raise ValueError("workers, chunk sizes, Top-N, and pair cap must be positive")
    artifact_catalog = artifact_catalog.resolve(); query_path = query_path.resolve()
    candidate_schedule = candidate_schedule.resolve(); output_dir = output_dir.resolve()
    ids = load_candidate_ids(candidate_schedule)
    config = {
        "artifact_catalog": str(artifact_catalog),
        "artifact_catalog_sha256": _sha256(artifact_catalog),
        "query": str(query_path), "query_sha256": _sha256(query_path),
        "candidate_schedule": str(candidate_schedule),
        "candidate_schedule_sha256": _sha256(candidate_schedule),
        "candidate_count": len(ids), "candidate_ids_sha256": _ids_sha256(ids),
        "coarse_chunk_size": coarse_chunk_size,
        "refine_chunk_size": refine_chunk_size,
        "top_n_per_objective": top_n_per_objective,
        "sigma_angstrom": sigma, "cutoff_angstrom": cutoff,
        "pair_tolerance_angstrom": pair_tolerance,
        "axial_samples": axial_samples, "max_pair_seeds": max_pair_seeds,
        "coarse_storage": "slim-sharded", "objectives": list(OBJECTIVES),
    }
    canonical = json.dumps(config, sort_keys=True, separators=(",", ":")).encode()
    config_hash = hashlib.sha256(canonical).hexdigest()
    output_dir.mkdir(parents=True, exist_ok=True)
    run_manifest_path = output_dir / "run-manifest.json"
    if run_manifest_path.is_file():
        run_manifest = json.loads(run_manifest_path.read_text(encoding="utf-8"))
        if (run_manifest.get("format") != SCALED_FORMAT
                or run_manifest.get("config_hash") != config_hash):
            raise ValueError("existing scaled Gaussian run has different inputs or parameters")
    else:
        run_manifest = {
            "format": SCALED_FORMAT, "version": SCALED_VERSION,
            "config_hash": config_hash, "config": config,
            "status": "initialized", "stages": {},
            "scale_claim": "architecture only; physical billion-scale validation pending",
        }
        _atomic_json(run_manifest_path, run_manifest)
    common = {"sigma": sigma, "cutoff": cutoff,
              "pair_tolerance": pair_tolerance, "axial_samples": axial_samples}
    coarse_dir = output_dir / "coarse"
    if stage in {"coarse", "all"}:
        _, coarse_manifest = _execute_stage(
            coarse_dir, "coarse", ids, artifact_catalog=artifact_catalog,
            query_path=query_path, parameters=dict(common, max_pair_seeds=0),
            config_hash=config_hash, workers=workers, chunk_size=coarse_chunk_size,
            resume=resume, progress_every=progress_every, storage="slim", merge=False)
        run_manifest["stages"]["coarse"] = coarse_manifest
        run_manifest["status"] = "coarse_complete"
        _atomic_json(run_manifest_path, run_manifest)
    elif not (coarse_dir / "stage-manifest.json").is_file():
        raise ValueError("refine stage requires completed slim coarse chunks")
    if stage == "coarse":
        return run_manifest
    for chunk_index, start in enumerate(range(0, len(ids), coarse_chunk_size)):
        stop = min(len(ids), start + coarse_chunk_size)
        if _validate_chunk(
                coarse_dir, "coarse", chunk_index, start, stop, ids[start:stop],
                config_hash, "slim") is None:
            raise ValueError(
                "invalid slim coarse chunk; rerun with --stage coarse to repair it")
    selection_path = output_dir / "refine" / "selected-candidates.npz"
    selection_path.parent.mkdir(parents=True, exist_ok=True)
    selected_ids, selection_manifest = _select_refinement_candidates_streaming(
        coarse_dir, ids, coarse_chunk_size, selection_path,
        top_n_per_objective, config_hash)
    run_manifest["stages"]["selection"] = selection_manifest
    _atomic_json(run_manifest_path, run_manifest)
    refine_path, refine_manifest = _execute_stage(
        output_dir / "refine", "refine", selected_ids,
        artifact_catalog=artifact_catalog, query_path=query_path,
        parameters=dict(common, max_pair_seeds=max_pair_seeds),
        config_hash=config_hash, workers=workers, chunk_size=refine_chunk_size,
        resume=resume, progress_every=progress_every)
    run_manifest["stages"]["refine"] = refine_manifest
    run_manifest["status"] = "complete"
    run_manifest["coarse_result"] = str(coarse_dir.resolve())
    run_manifest["final_result"] = str(refine_path)
    run_manifest["final_result_sha256"] = _sha256(refine_path)
    _atomic_json(run_manifest_path, run_manifest)
    return run_manifest
