from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
from pathlib import Path
from typing import Iterable, Mapping, Sequence

import numpy as np

from .conformer_artifacts import FEATURE_DTYPE, FEATURE_TYPES, META_DTYPE
from .gaussian_overlay import pair_alignment_seeds, principal_axis_seeds, score_overlay


QUERY_FORMAT = "aidd-gaussian-query"
RESULT_FORMAT = "aidd-gaussian-batch-scores"
SCHEMA_VERSION = 1
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


def score_gaussian_candidates(artifact_catalog: Path, query_path: Path,
                              candidate_ids_path: Path, output_path: Path, *,
                              sigma: float = 1.0, cutoff: float | None = 4.5,
                              pair_tolerance: float = 2.0, axial_samples: int = 6,
                              max_pair_seeds: int = 512) -> dict:
    if max_pair_seeds < 0:
        raise ValueError("max_pair_seeds must be non-negative")
    reader = ArtifactCatalogReader(artifact_catalog)
    query_manifest, query = _load_query(query_path.resolve())
    ids = load_candidate_ids(candidate_ids_path.resolve())
    query_shape = query["shape_points"]
    query_features = query["feature_points"]
    query_types = query["feature_types"]
    query_weights = query["anchored_weights"]
    anchor_indices = query["anchor_feature_indices"]
    seed_ids, molecule_ids, conformer_ids = [], [], []
    values = {name: {"objective": [], "shape_tanimoto": [], "shape_tversky": [],
                     "color_tanimoto": [], "color_tversky": [], "shape_cross": [],
                     "shape_query_self": [], "shape_candidate_self": [], "color_cross": [],
                     "color_query_self": [], "color_candidate_self": [], "transform": []}
              for name in OBJECTIVES}
    for gid in ids:
        candidate = reader.get(int(gid))
        molecule_ids.append(candidate.molecule_id); conformer_ids.append(candidate.conformer_id)
        seeds = list(principal_axis_seeds(candidate.shape_points, query_shape))
        if len(anchor_indices) >= 2 and len(candidate.feature_points) >= 2:
            pair_seeds = pair_alignment_seeds(
                candidate.feature_points, candidate.feature_types,
                query_features[anchor_indices], query_types[anchor_indices],
                tolerance=pair_tolerance, axial_samples=axial_samples)
            seeds.extend(pair_seeds[:max_pair_seeds])
        best = {name: None for name in OBJECTIVES}
        for seed_index, seed in enumerate(seeds):
            base_kwargs = dict(
                query_shape=query_shape, candidate_shape=candidate.shape_points,
                query_features=query_features, query_feature_types=query_types,
                candidate_features=candidate.feature_points,
                candidate_feature_types=candidate.feature_types,
                transform_matrix=seed.transform_matrix, sigma=sigma, cutoff=cutoff)
            unweighted = score_overlay(**base_kwargs)
            anchored = score_overlay(**base_kwargs, query_feature_weights=query_weights)
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
    arrays = {"global_ids": ids,
              "molecule_ids": np.asarray(molecule_ids, dtype="U16"),
              "conformer_ids": np.asarray(conformer_ids, dtype="U16"),
              "best_seed_ids": np.asarray(seed_ids, dtype="U96"),
              "objective_names": np.asarray(OBJECTIVES, dtype="U40")}
    for name, records in values.items():
        for metric, payload in records.items():
            arrays[f"{name}__{metric}"] = np.asarray(payload, dtype=np.float64)
    output_path = output_path.resolve()
    if output_path.suffix.lower() != ".npz":
        raise ValueError("Gaussian score output must use .npz")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(output_path, **arrays)
    manifest = {
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
        "query_manifest": query_manifest,
        "limitations": ["artifact-v1 has atom-centered features only",
                        "projected pharmacophore color was not scored",
                        "scores rerank and never alter L1 admission"],
    }
    output_path.with_suffix(".manifest.json").write_text(
        json.dumps(manifest, indent=2), encoding="utf-8")
    return manifest
