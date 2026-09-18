from __future__ import annotations

import hashlib
import json
import math
import os
from pathlib import Path
import time
from typing import Sequence

import numpy as np

from .chemical_companion import ChemicalCompanionReader
from .chemical_geometry import DIRECTION_AXIAL, DIRECTION_NONE, DIRECTION_SIGNED, transform_directions
from .gaussian_batch import ArtifactCatalogReader
from .gaussian_overlay import apply_transform


FORMAT = "aidd-key-interaction-matches"
VERSION = 1


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _average_ranks(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=np.float64)
    order = np.argsort(values, kind="stable")
    ranks = np.empty(len(values), dtype=np.float64)
    start = 0
    while start < len(values):
        stop = start + 1
        while stop < len(values) and values[order[stop]] == values[order[start]]:
            stop += 1
        ranks[order[start:stop]] = 0.5 * (start + stop - 1)
        start = stop
    return ranks


def _rank_summary(baseline: np.ndarray, interaction: np.ndarray,
                  assignments: np.ndarray) -> dict:
    if not len(interaction):
        return {
            "interaction_min_median_max": [None, None, None],
            "matched_anchors_median": None,
            "spearman_rho_vs_rigid_objective": None,
            "top_k_overlap": {},
        }
    rho = None
    if len(interaction) > 1 and np.ptp(baseline) > 0 and np.ptp(interaction) > 0:
        baseline_ranks = _average_ranks(baseline)
        interaction_ranks = _average_ranks(interaction)
        rho_value = float(np.corrcoef(baseline_ranks, interaction_ranks)[0, 1])
        rho = rho_value if math.isfinite(rho_value) else None
    overlaps = {}
    for requested in (100, 500, 1000):
        k = min(requested, len(interaction))
        left = set(np.argsort(-baseline, kind="stable")[:k].tolist())
        right = set(np.argsort(-interaction, kind="stable")[:k].tolist())
        count = len(left & right)
        overlaps[str(k)] = {"count": count, "fraction": count / k}
    return {
        "interaction_min_median_max": [
            float(np.min(interaction)), float(np.median(interaction)),
            float(np.max(interaction))],
        "matched_anchors_median": float(np.median(np.sum(assignments >= 0, axis=1))),
        "spearman_rho_vs_rigid_objective": rho,
        "top_k_overlap": overlaps,
    }


def _maximum_weight_assignment(values: np.ndarray, weights: np.ndarray) -> np.ndarray:
    """Deterministic rectangular Hungarian assignment; -1 denotes unmatched."""
    scores = np.asarray(values, dtype=np.float64)
    weights = np.asarray(weights, dtype=np.float64)
    if scores.ndim != 2 or weights.shape != (scores.shape[0],):
        raise ValueError("assignment arrays have inconsistent shapes")
    rows, real_columns = scores.shape
    columns = real_columns + rows
    cost = np.zeros((rows, columns), dtype=np.float64)
    cost[:, :real_columns] = -(scores * weights[:, None])
    u = np.zeros(rows + 1)
    v = np.zeros(columns + 1)
    p = np.zeros(columns + 1, dtype=np.int64)
    way = np.zeros(columns + 1, dtype=np.int64)
    for row in range(1, rows + 1):
        p[0] = row
        column0 = 0
        minimum = np.full(columns + 1, np.inf)
        used = np.zeros(columns + 1, dtype=bool)
        while True:
            used[column0] = True
            row0 = p[column0]
            delta = np.inf
            column1 = 0
            for column in range(1, columns + 1):
                if used[column]:
                    continue
                current = cost[row0 - 1, column - 1] - u[row0] - v[column]
                if current < minimum[column]:
                    minimum[column] = current
                    way[column] = column0
                if minimum[column] < delta:
                    delta = minimum[column]
                    column1 = column
            for column in range(columns + 1):
                if used[column]:
                    u[p[column]] += delta
                    v[column] -= delta
                else:
                    minimum[column] -= delta
            column0 = column1
            if p[column0] == 0:
                break
        while True:
            column1 = way[column0]
            p[column0] = p[column1]
            column0 = column1
            if column0 == 0:
                break
    result = np.full(rows, -1, dtype=np.int32)
    for column in range(1, columns + 1):
        row = int(p[column]) - 1
        candidate = column - 1
        if row >= 0 and candidate < real_columns and scores[row, candidate] > 0:
            result[row] = candidate
    return result


def interaction_match(
        query_points: Sequence[Sequence[float]], query_types: Sequence[int],
        query_directions: Sequence[Sequence[float]], query_kinds: Sequence[int],
        query_weights: Sequence[float], candidate_points: Sequence[Sequence[float]],
        candidate_types: Sequence[int], candidate_directions: Sequence[Sequence[float]],
        candidate_kinds: Sequence[int], *, sigma: float = 1.0,
        cutoff: float | None = 4.5, angular_power: float = 2.0) -> dict:
    qp = np.asarray(query_points, dtype=np.float64)
    cp = np.asarray(candidate_points, dtype=np.float64)
    qd = np.asarray(query_directions, dtype=np.float64)
    cd = np.asarray(candidate_directions, dtype=np.float64)
    qt = np.asarray(query_types)
    ct = np.asarray(candidate_types)
    qk = np.asarray(query_kinds, dtype=np.uint8)
    ck = np.asarray(candidate_kinds, dtype=np.uint8)
    weights = np.asarray(query_weights, dtype=np.float64)
    if (qp.ndim != 2 or cp.ndim != 2 or qp.shape[1:] != (3,) or cp.shape[1:] != (3,)
            or qd.shape != qp.shape or cd.shape != cp.shape
            or not (qt.shape == qk.shape == weights.shape == (len(qp),))
            or not (ct.shape == ck.shape == (len(cp),))):
        raise ValueError("interaction feature arrays have inconsistent shapes")
    known_kinds = {DIRECTION_NONE, DIRECTION_SIGNED, DIRECTION_AXIAL}
    if (not len(qp) or sigma <= 0 or not math.isfinite(sigma)
            or angular_power <= 0 or not math.isfinite(angular_power)
            or cutoff is not None and (cutoff <= 0 or not math.isfinite(cutoff))
            or not np.isfinite(qp).all() or not np.isfinite(cp).all()
            or not np.isfinite(qd).all() or not np.isfinite(cd).all()
            or not np.isfinite(weights).all() or np.any(weights < 0)
            or weights.sum() <= 0):
        raise ValueError("interaction matching requires anchors and positive finite parameters")
    if (not set(map(int, qk)).issubset(known_kinds)
            or not set(map(int, ck)).issubset(known_kinds)):
        raise ValueError("unknown direction kind")
    for directions, kinds in ((qd, qk), (cd, ck)):
        directional_rows = kinds != DIRECTION_NONE
        if np.any(directional_rows) and not np.allclose(
                np.linalg.norm(directions[directional_rows], axis=1),
                1.0, atol=5e-3):
            raise ValueError("valid feature directions must be unit vectors")
    delta = qp[:, None, :] - cp[None, :, :]
    squared = np.einsum("ijk,ijk->ij", delta, delta)
    spatial = np.exp(-squared / (2.0 * sigma * sigma))
    compatible = qt[:, None] == ct[None, :]
    if cutoff is not None:
        compatible &= squared <= cutoff * cutoff
    cosine = np.clip(qd @ cd.T, -1.0, 1.0)
    directional = qk[:, None] != DIRECTION_NONE
    same_kind = qk[:, None] == ck[None, :]
    agreement = np.ones_like(spatial)
    agreement = np.where(
        directional & (qk[:, None] == DIRECTION_SIGNED), np.maximum(cosine, 0.0), agreement)
    agreement = np.where(
        directional & (qk[:, None] == DIRECTION_AXIAL), np.abs(cosine), agreement)
    compatible &= (~directional) | same_kind
    pair_scores = np.where(compatible, spatial * agreement ** angular_power, 0.0)
    assignments = _maximum_weight_assignment(pair_scores, weights)
    anchor_scores = np.asarray([
        0.0 if candidate < 0 else pair_scores[index, candidate]
        for index, candidate in enumerate(assignments)], dtype=np.float64)
    score = float(np.dot(weights, anchor_scores) / weights.sum())
    if not math.isfinite(score):
        raise ValueError("interaction match score is not finite")
    return {"interaction_match_score": score, "assignments": assignments,
            "anchor_scores": anchor_scores, "pair_scores": pair_scores}


def score_interaction_matches(
        artifact_catalog: Path, chemical_companion: Path, query_path: Path,
        rigid_result: Path, output_path: Path, *, sigma: float = 1.0,
        cutoff: float | None = 4.5, angular_power: float = 2.0,
        engine: str = "reference") -> dict:
    started = time.perf_counter()
    if engine not in {"reference", "batched"}:
        raise ValueError("interaction engine must be reference or batched")
    paths = [Path(value).resolve() for value in
             (artifact_catalog, chemical_companion, query_path, rigid_result)]
    artifact_catalog, chemical_companion, query_path, rigid_result = paths
    output_path = Path(output_path).resolve()
    if output_path.suffix.lower() != ".npz":
        raise ValueError("interaction match output must use .npz")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    artifact_sha256 = _sha256(artifact_catalog)
    companion_document = json.loads(chemical_companion.read_text(encoding="utf-8"))
    if companion_document.get("artifact_v1_catalog_sha256") != artifact_sha256:
        raise ValueError("chemical companion does not match artifact catalog")
    query_manifest_path = query_path.with_suffix(".manifest.json")
    if query_manifest_path.is_file():
        query_manifest = json.loads(query_manifest_path.read_text(encoding="utf-8"))
        if query_manifest.get("query_npz_sha256") != _sha256(query_path):
            raise ValueError("Gaussian query checksum mismatch")
    rigid_manifest_path = rigid_result.with_suffix(".manifest.json")
    if rigid_manifest_path.is_file():
        rigid_manifest = json.loads(rigid_manifest_path.read_text(encoding="utf-8"))
        reported_hash = rigid_manifest.get(
            "result_npz_sha256", rigid_manifest.get("result_sha256"))
        if reported_hash != _sha256(rigid_result):
            raise ValueError("rigid result checksum mismatch")
    with np.load(query_path, allow_pickle=False) as archive:
        query = {name: archive[name].copy() for name in archive.files}
    required_query = {
        "feature_points", "feature_types", "feature_directions",
        "feature_direction_kinds", "anchored_weights", "anchor_feature_indices",
    }
    if not required_query.issubset(query):
        raise ValueError("query package lacks interaction matching arrays")
    anchors = np.asarray(query["anchor_feature_indices"], dtype=np.int64)
    feature_count = len(query["feature_points"])
    if (not len(anchors) or np.any(anchors < 0)
            or np.any(anchors >= feature_count) or len(np.unique(anchors)) != len(anchors)):
        raise ValueError("query key-interaction anchor indices are invalid")
    with np.load(rigid_result, allow_pickle=False) as archive:
        ids = np.asarray(archive["global_ids"], dtype=np.int64)
        objectives = tuple(map(str, archive["objective_names"]))
        transforms = {name: np.asarray(archive[f"{name}__transform"], dtype=np.float64)
                      for name in objectives}
        rigid_scores = {name: np.asarray(
            archive[f"{name}__objective"], dtype=np.float64) for name in objectives}
        molecule_ids = np.asarray(archive["molecule_ids"]).astype("U16")
        conformer_ids = np.asarray(archive["conformer_ids"]).astype("U16")
    if (not objectives or molecule_ids.shape != ids.shape or conformer_ids.shape != ids.shape
            or len(np.unique(ids)) != len(ids)
            or any(value.shape != (len(ids), 16) for value in transforms.values())
            or any(value.shape != ids.shape or not np.isfinite(value).all()
                   for value in rigid_scores.values())):
        raise ValueError("rigid result arrays have inconsistent shapes")
    reader = ArtifactCatalogReader(artifact_catalog) if engine == "reference" else None
    chemistry_reader = ChemicalCompanionReader(chemical_companion) if engine == "reference" else None
    scores = {name: [] for name in objectives}
    assignments = {name: [] for name in objectives}
    anchor_scores = {name: [] for name in objectives}
    for index, gid in (enumerate(ids) if engine == "reference" else ()):
        candidate = reader.get(int(gid)); chemistry = chemistry_reader.get(int(gid))
        if (candidate.molecule_id != molecule_ids[index] or candidate.conformer_id != conformer_ids[index]
                or chemistry.molecule_id != candidate.molecule_id
                or chemistry.conformer_id != candidate.conformer_id
                or len(chemistry.feature_directions) != len(candidate.feature_points)
                or len(chemistry.feature_kinds) != len(candidate.feature_points)):
            raise ValueError(f"artifact/companion/rigid identity mismatch at global ID {gid}")
        for objective in objectives:
            transform = transforms[objective][index]
            result = interaction_match(
                query["feature_points"][anchors], query["feature_types"][anchors],
                query["feature_directions"][anchors], query["feature_direction_kinds"][anchors],
                query["anchored_weights"][anchors],
                apply_transform(candidate.feature_points, transform), candidate.feature_types,
                transform_directions(chemistry.feature_directions, transform),
                chemistry.feature_kinds, sigma=sigma, cutoff=cutoff,
                angular_power=angular_power)
            scores[objective].append(result["interaction_match_score"])
            assignments[objective].append(result["assignments"])
            anchor_scores[objective].append(result["anchor_scores"])
    if engine == "batched":
        from .interaction_fast import score_batched
        scores, assignments, anchor_scores = score_batched(
            artifact_catalog, chemical_companion, query, ids, molecule_ids, conformer_ids,
            objectives, transforms, sigma=sigma, cutoff=cutoff, angular_power=angular_power)
    arrays: dict[str, np.ndarray] = {
        "global_ids": ids, "molecule_ids": molecule_ids, "conformer_ids": conformer_ids,
        "objective_names": np.asarray(objectives, dtype="U40"),
        "query_anchor_feature_indices": anchors,
    }
    for objective in objectives:
        arrays[f"{objective}__interaction_match_score"] = np.asarray(
            scores[objective], dtype=np.float32)
        arrays[f"{objective}__anchor_assignments"] = np.asarray(
            assignments[objective], dtype=np.int32).reshape(-1, len(anchors))
        arrays[f"{objective}__anchor_scores"] = np.asarray(
            anchor_scores[objective], dtype=np.float32).reshape(-1, len(anchors))
    temporary = output_path.with_name(output_path.name + ".partial")
    with temporary.open("wb") as stream:
        np.savez_compressed(stream, **arrays)
        stream.flush()
        os.fsync(stream.fileno())
    os.replace(temporary, output_path)
    analysis = {
        objective: _rank_summary(
            rigid_scores[objective], arrays[f"{objective}__interaction_match_score"],
            arrays[f"{objective}__anchor_assignments"])
        for objective in objectives
    }
    manifest = {
        "format": FORMAT, "version": VERSION, "status": "complete",
        "engine": engine,
        "output": str(output_path), "output_sha256": _sha256(output_path),
        "inputs": {"artifact_catalog": {"path": str(artifact_catalog), "sha256": artifact_sha256},
                   "chemical_companion": {"path": str(chemical_companion), "sha256": _sha256(chemical_companion)},
                   "query": {"path": str(query_path), "sha256": _sha256(query_path)},
                   "rigid_result": {"path": str(rigid_result), "sha256": _sha256(rigid_result)}},
        "candidates": len(ids), "objectives": list(objectives), "anchors": len(anchors),
        "parameters": {"sigma_angstrom": sigma, "cutoff_angstrom": cutoff,
                       "angular_power": angular_power},
        "analysis": analysis,
        "wall_seconds": time.perf_counter() - started,
        "policy": "single stored rigid pose; deterministic one-to-one query-anchor coverage; no admission or rank mutation",
    }
    manifest_path = output_path.with_suffix(".manifest.json")
    manifest_temporary = manifest_path.with_name(manifest_path.name + ".partial")
    manifest_temporary.write_text(
        json.dumps(manifest, indent=2), encoding="utf-8")
    os.replace(manifest_temporary, manifest_path)
    return manifest
