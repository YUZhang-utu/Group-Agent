from __future__ import annotations

from dataclasses import asdict, dataclass
import json
from pathlib import Path
from typing import Iterable, Sequence

import numpy as np


@dataclass(frozen=True)
class AnchorDefinition:
    anchor_id: str
    ligand_atom_indices: tuple[int, ...]
    feature_type: str
    atom_center: tuple[float, float, float]
    projected_points: tuple[tuple[float, float, float], ...]
    evidence: dict
    weight: float

    def __post_init__(self) -> None:
        if not self.anchor_id or not self.ligand_atom_indices:
            raise ValueError("anchor ID and ligand atom indices are required")
        required = {"interaction", "distance_angstrom", "angle_degrees", "burial"}
        if not required.issubset(self.evidence):
            raise ValueError(f"anchor evidence is missing {sorted(required - set(self.evidence))}")


@dataclass(frozen=True)
class QueryManifest:
    query_id: str
    feature_definition_version: str
    weights: dict[str, float]
    anchors: tuple[AnchorDefinition, ...]

    def __post_init__(self) -> None:
        ids = [anchor.anchor_id for anchor in self.anchors]
        if len(ids) != len(set(ids)):
            raise ValueError("anchor IDs must be unique")

    def write(self, query_directory: Path) -> Path:
        query_directory.mkdir(parents=True, exist_ok=True)
        path = query_directory / "query_manifest.json"
        path.write_text(json.dumps(asdict(self), indent=2), encoding="utf-8")
        return path


@dataclass(frozen=True)
class PoseScores:
    """All values are evaluated under this pose's named optimized transform."""

    objective_id: str
    transform_matrix: tuple[float, ...]
    shape_tanimoto: float
    shape_tversky_q: float
    color_tanimoto_unweighted: float
    color_tversky_unweighted: float
    color_tversky_anchored: float
    color_atomcentered: float
    color_projected: float
    per_anchor_overlap_raw: tuple[float, ...]
    per_anchor_self_overlap_raw: tuple[float, ...]

    def __post_init__(self) -> None:
        if not self.objective_id or len(self.transform_matrix) != 16:
            raise ValueError("a named objective and row-major 4x4 transform are required")
        if len(self.per_anchor_overlap_raw) != len(self.per_anchor_self_overlap_raw):
            raise ValueError("anchor cross/self overlap arrays must have equal lengths")


@dataclass(frozen=True)
class CandidateScores:
    """Lossless L2 record; multiple objectives have separate transforms."""

    mol_id: str
    conf_id: int
    poses: tuple[PoseScores, ...]

    def __post_init__(self) -> None:
        names = [pose.objective_id for pose in self.poses]
        if not names or len(names) != len(set(names)):
            raise ValueError("candidate pose objective IDs must be non-empty and unique")

    def pose(self, objective_id: str) -> PoseScores:
        for pose in self.poses:
            if pose.objective_id == objective_id:
                return pose
        raise KeyError(f"Candidate {self.conf_id} has no pose objective {objective_id}")

    def labelled_anchor_overlaps(self, manifest: QueryManifest,
                                 objective_id: str) -> dict[str, dict[str, float]]:
        pose = self.pose(objective_id)
        if len(pose.per_anchor_overlap_raw) != len(manifest.anchors):
            raise ValueError("candidate anchor arrays do not match query manifest")
        return {anchor.anchor_id: {"cross_overlap_raw": float(cross),
                                   "query_self_overlap_raw": float(self_overlap)}
                for anchor, cross, self_overlap in zip(
                    manifest.anchors, pose.per_anchor_overlap_raw,
                    pose.per_anchor_self_overlap_raw)}


def retain_l1_candidates(ids: Sequence[int], distances: Sequence[float], retain_n: int) -> np.ndarray:
    """Retain a broad L1 pool without consulting anchors or feature counts."""
    if retain_n <= 0:
        raise ValueError("retain_n must be positive")
    ids_array = np.asarray(ids, dtype=np.int64)
    distances_array = np.asarray(distances, dtype=np.float64)
    if ids_array.shape != distances_array.shape:
        raise ValueError("ids and distances must have equal shapes")
    valid = ids_array >= 0
    order = np.argsort(distances_array[valid], kind="stable")
    return ids_array[valid][order[:retain_n]]


def rerank(records: Iterable[CandidateScores], objective_id: str,
           score: str = "atomcentered") -> list[CandidateScores]:
    """Reorder an unchanged set using scores evaluated at one explicit pose."""
    rows = list(records)
    scorers = {
        "shape": lambda p: p.shape_tanimoto,
        "unweighted": lambda p: p.shape_tanimoto + p.color_tversky_unweighted,
        "atomcentered": lambda p: p.shape_tanimoto + p.color_atomcentered,
        "projected": lambda p: p.shape_tanimoto + p.color_projected,
        "max_color": lambda p: p.shape_tanimoto + max(p.color_atomcentered, p.color_projected),
        "anchored": lambda p: p.shape_tanimoto + p.color_tversky_anchored,
    }
    if score not in scorers:
        raise ValueError(f"Unknown reranking score: {score}")
    return sorted(rows, key=lambda row: (-scorers[score](row.pose(objective_id)), row.conf_id))


def ranking_diagnostic(unweighted_ids: Sequence[int], anchored_ids: Sequence[int], top_n: int) -> dict:
    left, right = list(unweighted_ids), list(anchored_ids)
    if len(left) != len(right) or set(left) != set(right):
        raise ValueError("rankings must contain the same candidate IDs")
    if top_n <= 0:
        raise ValueError("top_n must be positive")
    right_rank = {value: rank for rank, value in enumerate(right)}
    x = np.arange(len(left), dtype=np.float64)
    y = np.asarray([right_rank[value] for value in left], dtype=np.float64)
    rho = 1.0 if len(left) < 2 else float(np.corrcoef(x, y)[0, 1])
    n = min(top_n, len(left))
    overlap = len(set(left[:n]) & set(right[:n])) / n if n else 1.0
    return {"spearman_rho": rho, f"top_{n}_overlap": float(overlap)}
