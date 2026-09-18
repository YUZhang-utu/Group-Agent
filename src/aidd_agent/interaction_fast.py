"""Opt-in batched geometry with the unchanged deterministic assignment solver."""
from __future__ import annotations

from bisect import bisect_right
from pathlib import Path
from types import SimpleNamespace
import json

import numpy as np

from .chemical_companion import CHEM_META_DTYPE, FEATURE_DIRECTION_DTYPE
from .conformer_artifacts import META_DTYPE, FEATURE_DTYPE


class InteractionFeatureReader:
    """Map only the fields E031 consumes, without topology/torsion reconstruction."""

    def __init__(self, artifact, chemical):
        self.artifact = Path(artifact).resolve(); self.chemical = Path(chemical).resolve()
        load = lambda p: json.loads(p.read_text(encoding="utf-8"))
        a, c = load(self.artifact), load(self.chemical)
        self.rows = sorted(a["shards"], key=lambda r: r["global_id_start"])
        self.peers = {r["name"]: r for r in c["shards"]}
        if len(self.peers) != len(c["shards"]) or set(self.peers) != {r["name"] for r in self.rows}:
            raise ValueError("interaction catalog shard mismatch")
        self.starts = [r["global_id_start"] for r in self.rows]
        if len(set(self.starts)) != len(self.starts):
            raise ValueError("duplicate interaction shard ranges")
        self.cache = {}; self.maps = []

    def close(self):
        self.cache.clear()
        for array in self.maps:
            array._mmap.close()
        self.maps.clear()

    def _map(self, path, dtype):
        if path.stat().st_size == 0:
            return np.empty(0, dtype=dtype)
        array = np.memmap(path, mode="r", dtype=dtype)
        self.maps.append(array)
        return array

    def _open(self, row):
        name = row["name"]
        if name in self.cache:
            return self.cache[name]
        peer = self.peers[name]
        a, c = Path(row["path"]), Path(peer["path"])
        if not a.is_dir(): a = self.artifact.parent / name
        if not c.is_dir(): c = self.chemical.parent / name
        am = json.loads((a / "manifest.json").read_text(encoding="utf-8"))
        cm = json.loads((c / "manifest.json").read_text(encoding="utf-8"))
        expected = (row["global_id_start"], row["conformers"])
        if any((m["global_id_start"], m["conformers"]) != expected for m in (am, cm)):
            raise ValueError("interaction shard range mismatch")
        scale = float(am["coordinate_scale"])
        if not np.isfinite(scale) or scale <= 0:
            raise ValueError("invalid interaction coordinate scale")
        arrays = dict(meta=self._map(a / "meta.bin", META_DTYPE), features=self._map(a / "feats.bin", FEATURE_DTYPE),
                      chemical_meta=self._map(c / "chem-meta.bin", CHEM_META_DTYPE),
                      directions=self._map(c / "feature-directions.bin", FEATURE_DIRECTION_DTYPE),
                      mol=self._map(a / "molecule_ids.bin", "S16"), conf=self._map(a / "conformer_ids.bin", "S16"),
                      cmol=self._map(c / "molecule_ids.bin", "S16"), cconf=self._map(c / "conformer_ids.bin", "S16"))
        if any(len(arrays[k]) != expected[1] for k in ("meta", "chemical_meta", "mol", "conf", "cmol", "cconf")):
            raise ValueError("interaction shard row count mismatch")
        self.cache[name] = (scale, arrays)
        return scale, arrays

    def get(self, gid):
        if isinstance(gid, (bool, np.bool_)) or int(gid) != gid:
            raise ValueError("global ID must be an integer")
        gid = int(gid); index = bisect_right(self.starts, gid) - 1
        if index < 0 or gid >= self.starts[index] + self.rows[index]["conformers"]:
            raise IndexError("global ID outside interaction catalog")
        scale, a = self._open(self.rows[index]); local = gid - self.starts[index]
        m, c = a["meta"][local], a["chemical_meta"][local]
        if (m["global_id"] != gid or c["global_id"] != gid or m["features"] != c["features"]
                or a["mol"][local] != a["cmol"][local] or a["conf"][local] != a["cconf"][local]):
            raise ValueError(f"interaction feature identity mismatch at {gid}")
        fo, co, n = int(m["feature_offset"]), int(c["feature_offset"]), int(m["features"])
        if fo + n > len(a["features"]) or co + n > len(a["directions"]):
            raise ValueError("interaction feature offset out of bounds")
        f, d = a["features"][fo:fo+n], a["directions"][co:co+n]
        return SimpleNamespace(global_id=gid, molecule_id=bytes(a["mol"][local]).decode(),
            conformer_id=bytes(a["conf"][local]).decode(),
            feature_points=np.asarray(m["origin"], dtype=np.float64) + np.asarray(f["xyz"], dtype=np.float64) / scale,
            feature_types=np.array(f["type"], dtype=np.uint8),
            feature_directions=np.array(d["direction"], dtype=np.float64), feature_kinds=np.array(d["kind"], dtype=np.uint8))


def match_batch(query, points, types, directions, kinds, transforms, *, sigma=1., cutoff=4.5, angular_power=2.):
    """Same geometry/assignment semantics as interaction_match, batched by F."""
    from .interaction_matching import interaction_match, _maximum_weight_assignment
    qp, qt, qd, qk, weights = query
    # Keep reference query and parameter validation once per bounded group.
    interaction_match(qp, qt, qd, qk, weights, np.empty((0, 3)), [], np.empty((0, 3)), [],
                      sigma=sigma, cutoff=cutoff, angular_power=angular_power)
    qp, qd, weights = [np.asarray(v, dtype=np.float64) for v in (qp, qd, weights)]
    qt, qk = np.asarray(qt), np.asarray(qk, dtype=np.uint8)
    points, directions = np.asarray(points, dtype=np.float64), np.asarray(directions, dtype=np.float64)
    types, kinds = np.asarray(types), np.asarray(kinds, dtype=np.uint8)
    matrices = np.asarray(transforms, dtype=np.float64)
    if (points.ndim != 3 or points.shape[-1] != 3 or directions.shape != points.shape
            or types.shape != points.shape[:2] or kinds.shape != types.shape
            or matrices.shape not in ((len(points), 16), (len(points), 4, 4))):
        raise ValueError("invalid batched feature/transform shape")
    matrices = matrices.reshape(-1, 4, 4)
    if (not np.isfinite(points).all() or not np.isfinite(directions).all()
            or not np.isfinite(matrices).all() or not np.allclose(matrices[:, 3], [0, 0, 0, 1])
            or not np.isin(kinds, [0, 1, 2]).all()):
        raise ValueError("invalid batched feature/transform values")
    rotation = matrices[:, :3, :3].transpose(0, 2, 1)
    cp = points @ rotation + matrices[:, None, :3, 3]
    cd = directions @ rotation
    lengths = np.linalg.norm(cd, axis=2); valid = lengths > 1e-12
    cd[valid] /= lengths[valid, None]; cd[~valid] = 0
    if not np.allclose(np.linalg.norm(cd[kinds != 0], axis=1), 1., atol=5e-3):
        raise ValueError("valid feature directions must be unit vectors")
    delta = qp[None, :, None, :] - cp[:, None, :, :]
    squared = np.einsum("bijk,bijk->bij", delta, delta)
    spatial = np.exp(-squared / (2. * sigma * sigma))
    compatible = qt[None, :, None] == types[:, None, :]
    if cutoff is not None: compatible &= squared <= cutoff * cutoff
    cosine = np.clip(qd[None] @ cd.transpose(0, 2, 1), -1., 1.)
    directional = qk[None, :, None] != 0
    same_kind = qk[None, :, None] == kinds[:, None, :]
    agreement = np.ones_like(spatial)
    agreement = np.where(directional & (qk[None, :, None] == 1), np.maximum(cosine, 0.), agreement)
    agreement = np.where(directional & (qk[None, :, None] == 2), np.abs(cosine), agreement)
    compatible &= (~directional) | same_kind
    pairs = np.where(compatible, spatial * agreement ** angular_power, 0.)
    assignment = np.empty((len(points), len(qp)), dtype=np.int32)
    anchor_scores = np.zeros((len(points), len(qp)), dtype=np.float64)
    scores = np.empty(len(points), dtype=np.float64)
    for i, values in enumerate(pairs):
        assignment[i] = _maximum_weight_assignment(values, weights)
        for j, k in enumerate(assignment[i]):
            if k >= 0: anchor_scores[i, j] = values[j, k]
        scores[i] = np.dot(weights, anchor_scores[i]) / weights.sum()
    if not np.isfinite(scores).all(): raise ValueError("nonfinite interaction scores")
    return scores, assignment, anchor_scores


def score_batched(artifact, chemical, query, ids, molecule_ids, conformer_ids, objectives, transforms,
                  *, sigma, cutoff, angular_power, chunk_size=256):
    anchors = query["anchor_feature_indices"]
    q = tuple(query[key][anchors] for key in ("feature_points", "feature_types", "feature_directions",
                                             "feature_direction_kinds", "anchored_weights"))
    scores = {k: np.empty(len(ids), dtype=np.float64) for k in objectives}
    assignments = {k: np.empty((len(ids), len(anchors)), dtype=np.int32) for k in objectives}
    anchor_scores = {k: np.empty((len(ids), len(anchors)), dtype=np.float64) for k in objectives}
    reader = InteractionFeatureReader(artifact, chemical)
    try:
        for start in range(0, len(ids), chunk_size):
            groups = {}
            for index in range(start, min(len(ids), start + chunk_size)):
                c = reader.get(ids[index])
                if c.molecule_id != molecule_ids[index] or c.conformer_id != conformer_ids[index]:
                    raise ValueError(f"rigid/interaction identity mismatch at {ids[index]}")
                groups.setdefault(len(c.feature_points), []).append((index, c))
            for rows in groups.values():
                indices = np.asarray([i for i, _ in rows]); records = [c for _, c in rows]
                arrays = [np.asarray([getattr(c, name) for c in records]) for name in
                          ("feature_points", "feature_types", "feature_directions", "feature_kinds")]
                for objective in objectives:
                    result = match_batch(q, *arrays, transforms[objective][indices],
                                         sigma=sigma, cutoff=cutoff, angular_power=angular_power)
                    scores[objective][indices], assignments[objective][indices], anchor_scores[objective][indices] = result
    finally:
        reader.close()
    return scores, assignments, anchor_scores
