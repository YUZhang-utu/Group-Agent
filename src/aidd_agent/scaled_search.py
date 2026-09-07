from __future__ import annotations

from dataclasses import dataclass
import hashlib
import heapq
import json
import os
from pathlib import Path
from typing import Mapping, Sequence

import numpy as np


RANKED_FORMAT = "aidd-ranked-shard-merge"
SCHEDULE_FORMAT = "aidd-tiered-candidate-schedule"
CHANNELS = ("baseline", "strict", "balanced", "loose")
CHANNEL_FLAGS = {"baseline": 1, "strict": 2, "balanced": 4, "loose": 8}


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _ranked_hash(path: Path) -> str:
    if path.is_file():
        return _sha256(path)
    digest = hashlib.sha256()
    for name in ("global_ids.npy", "scores.npy"):
        child = path / name
        digest.update(name.encode())
        digest.update(_sha256(child).encode())
    return digest.hexdigest()


def _atomic_savez(path: Path, arrays: Mapping[str, np.ndarray]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    partial = path.with_name(path.name + ".partial")
    with partial.open("wb") as stream:
        np.savez(stream, **arrays)
        stream.flush()
        os.fsync(stream.fileno())
    partial.replace(path)


def _atomic_json(path: Path, value: Mapping) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    partial = path.with_name(path.name + ".partial")
    partial.write_text(json.dumps(value, indent=2, sort_keys=True), encoding="utf-8")
    partial.replace(path)


def _load_ranked(path: Path, limit: int | None = None) -> tuple[np.ndarray, np.ndarray]:
    if path.is_dir():
        ids = np.load(path / "global_ids.npy", allow_pickle=False, mmap_mode="r")
        scores = np.load(path / "scores.npy", allow_pickle=False, mmap_mode="r")
        if ids.dtype != np.dtype("<i8") or scores.dtype not in (np.dtype("<f4"), np.dtype("<f8")):
            raise ValueError(f"ranked directory requires int64 IDs and float scores: {path}")
    else:
        with np.load(path, allow_pickle=False) as archive:
            ids = np.asarray(archive["global_ids"], dtype=np.int64)
            score_name = "scores" if "scores" in archive else "retrieval_scores"
            scores = np.asarray(archive[score_name], dtype=np.float64)
    if limit is not None:
        ids = ids[:limit]
        scores = scores[:limit]
    if ids.ndim != 1 or scores.shape != ids.shape:
        raise ValueError(f"ranked arrays have inconsistent shapes: {path}")
    if len(np.unique(ids)) != len(ids):
        raise ValueError(f"ranked input contains duplicate IDs: {path}")
    if np.any(ids < 0) or not np.all(np.isfinite(scores)):
        raise ValueError(f"ranked input contains invalid IDs or scores: {path}")
    expected = np.lexsort((ids, -scores))
    if not np.array_equal(expected, np.arange(len(ids))):
        raise ValueError(f"ranked input must be sorted by score descending then ID: {path}")
    return ids, scores


def merge_ranked_shards(shard_results: Sequence[Path], output_path: Path, *,
                        top_k: int) -> dict:
    """Exactly merge pre-ranked shard prefixes using O(top_k + shards) state."""
    if top_k <= 0 or not shard_results:
        raise ValueError("top_k and at least one shard result are required")
    paths = [Path(path).resolve() for path in shard_results]
    ranked = [_load_ranked(path) for path in paths]
    heap: list[tuple[float, int, int, int]] = []
    for shard, (ids, scores) in enumerate(ranked):
        if len(ids):
            heapq.heappush(heap, (-float(scores[0]), int(ids[0]), shard, 0))
    out_ids, out_scores, out_shards, out_ranks = [], [], [], []
    seen: set[int] = set(); duplicates = 0
    while heap and len(out_ids) < top_k:
        negative, gid, shard, rank = heapq.heappop(heap)
        if gid not in seen:
            seen.add(gid)
            out_ids.append(gid); out_scores.append(-negative)
            out_shards.append(shard); out_ranks.append(rank)
        else:
            duplicates += 1
        next_rank = rank + 1
        ids, scores = ranked[shard]
        if next_rank < len(ids):
            heapq.heappush(
                heap, (-float(scores[next_rank]), int(ids[next_rank]), shard, next_rank))
    arrays = {
        "global_ids": np.asarray(out_ids, dtype=np.int64),
        "retrieval_scores": np.asarray(out_scores, dtype=np.float64),
        "source_shard": np.asarray(out_shards, dtype=np.int32),
        "source_rank": np.asarray(out_ranks, dtype=np.int32),
    }
    output_path = Path(output_path).resolve()
    _atomic_savez(output_path, arrays)
    manifest = {
        "format": RANKED_FORMAT, "version": 1, "top_k": top_k,
        "results": len(out_ids), "duplicate_ids_merged": duplicates,
        "state_bound_rows": top_k + max(map(lambda item: len(item[0]), ranked)) + len(paths),
        "shards": [{"path": str(path), "sha256": _ranked_hash(path),
                    "available_rows": len(ranked[index][0])}
                   for index, path in enumerate(paths)],
        "result": str(output_path), "result_sha256": _sha256(output_path),
        "tie_policy": "score descending, stable global ID ascending",
        "memory_policy": "mmap shard directories; O(global_top_k + max_local_top_k + shards)",
    }
    _atomic_json(output_path.with_suffix(".manifest.json"), manifest)
    return manifest


@dataclass(frozen=True)
class CandidateBudgets:
    baseline: int = 100_000
    strict: int = 100_000
    balanced: int = 50_000
    loose: int = 10_000

    def as_dict(self) -> dict[str, int]:
        result = {name: int(getattr(self, name)) for name in CHANNELS}
        if any(value < 0 for value in result.values()) or result["baseline"] <= 0:
            raise ValueError("baseline budget must be positive and other budgets non-negative")
        return result


def build_candidate_schedule(channel_results: Mapping[str, Path], output_path: Path,
                             *, budgets: CandidateBudgets = CandidateBudgets()) -> dict:
    """Build a deterministic fixed-compute schedule without removing baseline IDs."""
    unknown = set(channel_results) - set(CHANNELS)
    if unknown or "baseline" not in channel_results:
        raise ValueError(f"baseline is required and unknown channels are invalid: {sorted(unknown)}")
    limits = budgets.as_dict()
    loaded = {name: _load_ranked(Path(path).resolve())
              for name, path in channel_results.items()}
    admitted: list[int] = []
    positions: dict[int, int] = {}
    flags: list[int] = []
    first_channel: list[int] = []
    scores: list[list[float]] = []
    ranks: list[list[int]] = []
    channel_stats = {}
    for channel_index, name in enumerate(CHANNELS):
        if name not in loaded:
            channel_stats[name] = {"available": 0, "new_admitted": 0, "duplicates_annotated": 0}
            continue
        ids, channel_scores = loaded[name]
        new_count = duplicate_count = scanned = 0
        for rank, (gid_value, score_value) in enumerate(zip(ids, channel_scores)):
            if new_count >= limits[name]:
                break
            gid = int(gid_value); scanned += 1
            if gid in positions:
                row = positions[gid]
                flags[row] |= CHANNEL_FLAGS[name]
                scores[row][channel_index] = float(score_value)
                ranks[row][channel_index] = rank
                duplicate_count += 1
                continue
            positions[gid] = len(admitted)
            admitted.append(gid); flags.append(CHANNEL_FLAGS[name])
            first_channel.append(channel_index)
            row_scores = [float("nan")] * len(CHANNELS)
            row_ranks = [-1] * len(CHANNELS)
            row_scores[channel_index] = float(score_value); row_ranks[channel_index] = rank
            scores.append(row_scores); ranks.append(row_ranks)
            new_count += 1
        channel_stats[name] = {
            "available": len(ids), "rows_scanned": scanned,
            "new_admitted": new_count, "duplicates_annotated": duplicate_count,
            "budget": limits[name],
        }
    arrays = {
        "global_ids": np.asarray(admitted, dtype=np.int64),
        "source_flags": np.asarray(flags, dtype=np.uint8),
        "first_channel": np.asarray(first_channel, dtype=np.uint8),
        "channel_scores": np.asarray(scores, dtype=np.float32).reshape(-1, len(CHANNELS)),
        "channel_ranks": np.asarray(ranks, dtype=np.int64).reshape(-1, len(CHANNELS)),
        "channel_names": np.asarray(CHANNELS, dtype="U16"),
    }
    output_path = Path(output_path).resolve()
    _atomic_savez(output_path, arrays)
    baseline_prefix = loaded["baseline"][0][:limits["baseline"]]
    if not np.array_equal(arrays["global_ids"][:len(baseline_prefix)], baseline_prefix):
        raise AssertionError("baseline prefix was not preserved")
    manifest = {
        "format": SCHEDULE_FORMAT, "version": 1,
        "budgets": limits, "maximum_admission": sum(limits.values()),
        "admitted": len(admitted), "baseline_preserved": True,
        "channels": channel_stats,
        "evidence": {name: {"path": str(Path(path).resolve()),
                            "sha256": _ranked_hash(Path(path).resolve())}
                     for name, path in channel_results.items()},
        "result": str(output_path), "result_sha256": _sha256(output_path),
        "policy": "baseline first; strict, balanced, loose add only; duplicates merge flags",
    }
    _atomic_json(output_path.with_suffix(".manifest.json"), manifest)
    return manifest
