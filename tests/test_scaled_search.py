from pathlib import Path

import numpy as np
import pytest

from aidd_agent.scaled_search import (
    CandidateBudgets, build_candidate_schedule, merge_ranked_shards,
)


def _ranked(path: Path, rows):
    np.savez(path, global_ids=np.asarray([row[0] for row in rows], dtype=np.int64),
             scores=np.asarray([row[1] for row in rows], dtype=np.float32))
    return path


def _ranked_directory(path: Path, rows):
    path.mkdir()
    np.save(path / "global_ids.npy",
            np.asarray([row[0] for row in rows], dtype="<i8"))
    np.save(path / "scores.npy",
            np.asarray([row[1] for row in rows], dtype="<f4"))
    return path


def test_heap_shard_merge_matches_exhaustive_and_resolves_ties(tmp_path: Path):
    shards = [
        _ranked(tmp_path / "s0.npz", [(9, .9), (1, .8), (7, .5)]),
        _ranked(tmp_path / "s1.npz", [(3, .95), (2, .8), (8, .4)]),
        _ranked(tmp_path / "s2.npz", [(4, .8), (6, .7)]),
    ]
    output = tmp_path / "merged.npz"
    manifest = merge_ranked_shards(shards, output, top_k=6)
    with np.load(output) as result:
        assert result["global_ids"].tolist() == [3, 9, 1, 2, 4, 6]
    assert manifest["state_bound_rows"] == 12
    assert manifest["results"] == 6


def test_schedule_preserves_baseline_and_keeps_tier_provenance(tmp_path: Path):
    inputs = {
        "baseline": _ranked(tmp_path / "base.npz", [(10, .9), (20, .8), (30, .7)]),
        "strict": _ranked(tmp_path / "strict.npz", [(20, .99), (40, .8), (50, .7)]),
        "balanced": _ranked(tmp_path / "balanced.npz", [(10, .95), (60, .7)]),
        "loose": _ranked(tmp_path / "loose.npz", [(70, .6), (80, .5)]),
    }
    output = tmp_path / "schedule.npz"
    manifest = build_candidate_schedule(
        inputs, output,
        budgets=CandidateBudgets(baseline=2, strict=2, balanced=1, loose=1))
    with np.load(output) as result:
        assert result["global_ids"].tolist() == [10, 20, 40, 50, 60, 70]
        assert result["source_flags"].tolist() == [5, 3, 2, 2, 4, 8]
        assert result["channel_names"].tolist() == [
            "baseline", "strict", "balanced", "loose"]
    assert manifest["baseline_preserved"] is True
    assert manifest["maximum_admission"] == 6


def test_ranked_inputs_must_be_deterministically_sorted(tmp_path: Path):
    bad = _ranked(tmp_path / "bad.npz", [(2, .5), (1, .7)])
    with pytest.raises(ValueError, match="must be sorted"):
        merge_ranked_shards([bad], tmp_path / "out.npz", top_k=1)


def test_memory_mapped_ranked_directories_are_supported(tmp_path: Path):
    left = _ranked_directory(tmp_path / "left", [(1, .9), (3, .6)])
    right = _ranked_directory(tmp_path / "right", [(2, .8), (4, .5)])
    output = tmp_path / "merged-directory-inputs.npz"
    manifest = merge_ranked_shards([left, right], output, top_k=3)
    with np.load(output) as result:
        assert result["global_ids"].tolist() == [1, 2, 3]
    assert manifest["memory_policy"].startswith("mmap shard directories")
