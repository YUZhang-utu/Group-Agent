from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import shutil
import sys
from types import SimpleNamespace

import numpy as np
import pytest

from aidd_agent.chunk_execution import bounded_results
from aidd_agent.e035_validation import scale_validation
from aidd_agent.fast_3d_search import compare_archives, indexed_candidates
from aidd_agent.gaussian_batch import run_scaled_gaussian_reranking
from aidd_agent import library_acceptance as ev
from test_expanded_wee1 import real_pose_fixture


def test_bounded_submission_and_failure():
    produced = []
    def tasks():
        for i in range(15):
            produced.append(i)
            yield i
    with ThreadPoolExecutor(2) as pool:
        stream = bounded_results(pool, lambda x: x*x, tasks(), 4)
        first = next(stream)
        assert len(produced) == 4
        assert sorted([first, *stream]) == [(i, i*i) for i in range(15)]
        def fail(x):
            raise ValueError("worker error")
        with pytest.raises(ValueError, match="worker error"):
            list(bounded_results(pool, fail, range(10), 2))


def test_spawn_stress_matches_serial_and_validates_resume(tmp_path):
    root, batch, query, _ = real_pose_fixture(tmp_path)
    source = tmp_path / "e034"
    for label in ("8bju", "1x8b"):
        shutil.copytree(query, source / label / "query")
    protocol = tmp_path / "protocol.json"
    ev.write(protocol, dict(test="parallel stress"))
    serial = scale_validation(batch, source, tmp_path / "serial", 3, protocol, workers=1, chunk_size=1)
    parallel = scale_validation(batch, source, tmp_path / "parallel", 3, protocol, workers=2, chunk_size=1)
    for key in ("global_rank_checks", "max_absolute_error", "assignment_mismatches", "sampled_ids_sha256"):
        assert serial[key] == parallel[key]
    for i in range(3):
        assert compare_archives(tmp_path / "serial" / f"chunk-{i:06d}.npz",
                                tmp_path / "parallel" / f"chunk-{i:06d}.npz")["passed"]
    assert parallel["execution"]["computed_chunks"] == 3
    # Interrupt-like partial resume: one receipt/data removed, others retained.
    (tmp_path / "parallel/chunk-000001.stage.json").unlink()
    (tmp_path / "parallel/chunk-000001.npz").unlink()
    resumed = scale_validation(batch, source, tmp_path / "parallel", 3, protocol, workers=2, chunk_size=1)
    assert resumed["execution"]["reused_chunks"] == 2
    assert resumed["global_rank_checks"] == serial["global_rank_checks"]
    (tmp_path / "parallel/chunk-000000.npz").write_bytes(b"corrupt")
    with pytest.raises(ValueError, match="changed completed output"):
        scale_validation(batch, source, tmp_path / "parallel", 3, protocol, workers=2, chunk_size=1)


def test_gaussian_small_parallel_chunks_preserve_all_arrays(tmp_path):
    root, batch, query, candidates = real_pose_fixture(tmp_path)
    args = (batch / "artifacts/catalog.json", query / "gaussian-query.npz", candidates)
    serial = run_scaled_gaussian_reranking(*args, root / "serial", workers=1, coarse_chunk_size=3, refine_chunk_size=3)
    parallel = run_scaled_gaussian_reranking(*args, root / "parallel", workers=2, coarse_chunk_size=1, refine_chunk_size=1)
    assert compare_archives(Path(serial["final_result"]), Path(parallel["final_result"]))["passed"]
    assert parallel["stages"]["refine"]["chunks_computed"] == 3


def test_indexed_candidates_fetches_only_found_ids_and_reranks(tmp_path):
    raw = np.zeros((3, 60), dtype=np.float32)
    raw[:, 0] = [2, 0, 1]
    class Index:
        def search(self, query, k):
            assert k == 3
            return np.zeros((1, 3)), np.asarray([[0, 2, 1]])
    def fetch(ids, vectors):
        assert vectors and ids.tolist() == [0, 2, 1]
        return np.asarray(["m0", "m2", "m1"]), raw[ids]
    output = tmp_path / "candidates.npz"
    result = indexed_candidates(Index(), SimpleNamespace(total=3, fetch=fetch),
                                np.zeros(60, dtype=np.float32), np.zeros(60, dtype=np.float32),
                                np.ones(60, dtype=np.float32), output)
    assert result["exact_reference_scan"] is False
    with np.load(output) as a:
        assert a["global_ids"].tolist() == [1, 2, 0]
        assert a["squared_l2"].tolist() == [0, 1, 4]


def test_e036_full_indexed_pipeline_and_completed_resume(tmp_path, monkeypatch):
    from aidd_agent import fast_3d_search as fast
    from aidd_agent.expanded_wee1 import run_pose_stages, fingerprint
    from aidd_agent.gaussian_batch import _load_query, write_gaussian_query
    root, batch, original_query, _ = real_pose_fixture(tmp_path)
    source = tmp_path / "e034"; source.mkdir()
    (source / "retrieval").mkdir()
    (batch / "faiss").mkdir()
    (batch / "run.lock").touch()
    index_file = batch / "faiss/index.faiss"; index_file.write_bytes(b"fixture index")
    ev.write(batch / "faiss/manifest.json", dict(index_sha256=ev.sha(index_file)))
    np.savez(batch / "faiss/transform.npz", mean=np.zeros(60, dtype=np.float32), std=np.ones(60, dtype=np.float32))
    class Index:
        ntotal = 3
        d = 60
        nlist = 128
        def search(self, query, k):
            return np.zeros((1, k)), np.arange(k)[None]
    class Corpus:
        total = 3
        def __init__(self, catalog): pass
        def fetch(self, ids, vectors):
            return np.asarray([f"MOL-{i}" for i in ids]), np.zeros((len(ids), 60), dtype=np.float32)
    monkeypatch.setattr(ev, "Corpus", Corpus)
    monkeypatch.setitem(sys.modules, "faiss", SimpleNamespace(read_index=lambda _: Index(), omp_set_num_threads=lambda _: None))
    monkeypatch.setitem(sys.modules, "fcntl", SimpleNamespace(LOCK_EX=1, LOCK_NB=2, flock=lambda *_: None))
    monkeypatch.setattr(fast, "ensure_file_descriptor_limit", lambda _: None)
    monkeypatch.setattr(fast, "runtime_versions", lambda: dict(test="fixture"))
    report = dict(status="complete", poses={})
    for label, qid in fast.QUERIES:
        qdir = source / label / "query"; qdir.mkdir(parents=True)
        _, qa = _load_query(original_query / "gaussian-query.npz")
        write_gaussian_query(qdir / "gaussian-query.npz", **qa, source=dict(query_id=qid))
        (qdir / (qid.split(":")[0]+".cif")).write_text("data_test\n")
        descriptor = qdir / "usrcat.npy"; np.save(descriptor, np.zeros(60, dtype=np.float32))
        ev.write(source / f"{label}-query.stage.json", dict(outputs=fingerprint([descriptor])))
        schedule = source / "retrieval" / f"{label}-np128-candidates.npz"
        indexed_candidates(Index(), Corpus(None), np.zeros(60, dtype=np.float32), np.zeros(60, dtype=np.float32),
                           np.ones(60, dtype=np.float32), schedule)
        report["poses"][qid] = run_pose_stages(source / label, batch, qdir, schedule, workers=1)
    ev.write(source / "report.json", report)
    ev.write(source / "EXECUTION_COMPLETE.json", dict(status="complete", report_sha256=ev.sha(source / "report.json")))
    ev.write(source / "protocol.json", dict(sources=fingerprint([batch / n for n in
        ("artifacts/catalog.json", "faiss/manifest.json", "faiss/transform.npz")])))
    before = fingerprint(p for p in source.rglob("*") if p.is_file())
    args = SimpleNamespace(batch=batch, e034=source, output=tmp_path / "e036", resume=False,
                           workers=2, coarse_chunk=1, refine_chunk=1, fresh_retrieval=True)
    fast.run(args)
    result = ev.read(args.output / "report.json")
    assert result["status"] == "complete"
    assert all(r["fresh_compute"] and r["gaussian_equivalence"]["passed"] for r in result["queries"])
    assert all(r["annotation_equivalence"]["passed"] for r in result["queries"])
    args.resume = True
    def never(*_, **__): pytest.fail("Completed query recomputed")
    monkeypatch.setattr(fast, "run_scaled_gaussian_reranking", never)
    fast.run(args)
    assert before == fingerprint(p for p in source.rglob("*") if p.is_file())
    (args.output / "8bju/interaction-matches.npz").write_bytes(b"bad")
    with pytest.raises(ValueError, match="changed completed output"):
        fast.run(args)
