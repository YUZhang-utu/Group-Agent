"""E036: calibrated whole-index search and equivalent, finely scheduled 3D refinement."""
from __future__ import annotations

import argparse
import os
from pathlib import Path
import platform
import time

import numpy as np

from . import library_acceptance as ev
from .e035_validation import verify_e034, compare_sidecars
from .expanded_wee1 import fingerprint, checked_stage, ensure_file_descriptor_limit, runtime_versions
from .gaussian_batch import _atomic_json, run_scaled_gaussian_reranking
from .interaction_matching import score_interaction_matches

QUERIES = (("8bju", "8BJU:QT9:A:601"), ("1x8b", "1X8B:824:A:901"))


def compare_archives(reference, candidate):
    """Scheduling must preserve every saved score, transform, ID and row order."""
    with np.load(reference, allow_pickle=False) as left, np.load(candidate, allow_pickle=False) as right:
        ev.require(set(left.files) == set(right.files), "Gaussian archive fields changed")
        differing = [key for key in left.files if not np.array_equal(left[key], right[key])]
    return dict(passed=not differing, differing_arrays=differing, policy="all arrays exact")


def indexed_candidates(index, corpus, raw_query, mean, std, output):
    ev.require(raw_query.shape == mean.shape == std.shape == (60,), "Invalid USRCAT shape")
    ev.require(np.isfinite(raw_query).all() and np.isfinite(mean).all()
               and np.isfinite(std).all() and (std > 0).all(), "Invalid USRCAT transform")
    query = np.ascontiguousarray((raw_query-mean)/std, dtype=np.float32)
    index.nprobe = 128
    started = time.perf_counter()
    _, found = index.search(query[None], min(10000, corpus.total))
    search_seconds = time.perf_counter()-started
    ids = found[0][found[0] >= 0]
    ev.require(len(ids) == min(10000, corpus.total) and len(np.unique(ids)) == len(ids),
               "Incomplete or duplicate retrieval")
    molecules, raw = corpus.fetch(ids, vectors=True)
    delta = (raw-mean)/std-query
    distances = np.einsum("ij,ij->i", delta, delta)
    ev.require(np.isfinite(distances).all(), "Nonfinite descriptor distances")
    order = np.lexsort((ids, distances))
    np.savez(output, global_ids=ids[order], molecule_ids=molecules[order], squared_l2=distances[order])
    return dict(search_seconds=search_seconds, search_fetch_rerank_write_seconds=time.perf_counter()-started,
                candidates=len(ids), nprobe=128, exact_reference_scan=False)


def run(args):
    import fcntl
    bounded_seeds = getattr(args, "bounded_pair_seeds", False)
    profile = getattr(args, "profile_first_chunk", False)
    batch, source, output = (p.resolve() for p in (args.batch, args.e034, args.output))
    for protected in (batch, source):
        ev.require(not output.is_relative_to(protected) and not protected.is_relative_to(output),
                   "Output overlaps protected inputs")
    ev.require(min(args.workers, args.coarse_chunk, args.refine_chunk) > 0, "Positive scheduling parameters required")
    ev.require(not output.exists() or args.resume, "Output exists; use --resume or a fresh directory")
    with (batch / "run.lock").open("r") as lock:
        fcntl.flock(lock.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        setup_start = time.perf_counter()
        baseline, inputs = verify_e034(batch, source)
        catalog = ev.read(batch / "artifacts/catalog.json")
        ensure_file_descriptor_limit(len(catalog["shards"]))
        schedules = {}
        for label, qid in QUERIES:
            cfg = baseline["poses"][qid]["gaussian"]["config"]
            schedule = source / "retrieval" / f"{label}-np128-candidates.npz"
            ev.require(ev.sha(schedule) == cfg["candidate_schedule_sha256"], "Candidate lineage mismatch")
            schedules[label] = schedule
            inputs.append(schedule)
        index = corpus = mean = std = None
        if args.fresh_retrieval:
            import faiss
            original = ev.read(source / "protocol.json")["sources"]
            for name in ("artifacts/catalog.json", "faiss/manifest.json", "faiss/transform.npz"):
                path = batch / name
                ev.require(original.get(str(path)) == ev.sha(path), f"Changed E034 library evidence: {name}")
                inputs.append(path)
            ev.require(ev.sha(batch / "faiss/index.faiss") == ev.read(batch / "faiss/manifest.json")["index_sha256"],
                       "FAISS index checksum mismatch")
            faiss.omp_set_num_threads(1)
            index = faiss.read_index(str(batch / "faiss/index.faiss"))
            corpus = ev.Corpus(catalog)
            ev.require(index.ntotal == corpus.total and index.d == 60 and index.nlist >= 128,
                       "Unexpected index dimensions/count")
            with np.load(batch / "faiss/transform.npz", allow_pickle=False) as tr:
                mean, std = tr["mean"].astype(np.float32), tr["std"].astype(np.float32)
            for label, _ in QUERIES:
                p = source / label / "query/usrcat.npy"
                # The original query receipt binds this descriptor to its crystal query.
                receipt = ev.read(source / f"{label}-query.stage.json")
                ev.require(receipt["outputs"].get(str(p)) == ev.sha(p), "USRCAT query checksum mismatch")
                inputs.extend([p, source / f"{label}-query.stage.json"])
        protocol = dict(format="aidd-e036-fast-3d", sources=fingerprint(inputs),
                        code=fingerprint(sorted(Path(__file__).parent.glob("*.py"))),
                        versions=runtime_versions(), host=platform.node(), cpu_count=os.cpu_count(),
                        workers=args.workers, coarse_chunk=args.coarse_chunk, refine_chunk=args.refine_chunk,
                        fresh_retrieval=args.fresh_retrieval, bounded_pair_seeds=bounded_seeds,
                        profile_first_chunk=profile,
                        threads={k: os.environ.get(k) for k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS")})
        if output.exists():
            ev.require((output / "protocol.json").is_file()
                       and ev.read(output / "protocol.json") == protocol, "Changed E036 protocol; use fresh output")
        output.mkdir(parents=True, exist_ok=True)
        _atomic_json(output / "protocol.json", protocol)
        preflight_seconds = time.perf_counter()-setup_start
        _atomic_json(output / "RUN_STATUS.json", dict(status="running"))
        try:
            rows = []
            for label, qid in QUERIES:
                root = output / label
                root.mkdir(exist_ok=True)
                cfg = baseline["poses"][qid]["gaussian"]["config"]
                def operation():
                    started = time.perf_counter()
                    schedule = schedules[label]
                    retrieval = dict(mode="reused_E034_candidates")
                    if args.fresh_retrieval:
                        schedule = root / "candidates.npz"
                        retrieval = indexed_candidates(index, corpus, np.load(source / label / "query/usrcat.npy",
                                                       allow_pickle=False), mean, std, schedule)
                        check = compare_archives(schedules[label], schedule)
                        ev.require(check["passed"], "Fresh retrieval differs from locked E034 candidates")
                    query = source / label / "query/gaussian-query.npz"
                    gaussian_start = time.perf_counter()
                    g = run_scaled_gaussian_reranking(batch / "artifacts/catalog.json", query, schedule,
                        root / "gaussian", workers=args.workers, coarse_chunk_size=args.coarse_chunk,
                        refine_chunk_size=args.refine_chunk, top_n_per_objective=cfg["top_n_per_objective"],
                        sigma=cfg["sigma_angstrom"], cutoff=cfg["cutoff_angstrom"],
                        pair_tolerance=cfg["pair_tolerance_angstrom"], axial_samples=cfg["axial_samples"],
                        max_pair_seeds=cfg["max_pair_seeds"], resume=True,
                        bounded_pair_seeds=bounded_seeds, profile_first_chunk=profile)
                    gaussian_seconds = time.perf_counter()-gaussian_start
                    rigid = Path(g["final_result"])
                    side = root / "interaction-matches.npz"
                    annotation_start = time.perf_counter()
                    score_interaction_matches(batch / "artifacts/catalog.json", batch / "chemical/catalog.json",
                                              query, rigid, side, engine="batched")
                    annotation_seconds = time.perf_counter()-annotation_start
                    execution_seconds = time.perf_counter()-started
                    check = compare_archives(source / label / "gaussian/refine/merged-scores.npz", rigid)
                    annotation_check = compare_sidecars(source / label / "interaction-matches.npz", side)
                    ev.require(check["passed"] and annotation_check["passed"], "E036 output equivalence failed")
                    fresh = all(g["stages"][s]["chunks_reused"] == 0 for s in ("coarse", "refine"))
                    row = dict(query_id=qid, retrieval=retrieval, gaussian_seconds=gaussian_seconds,
                               annotation_seconds=annotation_seconds, execution_seconds=execution_seconds,
                               fresh_compute=fresh, latency_eligible=fresh and not profile,
                               profiling_enabled=profile, gaussian_equivalence=check, annotation_equivalence=annotation_check,
                               gaussian=g, timing_scope="query execution incl serialization and candidate check; excludes preflight and final equivalence checks",
                               service_latency_claim=False)
                    return row, sorted(p for p in root.rglob("*") if p.is_file())
                rows.append(checked_stage(output, label, [output / "protocol.json"], operation))
            ev.require(fingerprint(inputs) == protocol["sources"], "Source mutation detected")
            report = dict(status="complete", preflight_index_load_seconds=preflight_seconds, queries=rows,
                          e031_changes_ranking=False, biological_quality="not_evaluated",
                          timing_note="Completed receipts retain historical timing; partial reuse is not a fresh latency measurement")
            _atomic_json(output / "report.json", report)
            lines = ["# E036 calibrated 3D search", "", "| Query | Gaussian (s) | E031 (s) | Execution (s) | Fresh compute |", "|---|---:|---:|---:|---|"]
            lines += [f"| {r['query_id']} | {r['gaussian_seconds']:.3f} | {r['annotation_seconds']:.3f} | {r['execution_seconds']:.3f} | {r['fresh_compute']} |" for r in rows]
            lines += ["", "All Gaussian arrays and E031 scores/assignments/ranks passed equivalence checks.",
                      "Preflight/index loading and final equivalence checks are separate. No exact full-library scan.",
                      "Two calibrated WEE1 queries; annotation-only E031; no biology or service SLA claim."]
            (output / "report.md").write_text("\n".join(lines)+"\n", encoding="utf-8")
            _atomic_json(output / "RUN_STATUS.json", dict(status="complete", report_sha256=ev.sha(output / "report.json")))
        except Exception as exc:
            _atomic_json(output / "RUN_STATUS.json", dict(status="failed", error=str(exc)))
            raise


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for flag in ("batch", "e034", "output"):
        p.add_argument("--"+flag, type=Path, required=True)
    p.add_argument("--workers", type=int, default=16)
    p.add_argument("--coarse-chunk", type=int, default=500)
    p.add_argument("--refine-chunk", type=int, default=64)
    p.add_argument("--fresh-retrieval", action="store_true")
    p.add_argument("--resume", action="store_true")
    p.add_argument("--bounded-pair-seeds", action="store_true")
    p.add_argument("--profile-first-chunk", action="store_true")
    run(p.parse_args())


if __name__ == "__main__":
    main()
