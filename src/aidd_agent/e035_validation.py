"""E035: immutable pose review and reference-equivalent interaction benchmarking."""
from __future__ import annotations

import argparse
import cProfile
import io
import json
import os
from pathlib import Path
import platform
import pstats
import sys
import time
import multiprocessing
from concurrent.futures import ProcessPoolExecutor
from types import SimpleNamespace

import numpy as np

from . import library_acceptance as ev
from .expanded_wee1 import checked_stage, fingerprint, ensure_file_descriptor_limit
from .gaussian_batch import _atomic_json, _load_query, ArtifactCatalogReader
from .chemical_companion import ChemicalCompanionReader
from .gaussian_overlay import apply_transform
from .chemical_geometry import transform_directions
from .interaction_matching import interaction_match, score_interaction_matches
from .interaction_fast import InteractionFeatureReader, match_batch
from .interaction_review import write_review
from .chunk_execution import bounded_results


def molecule_winners(ids, molecules, scores):
    order = np.lexsort((ids, -scores))
    _, first = np.unique(molecules[order], return_index=True)
    return ids[order[np.sort(first)]]


def compare_score_arrays(ids, molecules, reference, candidate):
    if reference.shape != candidate.shape or not np.isfinite(reference).all() or not np.isfinite(candidate).all():
        return dict(passed=False, reason="shape_or_nonfinite")
    error = float(np.max(np.abs(reference.astype(float) - candidate.astype(float)), initial=0.))
    ranks = np.array_equal(np.lexsort((ids, -reference)), np.lexsort((ids, -candidate)))
    molecule_ranks = np.array_equal(molecule_winners(ids, molecules, reference), molecule_winners(ids, molecules, candidate))
    return dict(passed=error <= 1e-6 and ranks and molecule_ranks, max_absolute_error=error,
                conformer_ranking_exact=bool(ranks), molecule_representatives_and_ranking_exact=bool(molecule_ranks))


def compare_sidecars(reference, candidate):
    with np.load(reference, allow_pickle=False) as left, np.load(candidate, allow_pickle=False) as right:
        ev.require(set(left.files) == set(right.files), "Sidecar array key mismatch")
        for key in ("global_ids", "molecule_ids", "conformer_ids", "objective_names", "query_anchor_feature_indices"):
            ev.require(np.array_equal(left[key], right[key]), f"Sidecar identity/order mismatch: {key}")
        results = {}
        for objective in left["objective_names"].astype(str):
            prefix = objective + "__"
            score = compare_score_arrays(left["global_ids"], left["molecule_ids"],
                                         left[prefix + "interaction_match_score"], right[prefix + "interaction_match_score"])
            a, b = left[prefix + "anchor_scores"], right[prefix + "anchor_scores"]
            valid = a.shape == b.shape and np.isfinite(a).all() and np.isfinite(b).all()
            error = float(np.max(np.abs(a.astype(float) - b.astype(float)), initial=0.)) if valid else None
            exact_assignment = np.array_equal(left[prefix + "anchor_assignments"], right[prefix + "anchor_assignments"])
            score.update(anchor_scores_max_absolute_error=error, assignments_exact=bool(exact_assignment))
            score["passed"] = bool(score["passed"] and exact_assignment and valid and error <= 1e-6)
            results[objective] = score
        return dict(passed=all(r["passed"] for r in results.values()), objectives=results)


def verify_e034(batch, root):
    required = [root / name for name in ("report.json", "EXECUTION_COMPLETE.json", "protocol.json")]
    ev.require(all(p.is_file() for p in required), "Need original complete E034 directory, not copied report text")
    marker, report = ev.read(required[1]), ev.read(required[0])
    ev.require(marker["report_sha256"] == ev.sha(required[0]) and marker["status"] == report["status"] == "complete",
               "E034 report/completion mismatch")
    artifact, chemical = batch / "artifacts/catalog.json", batch / "chemical/catalog.json"
    inputs = list(required) + [artifact, chemical]
    for label, query_id in (("8bju", "8BJU:QT9:A:601"), ("1x8b", "1X8B:824:A:901")):
        row = report["poses"][query_id]
        side = root / label / "interaction-matches.npz"
        rigid = root / label / "gaussian/refine/merged-scores.npz"
        query = root / label / "query/gaussian-query.npz"
        qm, _ = _load_query(query)
        ev.require(qm["source"]["query_id"] == query_id, "E034 query identity mismatch")
        manifest = ev.read(side.with_suffix(".manifest.json"))
        ev.require(manifest == row["interaction"]["manifest"], "E034 embedded/sidecar manifest mismatch")
        ev.require(manifest["output_sha256"] == ev.sha(side), "E034 sidecar checksum mismatch")
        for key, path in (("artifact_catalog", artifact), ("chemical_companion", chemical), ("query", query), ("rigid_result", rigid)):
            ev.require(manifest["inputs"][key]["sha256"] == ev.sha(path), f"E034 {key} checksum mismatch")
        ev.require(row["gaussian"]["config"]["query_sha256"] == ev.sha(query)
                   and row["gaussian"]["config"]["artifact_catalog_sha256"] == ev.sha(artifact)
                   and row["gaussian"]["final_result_sha256"] == ev.sha(rigid), "E034 Gaussian lineage mismatch")
        inputs.extend([side, side.with_suffix(".manifest.json"), rigid, rigid.with_suffix(".manifest.json"),
                       query, query.with_suffix(".manifest.json"), root / label / "query" / (query_id.split(":")[0] + ".cif")])
    ev.require(ev.read(chemical)["artifact_v1_catalog_sha256"] == ev.sha(artifact), "chemical/artifact lineage mismatch")
    return report, inputs


def benchmark_query(batch, e034, output, label, query_id, report, repeats=3):
    output.mkdir(parents=True, exist_ok=True)
    artifact, chemical = batch / "artifacts/catalog.json", batch / "chemical/catalog.json"
    query_dir = e034 / label / "query"; query = query_dir / "gaussian-query.npz"
    rigid = e034 / label / "gaussian/refine/merged-scores.npz"
    original = e034 / label / "interaction-matches.npz"
    args = (artifact, chemical, query, rigid)
    before = fingerprint([*args, original])
    profiler = cProfile.Profile()
    profiler.runcall(score_interaction_matches, *args, output / "profile-reference.npz")
    profiler.dump_stats(str(output / "reference.pstats"))
    stream = io.StringIO(); pstats.Stats(profiler, stream=stream).sort_stats("cumulative").print_stats(45)
    (output / "profile.txt").write_text(stream.getvalue(), encoding="utf-8")
    timings = {"reference": [], "batched": []}; checks = []
    for repeat in range(repeats):
        for engine in (("reference", "batched") if repeat % 2 == 0 else ("batched", "reference")):
            path = output / f"{engine}-{repeat}.npz"
            start = time.perf_counter()
            score_interaction_matches(*args, path, engine=engine)
            timings[engine].append(time.perf_counter() - start)
            check = compare_sidecars(original, path); checks.append(dict(engine=engine, repeat=repeat, **check))
            ev.require(check["passed"], f"{label} {engine}: score/assignment/rank equivalence failed")
    ev.require(before == fingerprint([*args, original]), "Benchmark mutated source artifacts")
    prior = ev.read(e034 / "protocol.json")
    baseline = report["poses"][query_id]["interaction"]["timing"]
    matching_host = prior.get("host") == platform.node() and prior.get("cpu_count") == os.cpu_count()
    eligible = matching_host and baseline["status"] == "measured_same_run"
    optimized = float(np.median(timings["batched"])); reference = float(np.median(timings["reference"]))
    denominator = float(baseline["refine_wall_seconds"])
    review = write_review(artifact, chemical, query_dir, rigid, original, output / "review", query_id)
    return dict(query_id=query_id, timings=timings, reference_median_seconds=reference, optimized_median_seconds=optimized,
                speedup=reference / optimized, equivalence_passed=True, checks=checks,
                historical_refine_seconds=denominator, host_and_cpu_count_match=matching_host,
                overhead_fraction=optimized / denominator if eligible else None,
                latency_gate=optimized / denominator <= .10 if eligible else None,
                timing_caveat="historical E034 Gaussian baseline, not rerun; host/CPU count checked, load/clocks not controlled",
                review=review)


def stratified_ids(catalog, count, seed):
    rows = sorted(catalog["shards"], key=lambda r: r["global_id_start"])
    sizes = np.asarray([int(r["conformers"]) for r in rows], dtype=np.int64)
    ev.require(len(rows) <= count <= int(sizes.sum()), "Scale size must cover all shards and cannot exceed library size")
    ev.require(np.all(sizes > 0), "Empty shard")
    capacity = sizes - 1; remaining = count - len(rows)
    ideal = capacity * (remaining / max(1, int(capacity.sum())))
    take = np.floor(ideal).astype(np.int64) + 1
    missing = count - int(take.sum())
    order = np.argsort(-(ideal - np.floor(ideal)), kind="stable")
    for i in order:
        if missing and take[i] < sizes[i]: take[i] += 1; missing -= 1
    ev.require(missing == 0, "Unable to allocate stratified sample")
    rng = np.random.default_rng(seed)
    ids = np.concatenate([int(row["global_id_start"]) + np.sort(rng.choice(int(n), size=int(k), replace=False))
                          for row, n, k in zip(rows, sizes, take)])
    ev.require(len(ids) == count and len(np.unique(ids)) == count, "Nonunique scale sample")
    return ids, {r["name"]: int(k) for r, k in zip(rows, take)}


def stress_chunk(reader, ids, queries, output, *, reference_readers):
    """Real library features with constructed orientations, never called refined poses."""
    records = [reader.get(gid) for gid in ids]
    reference_records = []
    for gid, fast in zip(ids, records):
        candidate, chemistry = reference_readers[0].get(gid), reference_readers[1].get(gid)
        ev.require(candidate.molecule_id == chemistry.molecule_id == fast.molecule_id
                   and candidate.conformer_id == chemistry.conformer_id == fast.conformer_id,
                   "Scale reference/optimized reader identity mismatch")
        reference = SimpleNamespace(feature_points=candidate.feature_points, feature_types=candidate.feature_types,
                                    feature_directions=chemistry.feature_directions, feature_kinds=chemistry.feature_kinds)
        for key in ("feature_points", "feature_types", "feature_directions", "feature_kinds"):
            ev.require(np.array_equal(getattr(reference, key), getattr(fast, key)),
                       f"Scale reference/optimized feature reader mismatch: {key}")
        reference_records.append(reference)
    arrays = dict(global_ids=ids, molecule_ids=np.asarray([c.molecule_id for c in records], dtype="U16"))
    groups = {}
    for i, c in enumerate(records): groups.setdefault(len(c.feature_points), []).append(i)
    rotations = (np.eye(3), np.asarray([[0., -1., 0.], [1., 0., 0.], [0., 0., 1.]]),
                 np.asarray([[1., 0., 0.], [0., 0., -1.], [0., 1., 0.]]))
    elapsed = dict(reference=0., batched=0.); error = 0.; mismatches = 0
    for qi, q in enumerate(queries):
        center = np.mean(q[0], axis=0)
        for ri, rotation in enumerate(rotations):
            ref_scores = np.empty(len(ids), dtype=np.float32); fast_scores = np.empty(len(ids), dtype=np.float32)
            for indices_list in groups.values():
                indices = np.asarray(indices_list)
                block = [records[i] for i in indices]
                points, types, directions, kinds = [np.asarray([getattr(c, field) for c in block]) for field in
                                                  ("feature_points", "feature_types", "feature_directions", "feature_kinds")]
                matrices = np.tile(np.eye(4), (len(block), 1, 1)); matrices[:, :3, :3] = rotation
                centers = points.mean(axis=1) if points.shape[1] else np.zeros((len(block), 3))
                matrices[:, :3, 3] = center - centers @ rotation.T
                start = time.perf_counter()
                refs = [interaction_match(*q, apply_transform(c.feature_points, matrix), c.feature_types,
                                          transform_directions(c.feature_directions, matrix), c.feature_kinds)
                        for c, matrix in zip([reference_records[i] for i in indices], matrices)]
                elapsed["reference"] += time.perf_counter() - start
                start = time.perf_counter()
                scores, assignments, anchors = match_batch(q, points, types, directions, kinds, matrices)
                elapsed["batched"] += time.perf_counter() - start
                original = np.asarray([r["interaction_match_score"] for r in refs])
                original_anchors = np.asarray([r["anchor_scores"] for r in refs])
                original_assignments = np.asarray([r["assignments"] for r in refs])
                error = max(error, float(np.max(np.abs(original - scores), initial=0.)),
                            float(np.max(np.abs(original_anchors - anchors), initial=0.)))
                mismatches += int(np.sum(np.any(original_assignments != assignments, axis=1)))
                ref_scores[indices] = original; fast_scores[indices] = scores
            arrays[f"q{qi}_r{ri}_reference"] = ref_scores; arrays[f"q{qi}_r{ri}_batched"] = fast_scores
    ev.require(error <= 1e-6 and mismatches == 0, f"Stress equivalence failed: max_error={error}, assignment_mismatches={mismatches}")
    np.savez_compressed(output, **arrays)
    return dict(conformers=len(ids), comparisons=len(ids) * len(queries) * 3,
                max_absolute_error=error, assignment_mismatches=mismatches,
                scoring_seconds=elapsed, feature_count_min=min(len(c.feature_points) for c in records),
                feature_count_max=max(len(c.feature_points) for c in records))


_STRESS_WORKER = None


def _stress_initialize(artifact, chemical, queries):
    global _STRESS_WORKER
    if os.name != "nt":
        ensure_file_descriptor_limit(len(ev.read(artifact)["shards"]))
    _STRESS_WORKER = (InteractionFeatureReader(artifact, chemical), queries,
                      (ArtifactCatalogReader(artifact), ChemicalCompanionReader(chemical)))


def _stress_task(task):
    output, number, ids, protocol_path, sample = task
    reader, queries, reference_readers = _STRESS_WORKER
    path = output / f"chunk-{number:06d}.npz"
    reused = (output / f"chunk-{number:06d}.stage.json").exists()
    def operation():
        return stress_chunk(reader, ids, queries, path, reference_readers=reference_readers), [path]
    result = checked_stage(output, f"chunk-{number:06d}", [protocol_path, sample], operation)
    return result, reused


def scale_validation(batch, e034, output, count, protocol_path, workers=1, chunk_size=512):
    ev.require(workers > 0 and chunk_size > 0, "Scale workers and chunk size must be positive")
    started = time.perf_counter()
    output.mkdir(parents=True, exist_ok=True)
    artifact, chemical = batch / "artifacts/catalog.json", batch / "chemical/catalog.json"
    ids, coverage = stratified_ids(ev.read(artifact), count, 20260918)
    query_paths = [e034 / label / "query/gaussian-query.npz" for label in ("8bju", "1x8b")]
    queries = []
    for path in query_paths:
        _, query = _load_query(path); anchors = query["anchor_feature_indices"]
        queries.append(tuple(query[key][anchors] for key in ("feature_points", "feature_types", "feature_directions",
                                                             "feature_direction_kinds", "anchored_weights")))
    sample = output / "sampled-global-ids.npy"
    if sample.exists(): ev.require(np.array_equal(np.load(sample, allow_pickle=False), ids), "Changed scale sample")
    else: np.save(sample, ids)
    total = (len(ids) + chunk_size - 1) // chunk_size
    results = [None] * total
    paths = [output / f"chunk-{i:06d}.npz" for i in range(total)]
    tasks = ((output, i, ids[start:start+chunk_size], protocol_path, sample)
             for i, start in enumerate(range(0, len(ids), chunk_size)))
    completed = reused_count = 0
    checkpoint_reused = sum((output / f"chunk-{i:06d}.stage.json").exists() for i in range(total))
    compute_started = time.perf_counter()
    def record(index, value):
        nonlocal completed, reused_count
        result, reused = value
        results[index] = result
        completed += 1
        reused_count += int(reused)
        elapsed = time.perf_counter() - compute_started
        computed = completed - reused_count
        remaining = max(0, total - checkpoint_reused - computed)
        eta = elapsed / computed * remaining if computed else None
        progress = dict(completed_chunks=completed, total_chunks=total, reused_chunks=reused_count,
                        workers=workers, elapsed_seconds=elapsed, eta_seconds=eta,
                        eta_scope="chunk execution only; excludes final global rank checks")
        _atomic_json(output / "progress.json", progress)
        if completed == 1 or completed % 10 == 0 or completed == total:
            estimate = "pending" if eta is None else f"{eta / 60:.1f} min"
            ev.log(f"Scale {completed}/{total} ({100*completed/total:.1f}%); workers={workers}; "
                   f"reused={reused_count}; elapsed={elapsed/60:.1f} min; ETA={estimate}")
    if workers == 1:
        _stress_initialize(artifact, chemical, queries)
        try:
            for index, task in enumerate(tasks):
                record(index, _stress_task(task))
        finally:
            _STRESS_WORKER[0].close()
    else:
        for key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
            os.environ[key] = "1"
        with ProcessPoolExecutor(max_workers=workers, mp_context=multiprocessing.get_context("spawn"),
                                 initializer=_stress_initialize, initargs=(artifact, chemical, queries)) as pool:
            for index, value in bounded_results(pool, _stress_task, tasks, workers * 2):
                record(index, value)
    chunk_wall = time.perf_counter() - compute_started
    ev.log("Scale chunks complete; checking global conformer and molecule rankings")
    all_arrays = {}
    for path in paths:
        with np.load(path, allow_pickle=False) as a:
            for key in a.files: all_arrays.setdefault(key, []).append(a[key])
    all_arrays = {k: np.concatenate(v) for k, v in all_arrays.items()}
    ranks = {}
    for qi in range(len(queries)):
        for ri in range(3):
            prefix = f"q{qi}_r{ri}"
            ranks[prefix] = compare_score_arrays(all_arrays["global_ids"], all_arrays["molecule_ids"],
                                                all_arrays[prefix + "_reference"], all_arrays[prefix + "_batched"])
    ev.require(all(r["passed"] for r in ranks.values()), "Scale global score/rank equivalence failed")
    return dict(status="passed", unique_conformers=count, unique_molecules=len(np.unique(all_arrays["molecule_ids"])),
                execution=dict(workers=workers, chunk_size=chunk_size, chunks=total,
                               reused_chunks=reused_count, computed_chunks=total-reused_count,
                               chunk_wall_seconds_this_invocation=chunk_wall,
                               scale_wall_seconds_this_invocation=time.perf_counter()-started,
                               scoring_seconds_scope="sum of per-chunk compute; includes historical reused chunks; not elapsed wall time"),
                comparisons=sum(r["comparisons"] for r in results), sampled_ids_sha256=ev.sha(sample),
                shard_coverage=coverage, max_absolute_error=max(r["max_absolute_error"] for r in results),
                assignment_mismatches=sum(r["assignment_mismatches"] for r in results), global_rank_checks=ranks,
                reference_feature_reads_independent=True,
                scoring_seconds={engine: sum(r["scoring_seconds"][engine] for r in results) for engine in ("reference", "batched")},
                scope="real distinct library conformers, two queries x three constructed stress orientations; not refined/docking validation",
                timing_scope="geometry/assignment only; excludes common input reads, checkpoints and global ranking checks")


def render_report(report, output):
    lines = ["# E035 molecule review and equivalent E031 validation", "", f"Status: {report['status']}",
             "", "E031 remains annotation-only. No library admission or ranking was changed.", "",
             "| Query | Molecules | Top100 shared | Ref median (s) | Optimized median (s) | Speedup | Historical <=10% gate |",
             "|---|---:|---:|---:|---:|---:|---|"]
    for q in report["queries"]:
        lines.append(f"| {q['query_id']} | {q['review']['molecules']} | {q['review']['top_k_overlap']['100']['count']} | "
                     f"{q['reference_median_seconds']:.4f} | {q['optimized_median_seconds']:.4f} | {q['speedup']:.2f} | {q['latency_gate']} |")
    scale = report["scale"]
    lines += ["", f"Scale equivalence: {scale['unique_conformers']:,} distinct conformers; {scale['comparisons']:,} comparisons.",
              f"Max absolute error: {scale['max_absolute_error']}; assignment mismatches: {scale['assignment_mismatches']}.",
              "", "Scale uses real library features and constructed stress orientations, not Gaussian-refined poses.",
              "Latency includes full output serialization. Its denominator is the historical E034 refine measurement on a matching host/CPU count.",
              "Inspect query/review/molecule-ranks.csv, poses.json, SDFs and review.py. Match lines are not validated hydrogen bonds.",
              "No activity enrichment, pose-quality acceptance, or service latency guarantee is implied.", ""]
    (output / "report.md").write_text("\n".join(lines), encoding="utf-8")


def run(args):
    import fcntl
    batch, e034, output = [p.resolve() for p in (args.batch, args.e034, args.output)]
    for source in (batch, e034):
        ev.require(not output.is_relative_to(source) and not source.is_relative_to(output), "Output overlaps protected inputs")
    ev.require(args.repeats >= 3, "Need at least three alternating timing repeats")
    workers = getattr(args, "scale_workers", 1)
    chunk_size = getattr(args, "scale_chunk_size", 512)
    ev.require(workers > 0 and chunk_size > 0, "Scale workers and chunk size must be positive")
    ev.require(not output.exists() or args.resume, "Output exists; use --resume or a fresh directory")
    resume = args.resume and output.exists()
    if args.resume and not resume: ev.log("No E035 output yet; starting a fresh run")
    with (batch / "run.lock").open("r") as lock:
        fcntl.flock(lock.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        report, inputs = verify_e034(batch, e034)
        ensure_file_descriptor_limit(len(ev.read(batch / "artifacts/catalog.json")["shards"]))
        protocol = dict(format="aidd-e035", version=1, sources=fingerprint(inputs),
                        code=fingerprint(sorted(Path(__file__).parent.glob("*.py"))), numpy=np.__version__,
                        python=sys.version, platform=platform.platform(),
                        host=platform.node(), cpu_count=os.cpu_count(), repeats=args.repeats, scale=args.scale,
                        scale_workers=workers, scale_chunk_size=chunk_size,
                        threads={key: os.environ.get(key) for key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS")})
        if resume:
            ev.require((output / "protocol.json").is_file(), "E035 resume requires its protocol.json")
            ev.require(ev.read(output / "protocol.json") == protocol, "Changed E035 inputs/code/parameters/environment")
        output.mkdir(parents=True, exist_ok=resume)
        _atomic_json(output / "protocol.json", protocol)
        _atomic_json(output / "RUN_STATUS.json", dict(status="running"))
        try:
            results = []
            for label, query_id in (("8bju", "8BJU:QT9:A:601"), ("1x8b", "1X8B:824:A:901")):
                def operation(label=label, query_id=query_id):
                    result = benchmark_query(batch, e034, output / label, label, query_id, report, args.repeats)
                    return result, sorted(p for p in (output / label).rglob("*") if p.is_file())
                results.append(checked_stage(output, label, [output / "protocol.json"], operation))
            scale = scale_validation(batch, e034, output / "scale", args.scale, output / "protocol.json", workers, chunk_size)
            ev.require(protocol["sources"] == fingerprint(inputs), "E035 source mutation detected")
            result = dict(status="complete", queries=results, scale=scale, e031_changes_ranking=False,
                          biological_quality="not_evaluated", pose_review="exported_pending_human_review")
            _atomic_json(output / "report.json", result); render_report(result, output)
            _atomic_json(output / "RUN_STATUS.json", dict(status="complete", report_sha256=ev.sha(output / "report.json")))
            ev.log(f"E035 complete: {output / 'report.md'}")
        except Exception as exc:
            _atomic_json(output / "RUN_STATUS.json", dict(status="failed", error=str(exc)))
            raise


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--batch", type=Path, required=True)
    parser.add_argument("--e034", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--scale", type=int, default=100000)
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--scale-workers", type=int, default=1)
    parser.add_argument("--scale-chunk-size", type=int, default=512)
    parser.add_argument("--resume", action="store_true")
    run(parser.parse_args())


if __name__ == "__main__":
    main()
