from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import sys
import tempfile
import time

import numpy as np


def _faiss():
    try:
        import faiss
    except ImportError as exc:
        raise RuntimeError("faiss-cpu is required for the FAISS benchmark") from exc
    return faiss


def synthetic_usrcat(count: int, seed: int) -> np.ndarray:
    """Generate scale-heterogeneous, correlated 5x12 descriptor blocks."""
    rng = np.random.default_rng(seed)
    clusters = rng.normal(size=(256, 5, 4)).astype(np.float32)
    assignments = rng.integers(0, len(clusters), size=count)
    latent = clusters[assignments] + rng.normal(0, 0.35, size=(count, 5, 4))
    projection = rng.normal(size=(5, 4, 12)).astype(np.float32)
    values = np.einsum("nbi,bij->nbj", latent, projection)
    scales = np.asarray([5.0, 2.0, 1.0] * 4, dtype=np.float32)
    values *= scales[None, None, :]
    return np.ascontiguousarray(values.reshape(count, 60), dtype=np.float32)


def exact_neighbors(database: np.ndarray, queries: np.ndarray, k: int) -> np.ndarray:
    # Query batches keep the exact baseline bounded in RAM.
    result = []
    for query in queries:
        distances = np.einsum("ij,ij->i", database - query, database - query)
        selected = np.argpartition(distances, k - 1)[:k]
        result.append(selected[np.argsort(distances[selected], kind="stable")])
    return np.asarray(result)


def candidate_recall(truth: np.ndarray, candidates: np.ndarray) -> float:
    return float(np.mean([
        len(set(expected.tolist()).intersection(found.tolist())) / len(expected)
        for expected, found in zip(truth, candidates)
    ]))


def run_benchmark(count: int, queries: int, train_size: int, truth_k: int,
                  candidate_k: int, nlist: int, nprobe: int, seed: int,
                  output: Path) -> dict:
    faiss = _faiss()
    try:
        import psutil
    except ImportError as exc:
        raise RuntimeError("psutil is required to record benchmark memory usage") from exc
    process = psutil.Process(os.getpid())
    started = time.perf_counter()
    raw = synthetic_usrcat(count + queries, seed)
    database, query_raw = raw[:count], raw[count:]
    sample = database[np.random.default_rng(seed + 1).choice(
        count, size=min(train_size, count), replace=False)]
    mean, std = sample.mean(axis=0), sample.std(axis=0)
    std[std < 1e-6] = 1.0
    database = np.ascontiguousarray((database - mean) / std, dtype=np.float32)
    query = np.ascontiguousarray((query_raw - mean) / std, dtype=np.float32)
    exact_started = time.perf_counter()
    truth = exact_neighbors(database, query, truth_k)
    exact_seconds = time.perf_counter() - exact_started
    configurations = [(60, 20, "d60_m20"), (60, 30, "d60_m30"),
                      (64, 32, "d64_m32_zero_pad")]
    results = []
    for dimension, m, name in configurations:
        if dimension == 64:
            db = np.pad(database, ((0, 0), (0, 4)))
            q = np.pad(query, ((0, 0), (0, 4)))
            training = np.pad(sample.astype(np.float32), ((0, 0), (0, 4)))
            training[:, :60] = (training[:, :60] - mean) / std
        else:
            db, q = database, query
            training = np.ascontiguousarray((sample - mean) / std, dtype=np.float32)
        quantizer = faiss.IndexFlatL2(dimension)
        index = faiss.IndexIVFPQ(quantizer, dimension, nlist, m, 8)
        train_started = time.perf_counter(); index.train(training)
        train_seconds = time.perf_counter() - train_started
        add_started = time.perf_counter(); index.add(db)
        add_seconds = time.perf_counter() - add_started
        index.nprobe = nprobe
        query_started = time.perf_counter(); _, found = index.search(q, candidate_k)
        query_seconds = time.perf_counter() - query_started
        with tempfile.NamedTemporaryFile(suffix=".faiss", delete=False) as handle:
            index_path = Path(handle.name)
        try:
            faiss.write_index(index, str(index_path))
            index_bytes = index_path.stat().st_size
        finally:
            index_path.unlink(missing_ok=True)
        results.append({"name": name, "dimension": dimension, "m": m,
                        "code_bytes_per_vector": m, "train_seconds": train_seconds,
                        "add_seconds": add_seconds, "query_seconds": query_seconds,
                        "query_ms_each": 1000 * query_seconds / queries,
                        "candidate_recall": candidate_recall(truth, found),
                        "serialized_bytes": index_bytes,
                        "serialized_bytes_per_vector": index_bytes / count,
                        "rss_gb_after_query": process.memory_info().rss / 2**30})
        del index, db, q, training
    report = {"format": "aidd-faiss-usrcat-benchmark", "version": 1,
              "parameters": {"count": count, "queries": queries,
                             "train_size": min(train_size, count), "truth_k": truth_k,
                             "candidate_k": candidate_k, "nlist": nlist,
                             "nprobe": nprobe, "seed": seed},
              "standardization": {"method": "per-dimension-zscore",
                                  "mean": mean.tolist(), "std": std.tolist()},
              "exact_seconds": exact_seconds, "total_seconds": time.perf_counter() - started,
              "results": results}
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, indent=2), encoding="utf-8")
    return report


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--count", type=int, default=200_000)
    parser.add_argument("--queries", type=int, default=20)
    parser.add_argument("--train-size", type=int, default=50_000)
    parser.add_argument("--truth-k", type=int, default=100)
    parser.add_argument("--candidate-k", type=int, default=1_000)
    parser.add_argument("--nlist", type=int, default=512)
    parser.add_argument("--nprobe", type=int, default=32)
    parser.add_argument("--seed", type=int, default=20260902)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)
    report = run_benchmark(**vars(args))
    print(json.dumps({"parameters": report["parameters"],
                      "exact_seconds": report["exact_seconds"],
                      "total_seconds": report["total_seconds"],
                      "results": report["results"]}, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
